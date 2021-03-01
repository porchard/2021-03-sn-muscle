#!/usr/bin/env nextflow

IONICE = 'ionice -c2 -n7'
GENOTYPES = params.genotypes
DEMUXLET_MASK = params.demuxlet_mask
INITIAL_THRESHOLDS_ATAC = params.initial_thresholds_atac
INITIAL_THRESHOLDS_RNA = params.initial_thresholds_rna
LIBRARY_LABELS = params.library_labels
RNA_METRICS = params.rna_metrics

genotypes = Channel.fromPath(GENOTYPES)
demuxlet_mask = Channel.fromPath(DEMUXLET_MASK)
demuxlet_unmasked_vcf_in = Channel.fromPath(GENOTYPES).map({it -> ['unmasked', it]})
call_rna_metrics_in = Channel.fromPath(RNA_METRICS)

// Generic data
AUTOSOMAL_REFERENCES = ['hg19': (1..22).collect({it -> 'chr' + it}),
			'hg38': (1..22).collect({it -> 'chr' + it}),
			'rn5': (1..20).collect({it -> 'chr' + it}),
			'rn6': (1..20).collect({it -> 'chr' + it}),
			'mm9': (1..19).collect({it -> 'chr' + it}),
			'mm10': (1..19).collect({it -> 'chr' + it})
]

ORGANISMS = ['hg19': 'human', 
		'hg38': 'human',
		'rn5': 'rat',
		'rn6': 'rat',
		'mm9': 'mouse',
		'mm10': 'mouse']

get_autosomes = {
	genome ->
	AUTOSOMAL_REFERENCES[genome]
}

get_gtf = {
	genome ->
	params.gtf[genome]
}

get_macs2_genome_size = {
	genome ->
	return MACS2_GENOME_SIZE[genome]
}

get_rnaseq_qc = {
	library ->
	return params.libraries[library].qc
}

get_ensembl = {
	genome ->
	return params.ensembl[genome]
}

is_dual_modality = {
	library ->
	return params.dual_modality.contains(library)
}

get_primary_ataqv_json = {
	library ->
	return params.libraries[library].ataqv_json
}

get_primary_counts_matrix = {
	library ->
	return params.libraries[library].counts
}

make_excluded_regions_arg = {
	genome ->
	return params.blacklist[genome].collect({'--excluded-region-file ' + it}).join(' ')
}

get_genome_size = {
	genome ->
	MACS2_GENOME_SIZE[genome]
}

get_genome = {
	library ->
	params.libraries[library].genome
}

get_tss = {
	genome ->
	params.tss[genome]
}

get_organism = {
	genome ->
	ORGANISMS[genome]
}

get_chrom_sizes = {
	genome ->
	params.chrom_sizes[genome]
}

get_gene_bed = {
	genome ->
	params.gene_bed[genome]
}

get_samples = {
	library ->
	params.libraries[library].samples
}

get_pruned = {
	library ->
	params.libraries[library].pruned
}

get_modality = {
	library ->
	params['libraries'][library]['modality']
}

get_starsolo_counts = {
	library ->
	params['libraries'][library]['starsolo'] // features, barcodes, matrix.mtx
}

libraries = params.libraries.keySet()

ATAC_LIBRARIES = []
RNA_LIBRARIES = []

for(library in libraries) {
	if (get_modality(library) == 'ATAC') {
		ATAC_LIBRARIES << library
	}
	if (get_modality(library) == 'RNA') {
		RNA_LIBRARIES << library
	}
}


get_ataqv_metrics_in = []
call_nuclei_rna_in = []
counts_to_tpm_matrix_atac_in = []
make_rna_qc_in = []
doubletfinder_counts_in = []
sort_bam_in = []
droplet_utils_in = []

DEMUXLET_LIBRARIES = ['63_20-hg19', '63_40-hg19', '63_20_rna-hg19', '63_40_rna-hg19']

for (library in ATAC_LIBRARIES) {
	get_ataqv_metrics_in << [library, file(get_primary_ataqv_json(library))]
	sort_bam_in << [library, file(get_pruned(library))]
}

for (library in RNA_LIBRARIES) {
	call_nuclei_rna_in << [library, file(get_primary_counts_matrix(library))]
	doubletfinder_counts_in << [library, file(get_primary_counts_matrix(library))]
	make_rna_qc_in << [library, file(get_rnaseq_qc(library))]
	sort_bam_in << [library, file(get_pruned(library))]
	droplet_utils_in << [library, file(params['libraries'][library]['starsolo_counts_dir'])]
}

process sort_bam {

	publishDir "${params.results}/sort-bam", mode: 'rellink'
	memory '22 GB'
	container "${params.containers.general}"
	cache 'lenient'
	maxForks 5
	cpus 10

	input:
	set val(library), file(bam) from Channel.from(sort_bam_in)
	

	output:
	set val(library), file("${library}.sorted.bam"), file("${library}.sorted.bam.bai") into demuxlet_in
	set val(library), file("${library}.sorted.bam"), file("${library}.sorted.bam.bai") into demuxlet_in_bam

	"""
	samtools sort -@ 9 -m 2G -o ${library}.sorted.bam -O BAM $bam
	samtools index ${library}.sorted.bam
	"""

}


process get_ataqv_metrics {
	
	memory '40 GB'
	container "${params.containers.general}"
	cache 'lenient'
	maxForks 1

	input:
	set val(library), file(ataqv_json) from Channel.from(get_ataqv_metrics_in)

	output:
	set val(library), file("${library}.metrics.txt") into call_nuclei_atac_in

	"""
	extractAtaqvMetric.py --files $ataqv_json --metrics tss_enrichment percent_hqaa hqaa total_reads total_autosomal_reads percent_mitochondrial percent_autosomal_duplicate percent_duplicate max_fraction_reads_from_single_autosome | perl -pe 's@.*.ataqv.json.gz\t@${library}-@' > ${library}.metrics.txt
	"""

}

process get_barcodes_for_demuxlet {

	container "${params.containers.general}"
	memory '20 GB'
	cache 'lenient'
	
	input:
	set val(library), file(bam), file(bam_index) from demuxlet_in

	output:
	file("${library}.barcodes-consider.txt") into demuxlet_run
	file("${library}.barcode-batch.*.txt") into demuxlet_barcodes

	when:
	DEMUXLET_LIBRARIES.contains(library)

	script:
	min_total = get_modality(library) == 'ATAC' ? 5000 : 250

	"""
	count-reads-per-barcode.py $bam > reads-per-barcode.txt
	cat reads-per-barcode.txt | awk '\$2>=$min_total' | cut -f1 > ${library}.barcodes-consider.txt
	split --additional-suffix=.txt --lines 2000 --suffix-length 10 ${library}.barcodes-consider.txt ${library}.barcode-batch.
	"""

}

demuxlet_barcodes_in = demuxlet_barcodes.flatten().map({it -> [it.getName().tokenize('.')[0], it]})

process filter_bam_before_demuxlet {

	container "${params.containers.general}"
	cache 'lenient'
	maxForks 20
	
	input:
	set val(library), file(bam), file(bam_index), file(barcodes) from demuxlet_in_bam.combine(demuxlet_barcodes_in, by: 0)

	output:
	set val(library), file("${library}.filtered.bam"), file("${library}.filtered.bam.bai") into demuxlet_filtered_in

	when:
	DEMUXLET_LIBRARIES.contains(library)

	"""
	filter-bam-by-barcode.py $bam ${library}.filtered.bam $barcodes
	samtools index ${library}.filtered.bam 
	"""

}

// clip the BAM, if it's ATAC
no_clip = Channel.create()
clip = Channel.create()
demuxlet_filtered_in.choice(no_clip, clip) { it -> get_modality(it[0]) == 'ATAC' ? 1 : 0 }

process clip_bam {

	input:
	set val(library), file(bam), file(index) from clip

	output:
	set val(library), file("${library}.clipped.bam"), file("${library}.clipped.bam.bai") into clipped

	"""
	bam clipOverlap --poolSize 9000000 --in $bam --out ${library}.clipped.bam
	samtools index ${library}.clipped.bam
	"""

}

process mask_vcf {

	input:
	file(vcf) from genotypes
	file(mask) from demuxlet_mask

	output:
	set val('masked'), file('genotypes.masked.vcf.gz') into demuxlet_masked_vcf_in

	"""
	bcftools view --targets-file ^$mask -Ob -o genotypes.masked.vcf.gz $vcf
	"""

}

demuxlet_vcfs_in = demuxlet_unmasked_vcf_in.mix(demuxlet_masked_vcf_in)

process demuxlet {

	container "${params.containers.demuxlet}"
	memory '10 GB'
	cache 'lenient'
	
	input:
	set val(library), file(bam), file(bam_index), val(masking), file(vcf) from no_clip.mix(clipped).combine(demuxlet_vcfs_in)

	output:
	file("${library}.single")
	file("${library}.sing2")
	set val(library), val(masking), file("${library}.best") into demuxlet_out_best

	when:
	DEMUXLET_LIBRARIES.contains(library)

	"""
	demuxlet --sam $bam --vcf $vcf --field GP --out ${library}
	"""

}

process concat_demuxlet {
	
	publishDir "${params.results}/demuxlet/${masking}", mode: 'rellink'

	input:
	set val(library), val(masking), file("*.best") from demuxlet_out_best.groupTuple(by: [0, 1])

	output:
	set val(library), val(masking), file("${library}.best.txt") into call_atac_nuclei_demuxlet
	set val(library), val(masking), file("${library}.best.txt") into call_rna_nuclei_demuxlet

	"""
	cat *.best | grep BARCODE | sort | uniq > ${library}.best.txt
	cat *.best | grep -v BARCODE | sort | uniq >> ${library}.best.txt
	"""

}

call_atac_demuxlet_in = call_atac_nuclei_demuxlet.filter({it -> it[1] == 'unmasked'}).filter({it -> get_modality(it[0]) == 'ATAC'}).map({it -> it[2]}) toSortedList()
call_atac_metrics_in = call_nuclei_atac_in.filter({it -> !is_dual_modality(it[0])}).map({it -> it[1]}).toSortedList() 


call_rna_demuxlet_masked_in = Channel.create()
call_rna_demuxlet_unmasked_in = Channel.create()
call_rna_nuclei_demuxlet.filter({it -> get_modality(it[0]) == 'RNA'}).choice(call_rna_demuxlet_masked_in, call_rna_demuxlet_unmasked_in) {it -> it[1] == 'masked' ? 0 : 1}

process call_atac_nuclei {

	publishDir "${params.results}/call-nuclei-atac"
	publishDir "${params.results}/figures"

	input:
	file(demux) from call_atac_demuxlet_in
	file(thresholds) from Channel.fromPath(INITIAL_THRESHOLDS_ATAC)
	file(library_labels) from Channel.fromPath(LIBRARY_LABELS)
	file(metrics) from call_atac_metrics_in

	output:
	file("*.png")
	file("*.txt")
	file("*.tsv")

	"""
	cat ${metrics.join(' ')} > metrics.txt
	atac-qc.py --thresholds $thresholds --demuxlet ${demux.join(' ')} --library-labels $library_labels --metrics metrics.txt
	"""

}

process call_rna_nuclei {

	publishDir "${params.results}/call-nuclei-rna"

	input:
	file("demux_masked/*") from call_rna_demuxlet_masked_in.map({x -> x[2]}).toSortedList()
	file("demux_unmasked/*") from call_rna_demuxlet_unmasked_in.map({x -> x[2]}).toSortedList()
	file(thresholds) from Channel.fromPath(INITIAL_THRESHOLDS_RNA)
	file(metrics) from call_rna_metrics_in.toSortedList()

	output:
	file("*.tsv")
	file("*.png")

	"""
	rna-qc.py --thresholds $thresholds --demuxlet-masked demux_masked/* --demuxlet-unmasked demux_unmasked/* --metrics ${metrics.join(' ')}
	"""

}
/*
process call_dual_modality_nuclei {

	publishDir "${params.results}/call-nuclei"

	input:
	file(atac_metrics) from ...
	file(rna_metrics) from ...

}

process droplet_utils {

	publishDir "${params.results}/droplet-utils", mode: 'rellink'
	container "${params.containers.seuratv4}"
	memory '20 GB'
	cache 'lenient'

	input:
	set val(library), file(starsolo_dir) from Channel.fromPath(droplet_utils_in)

	output:
	file("*.png")
	file("*.txt")

	"""
	droplet-utils.R $starsolo_dir $library 
	"""

}
*/
