#!/usr/bin/env nextflow

IONICE = 'ionice -c2 -n7'
CLUSTER_NAMES = params.cluster_names
LIBRARY_LABELS = params.library_labels
LIGER_IN_MATRICES = params.liger_in_matrices_glob

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

MACS2_GENOME_SIZE = [
    'rn4': 'mm',
    'rn5': 'mm',
    'rn6': 'mm',
    'mm9': 'mm',
    'mm10': 'mm',
    'hg19': 'hs',
    'hg38': 'hs'
]

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

// want:
// nucleus --> cluster assignments for each library
// library --> [pruned bam, genome, modality] for each library that was clustered

// then:
// make ATAC bams for each library and cluster
// call ATAC peaks on each library, get peak counts
// make master peaks, get peak counts

libraries = params.libraries.keySet()

ATAC_LIBRARIES = []
RNA_LIBRARIES = []
per_library_per_cluster_bam_in = []

for(library in libraries) {
	if (get_modality(library) == 'ATAC') {
		ATAC_LIBRARIES << library
		per_library_per_cluster_bam_in << [library, file(get_pruned(library)), file(params.clusters)]
	}
	if (get_modality(library) == 'RNA') {
		RNA_LIBRARIES << library
	}
}


process aggregate_bam {

        publishDir "${params.results}/aggregate-of-quality-nuclei/bam", mode: 'rellink'
        container "${params.containers.general}"

        input:
	set val(library), file("to-subset.bam"), file(clusters) from Channel.from(per_library_per_cluster_bam_in)

        output:
        set val(library), file("${library}.bam"), file("${library}.bam.bai")

        """
        grep -w $library $clusters | cut -f2 > keep-nuclei.txt
        filter-bam-by-barcode.py to-subset.bam ${library}.bam keep-nuclei.txt
        samtools index ${library}.bam
        rm keep-nuclei.txt
        """

}


// Process by cluster
process per_library_per_cluster_bam {

	publishDir "${params.results}/bam-per-library-per-cluster", mode: 'rellink'
	container "${params.containers.general}"
	maxForks 10

	input:
	set val(library), file("bam-to-split"), file(clusters) from Channel.from(per_library_per_cluster_bam_in)

	output:
	file("*.bam") into per_library_per_cluster_bam_out_1
	file("*.bam") into per_library_per_cluster_bam_out_2
	file("*.bam") into per_library_per_cluster_bam_out_3
	file("*.bam") into per_library_per_cluster_bam_out_4

	"""
	grep -w $library $clusters | cut -f2,3 > barcode_to_cluster.txt
        split-bam-by-cluster.py bam-to-split ${library}. barcode_to_cluster.txt
	"""

}

per_library_per_cluster_bam_redirect = per_library_per_cluster_bam_out_1.flatten().map({it -> [it.getName().tokenize('.')[0], it.getName().tokenize('.')[1], it]}) // library, cluster, bam
per_library_per_cluster_bam_redirect.map({it -> [it[1], get_genome(it[0]), get_modality(it[0]), it[2]]}).into{merge_in; tmp_in}
peak_counts_bam_in = per_library_per_cluster_bam_out_3.flatten().map({it -> [it.getName().tokenize('.')[1], get_genome(it.getName().tokenize('.')[0]), it.getName().tokenize('.')[0], it]}) // cluster, genome, library, bam
peak_counts_bam_in_2 = per_library_per_cluster_bam_out_4.flatten().map({it -> [get_genome(it.getName().tokenize('.')[0]), it.getName().tokenize('.')[1], it.getName().tokenize('.')[0], it]}) // genome, cluster, library, bam

process merge_cluster {

	publishDir "${params.results}/bam-per-cluster/bam", mode: 'rellink'
	container "${params.containers.general}"
	maxForks 5

	input:
	set val(cluster), val(genome), val(modality), file(bams) from merge_in.groupTuple(by: [0, 1, 2])

	output:
	set val(cluster), val(genome), val(modality), file("${cluster}.${genome}.${modality}.bam") into bamtobed_in
	set val('all'), val(genome), val(modality), file("${cluster}.${genome}.${modality}.bam") into make_aggregate_in

	"""
	samtools merge ${cluster}.${genome}.${modality}.bam ${bams.join(' ')}
	"""

}

process make_aggregate {

	publishDir "${params.results}/bam-per-cluster/bam", mode: 'rellink'
	container "${params.containers.general}"
	maxForks 5

	input:
	set val(cluster), val(genome), val(modality), file(bams) from make_aggregate_in.groupTuple(by: [0, 1, 2])

	output:
	set val(cluster), val(genome), val(modality), file("${cluster}.${genome}.${modality}.bam") into bamtobed_aggregate_in

	"""
	samtools merge ${cluster}.${genome}.${modality}.bam ${bams.join(' ')}
	"""

}

process per_cluster_bamtobed {

	container "${params.containers.general}"
	tag "${cluster}"
	maxForks 5

	input:
	set val(cluster), val(genome), val(modality), file(bam) from bamtobed_in.mix(bamtobed_aggregate_in)

	output:
	set val(cluster), val(genome), file("reads.bed") into broadpeaks_in
	set val(cluster), val(genome), file("reads.bed") into narrowpeaks_in

	when:
	modality == 'ATAC'

	"""
	bedtools bamtobed -i $bam > reads.bed
	"""

}

process peaks {

	publishDir "${params.results}/peaks/broad", mode: 'rellink'
	container "${params.containers.macs2}"
	tag "${cluster}"

	input:
	set val(cluster), val(genome), file(reads) from broadpeaks_in

	output:
	set val(cluster), val(genome), file("${cluster}-${genome}_treat_pileup.bdg") into cluster_bigwigs_in
	set val(cluster), val(genome), file("${cluster}-${genome}_peaks.broadPeak.noblacklist") into master_peaks_in
	set val(cluster), val(genome), file("${cluster}-${genome}_peaks.broadPeak.noblacklist") into peak_counts_peaks_in
	set val(cluster), val(genome), file("${cluster}-${genome}_peaks.broadPeak.noblacklist") into cluster_chromatin_state_overlap_hg19_in
	set val(cluster), val(genome), file("${cluster}-${genome}_peaks.broadPeak.noblacklist") into peak_liftover_in

	"""
        macs2 callpeak -t $reads --outdir . -f BED -n ${cluster}-${genome} --SPMR -g ${get_genome_size(genome)} --nomodel --shift -100 --seed 762873 --extsize 200 -B --broad --keep-dup all
        bedtools intersect -a ${cluster}-${genome}_peaks.broadPeak -b ${params.blacklist[genome].collect().join(' ')} -v | sort -k1,1 -k2n,2 > ${cluster}-${genome}_peaks.broadPeak.noblacklist
	"""

}

process narrow_peaks {

	publishDir "${params.results}/peaks/narrow", mode: 'rellink'
	container "${params.containers.macs2}"
	tag "${cluster}"

	input:
	set val(cluster), val(genome), file(reads) from narrowpeaks_in

	output:
	set val(genome), val(cluster), file("${cluster}-${genome}_peaks.narrowPeak.noblacklist") into enrichment_features_in

	"""
        macs2 callpeak -t $reads --outdir . -f BED -n ${cluster}-${genome} --SPMR -g ${get_genome_size(genome)} --nomodel --shift -37 --seed 762873 --extsize 73 -B --keep-dup all
        bedtools intersect -a ${cluster}-${genome}_peaks.narrowPeak -b ${params.blacklist[genome].collect().join(' ')} -v | sort -k1,1 -k2n,2 > ${cluster}-${genome}_peaks.narrowPeak.noblacklist
	"""

}

process rat_peak_liftover {

	publishDir "${params.results}/rat-peak-liftover", mode: 'rellink'
	container "${params.containers.bnmapper}"

	input:
	set val(cluster), val(genome), file(peaks) from peak_liftover_in

	output:
	set val(cluster), val(genome), file(outfile) into cluster_chromatin_state_overlap_rn6_in

	when:
	genome == "rn6"

	script:
	outfile = cluster + '-' + genome + '_peaks_in_hg19.bed'

	"""
	bnMapper.py -o ${outfile}.tmp $peaks /lab/work/porchard/data/chain/rn6ToHg19.over.chain.gz
	cat ${outfile}.tmp | sort -k1,1 -k2n,2 | bedtools merge -i stdin -d 50 > $outfile
	"""	

}


process cluster_chromatin_state_overlap {

        publishDir "${params.results}/chromhmm-overlap", mode: 'rellink'
	container "${params.containers.general}"

        input:
        set val(cluster), val(genome), file(peaks) from cluster_chromatin_state_overlap_hg19_in.filter({x -> x[1].toString() == 'hg19'}).mix(cluster_chromatin_state_overlap_rn6_in)

        output:
        set file("${cluster}-${genome}.chromhmm_overlap.txt"), file("${cluster}-${genome}.chromhmm_overlap.pdf")

	"""
        chromhmm_overlap.py $peaks ${get_tss(genome)} ${get_chrom_sizes(genome)} ${params.chromatin_state_glob} > ${cluster}-${genome}.chromhmm_overlap.txt
        plot_chromatin_state_overlap.R --overlap ${cluster}-${genome}.chromhmm_overlap.txt --out ${cluster}-${genome}.chromhmm_overlap.pdf
        """

}

peak_counts_in = peak_counts_bam_in.filter({x -> get_modality(x[2].toString()) == 'ATAC'}).combine(peak_counts_peaks_in, by: [0, 1])
process peak_counts {

	memory '22 GB'
	cpus 5
	container "${params.containers.general}"
	maxForks 10

	input:
	set val(cluster), val(genome), val(library), file(bam), file(peaks) from peak_counts_in

	output:
	set val(cluster), val(genome), file("${library}.${cluster}.counts.txt") into peak_counts_out

	"""
	sort-bed-by-bam.py $peaks $bam > peaks.sorted.bed
	cut -f1-3 peaks.sorted.bed | perl -pe 's/\\t/:/g' > peak-name.txt
	cut -f1-3 peaks.sorted.bed | paste - peak-name.txt > peaks.for-overlap.bed
	bedtools intersect -wa -wb -bed -sorted -a $bam -b peaks.for-overlap.bed | cut -f4,16 | perl -pe 's@(.*_.*)/\\d+\\t(.*)@\$1\\t\$2@' | uniq | sort --parallel=5 -S 20G | uniq | perl -pe 's@.*_(.*)\\t(.*)@\$1\\t\$2@' > overlaps.txt
        atac-overlaps-to-feature-counts.py overlaps.txt ${library} > ${library}.${cluster}.counts.txt
	"""

}

process per_cluster_peak_counts {

	publishDir "${params.results}/peak-counts", mode: 'rellink'
	container "${params.containers.general}"
	memory '5 GB'

	input:
	set val(cluster), val(genome), file("counts*.txt") from peak_counts_out.groupTuple(by: [0, 1])

	output:
	set val(cluster), val(genome), file("${cluster}-${genome}.counts.txt") into per_cluster_peak_counts_out

	"""
	cat counts*.txt > ${cluster}-${genome}.counts.txt
	"""

}

all_cluster_peaks_for_genome = master_peaks_in.filter({it -> it[0] != 'all'}).map({it -> [it[1], it[2]]}).groupTuple()
process master_peaks {

	publishDir "${params.results}/process-by-cluster-round-1/master-peaks", mode: 'rellink'
	container "${params.containers.general}"

	input:
	set val(genome), file(peaks) from all_cluster_peaks_for_genome

	output:
	set val(genome), file("master-peaks.${genome}.bed") into master_peaks_out

	"""
	cat ${peaks.join(' ')} | sort -k1,1 -k2n,2 | bedtools merge -i stdin > master-peaks.${genome}.bed
	"""

}

peak_counts_in_2 = peak_counts_bam_in_2.filter({x -> get_modality(x[2].toString()) == 'ATAC'}).combine(master_peaks_out, by: 0)
process master_peak_counts {

	memory '22 GB'
	cpus 5
	container "${params.containers.general}"
	publishDir "${params.results}/master-peak-counts"
	maxForks 10

	input:
	set val(genome), val(cluster), val(library), file(bam), file(peaks) from peak_counts_in_2

	output:
	set val(cluster), val(genome), file("${library}.${cluster}.counts.txt")

	"""
	sort-bed-by-bam.py $peaks $bam > peaks.sorted.bed
	cut -f1-3 peaks.sorted.bed | perl -pe 's/\\t/:/g' > peak-name.txt
	cut -f1-3 peaks.sorted.bed | paste - peak-name.txt > peaks.for-overlap.bed
	bedtools intersect -wa -wb -bed -sorted -a $bam -b peaks.for-overlap.bed | cut -f4,16 | perl -pe 's@(.*_.*)/\\d+\\t(.*)@\$1\\t\$2@' | uniq | sort --parallel=5 -S 20G | uniq | perl -pe 's@.*_(.*)\\t(.*)@\$1\\t\$2@' > overlaps.txt
        atac-overlaps-to-feature-counts.py overlaps.txt ${library} > ${library}.${cluster}.counts.txt
	"""

}

process cluster_bigwigs {

	publishDir "${params.results}/bigwigs/atac", mode: 'rellink'
	container "${params.containers.general}"
	memory '50 GB'
	tag "${cluster}-${genome}"

	input:
	set val(cluster), val(genome), file(bedgraph) from cluster_bigwigs_in

	output:
	set val(genome), file("${cluster}-${genome}.bw") into bigwig_out

	"""
	bedClip $bedgraph ${get_chrom_sizes(genome)} bedgraph.clipped.bdg
        LC_COLLATE=C sort -k1,1 -k2n,2 bedgraph.clipped.bdg > bedgraph.sorted.bdg
        bedGraphToBigWig bedgraph.sorted.bdg ${get_chrom_sizes(genome)} ${cluster}-${genome}.bw
	"""
	
}

process plot_atac_marker_genes {

	publishDir "${params.results}/figures", mode: 'rellink'

	input:
	set val(genome), file(bigwigs) from bigwig_out.groupTuple()

	output:
	file('marker-genes-atac.png')

	when:
	genome == 'hg19'

	"""
	#plot-marker-genes-fig2.py --gene-bed /lab/work/porchard/sn-muscle-project/data/gencode-coding/gencode.v19.annotation.gtf.gz --out marker-genes-atac.png $CLUSTER_NAMES ${bigwigs.join(' ')}
	ln -s all-hg19.bw aggregate.bw
	plot-marker-genes-fig2.py --gene-bed /lab/work/porchard/sn-muscle-project/data/gencode-bed/gencode.hg19.genes.bed --out marker-genes-atac.png $CLUSTER_NAMES ${bigwigs.join(' ')} aggregate.bw
	"""

}

process get_long_matrices {

        memory '60 GB'
        publishDir "${params.results}/long-matrices"

        input:
        file(x) from Channel.fromPath(LIGER_IN_MATRICES)

        output:
        file(nm) into long_matrices
        file(nm) into long_matrices_2

        script:
        nm = x.getName().replaceAll('.hdf5', '.txt')

        """
        to-old-features-format-matrix-2.py $x | grep -v -w 0 > $nm
        """

}

process make_heatmaps {

        memory '10 GB'
        publishDir "${params.results}/figures"

        input:
        file(mats) from long_matrices.toSortedList()
        file(clusters) from Channel.fromPath(params.clusters)

        output:
        file("*.png")

        """
        make-heatmap-fig2.py $clusters marker-gene-heatmap.png ${mats.join(' ')}
        """

}

process plot_myh_counts {

        memory '40 GB'
        publishDir "${params.results}/figures"
	errorStrategy 'ignore'

        input:
        file(mats) from Channel.fromPath(LIGER_IN_MATRICES).toSortedList()
        file(clusters) from Channel.fromPath(params.clusters)
        file(cluster_names) from Channel.fromPath(CLUSTER_NAMES)

        output:
        file("*.png")

        """
        plot-myh-gene-counts.py $clusters $cluster_names ${mats.join(' ')}
        """

}

process plot_nucleus_counts {

        memory '40 GB'
        publishDir "${params.results}/figures"
	errorStrategy 'ignore'

        input:
        file(clusters) from Channel.fromPath(params.clusters)
        file(cluster_names) from Channel.fromPath(CLUSTER_NAMES)
        file(library_labels) from Channel.fromPath(LIBRARY_LABELS)

        output:
        file("*.png")

        """
        nucleus-count-barplots.py $clusters $cluster_names $library_labels
        """

}


process compare_rubenstein {

        memory '10 GB'
        publishDir "${params.results}/figures"

        input:
        file(mats) from long_matrices_2.toSortedList()
        file(clusters) from Channel.fromPath(params.clusters)

        output:
        file("rubenstein-vs-our-fiber-type-lfcs.pdf")

        """
        compare-rubenstein.py $clusters /lab/work/porchard/sn-muscle-project/data/rubenstein-fiber-type-differential-genes/fiber-type-differential-genes.txt rubenstein-vs-our-fiber-type-lfcs.pdf ${mats.join(' ')}
        """

}

human_rna_mats = Channel.fromPath(LIGER_IN_MATRICES).filter({it -> it.getName().contains('RNA') && it.getName().contains('KSM')})
process plot_gene_expression {

	memory '50 GB'
	publishDir "${params.results}/figures"

	input:
	file(mats) from human_rna_mats.toSortedList()
	file(clusters) from Channel.fromPath(params.clusters)
        file(cluster_names) from Channel.fromPath(CLUSTER_NAMES)

	output:
	file("*.png")

	"""
	make-gene-stripplot-and-heatmap.py --genes MEF2A MEF2B MEF2C MEF2D PITX2 PITX3 FOS JUN --matrices ${mats.join(' ')} --clusters $clusters --cluster-names $cluster_names
	"""

}

human_rna_mats_2 = Channel.fromPath(LIGER_IN_MATRICES).filter({it -> it.getName().contains('RNA') && it.getName().contains('KSM')})
rat_rna_mats = Channel.fromPath(LIGER_IN_MATRICES).filter({it -> it.getName().contains('RNA') && it.getName().contains('rat')})
process plot_gene_expression_near_snps {

	memory '50 GB'
	publishDir "${params.results}/gene-expression-near-snps"

	input:
	file(mats) from human_rna_mats_2.toSortedList()
	file(clusters) from Channel.fromPath(params.clusters)
        file(cluster_names) from Channel.fromPath(CLUSTER_NAMES)

	output:
	file("*.png")

	"""
	#zcat /home/porchard/github/ataqv/data/tss/hg19.tss.refseq.bed.gz | grep -v 'hap' > hg19.tss.bed
	cp /lab/work/porchard/sn-muscle-project/data/gencode-tss-for-connections/hg19.bed hg19.tss.bed
	plot-expression-of-genes-near-snps.py hg19.tss.bed $clusters $cluster_names ${mats.join(' ')}
	"""

}
/*
process plot_gene_expression_near_snps_rat {

	memory '50 GB'
	publishDir "${params.results}/gene-expression-near-snps"

	input:
	file(mats) from rat_rna_mats.toSortedList()
	file(clusters) from Channel.fromPath(params.clusters)
        file(cluster_names) from Channel.fromPath(CLUSTER_NAMES)

	output:
	file("*.png")

	"""
	zcat /home/porchard/github/ataqv/data/tss/rn6.tss.refseq.bed.gz | grep -v 'hap' | grep -v 'chrUn' | grep -v '_random'  > rn6.tss.bed
	plot-expression-of-genes-near-snps.py rn6.tss.bed $clusters $cluster_names ${mats.join(' ')}
	"""

}
*/
