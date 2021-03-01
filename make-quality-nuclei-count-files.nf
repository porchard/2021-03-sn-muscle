#!/usr/bin/env nextflow

STARSOLO_DIR_GLOB = params.starsolo_dir_glob
ATAC_GENE_COUNT_FILE_GLOB = params.atac_gene_count_glob
NUCLEI = params.nuclei
individual_assignments = params.individual_assignments

HUMAN_ATAC_COUNTS = Channel.fromPath(ATAC_GENE_COUNT_FILE_GLOB).filter({it -> it.getName().contains('hg19')})
RAT_ATAC_COUNTS = Channel.fromPath(ATAC_GENE_COUNT_FILE_GLOB).filter({it -> it.getName().contains('rn6')})

process get_atac_gene_counts {

	publishDir "${params.results}/atac"
	memory '20 GB'
	cache 'deep'

	input:
	file(nuclei) from Channel.fromPath(NUCLEI)
	each file(counts) from Channel.fromPath(ATAC_GENE_COUNT_FILE_GLOB)
	
	output:
	file("${library}.hdf5") into liger_atac_in
	file("${library}.hdf5") into liger_atac_in_2

	script:
	library = counts.getName().replaceAll('.counts.txt', '')

	"""
	grep -w $library $nuclei | cut -f2 > keep-barcodes.txt
	atac-gene-counts-to-hdf5.py $counts keep-barcodes.txt ${library}.hdf5
	"""

}

process get_rna_gene_counts {
	
	publishDir "${params.results}/rna"
	memory '20 GB'
	cache 'deep'

	input:
	file(nuclei) from Channel.fromPath(NUCLEI)
	each file(starsolo) from Channel.fromPath(STARSOLO_DIR_GLOB, type: 'dir')
	
	output:
	set val(library), file("${library}.hdf5") into seurat_in
	set val(library), file("${library}.hdf5") into seurat_in_2
	set val(library), file("${library}.hdf5") into decontx_in

	script:
	library = starsolo.getName()

	"""
	grep -w $library $nuclei | cut -f2 > keep-barcodes.txt
	starsolo-counts-to-hdf5-gene-name.py ${starsolo}/Solo.out/GeneFull/raw/matrix.mtx ${starsolo}/Solo.out/GeneFull/raw/barcodes.tsv ${starsolo}/Solo.out/GeneFull/raw/features.tsv keep-barcodes.txt ${library} ${library}.hdf5
	"""

}

human_seurat_in = seurat_in.filter({it -> it[0].contains('hg19')}).map({it -> it[1]})
rat_seurat_in = seurat_in_2.filter({it -> it[0].contains('rn6')}).map({it -> it[1]})

process seurat {

	publishDir "${params.results}/seurat/hg19"
	container "${params.containers.seuratv4}"
	memory '50 GB'
	cache 'deep'

	input:
	file(x) from human_seurat_in.toSortedList()

	output:
	set file('clusters.txt'), file('umap.txt') into decontx_clusters
	file("*.png")

	"""
	run-seurat-by-individual.R $individual_assignments '*.hdf5'
	"""

}

process seurat_rat {

	publishDir "${params.results}/seurat/rn6"
	container "${params.containers.seuratv4}"
	memory '50 GB'
	cache 'deep'

	input:
	file(x) from rat_seurat_in.toSortedList()

	output:
	set file('clusters.txt'), file('umap.txt') into decontx_clusters_rat
	file("*.png")

	"""
	run-seurat-rat.R $individual_assignments '*.hdf5'
	"""

}

human_decontx_in = Channel.create()
rat_decontx_in = Channel.create()
decontx_in.choice(human_decontx_in, rat_decontx_in) { it -> it[0].contains('hg19') ? 0 : 1}

process decontX {

	publishDir "${params.results}/decontX"
	container "${params.containers.decontX}"
	memory '20 GB'
	cache 'deep'

	input:
	set val(library), file(uncorrected), file(clusters), file(umap) from human_decontx_in.combine(decontx_clusters).mix(rat_decontx_in.combine(decontx_clusters_rat))

	output:
	set val(library), file("${library}.decontaminated-counts-rounded.txt") into decontx_out
	file("*.png")
	file("*.txt")

	"""
	decontX.R $uncorrected $clusters $umap ${library}.
	"""

}

process make_corrected_hdf5 {

	publishDir "${params.results}/rna-corrected"
	memory '40 GB'
	cache 'deep'

	input:
	set val(library), file(txt) from decontx_out

	output:
	file("${library}.hdf5") into liger_rna_in
	file("${library}.hdf5") into liger_rna_in_2

	"""
	tsv-to-hd5.py $txt ${library}.hdf5 ${library} nucleus
	"""

}

process make_liger_in {

	publishDir "${params.results}/liger-in"
	memory '150 GB'

	input:
	file(rna) from liger_rna_in.toSortedList()
	file(atac) from liger_atac_in.toSortedList()
	file(x) from Channel.fromPath(individual_assignments)

	output:
	file("*.hdf5")

	"""
	make-liger-in.py --individuals $x --rna ${rna.join(' ')} --atac ${atac.join(' ')}
	"""
}

process make_liger_in_by_species_and_modality {

	publishDir "${params.results}/liger-in-by-species-and-modality"
	memory '150 GB'

	input:
	file(rna) from liger_rna_in_2.toSortedList()
	file(atac) from liger_atac_in_2.toSortedList()
	file(x) from Channel.fromPath(individual_assignments)

	output:
	file("*.hdf5")

	"""
	make-liger-in-merge-samples.py --individuals $x --rna ${rna.join(' ')} --atac ${atac.join(' ')}
	"""
}
