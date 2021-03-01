#!/usr/bin/env nextflow

STARSOLO_DIR_GLOB = params.starsolo_dir_glob
starsolo_dirs = Channel.fromPath(STARSOLO_DIR_GLOB, type: 'dir')

process rna_qc_file {

	publishDir "${params.results}/rnaseq-qc"

	input:
	file(solo) from starsolo_dirs

	output:
	file("${library}.rnaseq-qc.txt")

	script:
	library = solo.getName()
	bam = solo.getName() + '/Aligned.sortedByCoord.out.bam'
	counts_matrix = solo.getName() + '/Solo.out/GeneFull/raw/matrix.mtx'
	barcodes_tsv = solo.getName() + '/Solo.out/GeneFull/raw/barcodes.tsv'

	"""
	qc-from-starsolo.py $bam $counts_matrix $barcodes_tsv > ${library}.rnaseq-qc.txt
	"""

}
