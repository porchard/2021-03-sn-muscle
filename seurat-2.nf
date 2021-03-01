#!/usr/bin/env nextflow

LIGER_IN = params.liger_in

process cluster {

	memory '75 GB'
	publishDir "${params.results}/seurat-out"
	container "${params.containers.seuratv4}"

	input:
	file(mats) from Channel.fromPath(LIGER_IN).toSortedList()

	output:
	file("*.png")
	file("*.pdf")
	file("*.Rda")
	file("*.txt")

	"""
	seurat-joint-2.R ${mats.join(' ')}
	"""

}
