#!/usr/bin/env nextflow

CREDIBLE_SET_FILE_GLOB = params.cs_file_glob
DSVM_FILE_GLOB = params.dsvm_file_glob
PEAK_FILE_GLOB = params.peak_file_glob

process get_dsvm_z_scores {

	publishDir "${params.results}/dsvm"
	memory '30 GB'

	input:
	file(cs) from Channel.fromPath(CREDIBLE_SET_FILE_GLOB).toSortedList()
	each file(x) from Channel.fromPath(DSVM_FILE_GLOB)

	output:
	file("${cluster}.dsvm.txt") into dsvm_z_scores

	script:
	cluster = x.getName().replaceAll('.txt', '').replaceAll('snp_scores_', '')

	"""
	cat ${cs.join(' ')} | grep -v IndexSNP | cut -f2-3 | sort | uniq | perl -pe 's/^/chr/' > pos-file.txt	
	get-deltaSVM-z-scores-at-pos.py pos-file.txt $x > ${cluster}.dsvm.txt
	"""

}

process plot {

	publishDir "${params.results}/plot"

	input:
	file(dsvm) from dsvm_z_scores.toSortedList()
	file(cs) from Channel.fromPath(CREDIBLE_SET_FILE_GLOB).toSortedList()
	file(peaks) from Channel.fromPath(PEAK_FILE_GLOB).toSortedList()

	output:
	file("*.png")

	"""
	plot-diamante-loci.py --cs ${cs.join(' ')} --dsvm ${dsvm.join(' ')} --peaks ${peaks.join(' ')}
	"""

}
