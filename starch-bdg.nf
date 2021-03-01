#!/usr/bin/env nextflow

BDG_GLOB = params.bdg_glob

process update_read_names_and_starch {

	publishDir "${params.results}/starch"
	memory '35 GB'

	input:
	file(bedgraph) from Channel.fromPath(BDG_GLOB)

	output:
	file("${pref}.starch")

	script:
	pref = bedgraph.getName().replaceAll('_treat_pileup.bdg', '')

	"""
	sort -k1,1 -k2n,2 -k3n,3 -S 30G $bedgraph | starch - > ${pref}.starch
	"""

}
