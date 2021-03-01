#!/usr/bin/env nextflow

LIGER_IN_MATRICES = params.liger_in
LAMBDA = params.lambda
KS = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 28, 30, 35, 40]
SELECT_GENES = ['RNA', 'KSM2_RNA']

process factorize_KSM2_and_rat {

	memory '30 GB'
	publishDir "${params.results}/factorize/${k}/${sg}"
	container "${params.containers.ligerv05}"

	input:
	file(x) from Channel.fromPath(LIGER_IN_MATRICES).toSortedList()
	each k from KS

	output:
	set val(k), val(sg), file("liger-factorized.Rda") into factorize_out

	script:
	sg = 'KSM2_rat_RNA'

	"""
	liger-factorize-use-rat.R --mats ${x.join(',')} --factorization_k $k --select_genes KSM2_RNA --factorization_lambda $LAMBDA --out liger-factorized.Rda
	"""

}

process post_factorization_3 {

	memory '30 GB'
	publishDir "${params.results}/post-factorization-3/${sg}/${k}"
	publishDir "${params.results}/clusters/${sg}", pattern: "*.txt"

	input:
	set val(k), val(sg), file(rda) from factorize_out

	output:
	file("*.pdf")
	file("*.png")
	set val(sg), file("*.txt") into clusters_out

	"""
	liger-post-factorization-both-species-3.R --rds $rda
	"""

}
