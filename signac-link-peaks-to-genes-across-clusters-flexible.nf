#!/usr/bin/env nextflow

CLUSTERS = params.clusters
ATAC_HDF5 = params.atac
RNA_HDF5 = params.rna
METHODS = ['pearson', 'spearman']

// generate lists of genes to test (100 at a time)
process make_lists_of_genes {

	input:
	file(rna) from Channel.fromPath(RNA_HDF5)

	output:
	file("genes-batch.*.txt") into gene_lists
	file("genes-batch.*.txt") into gene_lists_2

	"""
	list-features-matrix.py $rna | cut -f2 > genes.txt	
	split --additional-suffix=.txt --lines 500 --suffix-length 10 genes.txt genes-batch.
	"""

}


process find_links {

	publishDir "${params.results}/link-${method}"
	memory '30 GB'
	errorStrategy 'ignore'

	input:
	file(clusters) from Channel.fromPath(CLUSTERS)
	file(rna) from Channel.fromPath(RNA_HDF5)
	file(atac) from Channel.fromPath(ATAC_HDF5)
	each file(genes) from gene_lists.flatten()
	each method from METHODS

	output:
	file("all.batch_${batch}.txt")

	script:
	batch = genes.getName().tokenize('.')[1]

	"""
	signac-link-peaks-to-genes-flexible.R $clusters all $genes $atac $rna $method 1500000 5 all.batch_${batch}.txt
	"""

}
