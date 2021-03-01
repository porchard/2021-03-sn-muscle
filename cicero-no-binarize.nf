#!/usr/bin/env nextflow

COUNTS_GLOB = params.counts_glob
SNPS = params.snps
TSS = params.tss
//VCF = params.vcf
UMAP = params.umap
//NUCLEUS_TO_INDIVIDUAL = params.nucleus_to_individual
WINDOWS = [500000, 1000000, 1500000, 1750000, 2000000]

process reformat_counts {

	publishDir "${params.results}/cicero-counts-in"

	input:
	file(x) from Channel.fromPath(COUNTS_GLOB)

	output:
	set val(cluster), val(genome), file("${cluster}.cicero-in.txt") into counts_all

	script:
	cluster = x.getName().replaceAll('.counts.txt', '').tokenize('-')[0]
	genome = x.getName().replaceAll('.counts.txt', '').tokenize('-')[1]
	
	"""
	cut -f1-2 $x | perl -pe 's/\\t/-/' > nucleus.txt
	cut -f3 $x | perl -pe 's/:/_/g' > peak.txt
	cut -f4 $x | perl -pe 's/:/_/g' > count.txt
	paste peak.txt nucleus.txt count.txt > ${cluster}.cicero-in.txt
	"""

}

cicero_in = counts_all.filter({it -> it[1] == 'hg19'}).map({it -> it[2]})

process cicero {

	publishDir "${params.results}/cicero"
	container "${params.containers.cicero}"
	memory '70 GB'

	input:
	file(umap) from Channel.fromPath(UMAP)
	each file(counts) from cicero_in
	each window from WINDOWS

	output:
	file("${cluster}.${window}.cicero.txt") into cicero_out

	script:
	cluster = counts.getName().tokenize('.')[0]

	"""
	cicero-no-binarize.R $umap $counts ${params.chrom_sizes['hg19']} $window ${cluster}.${window}.cicero.txt
	"""

}

process annotate_cicero {

	publishDir "${params.results}/annotate-cicero"
	memory '20 GB'
	executor 'local'

	input:
	file(snps) from Channel.fromPath(SNPS)
	file(tss) from Channel.fromPath(TSS)
	each file(x) from cicero_out

	output:
	set val(window), val(cluster), val(prefix), file("${prefix}.cicero-annotated.txt") into annotate_cicero_out

	script:
	prefix = x.getName().replaceAll('.cicero.txt', '')
	(cluster, window) = prefix.tokenize('.')

	"""
	annotate-cicero-tss-and-snps.py --tss-distance 2000 --tss-upstream $x $tss $snps ${params.chrom_sizes['hg19']} > ${prefix}.cicero-annotated.txt
	"""

}
