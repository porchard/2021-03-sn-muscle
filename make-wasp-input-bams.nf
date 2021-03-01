#!/usr/bin/env nextflow

NUCLEI_WITH_INDIVIDUAL_ASSIGNMENTS = params.nuclei_with_individual_assignments
CLUSTERS = params.clusters
BAM_GLOB = params.bam_glob

// read in all cluster numbers
clusters = Channel.fromPath(CLUSTERS).splitText().map({it -> it.trim().tokenize('\t')[2]}).unique()

// read in all library --> individual mappings
library_and_individual = Channel.fromPath(NUCLEI_WITH_INDIVIDUAL_ASSIGNMENTS).splitText().map({it -> it.trim().tokenize('\t')}).map({it -> [it[0], it[2]]}).unique()

bam_to_library = {
	f ->
	m = f.getName().replaceAll('.bam', '')
	return m
}

bams = Channel.fromPath(BAM_GLOB).map({it -> [bam_to_library(it), it]})

//println clusters
//println library_and_individual
//println bams

process subset_bam_by_individual {
	
	publishDir "${params.results}/split-library-by-individual"
	maxForks 10
	tag "${library} ${individual}"

	input:
	set val(library), file(bam), val(individual), file(individual_assignments) from bams.combine(library_and_individual, by: 0).combine(Channel.fromPath(NUCLEI_WITH_INDIVIDUAL_ASSIGNMENTS))

	output:
	set val(individual), val(library), file("${library}.${individual}.bam") into subset_out

	"""
	grep -w $library $individual_assignments | grep -w $individual | cut -f2 > barcodes.txt
	filter-bam-by-barcode.py $bam ${library}.${individual}.bam barcodes.txt
	"""

}

/*
process subset_bam_by_individual_and_cluster {

	publishDir "${params.results}/split-library-by-individual-and-cluster"
	maxForks 15
	tag "${library} ${individual} ${cluster}"

	input:
	set val(individual), val(library), file(bam), val(cluster), file(cluster_assignments), file(individual_assignments) from subset_out.combine(clusters).combine(Channel.fromPath(CLUSTERS)).combine(Channel.fromPath(NUCLEI_WITH_INDIVIDUAL_ASSIGNMENTS))

	output:
	file("${library}.${individual}.${cluster}.bam") 

	"""
	grep -w $library $individual_assignments | grep -w $individual | cut -f2 | sort | uniq > barcodes-from-individual-and-library.txt
	grep -w $library $cluster_assignments | grep -w $cluster | cut -f2 | sort | uniq  > barcodes-from-library-and-cluster.txt
	cat barcodes-from-library-and-cluster.txt barcodes-from-individual-and-library.txt | sort | uniq -d > barcodes-from-library-and-cluster-and-individual.txt
	filter-bam-by-barcode.py $bam ${library}.${individual}.${cluster}.bam barcodes-from-library-and-cluster-and-individual.txt
	"""

}
*/

process subset_bam_by_individual_and_cluster {

	publishDir "${params.results}/split-library-by-individual-and-cluster"
	maxForks 15
	tag "${library} ${individual}"

	input:
	set val(individual), val(library), file(bam), file(cluster_assignments) from subset_out.combine(Channel.fromPath(CLUSTERS))

	output:
	file("${library}.${individual}.*.bam") 

	"""
	grep -w $library $cluster_assignments | cut -f2-3 | sort | uniq  > barcode-to-cluster.txt
	split-bam-by-cluster.py $bam ${library}.${individual}. barcode-to-cluster.txt
	"""

}
