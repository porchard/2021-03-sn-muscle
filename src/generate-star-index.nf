#!/usr/bin/env nextflow


index_in = Channel.from([[params.genome, file(params.fasta), file(params.gtf)]])

process make_index {

	cpus 10
	memory '50 GB'

	publishDir "${params.results}"

	input:
	set val(genome), file(fasta), file(gtf) from index_in

	output:
	file("*")	

	"""
	ln -s $gtf ${genome}.gtf
	STAR --runThreadN 10 --runMode genomeGenerate --genomeDir ${params.results} --genomeFastaFiles $fasta --sjdbGTFfile ${genome}.gtf
	ln -s ${genome}.gtf annotation.gtf
	"""

}

