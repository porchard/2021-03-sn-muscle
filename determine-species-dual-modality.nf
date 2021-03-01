#!/usr/bin/env nextflow

ATAC_WHITELIST = params.atac_whitelist
RNA_WHITELIST = params.rna_whitelist
ATAC_GLOB = params.atac_glob
RNA_GLOB = params.rna_glob

atac = Channel.fromPath(ATAC_GLOB)
rna = Channel.fromPath(RNA_GLOB)

process get_readname_barcode_flag_atac {

	input:
	file(whitelist) from Channel.fromPath(ATAC_WHITELIST)
	each file(bam) from atac

	output:
	set val(library), file("${library}.txt") into atac_raw

	script:
	library = bam.getName().replaceAll('.bam', '')

	"""
	get-readname-barcode-flag.py $whitelist $bam > ${library}.txt
	"""

}

process get_readname_barcode_flag_rna {

	input:
	file(whitelist) from Channel.fromPath(RNA_WHITELIST)
	each file(bam) from rna

	output:
	set val(library), file("${library}.txt") into rna_raw

	script:
	library = bam.getName().replaceAll('.bam', '')

	"""
	get-readname-barcode-umi-flag.py $whitelist $bam > ${library}.txt
	"""

}

process sort_atac_by_barcode {

	publishDir "${params.results}/atac/sorted"
	
	cpus 15
	memory '25 GB'

	input:
	set val(library), file(x) from atac_raw

	output:
	file("${library}.sorted.txt") into atac_sorted

	"""
	sort -k2,2 -k1,1 -S 20G --parallel=15 $x > ${library}.sorted.txt
	"""
}

process sort_rna_by_barcode {

	publishDir "${params.results}/rna/sorted"
	
	cpus 15
	memory '25 GB'

	input:
	set val(library), file(x) from rna_raw

	output:
	file("${library}.sorted.txt") into rna_sorted

	"""
	sort -k2,2 -k1,1 -S 20G --parallel=15 $x > ${library}.sorted.txt
	"""
}

process gather_atac_counts {
	
	publishDir "${params.results}/atac/summarized"
	executor 'local'

	input:
	file(x) from atac_sorted.toSortedList()

	output:
	file('atac-processed.txt') into atac_processed

	script:
	genomes = x.collect({tmp -> tmp.getName().replace('.sorted.txt', '')})

	"""
	count-species-unique-mapping-atac-2.py ${x.join(' ')} ${genomes.join(' ')} > atac-processed.txt
	"""

}

process gather_rna_counts {
	
	publishDir "${params.results}/rna/summarized"
	executor 'local'

	input:
	file(x) from rna_sorted.toSortedList()

	output:
	file('rna-processed.txt') into rna_processed

	script:
	genomes = x.collect({tmp -> tmp.getName().replace('.sorted.txt', '')})

	"""
	count-species-unique-mapping-rna.py ${x.join(' ')} ${genomes.join(' ')} > rna-processed.txt
	"""

}

// TODO: update paths used in this script
process plot {

	publishDir "${params.results}/plot"
	executor 'local'

	input:
	file(x) from rna_processed
	file(y) from atac_processed

	output:
	file("*.png")

	"""
	multiome-genome-comparisons.py
	"""

}
