#!/usr/bin/env nextflow

ATAC_BAM_GLOB = params.atac_bam_glob
RNA_BAM_GLOB = params.rna_bam_glob

process make_bed_atac {

	publishDir "${params.results}/beds/atac"

	input:
	file(bam) from Channel.fromPath(ATAC_BAM_GLOB)

	output:
	file("${library}.pruned.bed.gz") into atac_bed_out

	script:
	library = bam.getName().replaceAll('.pruned.bam', '')

	"""
	bedtools bamtobed -i $bam | awk '\$5="."' | perl -pe 's/ /\\t/g' | gzip -c > ${library}.pruned.bed.gz
	"""

}

process add_barcodes_to_readname {

	input:
	file(bam) from Channel.fromPath(RNA_BAM_GLOB)

	output:
	file("${library}.with-barcodes.bam") into rna_in

	script:
	library = bam.getName().replaceAll('.before-dedup.bam', '')

	"""
	add-barcode-and-umi-to-readname.py $bam ${library}.with-barcodes.bam
	"""

}

process make_bed_rna {

	publishDir "${params.results}/beds/rna"

	input:
	file(bam) from rna_in

	output:
	file("${library}.pruned-no-dedup.bed.gz") into rna_bed_out

	script:
	library = bam.getName().replaceAll('.with-barcodes.bam', '')

	"""
	bedtools bamtobed -i $bam | awk '\$5="."' | perl -pe 's/ /\\t/g' | gzip -c > ${library}.pruned-no-dedup.bed.gz
	"""

}

process update_read_names_and_starch {

	publishDir "${params.results}/starch"
	memory '150 GB'
	errorStrategy 'ignore'

	input:
	file(bed) from atac_bed_out.mix(rna_bed_out)

	output:
	file("${pref}.starch")

	script:
	pref = bed.getName().replaceAll('.bed.gz', '')

	"""
	change-read-names.py $bed | cut -f1-4,6 | sort -k1,1 -k2n,2 -k3n,3 -S 70G | starch - > ${pref}.starch
	"""

}
