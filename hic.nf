#!/usr/bin/env nextflow

FITHIC_GLOB = params.fithic_glob
PCHIC_GLOB = params.pchic_glob
YUE_GLOB = params.yue_glob
BINS = params.bins
PEAKS = params.peaks
SNP_BED = params.snp_bed
TSS = params.tss
// /lab/work/porchard/data/gencode-tss/hg19.bed 

fithic = Channel.fromPath(FITHIC_GLOB).map({it -> [it.getName().tokenize('.')[1].tokenize('_')[0], it]}).groupTuple().combine(Channel.fromPath(BINS)) // tissue, file

process fithic_to_bed_6 {

	publishDir "${params.results}/bed6"
	errorStrategy 'ignore'

	input:
	set val(tissue), file(x), file(bins) from fithic

	output:
	file("${tissue}.shmitt.bed") into bed6

	"""
	fithic-to-bed6.py $bins 0.000001 ${x.join(' ')} > ${tissue}.shmitt.bed
	"""

}

process filter_yue {

	input:
	file(bed) from Channel.fromPath(YUE_GLOB)

	output:
	file(out) into yue_filtered

	script:
	out = bed.getName().replaceAll('.loops', '')

	"""
	grep -v _ $bed > $out
	"""

}

process annotate {

	publishDir "${params.results}/annotate"

	input:
	file(peaks) from Channel.fromPath(PEAKS)
	file(tss) from Channel.fromPath(TSS)
	each file(bed) from bed6.mix(yue_filtered)

	output:
	file("${bn}.annotated.bed") into hic_annotated
	
	script:
	bn = bed.getName().replaceAll('.bed', '')

	"""
	annotate-hic.py $bed $tss $peaks > ${bn}.annotated.bed
	"""

}

process pchic_to_bed {

	publishDir "${params.results}/bed"

	input:
	file(x) from Channel.fromPath(PCHIC_GLOB)

	output:
	file(bn) into pchic_out

	script:
	bn = x.getName().replaceAll('.xlsx', '.bed')

	"""
	pchic-to-bed.py $x > $bn
	"""

}

process annotate_pchic {

	publishDir "${params.results}/annotate-pchic"

	input:
	file(peaks) from Channel.fromPath(PEAKS)
	each file(x) from pchic_out

	output:
	file(out)

	script:
	out = x.getName().replaceAll('.bed', '.annotated.bed')

	"""
	bedtools intersect -wa -wb -a $peaks -b $x > $out
	"""

}

process plot {

	publishDir "${params.results}/figures"

	input:
	file(x) from hic_annotated.toSortedList()
	file(tss) from Channel.fromPath(TSS)
	file(snp_bed) from Channel.fromPath(SNP_BED)

	output:
	file("*.png")

	"""
	make-hic-heatmap.py $snp_bed $tss ${x.join(' ')}
	"""

}
