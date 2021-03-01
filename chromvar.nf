#!/usr/bin/env nextflow

FIMO_GLOB = params.fimo_glob
BAM_GLOB = params.bam_glob
PEAK_GLOB = params.peak_glob
MOTIF_TO_TF = params.motif_to_tf
CLUSTERS = params.clusters
CLUSTER_NAMES = params.cluster_names

get_motifs_in = Channel.fromPath(FIMO_GLOB).map({it -> [it.getName().replaceAll('.fimo.txt', ''), it]}).filter({it -> !it[0].contains('disc')}).filter({it -> it[1].countLines() > 0})
bam_in = Channel.fromPath(BAM_GLOB).map({it -> [it.getName().replaceAll('.bam', ''), it]}).filter({it -> params.library_to_modality[it[0]] == 'ATAC'})

process get_motifs {

        maxForks 30
        container "${params.containers.general}"
	publishDir "${params.results}/motif-locations"

        input:
        set val(motif), file(fimo) from get_motifs_in

        output:
        set val(motif), file("${motif}.bed") into motifs_out

        """
        cat $fimo | cut -f2-4 | grep -v sequence | awk '\$3=\$3+1' | perl -pe 's/ /\\t/g' > ${motif}.bed
        """

}

process get_peaks {

	input:
	file(beds) from Channel.fromPath(PEAK_GLOB).toSortedList()

	output:
	file('peaks.bed') into get_peaks_out
	file('peaks.bed') into get_peaks_out_2

	"""
	cat ${beds.join(' ')} | sort -k1,1 -k2n,2 -k3n,3 | bedtools merge -i stdin > merged.bed
	get-bed-centers.py merged.bed | bedtools slop -b 250 -g ${params.chrom_sizes['hg19']} | sort -k1,1 -k2n,2 > peaks.bed
	"""

}

process get_motifs_in_peaks {

	maxForks 30
	publishDir "${params.results}/motifs-in-peaks"

	input:
	set val(motif), file(motifs), file(peaks) from motifs_out.combine(get_peaks_out_2)

	output:
	file("${motif}.in-peaks.bed") into motifs_in_peaks

	"""
	sort -k1,1 -k2n,2 $motifs | bedtools intersect -a stdin -b $peaks -u > ${motif}.in-peaks.bed
	"""
	

}

process get_fragment_counts {

	publishDir "${params.results}/fragment-counts"
	cpus 10
	memory '25 GB'
	maxForks 5

	input:
	set val(library), file(bam), file(peaks) from bam_in.combine(get_peaks_out)

	output:
	file("${library}.counts.bed") into fragment_counts

	"""
	atac-bam-to-fragment-ends.py $bam | sort -k1,1 -k2n,2 --parallel=10 -S 20G > fragment-ends.bed
	bedtools intersect -sorted -a $peaks -b fragment-ends.bed -wa -wb | cut -f1-3,7 | sort -k1,1 -k2n,2 -k4,4 --parallel=10 -S 20G > intersect.bed
	cat intersect.bed | uniq | perl -pe 's/(.*)\\t(.*)_(.*)\$/\$1\\t\$3/' | sort -k1,1 -k2n,2 -k4,4 --parallel=10 -S 20G | bedtools groupby -g 1,2,3,4 -c 4 -o count > ${library}.counts.bed
	"""

}

	//file(motifs) from motifs_in_peaks.toSortedList()
process run_chromvar {

        container "${params.containers.chromvar}"
	memory '20 GB'

	input:
	file(counts) from fragment_counts.toSortedList()
	file(motif_to_tf) from Channel.fromPath(MOTIF_TO_TF)
	file(clusters) from Channel.fromPath(CLUSTERS)
	each file(motifs) from motifs_in_peaks.collate(100, false)

	output:
	file('diff-deviation.txt') into diff_deviations
	file('diff-variation.txt') into diff_variations
	file('deviation-zscores.txt') into deviation_scores

	"""
	chromVAR.R '*.counts.bed' '*.in-peaks.bed' $clusters $motif_to_tf
	"""

}

process concat_chromvar {
	
	publishDir "${params.results}/chromvar"

	input:
	file("diff-deviation.*.txt") from diff_deviations.toSortedList()
	file("diff-variation.*.txt") from diff_variations.toSortedList()
	file("deviation-zscores.*.txt") from deviation_scores.toSortedList()

	output:
	file("*.txt")
	file('deviation-zscores.txt') into plot_chromvar_in

	"""
	head -1 diff-deviation.1.txt > diff-deviation.txt
	head -1 diff-variation.1.txt > diff-variation.txt
	head -1 deviation-zscores.1.txt > deviation-zscores.txt
	cat diff-deviation.*.txt | grep -v p_value_adjusted >> diff-deviation.txt
	cat diff-variation.*.txt | grep -v p_value_adjusted >> diff-variation.txt
	cat deviation-zscores.*.txt | grep -v deviation_z >> deviation-zscores.txt
	"""

}

process plot_chromvar {

	publishDir "${params.results}/plot"

	input:
	file(deviations) from plot_chromvar_in

	output:
	file("*.png")
	file("*.html")

	"""
	plot-chromvar-static-and-interactive.py $deviations $CLUSTERS $CLUSTER_NAMES $MOTIF_TO_TF
	"""

}
