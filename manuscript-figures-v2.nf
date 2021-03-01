#!/usr/bin/nextflow

ROOT = params.projroot
CLUSTER_NAMES = ROOT + '/cluster-names.txt'
FANS_AND_LOADING_POSITION = 'chr8:41,427,421-41,535,659'
HUMAN_MUSCLE_BW_GLOB = params.human_muscle_bw_glob
RAT_MUSCLE_BW_GLOB = params.rat_muscle_bw_glob
OTHER_BW_GLOB = params.other_bw_glob

get_tss = {
	genome ->
	params.tss[genome]
}

get_chrom_sizes = {
	genome ->
	params.chrom_sizes[genome]
}

process rat_gb_shots {

	publishDir "${params.results}/gb-python"
	executor 'local'

	input:
	file(b) from Channel.fromPath('/lab/work/porchard/sn-muscle-project/data/gencode-bed/rn6.chromnames.bed')
	file(muscle_bws) from Channel.fromPath(RAT_MUSCLE_BW_GLOB).toSortedList()

	output:
	file("*.pdf")

	"""
	# ARL15, big
	make-rat-gb-plots.py --out arl15-rat.pdf --highlight 45977453 --gtf_bed $b --region chr2:45967438-45987469 --muscle_bws ${muscle_bws.join(' ')}
	"""
}

process human_gb_shots {

	publishDir "${params.results}/gb-python"
	executor 'local'

	input:
	file(b) from Channel.fromPath('/lab/work/porchard/sn-muscle-project/data/gencode-bed/gencode.hg19.genes.bed')
	file(other_bws) from Channel.fromPath(OTHER_BW_GLOB).toSortedList()
	file(muscle_bws) from Channel.fromPath(HUMAN_MUSCLE_BW_GLOB).toSortedList()
	file(chromhmm) from Channel.fromPath('/lab/work/porchard/sn-muscle-project/data/chromhmm/*.dense.bed')

	output:
	file("*.pdf")

	"""
	# ARL15, big
	make-human-gb-plots.py --out arl15-human-big.pdf --highlight 53271420 --gtf_bed $b --credible_set_file /lab/work/porchard/sn-muscle-project/data/diamante-credible-sets/genetic_credible_sets/g_credible_set_Eur_ARL15_5_53271420.txt --region chr5:53,269,334-53,278,795 --muscle_bws ${muscle_bws.join(' ')} --other_bws ${other_bws.join(' ')} --chromhmm ${chromhmm.join(' ')}
	make-human-gb-plots.py --out arl15-human.pdf --highlight 53271420 --gtf_bed $b --credible_set_file /lab/work/porchard/sn-muscle-project/data/diamante-credible-sets/genetic_credible_sets/g_credible_set_Eur_ARL15_5_53271420.txt --region chr5:53,269,334-53,278,795 --muscle_bws ${muscle_bws.join(' ')}
	# ITPR2, big
	make-human-gb-plots.py --out itpr2-human-big.pdf --highlight 26472562 --gtf_bed $b --credible_set_file /lab/work/porchard/sn-muscle-project/data/diamante-credible-sets/genetic_credible_sets/g_credible_set_Eur_ITPR2_12_26453283.txt --region chr12:26,436,640-26,491,955 --muscle_bws ${muscle_bws.join(' ')} --other_bws ${other_bws.join(' ')} --chromhmm ${chromhmm.join(' ')}
	make-human-gb-plots.py --out itpr2-human.pdf --highlight 26472562 --gtf_bed $b --credible_set_file /lab/work/porchard/sn-muscle-project/data/diamante-credible-sets/genetic_credible_sets/g_credible_set_Eur_ITPR2_12_26453283.txt --region chr12:26,436,640-26,491,955 --muscle_bws ${muscle_bws.join(' ')}
	"""

}

process rnaseq_barcode_rank_plots {

	publishDir "${params.results}/barcode-rank-plots"
	memory '10 GB'
	
	input:
	file(starsolo_dir) from Channel.fromPath(ROOT + "/work/rnaseq/results/starsolo/63_*_rna-hg19", type: 'dir')

	output:
	file("*.png")

	script:
	lib = starsolo_dir.getName()
	nuclei = lib.replaceAll('_rna-hg19', '').replaceAll('63_', '')

	"""
	make-rnaseq-barcode-rank-plot.py ${starsolo_dir}/Solo.out/GeneFull/raw/barcodes.tsv ${starsolo_dir}/Solo.out/GeneFull/raw/matrix.mtx 1000 ${nuclei}000 ${lib}.png
	"""

}

process compute_allelic_bias {
	
	publishDir "${params.results}/allelic-bias"
	memory '40 GB'
	executor 'local'
	
	input:
	file(cluster_names) from Channel.fromPath(CLUSTER_NAMES)

	output:
	file("*.txt")
	file("*.tsv")

	"""
	compute-allelic-bias.py
	make-allelic-bias-supplementary-table.py allelic-bias.txt $cluster_names > allelic-bias-supplementary-table.tsv
	"""

}

process compare_deltasvm_to_allelic_bias {
	
	publishDir "${params.results}/compare-deltasvm-to-allelic-bias"
	memory '40 GB'
	executor 'local'

	output:
	file("*.png")

	"""
	deltasvm-vs-allelic-bias.py
	"""

}

process plot_seurat_links_spearman {

	publishDir "${params.results}/seurat-links/spearman"
	memory '10 GB'

	input:
	file(cic) from Channel.fromPath(ROOT + '/work/find-peak-gene-links-across-clusters-2/results/link-spearman/*').toSortedList()

	output:
	file("*.png")

	"""
	#zcat /home/porchard/github/ataqv/data/tss/hg19.tss.refseq.bed.gz | grep -v 'hap' > hg19.tss.bed
	cat /lab/work/porchard/sn-muscle-project/data/gencode-tss-for-connections/hg19.bed > hg19.tss.bed
	plot-seurat-snp-to-tss-heatmap.py hg19.tss.bed ${cic.join(' ')}
	"""

}

process plot_cicero_no_binarize_175Mb {

	publishDir "${params.results}/cicero-no-binarize-175Mb"
	memory '20 GB'

	input:
	file(cic) from Channel.fromPath(ROOT + '/work/cicero-no-binarize/results/cicero/*.1750000.cicero.txt').toSortedList()

	output:
	file("*.png")

	"""
	#zcat /home/porchard/github/ataqv/data/tss/hg19.tss.refseq.bed.gz | grep -v 'hap' > hg19.tss.bed
	cat /lab/work/porchard/sn-muscle-project/data/gencode-tss-for-connections/hg19.bed > hg19.tss.bed
	plot-cicero-snp-to-tss-heatmap-2.py hg19.tss.bed ${params.chrom_sizes['hg19']} ${cic.join(' ')}
	"""

}

process plot_cicero_no_binarize_rat_15Mb {

	publishDir "${params.results}/cicero-no-binarize-rat-15Mb"
	memory '20 GB'

	input:
	file(cic) from Channel.fromPath(ROOT + '/work/cicero-no-binarize-rat/results/cicero/*.1500000.cicero.txt').toSortedList()

	output:
	file("*.png")

	"""
	zcat /home/porchard/github/ataqv/data/tss/rn6.tss.refseq.bed.gz | grep -v 'hap' | grep -v 'chrUn' | grep -v '_random'  > rn6.tss.bed
	plot-cicero-snp-to-tss-heatmap-rat.py rn6.tss.bed ${params.chrom_sizes['hg19']} ${cic.join(' ')}
	"""

}

process fans_screenshot {

	publishDir "${params.results}/fans-vs-no-fans"
	executor 'local'
	cache false

	output:
	file("fans-gb.pdf")

	"""
	echo ${FANS_AND_LOADING_POSITION} | perl -pe 's/,//g; s/[:-]/\t/g' > locus.txt
	python3 /home/porchard/src/gb_screenshot_utility/gb_screenshot.py --user porchard --session 2020-sn-muscle-fans-no-nucleus-qc --bed locus.txt --prefix fans --viewer-width 700 --label-width 15 --extension 0 --text-size 24
	mv fans* fans-gb.pdf
	"""

}

process loading_screenshot {

	publishDir "${params.results}/20k-vs-40k"
	executor 'local'
	cache false

	output:
	file("loading-gb.pdf")

	"""
	echo ${FANS_AND_LOADING_POSITION} | perl -pe 's/,//g; s/[:-]/\t/g' > locus.txt
	python3 /home/porchard/src/gb_screenshot_utility/gb_screenshot.py --user porchard --session 2020-sn-muscle-loading-no-nucleus-qc --bed locus.txt --prefix loading --viewer-width 700 --label-width 15 --extension 0 --text-size 24
	mv loading* loading-gb.pdf
	"""

}

process logistic_regression_heatmap {

	publishDir "${params.results}/logistic-regression/both"
	cache 'deep'
	executor 'local'

	input:
	file('human.txt') from Channel.fromPath(ROOT + '/work/logistic-regression/hg19/results/enhancer-regression/model_results.txt')
	file('rat.txt') from Channel.fromPath(ROOT + '/work/logistic-regression/rn6/results/enhancer-regression/model_results.txt')
	file(cluster_names) from Channel.fromPath(CLUSTER_NAMES)
	
	output:
	file("enhancer-similarity-heatmap.pdf")
	file("*.html")

	"""
	plot-logistic-regression-results-both-species.py human.txt rat.txt $cluster_names ${ROOT}/data/roadmap-posteriors/roadmap_cell_types.txt
	"""

}

process plot_umap_by_modality_and_species {

	publishDir "${params.results}/umaps"
	cache 'deep'
	executor 'local'

	input:
	file(umap) from Channel.fromPath(ROOT + '/work/seurat-joint-2/results/seurat-out/seurat-umap.txt')
	file(clusters) from Channel.fromPath(ROOT + '/work/seurat-joint-2/results/seurat-out/seurat-clusters.txt')
	file(library_labels) from Channel.fromPath(ROOT + '/library-labels.txt')

	output:
	file("*.png")

	"""
	plot-umaps.py $umap $clusters $library_labels
	"""

}
/*
process marker_gene_screenshots_rn6 {

	publishDir "${params.results}/gb-screenshots/marker-genes"
	executor 'local'
	
	output:
	file("*.pdf")

	"""
	echo "chr10:53,588,847-53,909,741" >> regions.bed # MYH loci
	echo "chr10:760,194-769,669" >> regions.bed # MYH11
	echo "chr14:35,552,738-35,605,843" >> regions.bed # PDGFRA
	echo "chr4:158,080,951-158,094,049" >> regions.bed # VWF
	echo "chr4:156,742,565-156,775,939" >> regions.bed # CD163
	echo "chr5:158,278,488-158,327,893" >> regions.bed # PAX7
	cat regions.bed | perl -pe 's/,//g; s/[:-]/\t/g' > regions.tmp
	mv regions.tmp regions.bed
	python3 /home/porchard/src/gb_screenshot_utility/gb_screenshot.py --user porchard --session 2020-sn-muscle-cell-types-rn6 --bed regions.bed --prefix rn6 --viewer-width 1200 --label-width 28 --extension 0 --text-size 18
	"""

}

*/

process T2D_gwas_enrichments_hg19_no_bmiadj {

	publishDir "${params.results}/T2D-GWAS-enrichments-no-bmiadj/hg19"
	executor 'local'
	cache false
	errorStrategy 'ignore'
	
	input:
	file(single_human) from Channel.fromPath("${ROOT}/work/ldsc/T2D/hg19/one-model-per-cluster/results/partitioned-heritability/*.results").toSortedList()
	file(joint_human) from Channel.fromPath("${ROOT}/work/ldsc/T2D/hg19/joint-model/results/partitioned-heritability/*.results").toSortedList()

	output:
	file("*.pdf")

	"""
	plot-T2D-enrichment-no-bmiadj.R ${CLUSTER_NAMES} ${single_human.join(' ')} ${joint_human.join(' ')}
	"""

}

process T2D_gwas_enrichments_rn6_no_bmiadj {

	publishDir "${params.results}/T2D-GWAS-enrichments-no-bmiadj/rn6"
	executor 'local'
	cache false
	
	input:
	file(single_human) from Channel.fromPath("${ROOT}/work/ldsc/T2D/rn6/one-model-per-cluster/results/partitioned-heritability/*.results").toSortedList()
	file(joint_human) from Channel.fromPath("${ROOT}/work/ldsc/T2D/rn6/joint-model/results/partitioned-heritability/*.results").toSortedList()

	output:
	file("*.pdf")

	"""
	plot-T2D-enrichment-no-bmiadj.R ${CLUSTER_NAMES} ${single_human.join(' ')} ${joint_human.join(' ')}
	"""

}

process ukb_enrichment_table {

	publishDir "${params.results}/tables"
	executor 'local'

	input:
	file(cluster_names) from Channel.fromPath(CLUSTER_NAMES)
	file(phenotype_tsv) from Channel.fromPath("${ROOT}/data/ukb-summary-stats/ukb31063_h2_all.02Oct2019.tsv.gz")
	file(results) from Channel.fromPath("${ROOT}/work/ldsc/UKB-hesc/hg19/joint-model/results/partitioned-heritability/*.results").toSortedList()

	output:
	file('LDSC-UKB-Z-scores.tsv')

	"""
	make-ldsc-ukb-table.R $cluster_names $phenotype_tsv ${results.join(' ')}
	"""

}

process diamante_overlap_table {

	publishDir "${params.results}/tables"
	executor 'local'

	input:
	file(cluster_names) from Channel.fromPath(CLUSTER_NAMES)
	file(beta_peaks) from Channel.fromPath("${ROOT}/data/other-annotations/beta_ATAC.bed")
	file(adipose_peaks) from Channel.fromPath("${ROOT}/data/other-annotations/Adipose*.broadPeak").toSortedList()
	file(islet_peaks) from Channel.fromPath("${ROOT}/data/other-annotations/HP*").mix(Channel.fromPath("${ROOT}/data/other-annotations/AB*")).mix(Channel.fromPath("${ROOT}/data/other-annotations/AC*")).toSortedList()

	output:
	file('diamante-overlap-summary.csv')

	"""
	cat ${ROOT}/data/diamante-credible-sets/genetic_credible_sets/*.txt | grep -v IndexSNP | awk '{print(\$2, \$3, \$1)}' | perl -pe 's/ /\\t/g; s/^/chr/' > diamante-loci.txt
	mkdir bed
	cp ${ROOT}/data/chromhmm/Islets.all_promoter.bed bed/Islets.TSS_state.bed
	cp ${ROOT}/data/chromhmm/Adipose.all_promoter.bed bed/Adipose.TSS_state.bed
	cp ${ROOT}/data/chromhmm/Liver.all_promoter.bed bed/Liver.TSS_state.bed
	cp ${ROOT}/data/chromhmm/Islets.all_enhancer.bed bed/Islets.enhancer_state.bed
	cp ${ROOT}/data/chromhmm/Adipose.all_enhancer.bed bed/Adipose.enhancer_state.bed
	cp ${ROOT}/data/chromhmm/Liver.all_enhancer.bed bed/Liver.enhancer_state.bed
	cp ${ROOT}/data/gencode-coding/coding.bed bed/
	cat ${adipose_peaks.join(' ')} | sort -k1,1 -k2n,2 | bedtools merge -i stdin > bed/adipose_ATAC.bed
	cat ${islet_peaks.join(' ')} | sort -k1,1 -k2n,2 | bedtools merge -i stdin > bed/islet_ATAC.bed
	cat $beta_peaks > bed/beta_ATAC.bed
	cp ${ROOT}/work/process-by-cluster/results/peaks/broad/0-hg19_peaks.broadPeak.noblacklist bed/cluster_0_ATAC.bed
	cp ${ROOT}/work/process-by-cluster/results/peaks/broad/1-hg19_peaks.broadPeak.noblacklist bed/cluster_1_ATAC.bed
	cp ${ROOT}/work/process-by-cluster/results/peaks/broad/2-hg19_peaks.broadPeak.noblacklist bed/cluster_2_ATAC.bed
	cp ${ROOT}/work/process-by-cluster/results/peaks/broad/3-hg19_peaks.broadPeak.noblacklist bed/cluster_3_ATAC.bed
	cp ${ROOT}/work/process-by-cluster/results/peaks/broad/4-hg19_peaks.broadPeak.noblacklist bed/cluster_4_ATAC.bed
	cp ${ROOT}/work/process-by-cluster/results/peaks/broad/5-hg19_peaks.broadPeak.noblacklist bed/cluster_5_ATAC.bed
	cp ${ROOT}/work/process-by-cluster/results/peaks/broad/6-hg19_peaks.broadPeak.noblacklist bed/cluster_6_ATAC.bed
	annot-snp-bed.py --bed-files bed/* --pos-file diamante-loci.txt > annotated-loci.txt	
	make-diamante-table.R $cluster_names annotated-loci.txt
	"""

}


process ARL15_screenshots {

	publishDir "${params.results}/gb-screenshots/ARL15"
	executor 'local'
	cache false

	output:
	file("*.pdf")

	"""
	echo "chr5:53,269,334-53,278,795" >> regions.bed	
	cat regions.bed | perl -pe 's/,//g; s/[:-]/\t/g' > regions.tmp
	mv regions.tmp regions.bed
	python3 /home/porchard/src/gb_screenshot_utility/gb_screenshot.py --user porchard --session 2020-sn-muscle-GWAS-overlap-muscle-cell-types-only --bed regions.bed --prefix ARL15 --viewer-width 1200 --label-width 24 --extension 0 --text-size 18 --highlight-color FF0000 --highlight-region chr5:53271399-53271440 --genome hg19
	echo "chr5:53,269,196-53,277,711" > regions.bed	
	cat regions.bed | perl -pe 's/,//g; s/[:-]/\t/g' > regions.tmp
	mv regions.tmp regions.bed
	python3 /home/porchard/src/gb_screenshot_utility/gb_screenshot.py --user porchard --session 2020-sn-muscle-GWAS-overlap-ATAC --bed regions.bed --prefix ARL15-other-cell-types --viewer-width 1000 --label-width 24 --extension 0 --text-size 24 --highlight-color FF0000 --highlight-region chr5:53271399-53271440 --genome hg19
	echo "chr2:45977100-45977546" > regions.rn6.bed
	echo "chr2:45977438-45977469" >> regions.rn6.bed
	cat regions.rn6.bed | perl -pe 's/,//g; s/[:-]/\t/g' > regions.tmp
	mv regions.tmp regions.rn6.bed
	python3 /home/porchard/src/gb_screenshot_utility/gb_screenshot.py --user porchard --session 2020-sn-muscle-cell-types-rn6 --genome rn6 --bed regions.rn6.bed --prefix ARL15_rn6 --viewer-width 1200 --label-width 28 --extension 10000 --highlight --highlight-color FF0000 --text-size 18
	"""

}

process ITPR2_screenshots {

	publishDir "${params.results}/gb-screenshots/ITPR2"
	executor 'local'
	cache false

	output:
	file("*.pdf")

	"""
	echo "chr12:26,436,640-26,491,955" >> regions.bed
	echo "chr12:26,469,532-26,475,063" >> regions.bed
	cat regions.bed | perl -pe 's/,//g; s/[:-]/\t/g' > regions.tmp
	mv regions.tmp regions.bed
	python3 /home/porchard/src/gb_screenshot_utility/gb_screenshot.py --user porchard --session 2020-sn-muscle-GWAS-overlap-muscle-cell-types-only --bed regions.bed --prefix ITPR2 --viewer-width 1200 --label-width 24 --extension 0 --text-size 18 --highlight-color FF0000 --highlight-region chr12:26472541-26472582 --genome hg19
	python3 /home/porchard/src/gb_screenshot_utility/gb_screenshot.py --user porchard --session 2020-sn-muscle-GWAS-overlap-ATAC --bed regions.bed --prefix ITPR2-other-cell-types --viewer-width 1200 --label-width 24 --extension 0 --text-size 18 --highlight-color FF0000 --highlight-region chr12:26472541-26472582 --genome hg19
	echo "chr4:180399117-180400404" > regions.rn6.bed
	echo "chr4:180400064-180400103" >> regions.rn6.bed
	cat regions.rn6.bed | perl -pe 's/,//g; s/[:-]/\t/g' > regions.tmp
	mv regions.tmp regions.rn6.bed
	python3 /home/porchard/src/gb_screenshot_utility/gb_screenshot.py --user porchard --session 2020-sn-muscle-cell-types-rn6 --bed regions.rn6.bed --genome rn6 --prefix ITPR2_rn6 --viewer-width 1200 --label-width 28 --extension 10000 --highlight --highlight-color FF0000 --text-size 18
	"""

}


process chromatin_state_legend {

	publishDir "${params.results}/chromatin-state-legend"
	executor 'local'

	input:
	file(bed) from Channel.fromPath("${ROOT}/data/chromhmm/Adipose.dense.bed")

	output:
	file("*.pdf")

	"""
	make_chromatin_state_legend.R $bed
	"""	

}


process ataqv_metrics {

	publishDir "${params.results}/fans-vs-no-fans", pattern: "fans*"
	publishDir "${params.results}/20k-vs-40k", pattern: "loading*"
	executor 'local'

	input:
	file(ataqv) from Channel.fromPath("${ROOT}/work/process-as-bulk/results/ataqv/*.json.gz").mix(Channel.fromPath("${ROOT}/work/bulk-atacseq/results/ataqv/320-NM-*.ataqv.json.gz")).toSortedList()

	output:
	file("*.pdf")

	"""
	cp ${ROOT}/library-labels.txt .
	extractAtaqvMetric.py --files ${ataqv.join(' ')} --metrics total_peaks hqaa_overlapping_peaks_percent | perl -pe 's/.ataqv.json.gz//' | cut -f1,3,4 > metrics.txt
	extractFragmentLengths.py --files ${ataqv.join(' ')} | perl -pe 's/.ataqv.json.gz//' > fragment-lengths.txt
	extractTssCoverage.py --files ${ataqv.join(' ')} | perl -pe 's/.ataqv.json.gz//' > tss-coverage.txt
	plot-ataqv-metrics.R
	"""

}

process chromatin_state_overlap {
	
	container "${params.containers.general}"
	executor 'local'

	input:
	file(peaks) from Channel.fromPath("${ROOT}/work/process-as-bulk/results/macs2/*.noblacklist")

	output:
	file("${library}.chromhmm_overlap.txt") into plot_aggregate_chromatin_state_overlap_in

	script:
	genome = 'hg19'
	library = peaks.getName().replaceAll('_peaks.broadPeak.noblacklist', '')

	"""
	chromhmm_overlap.py --tss-extension 5000 $peaks ${get_tss(genome)} ${get_chrom_sizes(genome)} ${params.chromatin_state_glob} > ${library}.chromhmm_overlap.txt
	"""

}


process plot_aggregate_chromatin_state_overlap {
	
	publishDir "${params.results}/fans-vs-no-fans", pattern: "fans*"
	publishDir "${params.results}/20k-vs-40k", pattern: "loading*"
	container "${params.containers.rplot}"
	executor 'local'

	input:
	file(coverage) from plot_aggregate_chromatin_state_overlap_in.toSortedList()

	output:
	file("*.pdf")

	"""
	plot_chromatin_state_overlap_2.R ${ROOT}/library-labels.txt ${coverage.join(' ')}
	"""

}

process per_library_gene_counts {
	
	executor 'local'

	input:
	file(x) from Channel.fromPath("${ROOT}/work/counts/results/uncorrected-counts/*.features.txt").toSortedList()

	output:
	file("per-library-gene-counts.txt") into per_library_gene_counts_out

	"""
	aggregate-counts-per-library.py ${x.join(' ')} > per-library-gene-counts.txt
	"""

}

process rnaseq_aggregate_correlation {
	
	executor 'local'
	
	publishDir "${params.results}/fans-vs-no-fans", pattern: "fans*"
	publishDir "${params.results}/20k-vs-40k", pattern: "loading*"
       container "${params.containers.rplot}"

	input:
	file(x) from per_library_gene_counts_out

	output:
	file("*.png")

	"""
	cp ${ROOT}/library-labels.txt .
	rnaseq-aggregate-correlation.R $x
	"""

}

process atac_bulk_correlation_fans {

	publishDir "${params.results}/fans"
	container "${params.containers.rplot}"
	executor 'local'

	input:
	file(x) from Channel.fromPath("${ROOT}/work/atac-correlation-fans/results/counts/*.bed").toSortedList()

	output:
	file("fans-atac-correlation.png")

	"""
	atacseq-aggregate-correlation.R ${ROOT}/library-labels.txt ${x.join(' ')}
	ln -s all-atac-correlation.png fans-atac-correlation.png
	"""

}

process atac_bulk_correlation_loading {

	publishDir "${params.results}/loading"
	container "${params.containers.rplot}"
	executor 'local'

	input:
	file(x) from Channel.fromPath("${ROOT}/work/atac-correlation-loading/results/counts/*.bed").toSortedList()

	output:
	file("loading-atac-correlation.png")

	"""
	atacseq-aggregate-correlation.R ${ROOT}/library-labels.txt ${x.join(' ')}
	ln -s all-atac-correlation.png loading-atac-correlation.png
	"""

}

process plot_qc_thresholds_rna {

       publishDir "${params.results}/fans-vs-no-fans", pattern: "*-fans-vs-no-fans.png"
       publishDir "${params.results}/20k-vs-40k", pattern: "*-20k-vs-40k.png"
       publishDir "${params.results}/qc-for-downstream-libraries", pattern: "*-used-downstream.png"
       executor 'local'
       container "${params.containers.rplot}"

       input:
       file(qc) from Channel.fromPath("${ROOT}/work/rnaseq-qc/results/rnaseq-qc/*.txt").toSortedList()
       file(initial_thresholds) from Channel.fromPath("${ROOT}/initial-thresholds-rna.txt")
       file(library_labels) from Channel.fromPath("${ROOT}/library-labels.txt")

       output:
       file("*.png")

       """
       call-nuclei-rna-with-mitochondrial.R $initial_thresholds $library_labels ${qc.join(' ')}
       """

}
