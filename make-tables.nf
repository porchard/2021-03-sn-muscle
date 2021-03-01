#!/usr/bin/env nextflow

process make_atac_qc_table_and_figures {
	
	publishDir "${params.results}/figures"
	
	output:
	file("*.png")

	"""
	cp /lab/work/porchard/sn-muscle-project/work/qc/work/46/c20517dbc1b8a50b5c6d390a752c8a/metrics.txt .
	cp /lab/work/porchard/sn-muscle-project/work/qc/work/46/c20517dbc1b8a50b5c6d390a752c8a/*.best.txt .
	cp /lab/work/porchard/sn-muscle-project/work/qc/work/46/c20517dbc1b8a50b5c6d390a752c8a/initial-thresholds-atac.txt .
	cp /lab/work/porchard/sn-muscle-project/work/qc/work/46/c20517dbc1b8a50b5c6d390a752c8a/library-labels.txt .
	atac-qc-make-table.py --thresholds initial-thresholds-atac.txt --demuxlet 63_40-hg19.best.txt 63_20-hg19.best.txt --library-labels library-labels.txt --metrics metrics.txt
	"""

}

process make_atac_qc_table {

	publishDir "${params.results}/tables"

	output:
	file("*.csv")

	"""
	make-atac-qc-threshold-supplementary-table.py /lab/work/porchard/sn-muscle-project/initial-thresholds-atac.txt > atac-qc-thresholds.csv
	echo "Supplemental Table S1: snATAC-seq per-nucleus QC thresholds." >> atac-qc-thresholds.csv
	"""
}

process make_rna_qc_table {

	publishDir "${params.results}/tables"

	output:
	file("*.csv")

	"""
	make-rna-qc-threshold-supplementary-table.py /lab/work/porchard/sn-muscle-project/initial-thresholds-rna.txt > rna-qc-thresholds.csv
	echo "Supplemental Table S2: snRNA-seq per-nucleus QC thresholds." >> rna-qc-thresholds.csv
	"""
}

process make_tables_s3_and_s5 {

	publishDir "${params.results}/tables"

	output:
	file("*.tsv")

	"""
	make-table-s3-and-s5.py
	echo "Supplemental Table S3: Per library summary statistics (number nuclei per sample, mean and median UMIs/ATAC fragments per library)" >> per-library-summary-stats.tsv
	echo "Supplemental Table S5: Nucleus counts per cell type, species, and modality." >> per-species-per-modality-nucleus-counts.tsv
	"""

}

process make_deltasvm_performance_table {

	publishDir "${params.results}/tables"

	output:
	file("gkmsvm-model-results.csv")

	"""
	make-gkmsvm-model-performance-table.py /lab/work/porchard/sn-muscle-project/work/delta-svm/all-keep-40k/results/cluster_*/model_scores.tsv > gkmsvm-model-results.csv
	"""

}
