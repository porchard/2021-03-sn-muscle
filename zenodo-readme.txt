* atac-bedgraphs contains ATAC signal tracks for bulk data, single nucleus data processed as bulk (included in the 'bulk' subdirectory), and single nucleus data aggregated within cell type clusters (for human and rat). To reduce repo size, these files are 'starched' with bedops v. 2.4.26
* atac-peaks contains MACS2 ATAC broadpeak calls, organized in a similar manner as the atac-bigwigs directory
* clustering contains the per-nucleus cluster assignments (columns: library, barcode, cluster) and the umap coordinates (library, barcode, umap dim. 1, umap dim. 2)
* gene-counts includes per-gene counts used as input for LIGER clustering. There is one matrix per modality and biological sample.
* filtered-read-locations directory contain starched bed-like files. Note that we've removed all alignment CIGAR information for privacy purposes; this may cause read locations to appear much longer than actual read length, especially for RNA data (due to splicing):
** For ATAC: read mapping locations after filtering to high-quality alignments and duplicate removal. For single-nucleus libraries, corrected nucleus barcode is encoded in read name.
** For RNA: read mapping locations after filtering, but before duplicate removal. Corrected nucleus barcode and UMI are encoded (in that order) in read names.
* luciferase contains an Excel spreadsheet of luciferase results for SNP rs702634
