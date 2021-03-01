###### CHANGE THESE PATHS
ROOT=/lab/work/porchard/sn-muscle-project
HG19_FASTA=/lab/data/reference/human/hg19/hg19.fa
HG19_CHROM_SIZES=/lab/data/reference/human/hg19/hg19.chrom_sizes
HG19_BLACKLISTS=/lab/data/reference/human/hg19/annot/wgEncodeDacMapabilityConsensusExcludable.bed.gz /lab/data/reference/human/hg19/annot/wgEncodeDukeMapabilityRegionsExcludable.bed.gz
RN6_FASTA=/lab/work/porchard/data/fasta/rn6.fa
RN6_CHROM_SIZES=/lab/work/porchard/data/chrom_sizes/rn6.chrom_sizes
RN6_BLACKLISTS=/lab/work/porchard/data/mappability.new/rn6.blacklist.1.bed.gz
HG38_TO_HG19_CHAIN=/lab/work/porchard/data/chain/hg38ToHg19.over.chain.gz

##### 

WORK=$(ROOT)/work
DATA=$(ROOT)/data
FIGURES=$(ROOT)/figures
SRC=$(ROOT)/src
BIN=$(ROOT)/bin
CONTROL=$(ROOT)/control
SAMPLE_INFO=$(ROOT)/sample_info/sample_info.txt
NAME = sn-muscle-project
CLUSTER_NAMES=$(ROOT)/2019-01-04-cluster-names.txt
SHARED_OPEN_CHROMATIN=$(WORK)/meuleman-with-hesc/results/common/common_open_chromatin.bed
PATH_TO_BULK=$(WORK)/bulk-atacseq/results

.PHONY: data sample_info all

ANALYSIS = $(WORK)/$@
CONFIG = $(ANALYSIS)/config.json
PIPELINE = $(ANALYSIS)/pipeline

##### https://stackoverflow.com/questions/7039811/how-do-i-process-extremely-long-lists-of-files-in-a-make-recipe
define NL


endef
#####

### SET UP DATA FILES ###
data: repeats liger-features dbSNP-vcf orthologues roadmap-posteriors ukb-summary-stats chromhmm other-annotations reformat-ukbb ldsc-baseline diamante-summary-stats manning-summary-stats gencode-coding gencode-bed meuleman-data meuleman encode-cisregulatory-elements encode-cres 1000G-SNPs motifs

cisbp-motifs: ANALYSIS=$(DATA)/$@
cisbp-motifs:
	mkdir -p $(ANALYSIS)
	cat /lab/work/porchard/snp-aware-pwm-scanning/data/TF_Information.txt | cut -f4,7 | awk '$$1!="."' | grep -v -w Motif_ID > $(ANALYSIS)/motif-id-to-tf-name.txt


gencode-tss-for-connections: ANALYSIS=$(DATA)/$@
gencode-tss-for-connections:
	mkdir -p $(ANALYSIS)
	cp /lab/data/reference/human/hg19/annot/gencode.v19.annotation.gtf.gz $(ANALYSIS)/
	cd $(ANALYSIS) && python $(BIN)/gencodeGTFtoTSS.py gencode.v19.annotation.gtf.gz | sort | uniq | sort -k1,1 -k2n,2 > $(ANALYSIS)/hg19.bed

repeats: ANALYSIS=$(DATA)/$@
repeats:
	mkdir -p $(ANALYSIS)
	mysql --host=genome-mysql.cse.ucsc.edu --user=genome -D hg19 -e "SELECT * FROM simpleRepeat" > $(ANALYSIS)/trf_table.txt
	cut -f2-4 $(ANALYSIS)/trf_table.txt | grep -v chromStart | sort -k1,1 -k2n,2 | bedtools merge -i stdin > $(ANALYSIS)/trf.bed


meuleman-data: ANALYSIS=$(DATA)/meuleman
meuleman-data:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && wget https://www.meuleman.org/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz
	cd $(ANALYSIS) && wget https://www.meuleman.org/DHS_Index_and_Vocabulary_metadata.tsv

# DONE
meuleman-with-hesc:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --dhs_index $(DATA)/meuleman/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz --chain $(HG38_TO_HG19_CHAIN) --metadata $(DATA)/meuleman/DHS_Index_and_Vocabulary_metadata.tsv --dhs_matrix $(DATA)/meuleman/dat_bin_FDR01_hg38.RData $(ROOT)/prep-meuleman-with-hesc.nf &

# DONE
encode-cisregulatory-elements: ANALYSIS=$(DATA)/$@
encode-cisregulatory-elements:
	cd $(ANALYSIS) && bash commands

encode-cres:
	mkdir -p $(ANALYSIS)
	echo "#!/bin/bash" > $(PIPELINE)
	echo "#SBATCH --mem=10G" >> $(PIPELINE)
	cd $(ANALYSIS) && echo "python $(ROOT)/bin/sort_uniq_gzip.py $(DATA)/encode-cisregulatory-elements/*.bed.gz | sort -k1,1 -k2n,2 > encode-cres.bed" >> $(PIPELINE) && sbatch $(PIPELINE)

gencode-coding: ANALYSIS=$(DATA)/$@
gencode-coding:
	mkdir $(ANALYSIS) && cd $(ANALYSIS) && wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz && zcat gencode.v19.annotation.gtf.gz | awk '$$3=="CDS"' | cut -f1,4,5 | sort -k1,1 -k2n,2 -k3n,3 | bedtools merge -i stdin > coding.bed

gencode-bed: ANALYSIS=$(DATA)/$@
gencode-bed:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gff3.gz && gunzip gencode.v19.annotation.gff3.gz
	gff3ToGenePred $(ANALYSIS)/gencode.v19.annotation.gff3 $(ANALYSIS)/gencode.genePred
	genePredToBed $(ANALYSIS)/gencode.genePred $(ANALYSIS)/gencode.hg19.bed
	cd $(ANALYSIS) && wget ftp://ftp.ensembl.org/pub/release-101/gff3/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.101.gff3.gz && zcat Rattus_norvegicus.Rnor_6.0.101.gff3.gz > rn6.gff3
	gff3ToGenePred -useName $(ANALYSIS)/rn6.gff3 $(ANALYSIS)/rn6.genePred
	genePredToBed $(ANALYSIS)/rn6.genePred $(ANALYSIS)/rn6.bed
	cat $(ANALYSIS)/rn6.bed | perl -pe 's/^/chr/' > $(ANALYSIS)/rn6.chromnames.bed
	cd $(ANALYSIS) && python update-bed-to-gene-names.py gencode.v19.annotation.gff3 gencode.hg19.bed > gencode.hg19.genes.bed

motifs: ANALYSIS=$(DATA)/$@
motifs: MOTIFS=MEF2_known10 AP1_known5 PITX2_2
motifs:
	mkdir -p $(ANALYSIS)
	cp /lab/work/porchard/fimo/work/background/hg19.background.txt $(ANALYSIS)/
	$(foreach m,$(MOTIFS),cp /lab/work/porchard/fimo/work/ENCODE2013/data/$(m).meme $(ANALYSIS)/$(NL))
	$(foreach m,$(MOTIFS),meme2plain $(ANALYSIS)/$(m).meme | grep -v "^>" > $(ANALYSIS)/$(m).txt$(NL))

liger-features: ANALYSIS=$(DATA)/$@
liger-features:
	# For each genome:
	# get the GTF
	# get the annotation report
	# convert the chromosome names to ucsc style
	# then, using python script:
	# keep curated transcripts
	# for each gene, merge all the transcripts, and keep 3 kb upstream (even if it overlaps with other stuff)
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Rattus_norvegicus/annotation_releases/current/GCF_000001895.5_Rnor_6.0/GCF_000001895.5_Rnor_6.0_assembly_report.txt && wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Rattus_norvegicus/annotation_releases/current/GCF_000001895.5_Rnor_6.0/GCF_000001895.5_Rnor_6.0_genomic.gtf.gz
	cd $(ANALYSIS) && wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/105.20190906/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz && wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/105.20190906/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_assembly_report.txt
	python $(ROOT)/bin/refseq-annotation-and-report-to-feature-file.py $(ANALYSIS)/GCF_000001895.5_Rnor_6.0_genomic.gtf.gz $(ANALYSIS)/GCF_000001895.5_Rnor_6.0_assembly_report.txt | grep -v -e 'random' -e 'chrUn' | sort -k1,1 -k2n,2 -k3n,3 > $(ANALYSIS)/rn6.genes.bed
	python $(ROOT)/bin/make-liger-features.py $(ANALYSIS)/rn6.genes.bed $(RN6_CHROM_SIZES) | sort -k1,1 -k2n,2 -k3n,3 > $(ANALYSIS)/rn6.features.bed
	python $(ROOT)/bin/refseq-annotation-and-report-to-feature-file.py $(ANALYSIS)/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz $(ANALYSIS)/GCF_000001405.25_GRCh37.p13_assembly_report.txt | grep -v -e 'random' -e 'chrUn' | sort -k1,1 -k2n,2 -k3n,3 > $(ANALYSIS)/hg19.genes.bed
	python $(ROOT)/bin/make-liger-features.py $(ANALYSIS)/hg19.genes.bed $(HG19_CHROM_SIZES) | sort -k1,1 -k2n,2 -k3n,3 > $(ANALYSIS)/hg19.features.bed
	# blacklist filter
	bedtools intersect -a $(ANALYSIS)/rn6.features.bed -b $(RN6_BLACKLISTS) -v > $(ANALYSIS)/rn6.features.noblacklist 
	bedtools intersect -a $(ANALYSIS)/hg19.features.bed -b $(HG19_BLACKLISTS) -v > $(ANALYSIS)/hg19.features.noblacklist 

diamante-summary-stats: ANALYSIS=$(DATA)/$@
diamante-summary-stats:
	mkdir -p $(ANALYSIS)
	ln -sf /lab/work/porchard/data/gwas/diamante/Mahajan.NatGenet2018b.T2D.European.txt $(ANALYSIS)/ # TODO: drop from final file
	ln -sf /lab/work/porchard/data/gwas/diamante/Mahajan.NatGenet2018b.T2Dbmiadj.European.txt $(ANALYSIS)/ # TODO: drop from final file
	cd $(ANALYSIS) && nohup python $(ROOT)/bin/prep-diamante-for-ldsc.py --vcf $(DATA)/dbSNP-vcf/All_20170710.vcf.gz $(ANALYSIS)/Mahajan.NatGenet2018b.T2D.European.txt diamante.T2D.European.ldsc.txt.gz &
	cd $(ANALYSIS) && nohup python $(ROOT)/bin/prep-diamante-for-ldsc.py --vcf $(DATA)/dbSNP-vcf/All_20170710.vcf.gz $(ANALYSIS)/Mahajan.NatGenet2018b.T2Dbmiadj.European.txt diamante.T2Dbmiadj.European.ldsc.txt.gz &

manning-summary-stats: ANALYSIS=$(DATA)/$@
manning-summary-stats:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && wget ftp://ftp.sanger.ac.uk/pub/magic/MAGIC_Manning_et_al_lnFastingInsulin_MainEffect.txt.gz && wget ftp://ftp.sanger.ac.uk/pub/magic/MAGIC_Manning_et_al_FastingGlucose_MainEffect.txt.gz
	cd $(ANALYSIS) && Rscript $(ROOT)/bin/reformat-manning-summarystats-for-ldsc.R

ukb-summary-stats: ANALYSIS = $(DATA)/$@
ukb-summary-stats:
	cd $(ANALYSIS) && wget -O ukb31063_h2_all.02Oct2019.tsv.gz https://www.dropbox.com/s/ipeqyhrpdqav5uh/ukb31063_h2_all.02Oct2019.tsv.gz?dl=1
	cd $(ANALYSIS) && Rscript $(ROOT)/bin/choose-ukb-traits-new.R ukb31063_h2_all.02Oct2019.tsv.gz manifest.tsv && echo "# drmr:job" > wget-commands.sh && cut -f4 manifest-subset.tsv | perl -pe 's/.bgz$$/.gz/' >> wget-commands.sh && drmrarray -s 20 wget-commands.sh

ldsc-data: ANALYSIS=$(DATA)/$@
ldsc-data:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baselineLD_v2.2_ldscores.tgz
	cd $(ANALYSIS) && wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_frq.tgz && tar -xvzf 1000G_Phase3_frq.tgz
	cd $(ANALYSIS) && wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz && tar -xvzf 1000G_Phase3_plinkfiles.tgz
	cd $(ANALYSIS) && wget https://data.broadinstitute.org/alkesgroup/LDSCORE/hapmap3_snps.tgz && tar -xvzf hapmap3_snps.tgv
	cd $(ANALYSIS) && wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_weights_hm3_no_MHC.tgz && tar -xvzf 1000G_Phase3_weights_hm3_no_MHC.tgz
	cd $(ANALYSIS) && wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2 && bzip2 -d w_hm3.snplist.bz2
	cd $(ANALYSIS) && wget https://data.broadinstitute.org/alkesgroup/LDSCORE/weights_hm3_no_hla.tgz && tar -xvzf weights_hm3_no_hla.tgz
	cd $(ROOT)/bin && python make-snp-id-conversions.py

ldsc-baseline-min-maf:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume -with-trace -with-report --results $(ANALYSIS)/results --projroot $(ROOT) $(ROOT)/ldsc-update-baseline-model-min-maf.nf &

chromhmm:
	mkdir -p $(DATA)/$@
	cd $(DATA)/$@ && wget https://zenodo.org/record/3524578/files/islet-cage-zenodo.tar.gz && tar -xvzf islet-cage-zenodo.tar.gz
	cd $(DATA)/$@ && cp data/chromhmm/selected_annotated_states/files_by_state/* .
	cd $(DATA)/$@ && rename 's/cell4_11.//' *
	$(foreach t,Adipose SkeletalMuscle Liver Islets,cd $(DATA)/$@ && python $(BIN)/make-chromhmm-dense.py $(foreach s,Active_enhancer Active_TSS Bivalent_poised_TSS Flanking_TSS Genic_enhancer Quiescent_low_signal Repressed_polycomb Strong_transcription Weak_enhancer Weak_transcription Weak_TSS,$(t).$(s).bed) | sort -k1,1 -k2n,2 > $(t).dense.bed$(NL))
	cd $(DATA)/$@ && rm -r data

other-annotations:
	mkdir -p $(DATA)/$@
	ln -s /lab/work/vivekrai/2019_mohlke_adipose/work/macs2/Adipose*.broadPeak $(DATA)/$@/
	python $(ROOT)/bin/peak-sharing.py $(DATA)/$@/Adipose* | awk '$$4>=2' | cut -f1-3 | bedtools intersect -a stdin -b /lab/work/porchard/sn-muscle-project/data/mappability/hg19* -v > adipose.bed
	cat /home/vivekrai/analyses/2018_NIH_Islets.snatacseq.v2/work/2019-03-01_clustering-final/peaks/1_peaks.broadPeak.noblacklist | cut -f1-3 | sort -k1,1 -k2n,2 > $(DATA)/$@/beta_ATAC.bed
	cp /lab/work/vivekrai/2017_NIH_Islets.atacseq/work/2019-05-09_process-samples-subset/macs2/*.noblacklist $(DATA)/$@
	cd $(DATA)/$@ && rm EndoC*

roadmap-posteriors:
	cd $(DATA)/$@ && nohup bash $(ROOT)/src/download-roadmap-posteriors.sh &

reformat-ukbb:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --projroot $(ROOT) -with-singularity $(ROOT)/singularity/general/general.simg --results $(ANALYSIS) $(ROOT)/reformat-ukbb.nf &

rnaseq:
	mkdir -p $(ANALYSIS)
	python $(CONTROL)/rnaseq/make_config.py $(ROOT) $(ROOT)/sample_info/sample_info.txt > $(CONFIG)
	cd $(ANALYSIS) && nohup nextflow run -resume --chemistry V3 -with-trace -params-file $(CONFIG) -with-singularity /lab/work/porchard/singularity/archive/snRNA/2019-11-15/snRNA.simg -with-report -qs 300 --results $(ANALYSIS)/results /home/porchard/github/snRNAseq-NextFlow/main.nf &

rnaseq-dual-modality:
	mkdir -p $(ANALYSIS)
	ln -s /lab/work/porchard/Nova-315/work/rnaseq/results $(ANALYSIS)/

bulk-atacseq:
	mkdir -p $(ANALYSIS)
	python $(CONTROL)/$@/make_config.py $(ROOT) > $(ANALYSIS)/config.json
	cd $(ANALYSIS) && nohup nextflow run -resume -with-singularity /lab/work/porchard/singularity/archive/ATAC/2019-11-20/ATAC.simg -with-dag -with-timeline -with-trace -with-report -params-file $(ANALYSIS)/config.json --results $(ANALYSIS)/results /home/porchard/github/ATACseq-NextFlow/main.nf &

atacseq:
	mkdir -p $(ANALYSIS)
	python $(CONTROL)/$@/make_config.py $(ROOT) $(ROOT)/sample_info/sample_info.with_facs.txt > $(CONFIG)
	cp $(ROOT)/atac.nextflow.config $(ANALYSIS)/nextflow.config
	cd $(ANALYSIS) && nohup nextflow run -resume -with-trace -params-file $(CONFIG) -with-singularity /lab/work/porchard/singularity/archive/snATAC/2019-11-07/snATAC.simg --low_read_count_threshold 1000 -with-report -qs 300 --results $(ANALYSIS)/results /home/porchard/github/snATACseq-NextFlow/main.nf &

atacseq-dual-modality:
	mkdir -p $(ANALYSIS)
	ln -s /lab/work/porchard/Nova-303/work/atacseq/results $(ANALYSIS)/

process-as-bulk:
	mkdir -p $(ANALYSIS)/data/bams
	cp $(WORK)/atacseq/results/merge/133* $(ANALYSIS)/data/bams/
	cp $(WORK)/atacseq/results/merge/63_* $(ANALYSIS)/data/bams/
	cp $(ROOT)/atac.nextflow.config $(ANALYSIS)/nextflow.config
	cd $(ANALYSIS) && nohup nextflow run -resume --mapped_bam_glob '$(ANALYSIS)/data/bams/*.bam' -with-singularity $(ROOT)/singularity/ATAC.simg --results $(ANALYSIS)/results $(ROOT)/process-snATAC-as-bulk.nf &

counts-rna: ANALYSIS=$(WORK)/counts
counts-rna:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --projroot $(ROOT) $(ROOT)/counts.nf &

fragment-counts-atac:
	mkdir -p $(ANALYSIS)/data
	ln -s $(WORK)/atacseq/results/prune/*.bam $(ANALYSIS)/data/
	ln -s $(WORK)/atacseq-dual-modality/results/prune/*.bam $(ANALYSIS)/data/
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --projroot $(ROOT) --bam_glob '$(ANALYSIS)/data/*.bam' $(ROOT)/fragment-counts-atac.nf &

rnaseq-qc:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --starsolo_dir_glob '$(WORK)/rnaseq/results/starsolo/*' $(ROOT)/make-rna-qc.nf &

qc:
	mkdir -p $(ANALYSIS)
	python $(ROOT)/control/$@/make_config.py $(ROOT) > $(ANALYSIS)/config.json
	cd $(ANALYSIS) && nohup nextflow run -resume --rna_metrics '$(WORK)/rnaseq-qc/results/rnaseq-qc/*.txt' --initial_thresholds_rna $(ROOT)/initial-thresholds-rna.txt --initial_thresholds_atac $(ROOT)/initial-thresholds-atac.txt --library_labels $(ROOT)/library-labels.txt --results $(ANALYSIS)/results --genotypes $(DATA)/genotypes/genotypes.noblacklist.vcf.gz --demuxlet_mask $(WORK)/demuxlet-mask/demuxlet-mask.bed -params-file config.json $(ROOT)/$@.nf &


atac-correlation-fans: LIBRARIES=$(shell seq 133151 133154)
atac-correlation-fans: BULK_LIBRARIES=320-NM-1 320-NM-2
atac-correlation-fans:
	mkdir -p $(ANALYSIS)/data/peaks
	mkdir -p $(ANALYSIS)/data/bams
	$(foreach l,$(LIBRARIES),ln -sf $(WORK)/process-as-bulk/results/prune/$(l)-hg19.pruned.bam $(ANALYSIS)/data/bams/$(NL))
	$(foreach l,$(LIBRARIES),ln -sf $(WORK)/process-as-bulk/results/macs2/$(l)-hg19_peaks.broadPeak.noblacklist $(ANALYSIS)/data/peaks/$(NL))
	$(foreach l,$(BULK_LIBRARIES),ln -sf $(PATH_TO_BULK)/macs2/$(l)_peaks.broadPeak.noblacklist $(ANALYSIS)/data/peaks/$(NL))
	$(foreach l,$(BULK_LIBRARIES),ln -sf $(PATH_TO_BULK)/prune/$(l).pruned.bam $(ANALYSIS)/data/bams/$(NL))
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --bam_glob '$(ANALYSIS)/data/bams/*' --peak_glob '$(ANALYSIS)/data/peaks/*' $(ROOT)/atac-correlation.nf &

atac-correlation-loading: LIBRARIES=63_20 63_40
atac-correlation-loading: BULK_LIBRARIES=320-NM-1 320-NM-2 320-NM-3 320-NM-4
atac-correlation-loading:
	mkdir -p $(ANALYSIS)/data/peaks
	mkdir -p $(ANALYSIS)/data/bams
	$(foreach l,$(LIBRARIES),ln -sf $(WORK)/process-as-bulk/results/prune/$(l)-hg19.pruned.bam $(ANALYSIS)/data/bams/$(NL))
	$(foreach l,$(LIBRARIES),ln -sf $(WORK)/process-as-bulk/results/macs2/$(l)-hg19_peaks.broadPeak.noblacklist $(ANALYSIS)/data/peaks/$(NL))
	$(foreach l,$(BULK_LIBRARIES),ln -sf $(PATH_TO_BULK)/macs2/$(l)_peaks.broadPeak.noblacklist $(ANALYSIS)/data/peaks/$(NL))
	$(foreach l,$(BULK_LIBRARIES),ln -sf $(PATH_TO_BULK)/prune/$(l).pruned.bam $(ANALYSIS)/data/bams/$(NL))
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --bam_glob '$(ANALYSIS)/data/bams/*' --peak_glob '$(ANALYSIS)/data/peaks/*' $(ROOT)/atac-correlation.nf &

project-human-snps-to-rat:
	mkdir -p $(ANALYSIS)
	cat $(DATA)/diamante-credible-sets/genetic_credible_sets/* | grep -v Pos | cut -f2,3 | awk '{print($$1, $$2-20, $$2+20)}' | perl -pe 's/^/chr/; s/ /\t/g' > $(ANALYSIS)/snps-and-flanking.bed
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --bed_glob '$(ANALYSIS)/snps-and-flanking.bed' $(ROOT)/find-rat-snp.nf &

project-human-snps-to-rat-2:
	mkdir -p $(ANALYSIS)
	cat $(DATA)/diamante-credible-sets/genetic_credible_sets/* | grep -v Pos | cut -f2,3 | awk '{print($$1, $$2-1, $$2)}' | perl -pe 's/^/chr/; s/ /\t/g' > $(ANALYSIS)/snps-and-flanking.bed
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --bed_glob '$(ANALYSIS)/snps-and-flanking.bed' $(ROOT)/find-rat-snp.nf &

genotypes: ANALYSIS = $(DATA)/genotypes
genotypes:
	mkdir -p $(ANALYSIS)
	echo "7_KSm2 KSM2" > $(ANALYSIS)/rename.txt
	echo "28_Ksm1 KSM1" >> $(ANALYSIS)/rename.txt
	rm -rf rename-chroms.txt
	$(foreach c,$(shell seq 1 22 | sort),echo '$(c) chr$(c)' >> rename-chroms.txt$(NL))
	bcftools concat $(foreach c,$(shell seq 1 22 | sort),/lab/work/arushiv/muscle-sn/analyses/genotype_imputation/results/imputed_filtered_others/chr$(c).genotypes-rsid.vcf.gz) | bcftools view --samples 7_KSm2,28_Ksm1 | bcftools reheader --samples $(ANALYSIS)/rename.txt | bcftools annotate --rename-chrs rename-chroms.txt | bcftools view -o $(ANALYSIS)/genotypes.vcf.gz -Ob
	rm -rf rename-chroms.txt
	zcat $(DATA)/mappability/hg19.blacklist.* | sort -k1,1 -k2n,2 | bedtools merge > $(ANALYSIS)/blacklist.bed
	bcftools view --targets-file ^$(ANALYSIS)/blacklist.bed -Ob -o $(ANALYSIS)/genotypes.noblacklist.vcf.gz $(ANALYSIS)/genotypes.vcf.gz
	rm $(ANALYSIS)/blacklist.bed

diamante-locuszoom:
	mkdir -p $(ANALYSIS)
	printf "MarkerName\tP-value\n" > $(ANALYSIS)/gwas.txt
	cat /lab/work/porchard/data/gwas/diamante/Mahajan.NatGenet2018b.T2D.European.txt | cut -f1,9 | grep -v Pvalue | perl -pe 's/^5:53271420\t/rs702634\t/; s/^12:26472562\t/rs7132434\t/' >> $(ANALYSIS)/gwas.txt
	cd $(ANALYSIS) && locuszoom --metal $(ANALYSIS)/gwas.txt --refsnp rs702634 --flank 300kb --build hg19 --pop EUR --source 1000G_Nov2014 theme=publication height=6 hiStart=53269334 hiEnd=53278795
	cd $(ANALYSIS) && locuszoom --metal $(ANALYSIS)/gwas.txt --refsnp rs7132434 --flank 300kb --build hg19 --pop EUR --source 1000G_Nov2014 theme=publication height=6 hiStart=26436640 hiEnd=26491955

diamante-bmiadj-locuszoom:
	mkdir -p $(ANALYSIS)
	printf "MarkerName\tP-value\n" > $(ANALYSIS)/gwas.txt
	cat /lab/work/porchard/data/gwas/diamante/Mahajan.NatGenet2018b.T2Dbmiadj.European.txt | cut -f1,9 | grep -v Pvalue >> $(ANALYSIS)/gwas.txt
	cd $(ANALYSIS) && locuszoom --metal $(ANALYSIS)/gwas.txt --refsnp 5:53271420 --flank 300kb --build hg19 --pop EUR --source 1000G_Nov2014 theme=publication height=6 hiStart=53269334 hiEnd=53278795
	cd $(ANALYSIS) && locuszoom --metal $(ANALYSIS)/gwas.txt --refsnp 12:26472562 --flank 300kb --build hg19 --pop EUR --source 1000G_Nov2014 theme=publication height=6 hiStart=26436640 hiEnd=26491955

snp-overlap-with-encode-cres:
	cat $(WORK)/encode-cres/encode-cres.bed | cut -f1-4 | uniq > $@.bed
	bedtools closest -k 2 -d -a arl15-snp.bed -b $@.bed
	bedtools closest -k 2 -d -a itpr2-snp.bed -b $@.bed

compare-to-liger:
	mkdir -p $(ANALYSIS)/old-clustering
	mkdir -p $(ANALYSIS)/new-clustering-1
	cd $(ANALYSIS)/old-clustering && python $(BIN)/compare-to-old-clustering.py $(WORK)/downstream-new-features/results/liger/round-1/second-louvain-clusters.txt /lab/work/porchard/sn-muscle-project/work/process-by-cluster/clusters.txt $(ROOT)/library-labels.txt seurat-vs-old-liger 'LIGER cluster' 'Seurat cluster'
	cat /lab/work/porchard/sn-muscle-project/work/liger-human-and-rat-with-lambda-10/results/post-factorization-3/KSM2_rat_RNA/10/clusters-10-0.04.txt | perl -pe 's/hg19[A-Z]*/hg19/; s/rn6[A-Z]*/rn6/; s/ /\t/g' > $(ANALYSIS)/clusters-no-dual-modality.txt; python $(BIN)/add-dual-modality-ATAC-to-clusters.py $(ANALYSIS)/clusters-no-dual-modality.txt > $(ANALYSIS)/clusters-new-liger.txt
	cd $(ANALYSIS)/new-clustering-1 && python $(BIN)/compare-to-old-clustering.py $(ANALYSIS)/clusters-new-liger.txt /lab/work/porchard/sn-muscle-project/work/process-by-cluster/clusters.txt $(ROOT)/library-labels.txt seurat-vs-new-liger 'LIGER cluster' 'Seurat cluster'

tables:
	mkdir -p $(ANALYSIS)/tables
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results $(ROOT)/make-tables.nf &

knit: ANALYSIS=$(ROOT)/knit
knit: all
	# FANS vs no FANS figures
	ln -sf $(WORK)/manuscript-figures-v2/results/fans-vs-no-fans/fans-chromhmm-overlap.pdf $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures-v2/results/fans-vs-no-fans/umis-vs-mitochondrial-fans-vs-no-fans.png $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures-v2/results/fans-vs-no-fans/fans-rna-correlation.png $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures-v2/results/fans-vs-no-fans/fans-gb.pdf $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures-v2/results/fans/fans-atac-correlation.png $(ANALYSIS)/
	# 20k vs 40k figures
	ln -sf $(WORK)/tables/results/figures/hqaa-vs-tss-enrich-loading.png $(ANALYSIS)/hqaa-vs-tss-enrichment-20k-vs-40k.png
	ln -sf $(WORK)/tables/results/figures/hqaa-vs-max-frac-autosomal-loading.png $(ANALYSIS)/hqaa-vs-max-fraction-reads-from-single-autosome-20k-vs-40k.png
	ln -sf $(WORK)/manuscript-figures-v2/results/20k-vs-40k/loading-chromhmm-overlap.pdf $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures-v2/results/20k-vs-40k/loading-gb.pdf $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures-v2/results/20k-vs-40k/umis-vs-mitochondrial-20k-vs-40k.png $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures-v2/results/20k-vs-40k/loading-rna-correlation.png $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures-v2/results/loading/loading-atac-correlation.png $(ANALYSIS)/
	# dual modality QC figures
	ln -sf $(WORK)/qc-plots-dual-modality/post-species-calling-qc.png $(ANALYSIS)/dual-modality-qc.png
	ln -sf $(WORK)/qc-plots-dual-modality/human-vs-rat-atac-hqaa.png $(ANALYSIS)/dual-modality-human-vs-rat.png
	ln -sf $(WORK)/determine-species-dual-modality/results/plot/*.png $(ANALYSIS)/
	# Libraries used downstream
	ln -sf $(WORK)/tables/results/figures/hqaa-vs-tss-enrich-all-used-downstream.png $(ANALYSIS)/hqaa-vs-tss-enrichment-used-downstream.png
	ln -sf $(WORK)/tables/results/figures/hqaa-vs-max-frac-autosomal-all-used-downstream.png $(ANALYSIS)/hqaa-vs-max-fraction-reads-from-single-autosome-used-downstream.png
	ln -sf $(WORK)/manuscript-figures-v2/results/qc-for-downstream-libraries/umis-vs-mitochondrial-used-downstream.png $(ANALYSIS)/
	# clustering comparison
	ln -sf $(WORK)/compare-to-liger/new-clustering-1/seurat-vs-new-liger.alluvial-new-to-old.png $(ANALYSIS)/seurat-vs-liger-alluvial.png
	ln -sf $(WORK)/compare-to-liger/new-clustering-1/seurat-vs-new-liger.compare-to-old.png $(ANALYSIS)/seurat-vs-liger-heatmap.png
	# Figures for the biology part of the manuscript
	ln -sf $(WORK)/process-by-cluster/results/gene-expression-near-snps/chr12_26472562.png $(ANALYSIS)/gene-expression-near-ITPR2-locus.png
	ln -sf $(WORK)/process-by-cluster/results/gene-expression-near-snps/chr5_53271420.png $(ANALYSIS)/gene-expression-near-ARL15-locus.png
	ln -sf $(WORK)/process-by-cluster/results/figures/rubenstein-vs-our-fiber-type-lfcs.pdf $(ANALYSIS)/
	ln -sf $(WORK)/process-by-cluster/results/figures/fiber-types-per-species-MYH7-vs-MYH1_2_4.png $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures-v2/results/compare-deltasvm-to-allelic-bias/dsvm-vs-allelic-bias-type-2-fibers-fdr-1-deltasvmZ-2-min-coverage-15.png $(ANALYSIS)/dsvm-vs-allelic-bias.png
	ln -sf $(WORK)/manuscript-figures-v2/results/barcode-rank-plots/*.png $(ANALYSIS)/
	ln -sf $(WORK)/chromvar/results/plot/*.png $(ANALYSIS)/
	ln -sf $(WORK)/hic-overlap-gencode/results/figures/*.png $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures-v2/results/cicero-no-binarize-175Mb/chr5_53271420.png $(ANALYSIS)/ARL15-cicero.png
	ln -sf $(WORK)/manuscript-figures-v2/results/cicero-no-binarize-175Mb/chr12_26472562.png $(ANALYSIS)/ITPR2-cicero.png
	ln -sf $(WORK)/manuscript-figures-v2/results/seurat-links/spearman/chr12_26472562-peak-gene-links-heatmap.png $(ANALYSIS)/ITPR2-signac.png
	ln -sf $(WORK)/manuscript-figures-v2/results/seurat-links/spearman/chr5_53271420-peak-gene-links-heatmap.png $(ANALYSIS)/ARL15-signac.png
	$(foreach m,MEF2A MEF2B MEF2C MEF2D FOS JUN,ln -sf $(WORK)/process-by-cluster/results/figures/$(m).png $(ANALYSIS)/$(NL))
	ln -sf $(WORK)/luciferase/luciferase-all.png $(ANALYSIS)/
	ln -sf $(WORK)/diamante-deltaSVM/results/plot/5_53271420.png $(ANALYSIS)/
	ln -sf $(WORK)/diamante-deltaSVM/results/plot/12_26453283.png $(ANALYSIS)/
	ln -sf $(WORK)/ldsc/UKB-hesc/hg19/joint-model/results/heatmaps-and-tables/LDSC-UKB-across-all-cell-types.pdf $(ANALYSIS)/
	ln -sf $(WORK)/ldsc/UKB-hesc/rn6/joint-model/results/heatmaps-and-tables/LDSC-UKB-across-all-cell-types.pdf $(ANALYSIS)/UKB-LDSC-by-rn6.pdf
	# BMI unadjusted
	ln -sf $(WORK)/manuscript-figures-v2/results/T2D-GWAS-enrichments-no-bmiadj/rn6/T2D-FIns.pdf $(ANALYSIS)/T2D-FIns-rn6.pdf
	ln -sf $(WORK)/manuscript-figures-v2/results/gb-python/arl15-human-big.pdf $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures-v2/results/gb-python/itpr2-human-big.pdf $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures-v2/results/chromatin-state-legend/chrom_state_legend.pdf $(ANALYSIS)/
	cd $(ANALYSIS) && pdflatex main.tex && pdflatex main.tex

diamante-credible-sets: ANALYSIS=$(DATA)/$@
diamante-credible-sets:
	mkdir -p $(ANALYSIS)
	cp -r /lab/data/gwas/2018_02_DIAMANTE/credible_sets/genetic_credible_sets $(ANALYSIS)/

hiC:
	mkdir -p $(DATA)/$@
	cd $(DATA)/$@ && wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE87nnn/GSE87112/suppl/GSE87112_file.tar.gz
	cd $(DATA)/$@ && tar -xvzf GSE87112_file.tar.gz && tar -xvzf FitHiC_primary_cohort.tgz
	Rscript $(SRC)/format-hiC-data.R $(DATA)/$@/psoas-hic-connections.txt $(HG19_CHROM_SIZES) $(DATA)/$@/FitHiC_primary_cohort/FitHiC_output.PO_*

snatacseq-as-bulk:
	mkdir -p $(ANALYSIS)
	python $(CONTROL)/$@/make_config.py $(ROOT) $(ROOT)/sample_info/sample_info.with_facs.txt > $(ANALYSIS)/config.json
	cd $(ANALYSIS) && nohup nextflow run -resume -with-singularity /lab/work/porchard/singularity/archive/ATAC/2019-11-20/ATAC.simg -with-dag -with-timeline -with-trace -with-report -params-file $(ANALYSIS)/config.json --results $(ANALYSIS)/results /home/porchard/github/ATACseq-NextFlow/main.nf &

dual-modality-quality-nuclei:
	mkdir -p $(ANALYSIS)
	cat /lab/work/porchard/daily/2020-11-18/labs/human-nuclei.tsv | grep -v rna_barcode | cut -f2 | perl -pe 's/^/1846_RNA-hg19\t/; s/$$/\tKSM2/' > $(ANALYSIS)/nuclei-with-individual-assignments.txt
	cat /lab/work/porchard/daily/2020-11-18/labs/human-nuclei.tsv | grep -v rna_barcode | cut -f3 | perl -pe 's/^/1846_ATAC-hg19\t/; s/$$/\tKSM2/' >> $(ANALYSIS)/nuclei-with-individual-assignments.txt
	cat /lab/work/porchard/daily/2020-11-18/labs/rat-nuclei.tsv | grep -v rna_barcode | cut -f2 | perl -pe 's/^/1846_RNA-rn6\t/; s/$$/\trat1/'  >> $(ANALYSIS)/nuclei-with-individual-assignments.txt
	cat /lab/work/porchard/daily/2020-11-18/labs/rat-nuclei.tsv | grep -v rna_barcode | cut -f3 | perl -pe 's/^/1846_ATAC-rn6\t/; s/$$/\trat1/' >> $(ANALYSIS)/nuclei-with-individual-assignments.txt
	cat $(ANALYSIS)/nuclei-with-individual-assignments.txt | cut -f1-2 > $(ANALYSIS)/nuclei.txt

quality-nuclei-count-files:
	mkdir -p $(ANALYSIS)/data/atac
	mkdir -p $(ANALYSIS)/data/rna
	cat $(WORK)/qc/results/call-nuclei-atac/pass-qc.tsv $(WORK)/qc/results/call-nuclei-rna/pass-qc.tsv $(WORK)/dual-modality-quality-nuclei/nuclei-with-individual-assignments.txt | sort | uniq > $(ANALYSIS)/data/nuclei-with-individual-assignments.txt
	cut -f1-2 $(ANALYSIS)/data/nuclei-with-individual-assignments.txt > $(ANALYSIS)/data/nuclei.txt
	ln -s $(WORK)/fragment-counts-atac/results/uncorrected-counts/*.counts.txt $(ANALYSIS)/data/atac/
	ln -s $(WORK)/rnaseq/results/starsolo/* $(ANALYSIS)/data/rna
	ln -s $(WORK)/rnaseq-dual-modality/results/starsolo/* $(ANALYSIS)/data/rna
	# remove FANS atac libraries
	rm $(ANALYSIS)/data/atac/133152-hg19*
	rm $(ANALYSIS)/data/atac/133154-hg19*
	cd $(ANALYSIS) && nohup nextflow run -resume --starsolo_dir_glob '$(ANALYSIS)/data/rna/*' --results $(ANALYSIS)/results --nuclei $(ANALYSIS)/data/nuclei.txt --atac_gene_count_glob '$(ANALYSIS)/data/atac/*' --individual_assignments $(ANALYSIS)/data/nuclei-with-individual-assignments.txt $(ROOT)/make-quality-nuclei-count-files.nf &

demuxlet-mask:
	mkdir -p $(ANALYSIS)
	python $(BIN)/make-demuxlet-mask.py --mask-top 0.3 $(WORK)/rnaseq/results/starsolo/63_40_rna-hg19/Solo.out/GeneFull/raw/features.tsv $(WORK)/rnaseq/results/starsolo/63_40_rna-hg19/Solo.out/GeneFull/raw/matrix.mtx /lab/data/reference/human/hg19/index/STAR/current/gencode.v19.annotation.gtf > $(ANALYSIS)/demuxlet-mask.bed

wasp-input-bams:
	mkdir -p $(ANALYSIS)/data
	cp $(WORK)/process-by-cluster/clusters.txt $(ANALYSIS)/data/clusters.txt
	cat $(WORK)/quality-nuclei-count-files/data/nuclei-with-individual-assignments.txt | grep -v -w 133152-hg19 | grep -v -w 133154-hg19 > $(ANALYSIS)/data/nuclei-with-individual-assignments.txt
	ln -s $(WORK)/atacseq/results/merge/*.bam $(ANALYSIS)/data/
	ln -s $(WORK)/atacseq-dual-modality/results/merge/*.bam $(ANALYSIS)/data/
	# delete the FANS atac libraries
	rm -rf $(ANALYSIS)/data/133152-hg19*
	rm -rf $(ANALYSIS)/data/133154-hg19*
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --bam_glob '$(ANALYSIS)/data/*.bam' --clusters $(ANALYSIS)/data/clusters.txt --nuclei_with_individual_assignments $(ANALYSIS)/data/nuclei-with-individual-assignments.txt $(ROOT)/make-wasp-input-bams.nf &

atac-allelic-bias-by-cell-type:
	mkdir -p $(ANALYSIS)/data
	ln -s $(WORK)/wasp-input-bams/results/split-library-by-individual-and-cluster/*hg19*.bam $(ANALYSIS)/data/
	cd $(ANALYSIS) && nohup nextflow run -resume --vcf $(DATA)/genotypes/genotypes.noblacklist.vcf.gz --results $(ANALYSIS)/results --bam_glob '$(ANALYSIS)/data/*.bam' $(ROOT)/atac-allelic-bias-by-cell-type.nf &

seurat-joint-2:
	mkdir -p $(ANALYSIS)/data
	ln -s $(WORK)/quality-nuclei-count-files/results/liger-in/* $(ANALYSIS)/data/
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --liger_in '$(ANALYSIS)/data/*' $(ROOT)/seurat-2.nf &

liger-human-and-rat-with-lambda-10:
	mkdir -p $(ANALYSIS)/data
	ln -s $(WORK)/quality-nuclei-count-files/results/liger-in/* $(ANALYSIS)/data/
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --liger_in '$(ANALYSIS)/data/*' --lambda 10 $(ROOT)/liger-human-and-rat-with-lambda.nf &

process-by-cluster:
	mkdir -p $(ANALYSIS)/data
	ln -s $(WORK)/atacseq/results/prune/*.bam $(ANALYSIS)/data/
	ln -s $(WORK)/atacseq-dual-modality/results/prune/*.bam $(ANALYSIS)/data/
	# remove FANS atac libraries
	rm $(ANALYSIS)/data/133152-hg19*
	rm $(ANALYSIS)/data/133154-hg19*
	python $(CONTROL)/process-by-cluster/make_config.py $(ANALYSIS)/data/*.bam > $(CONFIG)
	cp /lab/work/porchard/sn-muscle-project/work/seurat-joint-2/results/seurat-out/seurat-clusters.txt $(ANALYSIS)/clusters-no-dual-modality.txt
	python $(BIN)/add-dual-modality-ATAC-to-clusters.py $(ANALYSIS)/clusters-no-dual-modality.txt > $(ANALYSIS)/clusters.txt
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results -params-file $(CONFIG) --library_labels $(ROOT)/library-labels.txt --liger_in_matrices_glob '$(WORK)/liger-human-and-rat-with-lambda-5/data/*' --cluster_names $(ROOT)/cluster-names.txt --clusters $(ANALYSIS)/clusters.txt $(ROOT)/process-by-cluster.nf &

### GWAS ENRICHMENTS ###
ldsc-T2D-one-model-per-cluster-hg19: ANALYSIS=$(WORK)/ldsc/T2D/hg19/one-model-per-cluster
ldsc-T2D-one-model-per-cluster-hg19:
	mkdir -p $(ANALYSIS)/data
	$(foreach c,$(shell seq 0 6),cat $(WORK)/process-by-cluster/results/peaks/broad/$(c)-hg19_peaks.broadPeak.noblacklist | sort -k1,1 -k2n,2 | cut -f1-3 > $(ANALYSIS)/data/cluster_$(c).bed$(NL))
	cp $(DATA)/other-annotations/beta_ATAC.bed $(ANALYSIS)/data/beta_cell.bed
	cp $(DATA)/other-annotations/adipose.bed $(ANALYSIS)/data/adipose.bed
	cp $(WORK)/meuleman-with-hesc/results/tissues/Liver.bed $(ANALYSIS)/data/liver.bed
	cd $(ANALYSIS) && nohup nextflow run -resume -with-trace -with-report --other_annotations '$(SHARED_OPEN_CHROMATIN)' --projroot $(ROOT) --results $(ANALYSIS)/results --cluster_names $(CLUSTER_NAMES) $(ROOT)/ldsc-T2D-one-model-per-cell-type.nf --peak_glob '$(ANALYSIS)/data/*.bed' &

ldsc-T2D-one-model-per-cluster-rn6: ANALYSIS=$(WORK)/ldsc/T2D/rn6/one-model-per-cluster
ldsc-T2D-one-model-per-cluster-rn6:
	mkdir -p $(ANALYSIS)/data
	$(foreach c,$(shell seq 0 6),cat $(WORK)/process-by-cluster/results/rat-peak-liftover/$(c)-rn6_peaks_in_hg19.bed | sort -k1,1 -k2n,2 | cut -f1-3 > $(ANALYSIS)/data/cluster_$(c).bed$(NL))
	cp $(DATA)/other-annotations/beta_ATAC.bed $(ANALYSIS)/data/beta_cell.bed
	cp $(DATA)/other-annotations/adipose.bed $(ANALYSIS)/data/adipose.bed
	cp $(WORK)/meuleman-with-hesc/results/tissues/Liver.bed $(ANALYSIS)/data/liver.bed
	cd $(ANALYSIS) && nohup nextflow run -resume -with-trace -with-report --other_annotations '$(SHARED_OPEN_CHROMATIN)' --projroot $(ROOT) --results $(ANALYSIS)/results --cluster_names $(CLUSTER_NAMES) $(ROOT)/ldsc-T2D-one-model-per-cell-type.nf --peak_glob '$(ANALYSIS)/data/*.bed' &

ldsc-T2D-joint-model-hg19: ANALYSIS=$(WORK)/ldsc/T2D/hg19/joint-model
ldsc-T2D-joint-model-hg19:
	mkdir -p $(ANALYSIS)/data
	$(foreach c,$(shell seq 0 6),cat $(WORK)/process-by-cluster/results/peaks/broad/$(c)-hg19_peaks.broadPeak.noblacklist | sort -k1,1 -k2n,2 | cut -f1-3 > $(ANALYSIS)/data/cluster_$(c).bed$(NL))
	cp $(DATA)/other-annotations/beta_ATAC.bed $(ANALYSIS)/data/beta_cell.bed
	cp $(DATA)/other-annotations/adipose.bed $(ANALYSIS)/data/adipose.bed
	cp $(WORK)/meuleman-with-hesc/results/tissues/Liver.bed $(ANALYSIS)/data/liver.bed
	cd $(ANALYSIS) && nohup nextflow run -resume -with-trace -with-report --other_annotations '$(SHARED_OPEN_CHROMATIN)' --projroot $(ROOT) --results $(ANALYSIS)/results --cluster_names $(CLUSTER_NAMES) $(ROOT)/ldsc-T2D-joint-model.nf --peak_glob '$(ANALYSIS)/data/*.bed' &

ldsc-T2D-joint-model-rn6: ANALYSIS=$(WORK)/ldsc/T2D/rn6/joint-model
ldsc-T2D-joint-model-rn6:
	mkdir -p $(ANALYSIS)/data
	$(foreach c,$(shell seq 0 6),cat $(WORK)/process-by-cluster/results/rat-peak-liftover/$(c)-rn6_peaks_in_hg19.bed | sort -k1,1 -k2n,2 | cut -f1-3 > $(ANALYSIS)/data/cluster_$(c).bed$(NL))
	cp $(DATA)/other-annotations/beta_ATAC.bed $(ANALYSIS)/data/beta_cell.bed
	cp $(DATA)/other-annotations/adipose.bed $(ANALYSIS)/data/adipose.bed
	cp $(WORK)/meuleman-with-hesc/results/tissues/Liver.bed $(ANALYSIS)/data/liver.bed
	cd $(ANALYSIS) && nohup nextflow run -resume -with-trace -with-report --other_annotations '$(SHARED_OPEN_CHROMATIN)' --projroot $(ROOT) --results $(ANALYSIS)/results --cluster_names $(CLUSTER_NAMES) $(ROOT)/ldsc-T2D-joint-model.nf --peak_glob '$(ANALYSIS)/data/*.bed' &

ldsc-ukb-hg19-hesc: OTHER_ANNOTATIONS = $(shell ls $(WORK)/meuleman-with-hesc/results/tissues/*.bed) $(shell ls $(DATA)/other-annotations/*.bed)
ldsc-ukb-hg19-hesc: ANALYSIS=$(WORK)/ldsc/UKB-hesc/hg19/joint-model
ldsc-ukb-hg19-hesc:
	mkdir -p $(ANALYSIS)/data
	$(foreach f,$(OTHER_ANNOTATIONS),#cp $(f) $(ANALYSIS)/data/$(NL))
	$(foreach c,$(shell seq 0 6),#cp $(WORK)/process-by-cluster/results/peaks/broad/$(c)-hg19_peaks.broadPeak.noblacklist $(ANALYSIS)/data/cluster_$(c).bed$(NL))
	#cd $(ANALYSIS)/data && rm -rf Amion.bed Esophagus.bed Periodontal_Ligament.bed
	cd $(ANALYSIS) && nohup nextflow run -resume -with-trace -with-report --projroot $(ROOT) --results $(ANALYSIS)/results $(ROOT)/ldsc-UKB-new-baseline.nf --other_annotations '$(SHARED_OPEN_CHROMATIN)' --trait_list $(ROOT)/ukb-traits.txt --peak_glob '$(ANALYSIS)/data/*.bed' &

ldsc-ukb-rn6-hesc: OTHER_ANNOTATIONS = $(shell ls $(WORK)/meuleman-with-hesc/results/tissues/*.bed) $(shell ls $(DATA)/other-annotations/*.bed)
ldsc-ukb-rn6-hesc: ANALYSIS=$(WORK)/ldsc/UKB-hesc/rn6/joint-model
ldsc-ukb-rn6-hesc:
	mkdir -p $(ANALYSIS)/data
	$(foreach f,$(OTHER_ANNOTATIONS),#cp $(f) $(ANALYSIS)/data/$(NL))
	$(foreach c,$(shell seq 0 6),#cp $(WORK)/process-by-cluster/results/rat-peak-liftover/$(c)-rn6_peaks_in_hg19.bed $(ANALYSIS)/data/cluster_$(c).bed$(NL))
	#cd $(ANALYSIS)/data && rm -rf Amion.bed Esophagus.bed Periodontal_Ligament.bed
	cd $(ANALYSIS) && nohup nextflow run -resume -with-trace -with-report --projroot $(ROOT) --results $(ANALYSIS)/results $(ROOT)/ldsc-UKB-new-baseline.nf --other_annotations '$(SHARED_OPEN_CHROMATIN)' --trait_list $(ROOT)/ukb-traits.txt --peak_glob '$(ANALYSIS)/data/*.bed' &

find-peak-gene-links-across-clusters-2:
	mkdir -p $(ANALYSIS)
	cat $(WORK)/process-by-cluster/results/master-peak-counts/1846_ATAC-hg19.*.counts.txt > $(ANALYSIS)/counts.txt
	cat $(ANALYSIS)/counts.txt | cut -f3 | perl -pe 's/:/\t/g' > $(ANALYSIS)/peaks.txt
	cat $(ANALYSIS)/counts.txt | cut -f2,4 > $(ANALYSIS)/barcodes.txt
	paste $(ANALYSIS)/peaks.txt $(ANALYSIS)/barcodes.txt > $(ANALYSIS)/counts.txt
	cut -f4 $(ANALYSIS)/counts.txt | uniq | sort | uniq > $(ANALYSIS)/include-barcodes.txt
	python $(BIN)/atac-peak-counts-to-hdf5.py $(ANALYSIS)/counts.txt $(ANALYSIS)/include-barcodes.txt $(ANALYSIS)/atac_counts.hdf5
	cp $(WORK)/quality-nuclei-count-files/results/rna-corrected/1846_RNA-hg19.hdf5 $(ANALYSIS)/rna_counts.hdf5
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --clusters $(WORK)/process-by-cluster/clusters.txt --atac $(ANALYSIS)/atac_counts.hdf5 --rna $(ANALYSIS)/rna_counts.hdf5 $(ROOT)/signac-link-peaks-to-genes-across-clusters-flexible.nf &

manuscript-figures-v2:
	mkdir -p $(ANALYSIS)/data/other-tissues
	mkdir -p $(ANALYSIS)/data/human-muscle
	mkdir -p $(ANALYSIS)/data/rat-muscle
	ln -sf $(WORK)/browser-sessions/2020-sn-muscle-GWAS-overlap/results/bigwigs/* $(ANALYSIS)/data/other-tissues/
	ln -sf $(WORK)/process-by-cluster/results/bigwigs/atac/*hg19* $(ANALYSIS)/data/human-muscle/
	ln -sf $(WORK)/process-by-cluster/results/bigwigs/atac/*rn6* $(ANALYSIS)/data/rat-muscle/
	rm $(ANALYSIS)/data/other-tissues/liver*
	rm $(ANALYSIS)/data/human-muscle/all*
	rm $(ANALYSIS)/data/rat-muscle/all*
	cd $(ANALYSIS) && nohup nextflow run --rat_muscle_bw_glob '$(ANALYSIS)/data/rat-muscle/*' --human_muscle_bw_glob '$(ANALYSIS)/data/human-muscle/*' --other_bw_glob '$(ANALYSIS)/data/other-tissues/*' -resume --results $(ANALYSIS)/results --projroot $(ROOT) $(ROOT)/manuscript-figures-v2.nf &

delta-svm: NUMBER_PEAKS_KEEP=40000
delta-svm: FASTA=$(HG19_FASTA)
delta-svm: CHROM_SIZES=$(HG19_CHROM_SIZES)
delta-svm: GENOME=hg19
delta-svm: BLACKLISTS=$(HG19_BLACKLISTS)
delta-svm: SNP_FILE=/lab/work/porchard/sn-muscle-project/work/1000G-SNP-vcf/1000G-snps.vcf.gz
delta-svm: ANALYSIS=$(WORK)/delta-svm/all-keep-40k
delta-svm:
	mkdir -p $(ANALYSIS)/data
	$(foreach c,$(shell seq 0 6),cat $(WORK)/process-by-cluster/results/peaks/narrow/$(c)-$(GENOME)_peaks.narrowPeak.noblacklist | sort -k9n,9 | tail -n $(NUMBER_PEAKS_KEEP) > $(ANALYSIS)/data/cluster_$(c).narrow.bed$(NL))
	$(foreach c,$(shell seq 0 6),cp $(WORK)/process-by-cluster/results/peaks/broad/$(c)-$(GENOME)_peaks.broadPeak.noblacklist $(ANALYSIS)/data/cluster_$(c).broad.bed$(NL))
	cd $(ANALYSIS) && cp $(ROOT)/deltasvm-config.nf nextflow.config && NXF_VER=19.10.0 nohup ~/github/snATAC_deltaSVM/main.nf -resume --repeat_bed $(DATA)/repeats/trf.bed --max_repeat_content 0.6 --results_dir $(ANALYSIS)/results --kernel 2 --peak_dir '$(ANALYSIS)/data/cluster_*.narrow.bed' --broadpeak_dir '$(ANALYSIS)/data/cluster_*.broad.bed' --exclude '$(BLACKLISTS)' --ref $(FASTA) --human_ref $(HG19_FASTA) --snp_file $(SNP_FILE) --chrom_sizes $(CHROM_SIZES) --workDir $(ANALYSIS)/work &

gkmexplain:
	mkdir -p $(ANALYSIS)/data
	$(foreach c,$(shell seq 0 6),cp $(WORK)/delta-svm/all-keep-40k/results/cluster_$(c)/model.model.txt $(ANALYSIS)/data/cluster_$(c).model.txt$(NL))
	cd $(ANALYSIS) && nohup nextflow run -resume --meme_glob '$(DATA)/motifs/*.meme' --fimo_background $(DATA)/motifs/hg19.background.txt --plain_motif_glob '$(DATA)/motifs/*.txt' --snp_file $(ROOT)/snps-for-manuscript.txt --results $(ANALYSIS)/results --model_glob '$(ANALYSIS)/data/*' $(ROOT)/gkmexplain-and-fimo-new-plots.nf &

find-disrupted-motifs:
	mkdir -p $(ANALYSIS)/data/motifs
	cp /lab/work/porchard/fimo/work/ENCODE2013/data/*.meme $(ANALYSIS)/data/motifs/
	rm $(ANALYSIS)/data/motifs/*_disc*
	cd $(ANALYSIS) && nohup nextflow run -resume --meme_glob '$(ANALYSIS)/data/motifs/*.meme' --fimo_background $(DATA)/motifs/hg19.background.txt --snp_file $(ROOT)/snps-for-manuscript.txt --results $(ANALYSIS)/results $(ROOT)/find-disrupted-motifs.nf &

get-deltasvm-z-scores:
	mkdir -p $(ANALYSIS)
	cat snps-for-manuscript.txt | perl -pe 's/ /\t/g' | cut -f1,2 | grep -v BP | perl -pe 's/^/chr/' > $(ANALYSIS)/pos-file.txt
	$(foreach c,$(shell seq 0 6),python $(ROOT)/bin/get-deltaSVM-z-scores.py --deltaSVM-file $(WORK)/delta-svm/all-keep-40k/results/cluster_$(c)/snp_scores_cluster_$(c).txt --pos-file $(ANALYSIS)/pos-file.txt > $(ANALYSIS)/cluster_$(c).z-scores.txt$(NL))


logistic-regression-hg19: ANALYSIS=$(WORK)/logistic-regression/hg19
logistic-regression-hg19:
	mkdir -p $(ANALYSIS)/data
	cp $(WORK)/process-by-cluster/results/peaks/broad/*hg19*_peaks.broadPeak.noblacklist $(ANALYSIS)/data
	rm $(ANALYSIS)/data/all*
	cd $(ANALYSIS) && nohup nextflow run -resume -with-trace -with-report --results $(ANALYSIS)/results --broadpeak_glob '$(ANALYSIS)/data/*' --genome hg19 --projroot $(ROOT) $(ROOT)/logistic-regression.nf &

logistic-regression-rn6: ANALYSIS=$(WORK)/logistic-regression/rn6
logistic-regression-rn6:
	mkdir -p $(ANALYSIS)/data
	cp $(WORK)/process-by-cluster/results/peaks/broad/*rn6*_peaks.broadPeak.noblacklist $(ANALYSIS)/data
	rm $(ANALYSIS)/data/all*
	cd $(ANALYSIS) && nohup nextflow run -resume -with-trace -with-report --results $(ANALYSIS)/results --broadpeak_glob '$(ANALYSIS)/data/*' --genome rn6 --projroot $(ROOT) $(ROOT)/logistic-regression.nf &

chromvar:
	mkdir -p $(ANALYSIS)/data/motifs
	mkdir -p $(ANALYSIS)/data/peaks
	mkdir -p $(ANALYSIS)/data/bams
	ln -s /lab/work/porchard/fimo/work/cisBP/* $(ANALYSIS)/data/motifs/
	ln -s $(WORK)/process-by-cluster/results/aggregate-of-quality-nuclei/bam/*hg19.bam $(ANALYSIS)/data/bams/
	cp $(WORK)/process-by-cluster/results/peaks/broad/*hg19_peaks.broadPeak.noblacklist $(ANALYSIS)/data/peaks/
	rm $(ANALYSIS)/data/peaks/all*
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --motif_to_tf $(DATA)/cisbp-motifs/motif-id-to-tf-name.txt --cluster_names $(ROOT)/cluster-names.txt --clusters $(WORK)/process-by-cluster/clusters.txt --fimo_glob '$(ANALYSIS)/data/motifs/*' --bam_glob '$(ANALYSIS)/data/bams/*' --peak_glob '$(ANALYSIS)/data/peaks/*' $(ROOT)/$@.nf &

luciferase:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && python $(BIN)/luciferase.py > out.txt

cicero-no-binarize:
	mkdir -p $(ANALYSIS)/data
	ln -s $(WORK)/process-by-cluster/results/peak-counts/* $(ANALYSIS)/data/
	zcat /home/porchard/github/ataqv/data/tss/hg19.tss.refseq.bed.gz | grep -v 'hap' > $(ANALYSIS)/hg19.tss.bed
	python $(BIN)/add-dual-modality-ATAC-to-umap.py $(WORK)/seurat-joint-2/results/seurat-out/seurat-umap.txt | perl -pe 's/\t/-/' > $(ANALYSIS)/data/umap.txt
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --tss $(ANALYSIS)/hg19.tss.bed --counts_glob '$(ANALYSIS)/data/*.counts.txt' --snps $(ROOT)/cicero-snps.bed --umap $(ANALYSIS)/data/umap.txt $(ROOT)/cicero-no-binarize.nf &

diamante-deltaSVM:
	mkdir -p $(ANALYSIS)/data
	cp /lab/work/porchard/sn-muscle-project/work/process-by-cluster/results/peaks/broad/*hg19_peaks.broadPeak.noblacklist $(ANALYSIS)/data
	rm $(ANALYSIS)/data/all*
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --peak_file_glob '$(ANALYSIS)/data/*' --cs_file_glob '$(DATA)/diamante-credible-sets/genetic_credible_sets/*' --dsvm_file_glob '$(WORK)/delta-svm/all-keep-40k/results/cluster_*/snp_scores_cluster_*.txt' $(ROOT)/delta-svm-diamante-loci.nf &

hic-overlap-gencode:
	mkdir -p $(ANALYSIS)
	bedtools intersect -wa -b arl15-snp.bed -a work/process-by-cluster/results/process-by-cluster-round-1/master-peaks/master-peaks.hg19.bed | perl -pe 's/$$/\tARL15_peak/'> $(ANALYSIS)/peaks.bed
	bedtools intersect -wa -b itpr2-snp.bed -a work/process-by-cluster/results/process-by-cluster-round-1/master-peaks/master-peaks.hg19.bed | perl -pe 's/$$/\tITPR2_peak/' >> $(ANALYSIS)/peaks.bed
	cp /lab/work/porchard/sn-muscle-project/data/gencode-tss-for-connections/hg19.bed $(ANALYSIS)/hg19.tss.bed
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --pchic_glob '$(DATA)/pc-hiC/po-*.xlsx' --yue_glob '$(DATA)/hic-yue-lab/hg19/*' --peaks $(ANALYSIS)/peaks.bed --fithic_glob '$(DATA)/hiC/FitHiC_primary_cohort/*' --bins $(DATA)/hiC/bins.bed --snp_bed $(ROOT)/cicero-snps.bed --tss $(ANALYSIS)/hg19.tss.bed $(ROOT)/hic.nf &

qc-plots-dual-modality:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && python $(BIN)/sn-dual-modality-qc.py

plot-rna-and-atac-rat-human-mixing:
	mkdir -p $(ANALYSIS)
	grep 1846_ATAC-rn6 $(WORK)/process-by-cluster/clusters.txt | shuf | head -n 5 | cut -f2 > $(ANALYSIS)/barcodes.txt
	grep 1846_ATAC-hg19 $(WORK)/process-by-cluster/clusters.txt | shuf | head -n 5 | cut -f2 >> $(ANALYSIS)/barcodes.txt
	python $(BIN)/filter-bam-by-barcode.py $(WORK)/atacseq-dual-modality/results/mark_duplicates/1846_ATAC-hg19.md.bam $(ANALYSIS)/hg19-atac.bam $(ANALYSIS)/barcodes.txt
	python $(BIN)/filter-bam-by-barcode.py $(WORK)/atacseq-dual-modality/results/mark_duplicates/1846_ATAC-rn6.md.bam $(ANALYSIS)/rn6-atac.bam $(ANALYSIS)/barcodes.txt
	samtools index $(ANALYSIS)/hg19-atac.bam
	samtools index $(ANALYSIS)/rn6-atac.bam
	cd $(ANALYSIS) && python $(BIN)/count-species-unique-mapping-atac.py ~/github/snATACseq-NextFlow/737K-arc-v1.txt $(WORK)/atacseq-dual-modality/results/mark_duplicates/1846_ATAC-hg19.md.bam $(WORK)/atacseq-dual-modality/results/mark_duplicates/1846_ATAC-rn6.md.bam hg19 rn6 > $(ANALYSIS)/1846_ATAC.txt
	cd $(ANALYSIS) && python $(BIN)/count-species-unique-mapping-atac.py ~/github/snATACseq-NextFlow/737K-arc-v1.txt $(ANALYSIS)/hg19-atac.bam $(ANALYSIS)/rn6-atac.bam hg19 rn6 > $(ANALYSIS)/1846_ATAC.txt

determine-species-dual-modality:
	mkdir -p $(ANALYSIS)/data
	ln -s $(WORK)/atacseq-dual-modality/results/mark_duplicates/1846_ATAC-hg19.md.bam $(ANALYSIS)/data/1846_ATAC-hg19.bam
	ln -s $(WORK)/atacseq-dual-modality/results/mark_duplicates/1846_ATAC-rn6.md.bam $(ANALYSIS)/data/1846_ATAC-rn6.bam
	ln -s $(WORK)/rnaseq-dual-modality/results/starsolo/1846_RNA-hg19/Aligned.sortedByCoord.out.bam $(ANALYSIS)/data/1846_RNA-hg19.bam
	ln -s $(WORK)/rnaseq-dual-modality/results/starsolo/1846_RNA-rn6/Aligned.sortedByCoord.out.bam $(ANALYSIS)/data/1846_RNA-rn6.bam
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --atac_glob '$(ANALYSIS)/data/*ATAC*' --rna_glob '$(ANALYSIS)/data/*RNA*' --rna_whitelist ~/github/snRNAseq-NextFlow/737K-arc-v1.txt --atac_whitelist ~/github/snATACseq-NextFlow/737K-arc-v1.txt $(ROOT)/$@.nf &

test-droplet-utils: SOLO_DIR=/lab/work/porchard/Nova-315/work/rnaseq-chimeric/results/starsolo/1846_RNA-hg19_rn6/Solo.out
test-droplet-utils: 
	mkdir -p $(ANALYSIS) 
	cat $(SOLO_DIR)/GeneFull/raw/features.tsv | perl -pe 's/$$/\tGene Expression/' > $(ANALYSIS)/features.tsv
	cp $(SOLO_DIR)/GeneFull/raw/matrix.mtx $(ANALYSIS)/matrix.mtx
	cp $(SOLO_DIR)/GeneFull/raw/barcodes.tsv $(ANALYSIS)/barcodes.tsv
	cd $(ANALYSIS) && gzip *.tsv *.mtx
	cd $(ANALYSIS) && Rscript $(BIN)/droplet-utils.R $(ANALYSIS)/ 1846_RNA

make-read-beds:
	mkdir -p $(ANALYSIS)/data/atac
	mkdir -p $(ANALYSIS)/data/rna
	ln -s /lab/work/porchard/sn-muscle-project/work/atacseq/results/prune/*pruned.bam $(ANALYSIS)/data/atac/
	ln -s /lab/work/porchard/sn-muscle-project/work/rnaseq/results/prune/*before-dedup.bam $(ANALYSIS)/data/rna/
	ln -s /lab/work/porchard/sn-muscle-project/work/atacseq-dual-modality/results/prune/*pruned.bam $(ANALYSIS)/data/atac/
	ln -s /lab/work/porchard/sn-muscle-project/work/rnaseq-dual-modality/results/prune/*before-dedup.bam $(ANALYSIS)/data/rna/
	$(foreach r,1 2 3 4,ln -s $(WORK)/bulk-atacseq/results/prune/320-NM-$(r).pruned.bam $(ANALYSIS)/data/atac/$(NL))
	cd $(ANALYSIS) && nohup nextflow run -resume --rna_bam_glob '$(ANALYSIS)/data/rna/*.bam' --atac_bam_glob '$(ANALYSIS)/data/atac/*.bam' --results $(ANALYSIS)/results $(ROOT)/$@.nf &

starch-bdg:
	mkdir -p $(ANALYSIS)/data
	cp $(WORK)/process-by-cluster/results/peaks/broad/*treat_pileup.bdg $(ANALYSIS)/data/
	rm -rf $(ANALYSIS)/data/all*
	cp $(WORK)/bulk-atacseq/results/macs2/*treat_pileup.bdg $(ANALYSIS)/data/
	cp $(WORK)/process-as-bulk/results/macs2/*treat_pileup.bdg $(ANALYSIS)/data/ 
	cd $(ANALYSIS) && nohup nextflow run -resume --bdg_glob '$(ANALYSIS)/data/*.bdg' --results $(ANALYSIS)/results $(ROOT)/$@.nf &

zenodo:
	mkdir -p $(ANALYSIS)/filtered-read-locations/atac/bulk
	mkdir -p $(ANALYSIS)/filtered-read-locations/atac/single-nucleus
	mkdir -p $(ANALYSIS)/filtered-read-locations/rna
	mkdir -p $(ANALYSIS)/atac-peaks/bulk
	mkdir -p $(ANALYSIS)/atac-peaks/single-nucleus
	mkdir -p $(ANALYSIS)/atac-bedgraphs/single-nucleus
	mkdir -p $(ANALYSIS)/atac-bedgraphs/bulk
	mkdir -p $(ANALYSIS)/gene-counts
	mkdir -p $(ANALYSIS)/clustering
	mkdir -p $(ANALYSIS)/luciferase
	cp $(DATA)/luciferase/ARL15-rs702634_Luciferase_Rani_12-21-20.xlsx $(ANALYSIS)/luciferase/luciferase.xlsx
	cp zenodo-sample-info.txt $(ANALYSIS)/sample-info.txt
	cp zenodo-readme.txt $(ANALYSIS)/README.txt
	# cluster assignments
	cp -L $(WORK)/process-by-cluster/clusters.txt $(ANALYSIS)/clustering/clusters.txt
	# umap
	cp -L $(WORK)/seurat-joint-2/results/seurat-out/seurat-umap.txt $(ANALYSIS)/clustering/umap.txt
	# bedgraphs (starched)
	# SN bedgraphs
	cp $(foreach g,hg19 rn6,$(foreach c,$(shell seq 0 6),$(WORK)/starch-bdg/results/starch/$(c)-$(g).starch)) $(ANALYSIS)/atac-bedgraphs/single-nucleus/
	# Bulk bedgraphs
	cp -L $(WORK)/starch-bdg/results/starch/133*.starch $(ANALYSIS)/atac-bedgraphs/bulk/
	cp -L $(WORK)/starch-bdg/results/starch/63_*.starch $(ANALYSIS)/atac-bedgraphs/bulk/
	cp -L $(WORK)/starch-bdg/results/starch/320-NM-1.starch $(ANALYSIS)/atac-bedgraphs/bulk/HSM1.bulk.rep1.starch
	cp -L $(WORK)/starch-bdg/results/starch/320-NM-2.starch $(ANALYSIS)/atac-bedgraphs/bulk/HSM1.bulk.rep2.starch
	cp -L $(WORK)/starch-bdg/results/starch/320-NM-3.starch $(ANALYSIS)/atac-bedgraphs/bulk/HSM2.bulk.rep1.starch
	cp -L $(WORK)/starch-bdg/results/starch/320-NM-4.starch $(ANALYSIS)/atac-bedgraphs/bulk/HSM2.bulk.rep2.starch
	# ATAC-seq peaks. Only use SN peaks.
	$(foreach c,$(shell seq 0 6),cp -rL $(WORK)/process-by-cluster/results/peaks/broad/$(c)-*.broadPeak.noblacklist $(ANALYSIS)/atac-peaks/single-nucleus/$(NL))
	cp -rL $(WORK)/process-as-bulk/results/macs2/*.broadPeak.noblacklist $(ANALYSIS)/atac-peaks/bulk/
	cp -L $(WORK)/bulk-atacseq/results/macs2/320-NM-1_peaks.broadPeak.noblacklist $(ANALYSIS)/atac-peaks/bulk/HSM1.bulk.rep1.broadPeak.noblacklist
	cp -L $(WORK)/bulk-atacseq/results/macs2/320-NM-2_peaks.broadPeak.noblacklist $(ANALYSIS)/atac-peaks/bulk/HSM1.bulk.rep2.broadPeak.noblacklist
	cp -L $(WORK)/bulk-atacseq/results/macs2/320-NM-3_peaks.broadPeak.noblacklist $(ANALYSIS)/atac-peaks/bulk/HSM2.bulk.rep1.broadPeak.noblacklist
	cp -L $(WORK)/bulk-atacseq/results/macs2/320-NM-4_peaks.broadPeak.noblacklist $(ANALYSIS)/atac-peaks/bulk/HSM2.bulk.rep2.broadPeak.noblacklist
	# Gene counts	
	$(foreach m,RNA ATAC,cp -L $(WORK)/quality-nuclei-count-files/results/liger-in/KSM1_$(m).hdf5 $(ANALYSIS)/gene-counts/$(m)-HSM1.counts.h5$(NL))
	$(foreach m,RNA ATAC,cp -L $(WORK)/quality-nuclei-count-files/results/liger-in/KSM2_$(m).hdf5 $(ANALYSIS)/gene-counts/$(m)-HSM2.counts.h5$(NL))
	$(foreach m,RNA ATAC,cp -L $(WORK)/quality-nuclei-count-files/results/liger-in/rat1_$(m).hdf5 $(ANALYSIS)/gene-counts/$(m)-rat.counts.h5$(NL))
	# Read bed files (starched)
	cp $(WORK)/make-read-beds/results/starch/*pruned-no-dedup.starch $(ANALYSIS)/filtered-read-locations/rna/
	cp $(WORK)/make-read-beds/results/starch/1*.pruned.starch $(ANALYSIS)/filtered-read-locations/atac/single-nucleus/
	cp $(WORK)/make-read-beds/results/starch/63*.pruned.starch $(ANALYSIS)/filtered-read-locations/atac/single-nucleus/
	cp $(WORK)/make-read-beds/results/starch/320-NM-1.pruned.starch $(ANALYSIS)/filtered-read-locations/atac/bulk/HSM1.bulk.rep1.pruned.starch
	cp $(WORK)/make-read-beds/results/starch/320-NM-2.pruned.starch $(ANALYSIS)/filtered-read-locations/atac/bulk/HSM1.bulk.rep2.pruned.starch
	cp $(WORK)/make-read-beds/results/starch/320-NM-3.pruned.starch $(ANALYSIS)/filtered-read-locations/atac/bulk/HSM2.bulk.rep1.pruned.starch
	cp $(WORK)/make-read-beds/results/starch/320-NM-4.pruned.starch $(ANALYSIS)/filtered-read-locations/atac/bulk/HSM2.bulk.rep2.pruned.starch
	# gzip things
	cd $(ANALYSIS)/atac-peaks/single-nucleus && gzip *.noblacklist
	cd $(ANALYSIS)/atac-peaks/bulk && gzip *.noblacklist
	cd $(ANALYSIS)/gene-counts && gzip *.counts.h5
	cd $(ANALYSIS)/filtered-read-locations/rna && gzip *
	cd $(ANALYSIS) && tar -cvf sn-muscle-processed.tar atac-bedgraphs atac-peaks clustering gene-counts filtered-read-locations luciferase sample-info.txt README.txt
	cp zenodo-noreads-readme.txt $(ANALYSIS)/README-no-reads.txt
	cp zenodo-reads-readme.txt $(ANALYSIS)/README-reads.txt
	cd $(ANALYSIS) && tar -cvf sn-muscle-processed-no-reads.tar atac-bedgraphs atac-peaks clustering gene-counts luciferase sample-info.txt README-no-reads.txt
	cd $(ANALYSIS) && tar -cvf sn-muscle-read-locations.tar filtered-read-locations sample-info.txt README-reads.txt
