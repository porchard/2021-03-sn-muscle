echo "ENCODE CREs:"
cat work/encode-cres/encode-cres.bed | cut -f1-4 | uniq | sort -k1,1 -k2n,2 | bedtools closest -d -a rs227731-hg19.bed -b stdin
echo "Meuleman CREs:"
zcat data/meuleman/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz | grep -v seqname | sort -k1,1 -k2n,2 | bedtools closest -d -a rs227731-hg38.bed -b stdin

