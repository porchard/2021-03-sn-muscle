#!/usr/bin/env nextflow
// start with merged bam files from atacseq
// get per library, per cell type, per individual bam
// encode as: {individual}.{library}.{cluster}.
//

BAM_GLOB = params.bam_glob
VCF = params.vcf

library_to_readgroup = {
	f ->
	f.tokenize('.')[0]
}

library_to_individual = {
	f ->
	f.tokenize('.')[1]
}

library_to_cluster = {
	f ->
	f.tokenize('.')[2]
}

get_sample = {
	f ->
	library_to_individual(f)
}


samples = ['KSM1', 'KSM2']
CHROMS = (1..22).collect({it -> "chr" + it})
GENOME = 'hg19'

bam_chan = Channel.fromPath(BAM_GLOB).map({it -> [it.getName().replaceAll('.bam', ''), it]})
libraries = Channel.fromPath(BAM_GLOB).map({it -> it.getName().replaceAll('.bam', '')})

get_bwa_index = {
	genome ->
	params.bwa_index[genome]
}

get_genome = {
	library ->
	GENOME
}

process make_heterozygous_positions {

	tag "${library}"
	memory '50 GB'
	maxForks 15

	input:
	file(vcf) from Channel.fromPath(VCF)
	each library from libraries

	output:
	set val(library), file("${library}.het-positions.bed") into het_positions
	set val(library), file("${library}.het-positions.bed") into het_positions_3

	"""
	zcat ${params.blacklist[GENOME].join(' ')} | sort -k1,1 -k2n,2 | gzip -c > blacklist.bed.gz
	bcftools view -Oz -s ${get_sample(library)} -T ^blacklist.bed.gz -o ${library}.genotypes.vcf.gz $vcf
	vcf2bed.pl -v ${library}.genotypes.vcf.gz -t -s -k > ${library}.het-positions.bed
	"""

}

process make_heterozygous_positions_sample {

	tag "${samp}"
	memory '50 GB'
	maxForks 15

	input:
	file(vcf) from Channel.fromPath(VCF)
	each samp from samples

	output:
	set val(samp), file("${samp}.het-positions.bed") into het_positions_2
	set val(samp), file("${samp}.het-positions.bed") into het_positions_4

	"""
	zcat ${params.blacklist[GENOME].join(' ')} | sort -k1,1 -k2n,2 | gzip -c > blacklist.bed.gz
	bcftools view -Oz -s $samp -T ^blacklist.bed.gz -o ${samp}.genotypes.vcf.gz $vcf
	vcf2bed.pl -v ${samp}.genotypes.vcf.gz -t -s -k > ${samp}.het-positions.bed
	"""

}


process split_vcf {

	tag "${chrom}"
	memory '50 GB'

	input:
	file(vcf) from Channel.fromPath(VCF)
	each chrom from CHROMS

	output:
	set val(chrom), file("${chrom}.genotypes.vcf.gz") into split_vcf_out

	"""
	bcftools view -Oz -o ${chrom}.genotypes.prefilter.vcf.gz -t $chrom $vcf
	zcat ${params.blacklist[GENOME].join(' ')} | sort -k1,1 -k2n,2 | gzip -c > blacklist.bed.gz
	bcftools view -Oz -o ${chrom}.genotypes.vcf.gz -T ^blacklist.bed.gz ${chrom}.genotypes.prefilter.vcf.gz
	"""

}

// TODO: use geno_prob if we can...
process make_snph5 {

	input:
	file(vcfs) from split_vcf_out.map({x -> x[1]}).toSortedList()

	output:
	set file("haplotypes.h5"), file("snp_index.h5"), file("snp_tab.h5") into make_snph5_out

	"""
	snp2h5 --format vcf --chrom ${params.chrom_sizes[GENOME]} ${vcfs.join(' ')} --geno_prob geno_probs.h5 --haplotype haplotypes.h5 --snp_index snp_index.h5 --snp_tab snp_tab.h5
	"""	

}


process prune_initial {

	tag "${library}"
	maxForks 10

	input:
	set val(library), file("${library}.unfiltered.bam") from bam_chan
	
	output:
	set val(library), file("${library}.bam"), file("${library}.bam.bai") into find_intersecting_in

	"""
	samtools index ${library}.unfiltered.bam
	samtools view -b -h -f 3 -F 4 -F 8 -F 256 -F 2048 -q 30 ${library}.unfiltered.bam ${CHROMS.join(' ')} | samtools sort -m 3g -O bam -T bam_sort -o ${library}.bam -
	samtools index ${library}.bam
	"""

}

process find_intersecting_snps {

	tag "${library}"
	maxForks 10

	input:
	set val(library), file(bam), file(bam_index), file(haplotypes), file(snp_index), file(snp_tab) from find_intersecting_in.combine(make_snph5_out)

	output:
	set val(library), file("${library}.keep.bam") into merge_no_remap_in
	set val(library), file("${library}.to.remap.bam") into filter_remapped_reads_keep_in
	set val(library), file("${library}.remap.fq1.gz"), file("${library}.remap.fq2.gz") into remap_in


	"""
	find_intersecting_snps.py --is_paired_end --is_sorted --output_dir . --haplotype $haplotypes --snp_tab $snp_tab --snp_index $snp_index --samples ${get_sample(library)} $bam
	"""

}

process remap {

	tag "${library}"
	maxForks 10

	input:
	set val(library), file(fastq_1), file(fastq_2) from remap_in

	output:
	set val(library), file("${library}.remapped.bam") into remapped

	"""
	bwa mem -I 200,200,5000 -M ${get_bwa_index(GENOME)} ${fastq_1} ${fastq_2} | samtools sort -m 1g -O bam -T sort_tmp -o ${library}.remapped.bam -
	"""

}


process prune_remap {

	tag "${library}"

	input:
	set val(library), file("${library}.unfiltered.bam") from remapped

	output:
	set val(library), file("${library}.bam") into filter_remapped_reads_remapped_in
	
	"""
	samtools index ${library}.unfiltered.bam
	samtools view -b -h -f 3 -F 4 -F 8 -F 256 -F 2048 -q 30 ${library}.unfiltered.bam ${CHROMS.join(' ')} | samtools sort -m 3g -O bam -T bam_sort -o ${library}.bam -
	samtools index ${library}.bam
	"""

}

process filter_remapped_reads {

	tag "${library}"
	memory '20 GB'

	input:
	set val(library), file(to_remap), file(remapped) from filter_remapped_reads_keep_in.combine(filter_remapped_reads_remapped_in, by: 0)

	output:
	set val(library), file("${library}.good_remappings.bam") into merge_remap_in

	"""
	filter_remapped_reads.py $to_remap $remapped ${library}.good_remappings.unsorted.bam
	samtools sort -m 3g -O bam -T bam_sort -o ${library}.good_remappings.bam ${library}.good_remappings.unsorted.bam
	"""

}


process merge {

	tag "${library}"

	input:
	set val(library), file(no_remap), file(remapped) from merge_no_remap_in.combine(merge_remap_in, by: 0)

	output:
	set val(library), file("${library}.bam"), file("${library}.bam.bai") into rmdup_in

	"""
	samtools merge ${library}.unsorted.bam $no_remap $remapped
	samtools sort -m 3g -O bam -T sort_tmp -o ${library}.bam ${library}.unsorted.bam
	samtools index ${library}.bam
	"""
}

process rmdup {

	tag "${library}"

	input:
	set val(library), file(bam), file(index) from rmdup_in

	output:
	set val(library), file("${library}.rmdup.bam") into rmdup_out

	"""
	rmdup_pe_sn.py $bam ${library}.rmdup.unsorted.bam
	samtools sort -m 3g -O bam -T sort_tmp -o ${library}.rmdup.bam ${library}.rmdup.unsorted.bam
	"""

}

process clip {

	maxForks 15

	input:
	set val(library), file(bam) from rmdup_out

	output:
	set val(library), file("${library}.clipped.bam") into clipped
	set val(library), file("${library}.clipped.bam") into clipped_2
	set val(library), file("${library}.clipped.bam") into clipped_3

	"""
	bam clipOverlap --poolSize 9000000 --in $bam --out ${library}.clipped.bam
	"""

}


merge_from_sample_in = clipped_2.map({x -> [get_sample(x[0].toString()), x[0], x[1]]}).groupTuple()

process merge_from_sample {

	input:
	set val(samp), val(library), file(bams) from merge_from_sample_in

	output:
	set val(samp), file("${samp}.merged-libraries.bam") into merged_libraries_out
	set val(samp), file("${samp}.merged-libraries.bam") into merged_libraries_out_2

	"""
	samtools merge ${samp}.merged-libraries.bam ${bams.join(' ')}
	"""

}


process nucleotide_counts {

	maxForks 30
	publishDir "${params.results}/nucleotide-counts"

	input:
	set val(library), file(bam), file(het) from clipped_3.combine(het_positions_3, by: 0).mix(merged_libraries_out_2.combine(het_positions_4, by: 0))

	output:
	set val(library), file(het), file("${library}.counts.txt") 

	"""
	cut -f1,3 $het > fetch.bed
	nucleotide-counts.py --min-base-quality 20 --min-mapping-quality 30 fetch.bed $bam > ${library}.counts.txt
	"""

}


process mpileup {

	maxForks 30
	publishDir "${params.results}/mpileup"

	input:
	set val(library), file(bam), file(het) from clipped.combine(het_positions, by: 0).mix(merged_libraries_out.combine(het_positions_2, by: 0))

	output:
	set val(library), file(het), file("${library}.mpileup.out") into post_mpileup_in

	"""
	perl /lab/work/porchard/endoC/src/allelic_bias/mpileup.pl -ref ${params.fasta[GENOME]} -bam $bam -bed $het -map 30 -base 20 -d 9999999 -v 6 -F 516 0 | cut -f1-7 > ${library}.mpileup.out
	"""

}

/*
process post_mpileup {

	maxForks 30
	publishDir "${params.results}/post-mpileup"

	input:
	set val(library), file(het), file(mpileup_out) from post_mpileup_in

	output:
	set val(library), file("${library}.post_mpileup.out"), file("${library}.post_mpileup.err")

	"""
	perl /lab/work/porchard/endoC/src/allelic_bias/post_mpileup.pl -b $het -m $mpileup_out -v 6 -a > ${library}.post_mpileup.out 2> ${library}.post_mpileup.err
	"""
}
*/
