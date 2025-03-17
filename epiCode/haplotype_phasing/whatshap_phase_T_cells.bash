# Haplotype phasing of T cells using whatshap
# Based on Oberon's code:
# This script is based on the script used to phase GM12878 cells in the epiCode repository
# The script is used to phase the T cells using the same pipeline as the GM12878 cells

# For K562 cells:
# Using previous data aligned to hg19:
# whatshap polyphase /home/michalula/ref_genomes/K562/ENCFF960SSF.vcf   /home/michalula/data/cas9_nanopore/data/2024927_Cas9_R9_promethion_K562/to_hg19/sort_aligned_hg19_trim_bam_pass_merged.bam --ploidy p --reference /home/michalula/ref_genomes/hg19/GRCh37.p13.genome.fa  -o /home/michalula/data/cas9_nanopore/data/2024927_Cas9_R9_promethion_K562/to_hg19/hyplotyped/whatshap_output.vcf
# whatshap polyphase /home/michalula/ref_genomes/K562/ENCFF960SSF.vcf \
#   /home/michalula/data/cas9_nanopore/data/2024927_Cas9_R9_promethion_K562/to_hg19/sort_aligned_hg19_trim_bam_pass_merged.bam \
#   --ploidy 2 \
#   --reference /home/michalula/ref_genomes/hg19/GRCh37.p13.genome.fa \
#   -o /home/michalula/data/cas9_nanopore/data/2024927_Cas9_R9_promethion_K562/to_hg19/hyplotyped/whatshap_output_K562.vcf



# Oberon installed whatshap, samtools, and picard from conda

# For GM12878 cells:
picard LiftoverVcf I=NA12878.vcf.gz O=NA12878.hs1.vcf.gz CHAIN=/clusterfs/nilah/oberon/genomes/hg38ToHs1.over.chain.gz REJECT=NA12878.hs1.rejected.vcg.gz R=/clusterfs/nilah/oberon/genomes/chm13v2.0.fa RECOVER_SWAPPED_REF_ALT=true

whatshap haplotag \
    --reference /clusterfs/nilah/oberon/genomes/chm13v2.0.fa \
    --ignore-read-groups \
    /clusterfs/nilah/oberon/datasets/na12878/NA12878.hs1.vcf.gz \
    /global/scratch/users/dixonluinenburg/dimelo-accessibility/repressive-histone/20250114_H3K27me3_access_realtime_pass.trim.align.mapq60.sorted.bam \
    | samtools view -@ 4 -b -o /global/scratch/users/dixonluinenburg/dimelo-accessibility/repressive-histone/20250114_H3K27me3_access_realtime_pass.trim.align.mapq60.haplotagged.sorted.bam

samtools view -h 20250114_H3K27me3_access_realtime_pass.trim.align.mapq60.haplotagged.sorted.bam -@ 32 | grep -E '^@|HP:i:1' | samtools view -@ 32 -b -o 20250114_H3K27me3_access_realtime_pass.trim.align.mapq60.hp1.sorted.bam -

samtools view -h 20250114_H3K27me3_access_realtime_pass.trim.align.mapq60.haplotagged.sorted.bam -@ 32 | grep -E '^@|HP:i:2' | samtools view -@ 32 -b -o 20250114_H3K27me3_access_realtime_pass.trim.align.mapq60.hp2.sorted.bam -


# Human Primary T cells from a donor to IGI ES

# Schema for 1000G Ph3 Vars - 1000 Genomes Phase 3 Integrated Variant Calls: SNVs, Indels, SVs
#  	Database: hg19    Primary Table: tgpPhase3
# VCF File Download: https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
# Format description: The fields of a Variant Call Format data line
# field	description
# chrom	An identifier from the reference genome
# pos	The reference position, with the 1st base having position 1
# id	Semi-colon separated list of unique identifiers where available
# ref	Reference base(s)
# alt	Comma separated list of alternate non-reference alleles called on at least one of the samples
# qual	Phred-scaled quality score for the assertion made in ALT. i.e. give -10log_10 prob(call in ALT is wrong)
# filter	PASS if this position has passed all filters. Otherwise, a semicolon-separated list of codes for filters that fail
# info	Additional information encoded as a semicolon-separated series of short keys with optional comma-separated values
# format	If genotype columns are specified in header, a semicolon-separated list of of short keys starting with GT
# genotypes	If genotype columns are specified in header, a tab-separated set of genotype column values; each value is a colon-separated list of values corresponding to keys in the format column

cd /home/michalula/data/ref_genomes/hg19/1000Genomes_Ph3_Vars
wget  https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
conda activate whatshap-env    
conda install -c bioconda picard

# Get T2T-CHM13v2.0 reference genome:
cd /home/michalula/data/ref_genomes/t2t_v2_0
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
gunzip chm13v2.0.fa.gz

# # Get VCF files:
# Variant calls
# 1000 Genomes Project, recalled on T2T-CHM13v2.0. Now available for all chromosomes, for the entire 3,202 samples or the unrelated 2504 samples. Reference sets, bam, and vcf files are also available on AnVIL_T2T_CHRY.
# * 1KGP variant calls for all chromosomes. Jan. 3 2023. Annotation update.
# * 1KGP and SGDP bam / vcf released publicly on [AnVIL_T2T_CHRY](https://anvil.terra.bio/#workspaces/anvil-datastorage/AnVIL_T2T_CHRY). May 23, 2023. Data Update.
# * 1KGP AF release. Jul 6 2023. Annotation update.
# From: https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/all_samples_3202/
cd /home/michalula/data/ref_genomes/t2t_v2_0/haplotype_vcf/1000Genomes
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/all_samples_3202/1KGP.CHM13v2.0.chr1.recalibrated.snp_indel.pass.vcf.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/all_samples_3202/1KGP.CHM13v2.0.chr1.recalibrated.snp_indel.pass.vcf.gz.tbi


# Haplotype phasing of T cells using whatshap
cd /home/michalula/data/cas9_nanopore/data/20241226_MR_nCATs_TcellsPrES_unedit_P2R9/passed_fast5/5mCG/to_t2t_v2_0/
mkdir haplotyped

# whatshap haplotag \
#     --reference /home/michalula/data/ref_genomes/t2t_v2_0/chm13v2.0.fa \
#     --ignore-read-groups \
#     /home/michalula/data/ref_genomes/t2t_v2_0/haplotype_vcf/1000Genomes/1KGP.CHM13v2.0.chr1.recalibrated.snp_indel.pass.vcf.gz \
#     ./sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed.dna_r9.4.1_e8_sup@v3.3.5mCG.bam \
#     | samtools view -@ 4 -b -o ./haplotyped/sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed.dna_r9.4.1_e8_sup@v3.3.5mCG.haplotagged.bam
# ERROR about the vcf file
conda install -c bioconda bcftools
bcftools query -l /home/michalula/data/ref_genomes/t2t_v2_0/haplotype_vcf/1000Genomes/1KGP.CHM13v2.0.chr1.recalibrated.snp_indel.pass.vcf.gz
# ls /home/michalula/data/ref_genomes/t2t_v2_0/haplotype_vcf/
samtools faidx /home/michalula/data/ref_genomes/t2t_v2_0/chm13v2.0.fa
ls -lh /home/michalula/data/ref_genomes/t2t_v2_0/chm13v2.0.fa.fai

# let's randomly select NA21116 #TODO - explore which sample to use!!! 
whatshap haplotag \
    --reference /home/michalula/data/ref_genomes/t2t_v2_0/chm13v2.0.fa \
    --ignore-read-groups \
    --sample NA21116 \
    /home/michalula/data/ref_genomes/t2t_v2_0/haplotype_vcf/1000Genomes/1KGP.CHM13v2.0.chr1.recalibrated.snp_indel.pass.vcf.gz \
    ./sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed.dna_r9.4.1_e8_sup@v3.3.5mCG.bam \
    | samtools view -@ 4 -b -o ./haplotyped/sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed.dna_r9.4.1_e8_sup@v3.3.5mCG.haplotagged.bam



# samtools view -h 20250114_H3K27me3_access_realtime_pass.trim.align.mapq60.haplotagged.sorted.bam -@ 32 | grep -E '^@|HP:i:1' | samtools view -@ 32 -b -o 20250114_H3K27me3_access_realtime_pass.trim.align.mapq60.hp1.sorted.bam -
samtools view -h ./haplotyped/sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed.dna_r9.4.1_e8_sup@v3.3.5mCG.haplotagged.bam  -@ 32 | grep -E '^@|HP:i:1' | samtools view -@ 32 -b -o  ./haplotyped/sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed.dna_r9.4.1_e8_sup@v3.3.5mCG.Haplotype_1.bam 


# samtools view -h 20250114_H3K27me3_access_realtime_pass.trim.align.mapq60.haplotagged.sorted.bam -@ 32 | grep -E '^@|HP:i:2' | samtools view -@ 32 -b -o 20250114_H3K27me3_access_realtime_pass.trim.align.mapq60.hp2.sorted.bam -
samtools view -h ./haplotyped/sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed.dna_r9.4.1_e8_sup@v3.3.5mCG.haplotagged.bam  -@ 32 | grep -E '^@|HP:i:2' | samtools view -@ 32 -b -o  ./haplotyped/sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed.dna_r9.4.1_e8_sup@v3.3.5mCG.Haplotype_2.bam 

