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


conda activate whatshap-env 
# let's randomly select NA21116 #TODO - explore which sample to use!!! 
# whatshap haplotag \
#     --reference /home/michalula/data/ref_genomes/t2t_v2_0/chm13v2.0.fa \
#     --ignore-read-groups \
#     --sample NA21116 \
#     /home/michalula/data/ref_genomes/t2t_v2_0/haplotype_vcf/1000Genomes/1KGP.CHM13v2.0.chr1.recalibrated.snp_indel.pass.vcf.gz \
#     ./sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed.dna_r9.4.1_e8_sup@v3.3.5mCG.bam \
#     | samtools view -@ 4 -b -o ./haplotyped/sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed.dna_r9.4.1_e8_sup@v3.3.5mCG.haplotagged.bam
# # Started running whatshap haplotag ~3:20pm 
# # Finished running whatshap haplotag ~
# EROOR
# ChatGPT fix:
# whatshap haplotag \
#     --reference /home/michalula/data/ref_genomes/t2t_v2_0/chm13v2.0.fa \
#     --ignore-read-groups \
#     --sample NA21116 \
#     /home/michalula/data/ref_genomes/t2t_v2_0/haplotype_vcf/1000Genomes/1KGP.CHM13v2.0.chr1.recalibrated.snp_indel.pass.vcf.gz \
#     ./sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed.dna_r9.4.1_e8_sup@v3.3.5mCG.bam \
#     -o ./haplotyped/tmp_haplotagged.bam

# Step 2: Verify the WhatsHap Output BAM File
# After running the command, check if the BAM file exists and has content:
ls -lh ./haplotyped/tmp_haplotagged.bam

# Then, check the BAM file header:
samtools view -H ./haplotyped/tmp_haplotagged.bam
# Possible Outcomes
# If the BAM file is empty or missing → WhatsHap did not generate a valid file.
# If samtools view -H works → WhatsHap produced a valid file, and the issue is in piping.

# whatshap haplotag \
#     --reference /home/michalula/data/ref_genomes/t2t_v2_0/chm13v2.0.fa \
#     --ignore-read-groups \
#     --sample NA21116 \
#     /home/michalula/data/ref_genomes/t2t_v2_0/haplotype_vcf/1000Genomes/1KGP.CHM13v2.0.chr1.recalibrated.snp_indel.pass.vcf.gz \
#     ./sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed.dna_r9.4.1_e8_sup@v3.3.5mCG.bam \
#     -o ./haplotyped/tmp_haplotagged.bam


#  Reduce BAM File Size (Process by Chromosome)
# If memory is still an issue, process one chromosome at a time instead of the whole BAM:

samtools view -b ./sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed.dna_r9.4.1_e8_sup@v3.3.5mCG.bam chr1 > ./chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam

# Have to resort and reindex the BAM subser of chr1
samtools sort \
  -o ./sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam  \
   ./chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam   

samtools index \
  ./sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam
  

# start here: in whatshap-env
cd /home/michalula/data/cas9_nanopore/data/20241226_MR_nCATs_TcellsPrES_unedit_P2R9/passed_fast5/5mCG/to_t2t_v2_0/
whatshap haplotag \
    --reference /home/michalula/data/ref_genomes/t2t_v2_0/chm13v2.0.fa \
    --ignore-read-groups \
    --sample NA21116 \
    /home/michalula/data/ref_genomes/t2t_v2_0/haplotype_vcf/1000Genomes/1KGP.CHM13v2.0.chr1.recalibrated.snp_indel.pass.vcf.gz \
    ./sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam \
    -o ./haplotyped/haplotagged.chr1.sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed.dna_r9.4.1_e8_sup@v3.3.5mCG.bam
# Started running whatshap haplotag 7:40am on March 18, 2025
# Finished running whatshap haplotag ~

# Run $terminal> top > to see whe memory usage:
# top - 22:44:56 up 19 days,  3:45,  1 user,  load average: 1.32, 1.09, 0.65
# Tasks: 439 total,   3 running, 436 sleeping,   0 stopped,   0 zombie
# %Cpu(s):  3.9 us,  1.0 sy,  0.0 ni, 95.1 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
# MiB Mem : 128554.4 total, 102319.3 free,  16319.1 used,   9915.9 buff/cache
# MiB Swap:   2048.0 total,    404.1 free,   1643.8 used. 111001.5 avail Mem 

#     PID USER      PR  NI    VIRT    RES    SHR S  %CPU  %MEM     TIME+ COMMAND                                                               
# 1409811 michalu+  20   0   17.9g  14.3g  11136 R      100.0  11.4   4:12.89 whatshap   
#     PID USER      PR  NI    VIRT    RES    SHR S  %CPU  %MEM     TIME+ COMMAND                                                               
# 1409811 michalu+  20   0   44.3g  35.5g  11136 R       99.7  28.3   6:37.17 whatshap   

# top - 22:51:39 up 19 days,  3:52,  1 user,  load average: 1.12, 1.12, 0.84
# Tasks: 439 total,   4 running, 435 sleeping,   0 stopped,   0 zombie
# %Cpu(s):  4.2 us,  0.5 sy,  0.0 ni, 95.2 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
# MiB Mem : 128554.4 total,  35590.3 free,  72099.2 used,  20864.8 buff/cache
# MiB Swap:   2048.0 total,    418.5 free,   1629.5 used.  55221.4 avail Mem 
#     PID USER      PR  NI    VIRT    RES    SHR S  %CPU  %MEM     TIME+ COMMAND                                                               
# 1409811 michalu+  20   0   85.3g  68.6g  11136 R  99.7  54.7  10:55.69 whatshap   


# top - 22:54:40 up 19 days,  3:55,  1 user,  load average: 1.12, 1.13, 0.91
# Tasks: 437 total,   3 running, 434 sleeping,   0 stopped,   0 zombie
# %Cpu(s):  4.8 us,  0.3 sy,  0.0 ni, 94.9 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
# MiB Mem : 128554.4 total,   7328.4 free,  95427.8 used,  25798.2 buff/cache
# MiB Swap:   2048.0 total,    424.3 free,   1623.7 used.  31892.8 avail Mem 

#     PID USER      PR  NI    VIRT    RES    SHR S  %CPU  %MEM     TIME+ COMMAND                                                               
# 1409811 michalu+  20   0  113.5g  91.3g  11136 R 100.0  72.7  13:56.05 whatshap   

# top - 22:55:52 up 19 days,  3:56,  1 user,  load average: 1.19, 1.16, 0.94
# Tasks: 437 total,   2 running, 435 sleeping,   0 stopped,   0 zombie
# %Cpu(s):  4.1 us,  0.2 sy,  0.0 ni, 95.7 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
# MiB Mem : 128554.4 total,    897.0 free, 104778.7 used,  22878.7 buff/cache
# MiB Swap:   2048.0 total,    417.8 free,   1630.2 used.  22541.9 avail Mem 
#     PID USER      PR  NI    VIRT    RES    SHR S  %CPU  %MEM     TIME+ COMMAND                                                               
# 1409811 michalu+  20   0  124.9g 100.5g  11136 R      99.7  80.0  15:08.06 whatshap     

# top - 22:57:31 up 19 days,  3:58,  1 user,  load average: 1.17, 1.15, 0.96
# Tasks: 438 total,   3 running, 435 sleeping,   0 stopped,   0 zombie
# %Cpu(s):  4.1 us,  0.5 sy,  0.0 ni, 95.4 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
# MiB Mem : 128554.4 total,    864.2 free, 117336.6 used,  10353.6 buff/cache
# MiB Swap:   2048.0 total,    414.4 free,   1633.6 used.   9984.3 avail Mem 

#     PID USER      PR  NI    VIRT    RES    SHR S  %CPU  %MEM     TIME+ COMMAND                                                               
# 1409811 michalu+  20   0  140.2g 112.7g  11136 R 100.0  89.8  16:47.26 whatshap               

# top - 22:58:46 up 19 days,  3:59,  1 user,  load average: 1.21, 1.18, 0.99
# Tasks: 439 total,   4 running, 435 sleeping,   0 stopped,   0 zombie
# %Cpu(s):  4.1 us,  0.5 sy,  0.0 ni, 95.4 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
# MiB Mem : 128554.4 total,    917.4 free, 126972.3 used,    664.6 buff/cache
# MiB Swap:   2048.0 total,    402.3 free,   1645.7 used.    442.6 avail Mem 
#     PID USER      PR  NI    VIRT    RES    SHR S  %CPU  %MEM     TIME+ COMMAND                                                               
# 1409811 michalu+  20   0  151.9g 122.1g   8448 R  99.7  97.2  18:02.37 whatshap          

# top - 22:59:07 up 19 days,  4:00,  1 user,  load average: 2.29, 1.41, 1.07
# Tasks: 438 total,   4 running, 434 sleeping,   0 stopped,   0 zombie
# %Cpu(s):  0.7 us, 14.5 sy,  0.0 ni, 71.9 id, 12.8 wa,  0.0 hi,  0.0 si,  0.0 st
# MiB Mem : 128554.4 total,    579.9 free, 127509.2 used,    465.3 buff/cache
# MiB Swap:   2048.0 total,      0.1 free,   2047.9 used.    862.7 avail Mem 

#     PID USER      PR  NI    VIRT    RES    SHR S  %CPU  %MEM     TIME+ COMMAND                                                               
# 1138068 root      20   0 2430752   8488      0 S 172.5   0.0   1:18.06 snapd                                                                 
#     191 root      20   0       0      0      0 R  99.3   0.0  11:07.88 kswapd0                                                               
# 1409811 michalu+  20   0  153.2g 122.7g      0 R  28.5  97.7  18:12.77 whatshap                                                              
# 1407530 michalu+  20   0   53.0g 541368    960 R  18.5   0.4   0:24.27 node             

# top - 23:02:48 up 19 days,  4:03,  1 user,  load average: 3.15, 7.76, 4.33
# Tasks: 438 total,   1 running, 437 sleeping,   0 stopped,   0 zombie
# %Cpu(s):  0.1 us,  0.1 sy,  0.0 ni, 99.7 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
# MiB Mem : 128554.4 total, 126199.1 free,   1667.0 used,    688.3 buff/cache
# MiB Swap:   2048.0 total,    431.9 free,   1616.1 used. 125731.0 avail Mem 

#     PID USER      PR  NI    VIRT    RES    SHR S  %CPU  %MEM     TIME+ COMMAND                                                               
#  421741 minknow   35  15 2374072   6448   4404 S   2.7   0.0 481:07.57 control_main                                                          
#  421742 minknow   35  15 6952512   7188   5388 S   2.3   0.0   1424:31 control_main                                                          
# 1407500 michalu+  20   0   11.3g  60444  22464 S   0.7   0.0   0:42.33 node                                                                  
# 1407530 michalu+  20   0   53.0g 537372  21120 S   0.7   0.4   1:09.46 node   
> Killed process.. so selecting a subset of interest 

samtools view -b ./sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam \
    chr1:206560169-206614236 > \
    chr1_206560169_206614236.sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam
ls -lh chr1_206560169_206614236.sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam
# calculate how many reads are saved 
samtools view -c chr1_206560169_206614236.sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam

samtools index chr1_206560169_206614236.sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam

ls -lah  /home/michalula/data/cas9_nanopore/data/20241226_MR_nCATs_TcellsPrES_unedit_P2R9/passed_fast5/5mCG/to_t2t_v2_0/

# total 13G
# drwxrwxr-x 3 michalula michalula 4.0K Mar 17 23:13 .
# drwxrwxrwx 4 michalula michalula 4.0K Mar 17 04:54 ..
# -rw-rw-r-- 1 michalula michalula 5.2G Mar 17 05:26 align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed.dna_r9.4.1_e8_sup@v3.3.5mCG.bam
# -rw-rw-r-- 1 michalula michalula  24M Mar 17 23:13 chr1_206560169_206614236.sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam
# -rw-rw-r-- 1 michalula michalula 104K Mar 17 23:13 chr1_206560169_206614236.sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam.bai
# -rw-r--r-- 1 michalula michalula 482M Mar 17 10:35 chr1.sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed.dna_r9.4.1_e8_sup@v3.3.5mCG.bam
# -rw-rw-r-- 1 michalula michalula 482M Mar 17 11:43 chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed.dna_r9.4.1_e8_sup@v3.3.5mCG.bam
# -rw-rw-r-- 1 michalula michalula 482M Mar 17 11:43 chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_sup@v3mCG.bam
# -rw-rw-r-- 1 michalula michalula 482M Mar 17 11:43 chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam
# drwxrwxr-x 2 michalula michalula 4.0K Mar 17 23:19 haplotyped
# -rw-rw-r-- 1 michalula michalula 2.0K Mar 17 10:31 script_t2t_v2_align.bash
# -rw-rw-r-- 1 michalula michalula 5.1G Mar 17 05:30 sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed.dna_r9.4.1_e8_sup@v3.3.5mCG.bam
# -rw-rw-r-- 1 michalula michalula 3.0M Mar 17 05:37 sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed.dna_r9.4.1_e8_sup@v3.3.5mCG.bam.bai
# -rw-rw-r-- 1 michalula michalula 482M Mar 17 11:49 sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam
# -rw-rw-r-- 1 michalula michalula 267K Mar 17 11:49 sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam.bai
# -rw-rw-r-- 1 michalula michalula 214M Mar 17 05:42 summary_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed.dna_r9.4.1_e8_sup@v3.3.5mCG.tsv


whatshap haplotag \
    --reference /home/michalula/data/ref_genomes/t2t_v2_0/chm13v2.0.fa \
    --ignore-read-groups \
    --sample NA21116 \
    /home/michalula/data/ref_genomes/t2t_v2_0/haplotype_vcf/1000Genomes/1KGP.CHM13v2.0.chr1.recalibrated.snp_indel.pass.vcf.gz \
    chr1_206560169_206614236.sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam \
    -o ./haplotyped/haplotagged.chr1_206560169_206614236.sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam


# top - 23:17:49 up 19 days,  4:18,  1 user,  load average: 0.79, 0.73, 1.79
# Tasks: 438 total,   2 running, 436 sleeping,   0 stopped,   0 zombie
# %Cpu(s):  4.2 us,  0.1 sy,  0.0 ni, 95.7 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
# MiB Mem : 128554.4 total,    825.0 free, 120410.2 used,   7319.2 buff/cache
# MiB Swap:   2048.0 total,    439.0 free,   1609.0 used.   6894.5 avail Mem 
#     PID USER      PR  NI    VIRT    RES    SHR S  %CPU  %MEM     TIME+ COMMAND                                                                                                  
# 1427898 michalu+  20   0  115.8g 115.7g   7296 R 100.0  92.2   1:18.11 whatshap  
# top - 23:20:14 up 19 days,  4:21,  1 user,  load average: 1.17, 0.93, 1.71
# Tasks: 439 total,   2 running, 437 sleeping,   0 stopped,   0 zombie
# %Cpu(s):  4.1 us,  0.2 sy,  0.0 ni, 95.6 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
# MiB Mem : 128554.4 total, 111118.8 free,   9040.5 used,   8395.0 buff/cache
# MiB Swap:   2048.0 total,    439.7 free,   1608.3 used. 118280.3 avail Mem 
#     PID USER      PR  NI    VIRT    RES    SHR S  %CPU  %MEM     TIME+ COMMAND                                                                                                  
# 1427898 michalu+  20   0 9348468   7.1g  11328 R 100.0   5.7   3:42.40 whatshap    

# top - 23:27:11 up 19 days,  4:28,  1 user,  load average: 1.08, 1.07, 1.49
# Tasks: 438 total,   3 running, 435 sleeping,   0 stopped,   0 zombie
# %Cpu(s):  4.1 us,  0.2 sy,  0.0 ni, 95.7 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
# MiB Mem : 128554.4 total,  38325.9 free,  70322.2 used,  19906.2 buff/cache
# MiB Swap:   2048.0 total,    447.0 free,   1601.0 used.  56998.5 avail Mem 
#     PID USER      PR  NI    VIRT    RES    SHR S  %CPU  %MEM     TIME+ COMMAND                                                               
# 1427898 michalu+  20   0   83.1g  66.9g  11328 R  99.3  53.3  10:40.30 whatshap                                                              
#  421742 minknow   35  15 6952512   7188   5388 S   1.3   0.0   1424:55 control_main                                                          
#  421741 minknow   35  15 2374072   6256   4404 S   1.0   0.0 481:31.90 control_main    

# top - 23:29:15 up 19 days,  4:30,  1 user,  load average: 1.23, 1.13, 1.45
# Tasks: 440 total,   2 running, 438 sleeping,   0 stopped,   0 zombie
# %Cpu(s):  4.1 us,  1.4 sy,  0.0 ni, 94.5 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
# MiB Mem : 128554.4 total,  18904.7 free,  86385.6 used,  23264.1 buff/cache
# MiB Swap:   2048.0 total,    452.7 free,   1595.3 used.  40935.2 avail Mem 
#     PID USER      PR  NI    VIRT    RES    SHR S  %CPU  %MEM     TIME+ COMMAND                                                               
# 1427898 michalu+  20   0  102.5g  82.5g  11328 R 100.0  65.7  12:43.56 whatshap                                                              
#     814 root      20   0       0      0      0 S  19.9   0.0   2165:10 nvidia-modeset/kthread_q                                              
#  938520 gdm       20   0   24.2g  13796  12672 S   7.3   0.0 438:38.93 Xorg      
# top - 23:30:21 up 19 days,  4:31,  1 user,  load average: 1.13, 1.12, 1.43
# Tasks: 442 total,   2 running, 440 sleeping,   0 stopped,   0 zombie
# %Cpu(s):  4.1 us,  0.2 sy,  0.0 ni, 95.7 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
# MiB Mem : 128554.4 total,   8637.7 free,  94830.3 used,  25086.4 buff/cache
# MiB Swap:   2048.0 total,    452.7 free,   1595.3 used.  32490.5 avail Mem 
#     PID USER      PR  NI    VIRT    RES    SHR S  %CPU  %MEM     TIME+ COMMAND                                                               
# 1427898 michalu+  20   0  112.8g  90.7g  11328 R  99.7  72.2  13:49.70 whatshap                                                              
#  421741 minknow   35  15 2374072   6256   4404 S   1.0   0.0 481:33.84 control_main                                                          
#  421742 minknow   35  15 6952512   7188   5388 S   1.0   0.0   1424:57 control_main                  
# top - 23:32:06 up 19 days,  4:33,  1 user,  load average: 1.14, 1.12, 1.39
# Tasks: 440 total,   3 running, 437 sleeping,   0 stopped,   0 zombie
# %Cpu(s):  4.1 us,  0.7 sy,  0.0 ni, 95.2 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
# MiB Mem : 128554.4 total,    852.8 free, 108210.2 used,  19491.3 buff/cache
# MiB Swap:   2048.0 total,    446.8 free,   1601.2 used.  19110.6 avail Mem 
#     PID USER      PR  NI    VIRT    RES    SHR S  %CPU  %MEM     TIME+ COMMAND                                                               
# 1427898 michalu+  20   0  129.1g 103.8g  11328 R 100.0  82.7  15:34.92 whatshap                                                              
#     814 root      20   0       0      0      0 R   7.0   0.0   2165:23 nvidia-modeset/kthread_q                                              
#  938520 gdm       20   0   24.2g  13796  12672 S   4.3   0.0 438:43.88 Xorg            
# top - 23:33:00 up 19 days,  4:33,  1 user,  load average: 1.19, 1.13, 1.38
# Tasks: 442 total,   2 running, 440 sleeping,   0 stopped,   0 zombie
# %Cpu(s):  4.1 us,  0.2 sy,  0.0 ni, 95.6 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
# MiB Mem : 128554.4 total,    773.8 free, 115228.3 used,  12552.2 buff/cache
# MiB Swap:   2048.0 total,    444.1 free,   1603.9 used.  12092.5 avail Mem 
#     PID USER      PR  NI    VIRT    RES    SHR S  %CPU  %MEM     TIME+ COMMAND                                                               
# 1427898 michalu+  20   0  137.6g 110.6g  11328 R 100.0  88.1  16:29.03 whatshap                                                              
#  421741 minknow   35  15 2374072   6256   4404 S   1.0   0.0 481:35.52 control_main                                                          
# 1407500 michalu+  20   0   11.3g  66820  23040 S   1.0   0.1   0:44.64 node            
# top - 23:33:42 up 19 days,  4:34,  1 user,  load average: 1.09, 1.11, 1.36
# Tasks: 440 total,   2 running, 438 sleeping,   0 stopped,   0 zombie
# %Cpu(s):  4.1 us,  1.5 sy,  0.0 ni, 94.4 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
# MiB Mem : 128554.4 total,    749.8 free, 120576.1 used,   7228.5 buff/cache
# MiB Swap:   2048.0 total,    444.1 free,   1603.9 used.   6744.7 avail Mem 
#     PID USER      PR  NI    VIRT    RES    SHR S  %CPU  %MEM     TIME+ COMMAND                                                               
# 1427898 michalu+  20   0  144.1g 115.8g  11328 R 100.0  92.2  17:11.12 whatshap                                                              
#     814 root      20   0       0      0      0 S  20.3   0.0   2165:30 nvidia-modeset/kthread_q                                              
#  938520 gdm       20   0   24.2g  13796  12672 S   7.3   0.0 438:46.81 Xorg            
# top - 23:34:15 up 19 days,  4:35,  1 user,  load average: 1.05, 1.10, 1.35
# Tasks: 443 total,   2 running, 441 sleeping,   0 stopped,   0 zombie
# %Cpu(s):  4.1 us,  1.2 sy,  0.0 ni, 94.7 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
# MiB Mem : 128554.4 total,    731.7 free, 124656.6 used,   3166.0 buff/cache
# MiB Swap:   2048.0 total,    436.8 free,   1611.2 used.   2664.2 avail Mem 
#     PID USER      PR  NI    VIRT    RES    SHR S  %CPU  %MEM     TIME+ COMMAND                                                               
# 1427898 michalu+  20   0  149.0g 119.8g  10752 R 100.0  95.4  17:44.18 whatshap                                                              
#     814 root      20   0       0      0      0 S  16.9   0.0   2165:33 nvidia-modeset/kthread_q                                              
#  938520 gdm       20   0   24.2g  13796  12672 S   5.6   0.0 438:47.74 Xorg            
# top - 23:34:33 up 19 days,  4:35,  1 user,  load average: 1.11, 1.11, 1.34
# Tasks: 439 total,   3 running, 436 sleeping,   0 stopped,   0 zombie
# %Cpu(s):  4.1 us,  0.5 sy,  0.0 ni, 95.4 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
# MiB Mem : 128554.4 total,    828.3 free, 126927.7 used,    798.4 buff/cache
# MiB Swap:   2048.0 total,    203.7 free,   1844.3 used.    415.4 avail Mem 
#     PID USER      PR  NI    VIRT    RES    SHR S  %CPU  %MEM     TIME+ COMMAND                                                               
# 1427898 michalu+  20   0  152.1g 122.0g   8832 R  99.7  97.2  18:02.21 whatshap                                                              
#     191 root      20   0       0      0      0 R   6.6   0.0  12:17.61 kswapd0                                                               
#  421742 minknow   35  15 6952512   7188   5388 S   1.3   0.0   1425:00 control_main    
# top - 23:34:55 up 19 days,  4:35,  1 user,  load average: 2.62, 1.44, 1.45
# Tasks: 440 total,   5 running, 435 sleeping,   0 stopped,   0 zombie
# %Cpu(s):  0.6 us, 20.8 sy,  0.0 ni, 65.2 id, 13.4 wa,  0.0 hi,  0.0 si,  0.0 st
# MiB Mem : 128554.4 total,    579.8 free, 127511.1 used,    463.5 buff/cache
# MiB Swap:   2048.0 total,      0.1 free,   2047.9 used.    861.0 avail Mem 
#     PID USER      PR  NI    VIRT    RES    SHR S  %CPU  %MEM     TIME+ COMMAND                                                               
# 1138068 root      20   0 2430752   9032      0 S 251.0   0.0   2:49.68 snapd                                                                 
#     191 root      20   0       0      0      0 R 100.0   0.0  12:36.01 kswapd0                                                               
# 1427898 michalu+  20   0  153.2g 122.9g    960 D  32.5  97.9  18:12.24 whatshap      
# top - 23:37:05 up 19 days,  4:38,  1 user,  load average: 11.89, 8.91, 4.45
# Tasks: 451 total,   2 running, 449 sleeping,   0 stopped,   0 zombie
# %Cpu(s):  5.2 us,  0.4 sy,  0.0 ni, 94.2 id,  0.2 wa,  0.0 hi,  0.0 si,  0.0 st
# MiB Mem : 128554.4 total, 123708.3 free,   3303.3 used,   1542.7 buff/cache
# MiB Swap:   2048.0 total,    476.0 free,   1572.0 used. 124017.5 avail Mem 

#     PID USER      PR  NI    VIRT    RES    SHR S  %CPU  %MEM     TIME+ COMMAND                                                               
# 1438724 michalu+  20   0   11.1g 169816  41088 S  96.3   0.1   0:14.61 node                                                                  
# 1438713 michalu+  20   0   53.2g 901620  51840 R  33.9   0.7   0:06.36 node   


 Reduce the Number of Samples in the VCF
# Your VCF contains 3202 samples. By default, WhatsHap loads all samples into memory, even if you're only phasing one sample. This causes extreme memory consumption.
# Instead of using the full 1000 Genomes VCF, extract only your sample (NA21116) and create a sample-specific VCF.
# Step 1: Extract Only NA21116 from VCF 
conda install -c bioconda bcftools

bcftools view -s NA21116 -Oz -o /home/michalula/data/ref_genomes/t2t_v2_0/haplotype_vcf/1000Genomes/1sample_NA21116.1KGP.CHM13v2.0.chr1.recalibrated.snp_indel.pass.vcf.gz \
    /home/michalula/data/ref_genomes/t2t_v2_0/haplotype_vcf/1000Genomes/1KGP.CHM13v2.0.chr1.recalibrated.snp_indel.pass.vcf.gz
# # THIS TAKES MORE THAN 106 mins
# cp /home/michalula/data/ref_genomes/t2t_v2_0/haplotype_vcf/1000Genomes/1sample_NA21116.1KGP.CHM13v2.0.chr1.recalibrated.snp_indel.pass.vcf.gz \
#    /home/michalula/data/ref_genomes/t2t_v2_0/haplotype_vcf/1000Genomes/1_sample.NA21116.1KGP.CHM13v2.0.chr1.recalibrated.snp_indel.pass.vcf.gz

# mv /home/michalula/data/ref_genomes/t2t_v2_0/haplotype_vcf/1000Genomes/NA21116.1KGP.CHM13v2.0.chr1.recalibrated.snp_indel.pass.vcf.gz \
#    /home/michalula/data/ref_genomes/t2t_v2_0/haplotype_vcf/1000Genomes/1sample_NA21116.1KGP.CHM13v2.0.chr1.recalibrated.snp_indel.pass.vcf.gz
# Verify the Extraction
# Check that the new VCF contains only NA21116:
bcftools query -l /home/michalula/data/ref_genomes/t2t_v2_0/haplotype_vcf/1000Genomes/1sample_NA21116.1KGP.CHM13v2.0.chr1.recalibrated.snp_indel.pass.vcf.gz
# Index the new VCF:
bcftools index /home/michalula/data/ref_genomes/t2t_v2_0/haplotype_vcf/1000Genomes/1sample_NA21116.1KGP.CHM13v2.0.chr1.recalibrated.snp_indel.pass.vcf.gz
# Step 2: Use the Smaller VCF with WhatsHap
# Now, run WhatsHap using the sample-specific VCF instead of the full 1000 Genomes VCF:
whatshap haplotag \
    --reference /home/michalula/data/ref_genomes/t2t_v2_0/chm13v2.0.fa \
    --ignore-read-groups \
    --sample NA21116 \
    /home/michalula/data/ref_genomes/t2t_v2_0/haplotype_vcf/1000Genomes/1sample_NA21116.1KGP.CHM13v2.0.chr1.recalibrated.snp_indel.pass.vcf.gz \
    chr1_206560169_206614236.sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam \
    -o ./haplotyped/haplotagged.chr1_206560169_206614236.bam




# top - 01:31:24 up 19 days,  6:32,  1 user,  load average: 0.46, 0.59, 0.88
# Tasks: 438 total,   3 running, 435 sleeping,   0 stopped,   0 zombie
# %Cpu(s):  3.9 us,  0.4 sy,  0.0 ni, 95.6 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
# MiB Mem : 128554.4 total,    858.3 free,   7899.9 used, 119796.1 buff/cache
# MiB Swap:   2048.0 total,   1036.8 free,   1011.2 used. 119420.8 avail Mem 

# Found 1 sample(s) in input VCF
# Keeping 1 sample(s) for haplo-tagging
# Found 1 sample(s) in BAM file

Number of supplementary alignments: 0
Number of non-singleton groups: 0
Skipped 0 groups
Found 0 reads covering 0 variants

== SUMMARY ==
Total alignments processed:                      3876
Alignments that could be tagged:                    0
Alignments spanning multiple phase sets:            0
Finished in 71.3 s 

>> EMPTY files... ? :(
ls -lah ./haplotyped/haplotagged.chr1_206560169_206614236.bam
)

#Phase all chr 1:
whatshap haplotag \
    --reference /home/michalula/data/ref_genomes/t2t_v2_0/chm13v2.0.fa \
    --ignore-read-groups \
    --sample NA21116 \
    /home/michalula/data/ref_genomes/t2t_v2_0/haplotype_vcf/1000Genomes/1sample_NA21116.1KGP.CHM13v2.0.chr1.recalibrated.snp_indel.pass.vcf.gz \
    sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam \
    -o ./haplotyped/haplotagged.chr1.sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam

mv ./haplotyped/haplotagged.chr1.sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCGbam ./haplotyped/haplotagged.chr1.sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam


# samtools view -h 20250114_H3K27me3_access_realtime_pass.trim.align.mapq60.haplotagged.sorted.bam -@ 32 | grep -E '^@|HP:i:1' | samtools view -@ 32 -b -o 20250114_H3K27me3_access_realtime_pass.trim.align.mapq60.hp1.sorted.bam -
samtools view -h ./haplotyped/haplotagged.chr1.sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam -@ 32 | grep -E '^@|HP:i:1' | samtools view -@ 32 -b -o  ./haplotyped/Haplotype_1.chr1.sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed.dna_r9.4.1_e8_sup@v3.3.5mCG.bam 
samtools index  ./haplotyped/Haplotype_1.chr1.sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed.dna_r9.4.1_e8_sup@v3.3.5mCG.bam 

# samtools view -h 20250114_H3K27me3_access_realtime_pass.trim.align.mapq60.haplotagged.sorted.bam -@ 32 | grep -E '^@|HP:i:2' | samtools view -@ 32 -b -o 20250114_H3K27me3_access_realtime_pass.trim.align.mapq60.hp2.sorted.bam -
samtools view -h ./haplotyped/haplotagged.chr1.sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam -@ 32 | grep -E '^@|HP:i:2' | samtools view -@ 32 -b -o  ./haplotyped/Haplotype_2.chr1.sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed.dna_r9.4.1_e8_sup@v3.3.5mCG.bam 
samtools index  ./haplotyped/Haplotype_2.chr1.sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed.dna_r9.4.1_e8_sup@v3.3.5mCG.bam 

# Create a VCF File for our BAM::
bcftools mpileup -Ou -f /home/michalula/data/ref_genomes/t2t_v2_0/chm13v2.0.fa  \
 chr1_206560169_206614236.sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam \
 | bcftools call -mv -Oz -o  \
 vcf.chr1_206560169_206614236.sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.vcf.gz

# Use own generated VCF File:

bcftools index  vcf.chr1_206560169_206614236.sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.vcf.gz

bcftools view vcf.chr1_206560169_206614236.sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.vcf.gz \
    -r chr1:206560169-206614236

#Phase using own VCF:
whatshap haplotag \
    --reference /home/michalula/data/ref_genomes/t2t_v2_0/chm13v2.0.fa \
    vcf.chr1_206560169_206614236.sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.vcf.gz \
    sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam \
    -o ./haplotyped/haplotagged.own_chr1_206560169_206614236.sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam
