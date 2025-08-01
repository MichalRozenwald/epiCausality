# Setect only broader ROI 
samtools view -b /home/michalula/data/cas9_nanopore/data/20241226_MR_nCATs_TPrES_CROFF90prD6_P2R9/pod5_converted_basecall/5mCG/to_t2t_v2_0/sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TPrES_CROFF90prD6_P2R9.dna_r9.4.1_e8_sup@v3.3.5mCG.bam \
    chr1:206560169-206614236 > \
    /home/michalula/data/cas9_nanopore/data/20241226_MR_nCATs_TPrES_CROFF90prD6_P2R9/pod5_converted_basecall/5mCG/to_t2t_v2_0/chr1_206560169_206614236.sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TPrES_CROFF90prD6_P2R9.dna_r9.4.1_e8_sup@v3.3.5mCG.bam 

ls -lh /home/michalula/data/cas9_nanopore/data/20241226_MR_nCATs_TPrES_CROFF90prD6_P2R9/pod5_converted_basecall/5mCG/to_t2t_v2_0
cd /home/michalula/data/cas9_nanopore/data/20241226_MR_nCATs_TPrES_CROFF90prD6_P2R9/pod5_converted_basecall/5mCG/to_t2t_v2_0
# calculate how many reads are saved 

samtools view -c chr1_206560169_206614236.sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TPrES_CROFF90prD6_P2R9.dna_r9.4.1_e8_sup@v3.3.5mCG.bam 
# 10999
samtools index chr1_206560169_206614236.sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TPrES_CROFF90prD6_P2R9.dna_r9.4.1_e8_sup@v3.3.5mCG.bam 


# bcftools mpileup -Ou -f reference.fasta sorted_Tcell.bam | bcftools call -mv -Oz -o Tcell_variants.vcf.gz

bcftools mpileup -Ou -f /home/michalula/data/ref_genomes/t2t_v2_0/chm13v2.0.fa  \
 chr1_206560169_206614236.sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TPrES_CROFF90prD6_P2R9.dna_r9.4.1_e8_sup@v3.3.5mCG.bam  \
 | bcftools call -mv -Oz -o  \
 vcf.chr1_206560169_206614236.chr1_206560169_206614236.sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TPrES_CROFF90prD6_P2R9.dna_r9.4.1_e8_sup@v3.3.5mCG.vcf.gz


bcftools index  vcf.chr1_206560169_206614236.chr1_206560169_206614236.sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TPrES_CROFF90prD6_P2R9.dna_r9.4.1_e8_sup@v3.3.5mCG.vcf.gz

bcftools view -r chr1:206560169-206614236 vcf.chr1_206560169_206614236.chr1_206560169_206614236.sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TPrES_CROFF90prD6_P2R9.dna_r9.4.1_e8_sup@v3.3.5mCG.vcf.gz | head -20

# top - 05:42:23 up 19 days, 10:43,  1 user,  load average: 1.17, 1.02, 0.77
# Tasks: 442 total,   2 running, 440 sleeping,   0 stopped,   0 zombie
# %Cpu(s):  4.2 us,  0.1 sy,  0.0 ni, 95.7 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
# MiB Mem : 128554.4 total,  47447.2 free,   2883.7 used,  78223.4 buff/cache
# MiB Swap:   2048.0 total,   1036.6 free,   1011.4 used. 124437.2 avail Mem 

#     PID USER      PR  NI    VIRT    RES    SHR S  %CPU  %MEM     TIME+ COMMAND                                                        
# 1593328 michalu+  20   0  362344 355044   5184 R 100.0   0.3   9:33.73 bcftools                                                       
#  421741 minknow   35  15 2374072   7216   5172 S   1.0   0.0 490:11.37 control_main                                                   
#  421742 minknow   35  15 6952512   6612   4812 S   0.7   0.0   1433:39 control_main                                                   
#     955 minknow   20   0 9962
# top - 05:42:23 up 19 days, 10:43,  1 user,  load average: 1.17, 1.02, 0.77
# Tasks: 442 total,   2 running, 440 sleeping,   0 stopped,   0 zombie
# %Cpu(s):  4.2 us,  0.1 sy,  0.0 ni, 95.7 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
# MiB Mem : 128554.4 total,  47447.2 free,   2883.7 used,  78223.4 buff/cache
# MiB Swap:   2048.0 total,   1036.6 free,   1011.4 used. 124437.2 avail Mem 

#     PID USER      PR  NI    VIRT    RES    SHR S  %CPU  %MEM     TIME+ COMMAND                                                        
# 1593328 michalu+  20   0  362344 355044   5184 R 100.0   0.3   9:33.73 bcftools                                                       
#  421741 minknow   35  15 2374072   7216   5172 S   1.0   0.0 490:11.37 control_main                                                   
#  421742 minknow   35  15 6952512   6612   4812 S   0.7   0.0   1433:39 control_main                                                   
#     955 minknow   20   0 9962


# Still no vcf file generated correc
samtools mpileup -r chr1:206560169-206614236 \
    -f /home/michalula/data/ref_genomes/t2t_v2_0/chm13v2.0.fa \
    sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TPrES_CROFF90prD6_P2R9.dna_r9.4.1_e8_sup@v3.3.5mCG.bam | head -20

samtools mpileup -r chr1:206560169-206614236 \
>     -f /home/michalula/data/ref_genomes/t2t_v2_0/chm13v2.0.fa \
>     sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TPrES_CROFF90prD6_P2R9.dna_r9.4.1_e8_sup@v3.3.5mCG.bam | head -20
# [mpileup] 1 samples in 1 input files

# chr1    206560169       t       19      ,,,,,,.,,,,,,,,,,,,     @>>>@>>B?A55@6??:>=
# chr1    206560170       g       19      ,,,,,,.,,,,,,,,,,,,     ?==D?>=A>B55B7>=:>?
# chr1    206560171       t       23      ,,,,,,,,.,,,,,,,,,,,,,, >8<=3=7===<1:./85:0;9>8
# chr1    206560172       g       21      ,,,,,,,,.,,,,,,,,,,,,   >7<<387=<=;1:78:/=8?=
# chr1    206560173       a       20      ,,,,,,,,.,,,,,,,,,,,    =7<<290<=<;19779<8>;
# chr1    206560174       g       19      ,,,,,,,.,,,,,,,,,,,     ?6=>29=<=;1777:<9=;
# chr1    206560175       c       18      ,,*,,,,.,,,,,,,,,,      ;=/>29=:=;1547:;52
# chr1    206560176       t       22      ,,,,,,,,.,*,,,,,,,,,,,  ;=/A.7:@:?0=4<76=@>072
# chr1    206560177       t       25      ,,,,,,,,.,,,,,,,,,,,,,,,,       <@0?.;=A;A0>5>8007/AA?872
# chr1    206560178       g       23      ,,,,,,,,.,,,,,,,,,,,,,, .>/<.9=A;@3?696204B9972
# chr1    206560179       t       22      ,,,,,,,.,,,,,,,,,,,,,,  >/;/:=A;@4=676404@<872
# chr1    206560180       t       28      ,,,,,,,,,.,,,,,,,,,,,,,,,,,,    7B.:/2?D??B7C@89>875?4?:6A<0
# chr1    206560181       g       26      ,,,,,,,,.,,,,,,,,,,,,,,,,,      6>=.6?BC>A7CC88=575>5?>9A@
# chr1    206560182       g       27      ,,,,,,,,.,,,,,,,,,,,,,,,,,,     B<@?=D0D?@5?@:9>5;677@<<D6C
# chr1    206560183       a       23      ,,,,,,,,.,,,,,,,,,,,,,, >C>=:>.>=?F?0D3@0>7;A1A
# chr1    206560184       t       21      ,,,,,,,.,,,,,,,,,,,,,   9?9<9<9659:>29/7589.8
# chr1    206560185       c       20      ,,,,,,,.,,,,,,,,,,,,    9=2.7;1848:=28/74888
# chr1    206560186       a       17      ,,,,.,,,,,,,-1t,,,,,    7=2..38:727/72778
# chr1    206560187       t       16      ,,,,,,,,,,*,,,,,        7<2.38:726/61378
# chr1    206560188       c       16      ,*,,,,,,,,,,,,,,        82=238:18/7/7037



✅ Reads Contain Variants – Variant Calling Should Work
Your samtools mpileup output confirms that there are sequence mismatches (potential variants) in the region chr1:206560169-206614236.

Now, we need to make sure bcftools properly detects these variants and includes them in the VCF.


🛠 Fix: Correct Ploidy Setting
1️⃣ Specify Ploidy Using GRCh38 (for T2T-CHM13)
Since your data is aligned to T2T-CHM13 v2.0, which is derived from GRCh38, use:

bcftools mpileup -a DP,AD -Q 0 -q 0 -C 50 -Ou -f /home/michalula/data/ref_genomes/t2t_v2_0/chm13v2.0.fa \
    -r chr1:206560169-206614236 \
    sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TPrES_CROFF90prD6_P2R9.dna_r9.4.1_e8_sup@v3.3.5mCG.bam | \
bcftools call -mv -A --ploidy GRCh38 --keep-alts -Oz -o fixed_variants.chr1_206560169_206614236.vcf.gz

# 🔹 Explanation:
# --ploidy GRCh38 → Uses human diploid settings from GRCh38, which should match T2T-CHM13.
# --keep-alts → Retains alternative alleles.
# -A → Outputs all sites, including reference sites.

# 2️⃣ Verify That Variants Were Called
# Once complete, check if the new VCF contains variants:

bcftools index fixed_variants.chr1_206560169_206614236.vcf.gz
bcftools view -r chr1:206560169-206614236 fixed_variants.chr1_206560169_206614236.vcf.gz | head -20

# If variants appear, we have successfully generated a valid VCF! 🎉
# Nopee...

# 3️⃣ If Still Failing, Manually Set Ploidy
# If GRCh38 does not work, force bcftools to assume diploidy:
bcftools call -mv -A --ploidy-file <(echo -e "chr1\t206560169\t206614236\t*\t2") --keep-alts -Oz \
    -o fixed_variants.chr1_206560169_206614236.vcf.gz

🚀 Final Steps: Run WhatsHap
If variants are now detected, rerun WhatsHap:
whatshap haplotag \
    --reference /home/michalula/data/ref_genomes/t2t_v2_0/chm13v2.0.fa \
    fixed_variants.chr1_206560169_206614236.vcf.gz \
    sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TPrES_CROFF90prD6_P2R9.dna_r9.4.1_e8_sup@v3.3.5mCG.bam \
    -o ./haplotyped/haplotagged.chr1_206560169_206614236.bam

3️⃣ Force bcftools to Call Variants in All Sites
If depth and reference are correct, force variant calling even for reference-matching sites:

bcftools mpileup -a DP,AD -Q 0 -q 0 -C 50 --output-tags DP,AD -Ou -f /home/michalula/data/ref_genomes/t2t_v2_0/chm13v2.0.fa \
    -r chr1:206560169-206614236 \
    /home/michalula/data/cas9_nanopore/data/20241226_MR_nCATs_TPrES_CROFF90prD6_P2R9/pod5_converted_basecall/5mCG/to_t2t_v2_0/chr1_206560169_206614236.sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TPrES_CROFF90prD6_P2R9.dna_r9.4.1_e8_sup@v3.3.5mCG.bam  | \
bcftools call -mv -A --ploidy GRCh38 --keep-alts --output-type z \
    -o forced_variants_with_refs.chr1_206560169_206614236.vcf.gz

# Look at a broader set: chr1:206265172-206698324
