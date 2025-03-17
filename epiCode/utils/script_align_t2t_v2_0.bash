# Aligning reads to t2t v2.0
# date: 2025-03-16

cd /home/michalula/data/cas9_nanopore/data/20241226_nCATS_K562_ZFPOFFpostSort_HIGH/pod5_converted_basecall/5mCG/to_t2t_v2_0
dorado aligner /home/michalula/data/ref_genomes/t2t_v2_0/chm13v2.0.fa \
 /home/michalula/data/cas9_nanopore/data/20241226_nCATS_K562_ZFPOFFpostSort_HIGH/pod5_converted_basecall/5mCG/trim_20241226_nCATS_K562_ZFPOFFpostSort_HIGH_R9min_converted_fast5.dna_r9.4.1_e8_sup@v3.3.5mCG.bam    \
 >./align_t2t_v2_0_trim_20241226_nCATS_K562_ZFPOFFpostSort_HIGH_R9min_converted_fast5.dna_r9.4.1_e8_sup@v3.3.5mCG.bam   

samtools sort \
  -o ./sort_align_t2t_v2_0_trim_20241226_nCATS_K562_ZFPOFFpostSort_HIGH_R9min_converted_fast5.dna_r9.4.1_e8_sup@v3.3.5mCG.bam  \
   ./align_t2t_v2_0_trim_20241226_nCATS_K562_ZFPOFFpostSort_HIGH_R9min_converted_fast5.dna_r9.4.1_e8_sup@v3.3.5mCG.bam   

samtools index \
  ./sort_align_t2t_v2_0_trim_20241226_nCATS_K562_ZFPOFFpostSort_HIGH_R9min_converted_fast5.dna_r9.4.1_e8_sup@v3.3.5mCG.bam
  
dorado summary \
  ./sort_align_t2t_v2_0_trim_20241226_nCATS_K562_ZFPOFFpostSort_HIGH_R9min_converted_fast5.dna_r9.4.1_e8_sup@v3.3.5mCG.bam \
  > ./summary_sort_align_t2t_v2_0_trim_20241226_nCATS_K562_ZFPOFFpostSort_HIGH_R9min_converted_fast5.dna_r9.4.1_e8_sup@v3.3.5mCG.tsv
  
grep -o 'chr' \
  ./summary_sort_align_t2t_v2_0_trim_20241226_nCATS_K562_ZFPOFFpostSort_HIGH_R9min_converted_fast5.dna_r9.4.1_e8_sup@v3.3.5mCG.tsv \
   | wc -l
# ZFP off high aligned to T2T version 2:
>> 226 459
# ZFP off LOW aligned to T2T version 2:
>>850 533

# T cells post CRISPROff:
# >> 2 706 166 >>> WOOOW 2 MLN READS

# Unedited T cells Day too afret basecalling quality filtering 
# >> 775 506

# To T2T v1.1 alignment:
# >> 226 441

# # ZFP OFF:
# >> 850 454 

# # >> 55376

# No modification calls:
# >> 55376 

# When basecall on the fly:
# >> 23981


# # Calculate in a bam file how many bases there are across all reads in that file
# Method 1: Using samtools stats
# Run the following command: 
# samtools stats your_file.bam | grep "bases mapped (cigar)"
# This extracts the total number of bases mapped according to the CIGAR string.
# If you want only the number: 
# samtools stats your_file.bam | grep "bases mapped (cigar)" | cut -f3
samtools stats sort_align_t2t_v2_0_trim_20241226_nCATS_K562_ZFPOFFpostSort_HIGH_R9min_converted_fast5.dna_r9.4.1_e8_sup@v3.3.5mCG.bam  | grep "bases mapped (cigar)" | cut -f3
> 2328983440
> 7.370183e-02
samtools stats sort_align_t2t_v2_0_trim_20241226_nCATS_K562_ZFPOFFpostSort_HIGH_R9min_converted_fast5.dna_r9.4.1_e8_sup@v3.3.5mCG.bam | grep "bases mapped (cigar)"

samtools stats sort_align_t2t_v2_0_trim_20241226_nCATS_K562_ZFPOFFpostSort_HIGH_R9min_converted_fast5.dna_r9.4.1_e8_sup@v3.3.5mCG.bam | grep "bases mapped (cigar)" | awk '{print $NF}'

  SN      bases mapped (cigar):   2328983440      # more accurate
  SN      error rate:     7.370183e-02    # mismatches / bases mapped (cigar)

#To filter only the bases mapped value correctly, modify the command:
# samtools stats sort_align_t2t_v2_0_trim_20241226_nCATS_K562_ZFPOFFpostSort_HIGH_R9min_converted_fast5.dna_r9.4.1_e8_sup@v3.3.5mCG.bam | grep "bases mapped (cigar):" | awk '{print $NF}'

samtools stats sort_align_t2t_v2_0_trim_20241226_nCATS_K562_ZFPOFFpostSort_HIGH_R9min_converted_fast5.dna_r9.4.1_e8_sup@v3.3.5mCG.bam | grep "bases mapped (cigar):"


# Method 2: Using samtools view and awk
# samtools view your_file.bam | awk '{sum += length($10)} END {print sum}'
samtools view sort_align_t2t_v2_0_trim_20241226_nCATS_K562_ZFPOFFpostSort_HIGH_R9min_converted_fast5.dna_r9.4.1_e8_sup@v3.3.5mCG.bam | awk '{sum += length($10)} END {print sum}'
# 2.53461e+09

Alternative: Counting Only Mapped Reads
# To exclude unmapped reads, filter using the SAM flag (0x4 for unmapped reads):
samtools view -F 4 sort_align_t2t_v2_0_trim_20241226_nCATS_K562_ZFPOFFpostSort_HIGH_R9min_converted_fast5.dna_r9.4.1_e8_sup@v3.3.5mCG.bam | awk '{sum += length($10)} END {print sum}'
2.46353e+09

