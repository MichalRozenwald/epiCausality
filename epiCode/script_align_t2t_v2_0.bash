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
