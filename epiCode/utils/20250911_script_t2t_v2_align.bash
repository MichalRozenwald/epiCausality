# Aligning reads to t2t v2.0
# date: 2025-07-31

0. choose folder 
# Data Path 
/mnt/faststorage/michalula/T_cells/20250908_Day_35_T_CRoff_nCATs/20250908_1810_MN31715_FBD90703_2c0c1ee5/pod5
# /mnt/faststorage/michalula/T_cells/Day_28_post_EP_07212025/CRISPR_OFF/20250721_Day28_CRoff_T_cells_joint/no_sample_id/20250722_0004_MN31715_FBD30913_92318e69/pod5/
# /mnt/faststorage/michalula/T_cells/Day_28_post_EP_07212025/NT/20250721_Day28_NT_T_cells_joint/no_sample_id/20250721_2353_MN36507_FBD30760_a93f7f57/pod5
# /mnt/faststorage/michalula/T_cells/Day_28_post_EP_07212025/NT/MR_NT_Day28_07212025/no_sample_id/20250721_2353_MN36507_FBD30760_a93f7f57/pod5/
# /var/lib/minknow/data/MR_rerun_same_flowcell_NT_Day28_07222025/no_sample_id/20250722_1343_MN36507_FBD30760_cf22e917/pod5/
# /mnt/faststorage/michalula/T_cells/Day_28_post_EP_07212025/NT/MR_NT_Day28_07212025/no_sample_id/20250721_2353_MN36507_FBD30760_a93f7f57/pod5
# /mnt/faststorage/michalula/T_cells/Day_28_post_EP_07212025/CRISPR_OFF/MR_T_cells_CRoff_Day28_07212025/no_sample_id/20250722_0004_MN31715_FBD30913_92318e69/pod5
# cd /home/michalula/data/cas9_nanopore/data/20250721_nCATs_Tcells_UNEDITED_Day28/run_day_1and2

cd /home/michalula/data/cas9_nanopore/data/20250908_nCATs_T_CRoff_Day_35

# 2. basecall:  

# $ /path/to/dorado-x.y.z-linux-x64/bin/dorado basecaller hac pod5s/ > calls.bam


mkdir ./5mCG
  dorado basecaller \
   /home/michalula/software/dna_r9.4.1_e8_sup@v3.3	  \
   /mnt/faststorage/michalula/T_cells/20250908_Day_35_T_CRoff_nCATs/20250908_1810_MN31715_FBD90703_2c0c1ee5/pod5 \
    --modified-bases 5mCG \
    > ./5mCG/20250908_Day35_CROFF_Tcells_2Libraries_Minion_R9.dna_r9.4.1_e8_sup@v3.3.5mCG.bam   
 
 
dorado trim \
   --sequencing-kit SQK-CS9109 \
  ./5mCG/20250908_Day35_CROFF_Tcells_2Libraries_Minion_R9.dna_r9.4.1_e8_sup@v3.3.5mCG.bam     \
  >  ./5mCG/trim_20250908_Day35_CROFF_Tcells_2Libraries_Minion_R9.dna_r9.4.1_e8_sup@v3.3.5mCG.bam       


mkdir ./5mCG/to_t2t_v2_0
dorado aligner /home/michalula/data/ref_genomes/t2t_v2_0/chm13v2.0.fa \
 ./5mCG/trim_20250908_Day35_CROFF_Tcells_2Libraries_Minion_R9.dna_r9.4.1_e8_sup@v3.3.5mCG.bam      \
 > ./5mCG/to_t2t_v2_0/align_t2t_v2_0_trim_20250908_Day35_CROFF_Tcells_2Libraries_Minion_R9.dna_r9.4.1_e8_sup@v3.3.5mCG.bam   
 
 
samtools sort \
  -o ./5mCG/to_t2t_v2_0/sort_align_t2t_v2_0_trim_20250908_Day35_CROFF_Tcells_2Libraries_Minion_R9.dna_r9.4.1_e8_sup@v3.3.5mCG.bam    \
   ./5mCG/to_t2t_v2_0/align_t2t_v2_0_trim_20250908_Day35_CROFF_Tcells_2Libraries_Minion_R9.dna_r9.4.1_e8_sup@v3.3.5mCG.bam   

 
samtools index \
  ./5mCG/to_t2t_v2_0/sort_align_t2t_v2_0_trim_20250908_Day35_CROFF_Tcells_2Libraries_Minion_R9.dna_r9.4.1_e8_sup@v3.3.5mCG.bam
  

dorado summary \
  ./5mCG/to_t2t_v2_0/sort_align_t2t_v2_0_trim_20250908_Day35_CROFF_Tcells_2Libraries_Minion_R9.dna_r9.4.1_e8_sup@v3.3.5mCG.bam \
  > ./5mCG/to_t2t_v2_0/summary_sort_align_t2t_v2_0_trim_20250908_Day35_CROFF_Tcells_2Libraries_Minion_R9.dna_r9.4.1_e8_sup@v3.3.5mCG.tsv
  

grep -o 'chr' \
  ./5mCG/to_t2t_v2_0/summary_sort_align_t2t_v2_0_trim_20250908_Day35_CROFF_Tcells_2Libraries_Minion_R9.dna_r9.4.1_e8_sup@v3.3.5mCG.tsv \
   | wc -l

#  Day 35: 1,165,957 reads  - overlapping chromosome 

# Day 28: 10713 reads   - overlapping chromosome 1 

# 11964
# 8821
# Day 1: 5399 total overlapping chromosome 1 reads

# >> 226441

# # ZFP OFF:
# >> 850 454 

# # >> 55376


# No modification calls:
# >> 55376 

# When basecall on the fly:
# >> 23981


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

cd /home/michalula/data/cas9_nanopore/data/20241226_MR_nCATs_TPrES_CROFF90prD6_P2R9/pod5_converted_basecall/5mCG/to_t2t_v2_0
samtools view sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TPrES_CROFF90prD6_P2R9.dna_r9.4.1_e8_sup@v3.3.5mCG.bam | awk '{sum += length($10)} END {print sum}'
> 2.06885e+10
