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


bcftools index forced_variants.chr1_206560169_206614236.vcf.gz

bcftools view -r chr1:206560169-206614236 forced_variants.chr1_206560169_206614236.vcf.gz | head -20


top - 05:42:23 up 19 days, 10:43,  1 user,  load average: 1.17, 1.02, 0.77
Tasks: 442 total,   2 running, 440 sleeping,   0 stopped,   0 zombie
%Cpu(s):  4.2 us,  0.1 sy,  0.0 ni, 95.7 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
MiB Mem : 128554.4 total,  47447.2 free,   2883.7 used,  78223.4 buff/cache
MiB Swap:   2048.0 total,   1036.6 free,   1011.4 used. 124437.2 avail Mem 

    PID USER      PR  NI    VIRT    RES    SHR S  %CPU  %MEM     TIME+ COMMAND                                                        
1593328 michalu+  20   0  362344 355044   5184 R 100.0   0.3   9:33.73 bcftools                                                       
 421741 minknow   35  15 2374072   7216   5172 S   1.0   0.0 490:11.37 control_main                                                   
 421742 minknow   35  15 6952512   6612   4812 S   0.7   0.0   1433:39 control_main                                                   
    955 minknow   20   0 9962
top - 05:42:23 up 19 days, 10:43,  1 user,  load average: 1.17, 1.02, 0.77
Tasks: 442 total,   2 running, 440 sleeping,   0 stopped,   0 zombie
%Cpu(s):  4.2 us,  0.1 sy,  0.0 ni, 95.7 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
MiB Mem : 128554.4 total,  47447.2 free,   2883.7 used,  78223.4 buff/cache
MiB Swap:   2048.0 total,   1036.6 free,   1011.4 used. 124437.2 avail Mem 

    PID USER      PR  NI    VIRT    RES    SHR S  %CPU  %MEM     TIME+ COMMAND                                                        
1593328 michalu+  20   0  362344 355044   5184 R 100.0   0.3   9:33.73 bcftools                                                       
 421741 minknow   35  15 2374072   7216   5172 S   1.0   0.0 490:11.37 control_main                                                   
 421742 minknow   35  15 6952512   6612   4812 S   0.7   0.0   1433:39 control_main                                                   
    955 minknow   20   0 9962

    
# Look at a broader set: chr1:206265172-206698324