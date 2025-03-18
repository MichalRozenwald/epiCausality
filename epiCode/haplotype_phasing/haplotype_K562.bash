

# For K562 cells:
# Using previous data aligned to hg19:
whatshap polyphase /home/michalula/ref_genomes/K562/ENCFF960SSF.vcf   /home/michalula/data/cas9_nanopore/data/2024927_Cas9_R9_promethion_K562/to_hg19/sort_aligned_hg19_trim_bam_pass_merged.bam --ploidy p --reference /home/michalula/ref_genomes/hg19/GRCh37.p13.genome.fa  -o /home/michalula/data/cas9_nanopore/data/2024927_Cas9_R9_promethion_K562/to_hg19/hyplotyped/whatshap_output.vcf
whatshap polyphase /home/michalula/ref_genomes/K562/ENCFF960SSF.vcf \
  /home/michalula/data/cas9_nanopore/data/2024927_Cas9_R9_promethion_K562/to_hg19/sort_aligned_hg19_trim_bam_pass_merged.bam \
  --ploidy 2 \
  --reference /home/michalula/ref_genomes/hg19/GRCh37.p13.genome.fa \
  -o /home/michalula/data/cas9_nanopore/data/2024927_Cas9_R9_promethion_K562/to_hg19/hyplotyped/whatshap_output_K562.vcf


cd /home/michalula/data/cas9_nanopore/data/2024927_Cas9_R9_promethion_K562/to_hg19/
mkdir /home/michalula/data/cas9_nanopore/data/2024927_Cas9_R9_promethion_K562/to_hg19/haplotyped
whatshap haplotag \
    --reference /home/michalula/ref_genomes/hg19/GRCh37.p13.genome.fa \
    /home/michalula/ref_genomes/K562/ENCFF960SSF.vcf \
    /home/michalula/data/cas9_nanopore/data/2024927_Cas9_R9_promethion_K562/to_hg19/sort_aligned_hg19_trim_bam_pass_merged.bam \
    -o /home/michalula/data/cas9_nanopore/data/2024927_Cas9_R9_promethion_K562/to_hg19/haplotyped/haplotagged.sort_aligned_hg19_trim_bam_pass_merged.bam

# \
#     --ignore-read-groups \
#     --sample NA21116 \