# Haplotype phasing of T cells using whatshap
# Based on Oberon's code:
# This script is based on the script used to phase GM12878 cells in the epiCode repository
# The script is used to phase the T cells using the same pipeline as the GM12878 cells

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

