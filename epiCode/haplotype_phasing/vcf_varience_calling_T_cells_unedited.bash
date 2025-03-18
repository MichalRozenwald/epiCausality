# # Generate a VCF from a BAM File

bcftools mpileup -Ou -f reference.fasta sorted_Tcell.bam | bcftools call -mv -Oz -o Tcell_variants.vcf.gz

bcftools mpileup -Ou -f /home/michalula/data/ref_genomes/t2t_v2_0/chm13v2.0.fa  \
 chr1_206560169_206614236.sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam \
 | bcftools call -mv -Oz -o  \
 vcf.chr1_206560169_206614236.sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.vcf.gz

bcftools index  vcf.chr1_206560169_206614236.sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.vcf.gz

bcftools view vcf.chr1_206560169_206614236.sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.vcf.gz \
    -r chr1:206560169-206614236


 Step 2: Verify Reference Genome Consistency
# If the VCF file was generated from a different reference genome than the BAM file, WhatsHap will not find matching positions.

bcftools view -h vcf.chr1_206560169_206614236.sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.vcf.gz | grep "##reference"
samtools view -H sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam | grep "@SQ"
bcftools query -f '%CHROM\n' vcf.chr1_206560169_206614236.sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.vcf.gz | sort | uniq

# 1️⃣ Check if the VCF Contains Variants in the Selected Region
bcftools view -r chr1:206560169-206614236 vcf.chr1_206560169_206614236.sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.vcf.gz | head -20

# 🚨 Issue: No Variants Found in the Selected Region
# Your bcftools view -r chr1:206560169-206614236 command only returned the VCF header and no variant records.
# This means your VCF file does not contain any variants in the region chr1:206560169-206614236.

bcftools mpileup -Ou -f /home/michalula/data/ref_genomes/t2t_v2_0/chm13v2.0.fa \
    -r chr1:206560169-206614236 \
    sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam | \
bcftools call -mv -Oz -o new_variants.chr1_206560169_206614236.vcf.gz

bcftools index new_variants.chr1_206560169_206614236.vcf.gz

bcftools view -r chr1:206560169-206614236 new_variants.chr1_206560169_206614236.vcf.gz | head -20

whatshap haplotag \
    --reference /home/michalula/data/ref_genomes/t2t_v2_0/chm13v2.0.fa \
    new_variants.chr1_206560169_206614236.vcf.gz \
    sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam \
    -o ./haplotyped/haplotagged.chr1_206560169_206614236.bam

# No variants detected


# Re-run variant calling, making sure that variant sites are detected.
# 1️⃣ Generate a Raw Variant Call File
# Instead of only using mpileup, add call -mv: 
bcftools mpileup -a DP,AD -Ou -f /home/michalula/data/ref_genomes/t2t_v2_0/chm13v2.0.fa \
    -r chr1:206560169-206614236 \
    sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam | \
bcftools call -mv -Oz -o forced_variants.chr1_206560169_206614236.vcf.gz

bcftools index forced_variants.chr1_206560169_206614236.vcf.gz

bcftools view -r chr1:206560169-206614236 forced_variants.chr1_206560169_206614236.vcf.gz | head -20

# If reads exist but no variants are called, try forcing bcftools call to output all sites, including reference sites (-A flag):
bcftools mpileup -a DP,AD -Ou -f /home/michalula/data/ref_genomes/t2t_v2_0/chm13v2.0.fa \
    -r chr1:206560169-206614236 \
    sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam | \
bcftools call -mv -A -Oz -o forced_variants_with_refs.chr1_206560169_206614236.vcf.gz

# Then index the new VCF:
bcftools index forced_variants_with_refs.chr1_206560169_206614236.vcf.gz
# Check if variants appear:
bcftools view -r chr1:206560169-206614236 forced_variants_with_refs.chr1_206560169_206614236.vcf.gz | head -20
 # Still nothing...
 samtools depth -r chr1:206560169-206614236 sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam
samtools view -q 30 -c sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam chr1:206560169-206614236


################### FROM HERE: ################
# If reads exist but no variants are called, try forcing bcftools call to output all sites, including reference sites (-A flag):
bcftools mpileup -a DP,AD -Ou -f /home/michalula/data/ref_genomes/t2t_v2_0/up_chm13v2.0.fasta \
    -r chr1:206560169-206614236 \
    /home/michalula/data/cas9_nanopore/data/20241226_MR_nCATs_TcellsPrES_unedit_P2R9/passed_fast5/5mCG/to_t2t_v2_0/sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam | \
bcftools call -mv -A -Oz -o /home/michalula/data/cas9_nanopore/data/20241226_MR_nCATs_TcellsPrES_unedit_P2R9/passed_fast5/5mCG/to_t2t_v2_0/haplotyped/uppercaseRef_forced_variants_with_refs.chr1_206560169_206614236.vcf.gz

bcftools index /home/michalula/data/cas9_nanopore/data/20241226_MR_nCATs_TcellsPrES_unedit_P2R9/passed_fast5/5mCG/to_t2t_v2_0/haplotyped/uppercaseRef_forced_variants_with_refs.chr1_206560169_206614236.vcf.gz
bcftools view -r chr1:206560169-206614236 /home/michalula/data/cas9_nanopore/data/20241226_MR_nCATs_TcellsPrES_unedit_P2R9/passed_fast5/5mCG/to_t2t_v2_0/haplotyped/uppercaseRef_forced_variants_with_refs.chr1_206560169_206614236.vcf.gz

bcftools view -r chr1:206560169-206614236 /home/michalula/data/cas9_nanopore/data/20241226_MR_nCATs_TcellsPrES_unedit_P2R9/passed_fast5/5mCG/to_t2t_v2_0/haplotyped/uppercaseRef_forced_variants_with_refs.chr1_206560169_206614236.vcf.gz | head -20

bcftools view -r chr1:206560169-206614236 /home/michalula/data/cas9_nanopore/data/20241226_MR_nCATs_TcellsPrES_unedit_P2R9/passed_fast5/5mCG/to_t2t_v2_0/haplotyped/uppercaseRef_forced_variants_with_refs.chr1_206560169_206614236.vcf.gz > /home/michalula/data/cas9_nanopore/data/20241226_MR_nCATs_TcellsPrES_unedit_P2R9/passed_fast5/5mCG/to_t2t_v2_0/haplotyped/uppercaseRef_forced_variants_with_refs.chr1_206560169_206614236.vcf
# Look at a broader set: chr1:206265172-206698324 - to check if the haplotyping work at the brigger region

bcftools index /home/michalula/data/cas9_nanopore/data/20241226_MR_nCATs_TcellsPrES_unedit_P2R9/passed_fast5/5mCG/to_t2t_v2_0/haplotyped/uppercaseRef_forced_variants_with_refs.chr1_206560169_206614236.vcf.gz

whatshap haplotag \
    --reference /home/michalula/data/ref_genomes/t2t_v2_0/up_chm13v2.0.fasta \
    /home/michalula/data/cas9_nanopore/data/20241226_MR_nCATs_TcellsPrES_unedit_P2R9/passed_fast5/5mCG/to_t2t_v2_0/haplotyped/uppercaseRef_forced_variants_with_refs.chr1_206560169_206614236.vcf.gz \
    /home/michalula/data/cas9_nanopore/data/20241226_MR_nCATs_TcellsPrES_unedit_P2R9/passed_fast5/5mCG/to_t2t_v2_0/sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam \
    -o /home/michalula/data/cas9_nanopore/data/20241226_MR_nCATs_TcellsPrES_unedit_P2R9/passed_fast5/5mCG/to_t2t_v2_0/haplotyped/haplotagged.chr1_206560169_206614236.sort_chr1_sort_align_t2t_v2_0_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed_dna_r9_e8_supv3mCG.bam
