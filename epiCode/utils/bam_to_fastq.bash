#!/bin/bash
# Script to convert all BAM files in the current directory to FASTQ files

conda activate dimelo_v2_modkit_parsing
cd /mnt/faststorage/michalula/T_cells/Day_28_post_EP_07212025/CRISPR_OFF/20250721_Day28_CRoff_T_cells_joint/no_sample_id/20250722_0004_MN31715_FBD30913_92318e69/fastq
fastq_dir="/mnt/faststorage/michalula/T_cells/Day_28_post_EP_07212025/CRISPR_OFF/20250721_Day28_CRoff_T_cells_joint/no_sample_id/20250722_0004_MN31715_FBD30913_92318e69/fastq"
bam_dir="/mnt/faststorage/michalula/T_cells/Day_28_post_EP_07212025/CRISPR_OFF/20250721_Day28_CRoff_T_cells_joint/no_sample_id/20250722_0004_MN31715_FBD30913_92318e69/bam_pass"

# Check for samtools
if ! command -v samtools &> /dev/null; then
    echo "samtools could not be found. Please install samtools."
    exit 1
fi

mkdir "${fastq_dir}/bam_to_fastq"


for bam in $bam_dir/*.bam; do
    [ -e "$bam" ] || continue
    base="${bam%.bam}"
    bam_filename=$(basename "$base")
    echo "Converting $bam to ${base}.fastq"
    samtools fastq "$bam" > "${fastq_dir}/bam_to_fastq/${bam_filename}.fastq"
done

# Merge all FASTQ files into a single file
merged_fastq="${fastq_dir}/merged_20250721_Day28_CRoff_T_cells_nCATs.fastq"

cat "${fastq_dir}"/bam_to_fastq/*.fastq > "$merged_fastq" 
echo "All FASTQ files have been merged into $merged_fastq"

