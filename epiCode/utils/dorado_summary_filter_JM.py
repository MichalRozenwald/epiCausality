import argparse
from pathlib import Path
import subprocess
import csv

def main(
        dorado_summary: Path,
        input_bam: Path,
        min_mean_qscore: float,
        threads: int
    ) -> None:
    print("Finding reads to filter...")
    outgroup_file = Path("TEMP_outgroup.csv")
    with dorado_summary.open("r", newline="") as in_fp, outgroup_file.open("w") as out_fp:
        reader = csv.DictReader(in_fp, fieldnames=None, delimiter="\t")
        for row in reader:
            if float(row["mean_qscore_template"]) < min_mean_qscore and float(row["alignment_mapq"]) != 60:
                out_fp.write(row["read_id"])
                out_fp.write("\n")
    print("Filtering reads...")
    filtered_bam = Path("TEMP_filtered.bam")
    subprocess.run(["samtools", "view", "-b", "-N" f"^{outgroup_file}", "-o", filtered_bam, f"-@{threads}", input_bam])
    print("Sorting and indexing output...")
    output_bam = input_bam.with_suffix(f".filtered.q{min_mean_qscore}.bam")
    subprocess.run(["samtools", "sort", "-o", output_bam, f"-@{threads}", filtered_bam])
    subprocess.run(["samtools", "index", f"-@{threads}", output_bam])
    print("Cleaning up...")
    outgroup_file.unlink()
    filtered_bam.unlink()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="filter a dorado-generated BAM file based on the provided `dorado summary` file")
    parser.add_argument("input_bam", type=Path)
    parser.add_argument("dorado_summary", type=Path),
    parser.add_argument("-q", "--min_mean_qscore", type=float, default=10)
    parser.add_argument("-t", "--threads", type=int, default=1)

    cmd_args = parser.parse_args()
    main(**vars(cmd_args))
