import pysam
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, Tuple, Dict, List
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

''' 
1. > "How are the quality scores and alignment scores calculated?"
    Quality scores in your code are the per-base Phred scores from the BAM file, 
    accessed via read.query_qualities. These scores represent the probability that 
    a base call is incorrect, with higher scores indicating higher confidence.

    Alignment scores in your heatmap are calculated for each base in the read that 
    aligns to a reference position within the specified region. The calculation 
    depends on the mode parameter:

    For each aligned base (where both rpos and refpos are not None), the code compares the read base (base = qseq[rpos].upper()) to the reference base (refb = ref_arr[col]).
    If the bases match and are A/C/G/T:
    'quality_signed': +q (positive Phred score)
    'quality': q (Phred score)
    'binary' or 'ternary': 1.0
    If the bases mismatch (including ambiguous bases):
    'quality_signed': -q (negative Phred score)
    'quality': 0.0
    'binary': -1.0
    'ternary': 0.0
    If there is an indel or soft clip (i.e., rpos or refpos is None), the score is left as np.nan.
    Summary:

    Quality scores: Directly from the BAM file, per-base Phred scores.
    Alignment scores: Calculated per reference position, based on 
    match/mismatch and the selected mode, using the quality score for matches/mismatches. 
    Indels/clips are set to NaN.

2. > "What are rpos and refpos?"
rpos and refpos are variables representing positions in the read and reference, 
respectively, as returned by the read.get_aligned_pairs() method from pysam.

rpos: The position (index) of a base in the read sequence (read.query_sequence).
refpos: The position (index) of the corresponding base in the reference genome.
Each tuple (rpos, refpos) describes how a base in the read aligns to a base in the reference.

If both are not None, the base is aligned (match or mismatch).
If rpos is None, it's a deletion in the read (relative to the reference).
If refpos is None, it's an insertion or soft clip in the read (relative to the reference).
'''



def _parse_region(region: Optional[str], bam: pysam.AlignmentFile) -> Tuple[Optional[str], Optional[int], Optional[int], Optional[int]]:
    """
    Parse 'chr:start-end' → (chr, start, end, length). If region is None, return Nones.
    """
    if region is None:
        return None, None, None, None
    chrom, rng = region.split(":")
    start, end = map(int, rng.replace(",", "").split("-"))
    if end <= start:
        raise ValueError("Region end must be > start.")
    # Ensure contig exists in BAM header
    if chrom not in bam.references:
        raise ValueError(f"Chromosome {chrom} not present in BAM header.")
    return chrom, start, end, end - start


def _fetch_reference_slice(reference_fasta: str, chrom: str, start: int, end: int) -> np.ndarray:
    """
    Fetch uppercase reference bases for chrom:start-end as a numpy array of dtype '<U1'.
    """
    fa = pysam.FastaFile(reference_fasta)
    seq = fa.fetch(chrom, start, end).upper()
    fa.close()
    return np.array(list(seq), dtype="<U1")


def get_read_info_by_name(bam_path, read_name):
    """
    Print information about a read in a BAM file given its query name.

    Parameters
    ----------
    bam_path : str or Path
        Path to the BAM file.
    read_name : str
        Query name (ID) of the read to search for.

    Output
    ------
    Prints the index number, read name, reference name, start/end coordinates,
    mapping quality, and strand orientation for the first matching read.

    Goal
    ----
    Quickly locate and display alignment details for a specific read by its name.
    """
    with pysam.AlignmentFile(bam_path, "rb") as bamfile:
        for i, read in enumerate(bamfile):
            if read.query_name == read_name:
                print(f"Read Index Number: {i}")
                print(f"Read name: {read.query_name}")
                print(f"Reference name: {bamfile.get_reference_name(read.reference_id)}")
                print(f"Start: {read.reference_start}")
                print(f"End: {read.reference_end}")
                print(f"Mapping quality:  {read.mapping_quality}")
                print(f"Is reverse: {read.is_reverse}")
                break


def get_read_info_by_index(bam_path, read_index):
    """
    Print information about a read in a BAM file given its index.

    Parameters
    ----------
    bam_path : str or Path
        Path to the BAM file.
    read_index : int
        Index (0-based) of the read to retrieve.

    Output
    ------
    Prints the read index, read name, start, and end coordinates for the specified read.

    Goal
    ----
    Quickly display alignment details for a read at a specific index in the BAM file.
    """
    with pysam.AlignmentFile(bam_path, "rb") as bamfile:
        for i, read in enumerate(bamfile):
            if i == read_index:
                print(f"Read Index Number: {i}")
                print(f"Read name: {read.query_name}")
                print(f"Reference name: {bamfile.get_reference_name(read.reference_id)}")
                print(f"Start: {read.reference_start}")
                print(f"End: {read.reference_end}")
                print(f"Mapping quality:  {read.mapping_quality}")
                print(f"Is reverse: {read.is_reverse}")
                break



def build_alignment_heatmap(
    bam_path: str | Path,
    region: str,
    reference_fasta: str,
    mode: str = "quality_signed",
    max_reads: Optional[int] = None,
    primary_only: bool = True,
    min_mapq: Optional[int] = None,
    include_unmapped: bool = False,
) -> Tuple[np.ndarray, Dict]:
    """
    Build a per-read x per-position matrix over a region.
    Each cell is a score for that read's base aligned to that reference position.

    Parameters
    ----------
    bam_path : str | Path
        Input BAM (indexed recommended).
    region : str
        Genomic window 'chr:start-end'. Keep it reasonably small (e.g., ≤100kb) to avoid huge matrices.
    reference_fasta : str
        Path to reference FASTA (indexed .fai). Used to compare matches/mismatches.
    mode : {'quality_signed','quality','binary','ternary'}
        Scoring scheme:
          - 'quality_signed' (default): +Q for match, -Q for mismatch, NaN for indels/clips (Q = per-base Phred).
          - 'quality'         : Q for match, 0 for mismatch, NaN for indels/clips.
          - 'binary'          : 1 for match, -1 for mismatch, NaN for indels/clips.
          - 'ternary'         : 1 for match, 0 for mismatch, NaN for indels/clips.
    max_reads : int | None
        Limit number of reads (first N encountered) to control memory; None = all reads in region.
    primary_only : bool
        If True, skip secondary/supplementary alignments.
    min_mapq : int | None
        If set, only include reads with mapping_quality >= min_mapq.
    include_unmapped : bool
        If True, keep unmapped reads as all-NaN rows (rarely useful).

    Returns
    -------
    matrix : np.ndarray
        Shape (num_reads, region_length); dtype float; NaN where no aligned base.
    meta : dict
        {'read_names': [...], 'region': region, 'mode': mode, 'chrom': str, 'start': int, 'end': int}
    """
    bam_path = Path(bam_path)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        chrom, start, end, region_len = _parse_region(region, bam)

        # Preload reference for this region
        ref_arr = _fetch_reference_slice(reference_fasta, chrom, start, end)

        rows: List[np.ndarray] = []
        read_names: List[str] = []

        # Iterate reads overlapping window
        count = 0
        for read in bam.fetch(chrom, start, end):
            if primary_only and (read.is_secondary or read.is_supplementary):
                continue
            if (not include_unmapped) and read.is_unmapped:
                continue
            if (min_mapq is not None) and (read.mapping_quality < min_mapq):
                continue

            # Initialize this read's row as NaNs
            row = np.full(region_len, np.nan, dtype=float)

            # Skip if no sequence/qualities (can happen)
            if read.query_sequence is None or read.query_qualities is None:
                rows.append(row)
                read_names.append(read.query_name)
                count += 1
                if max_reads and count >= max_reads:
                    break
                continue

            qseq = read.query_sequence
            qquals = read.query_qualities  # array of Phred scores (0..~60)

            # Map aligned pairs
            # (read_pos, ref_pos) with both not None => aligned base (M/=/X)
            for rpos, refpos in read.get_aligned_pairs(matches_only=False, with_seq=False):
                if refpos is None:
                    # insertion relative to reference (or soft clip) → no reference column
                    continue
                if rpos is None:
                    # deletion relative to reference; keep as NaN (or 0, but NaN is clearer)
                    continue

                # Only fill if within our window
                if start <= refpos < end:
                    col = refpos - start
                    base = qseq[rpos].upper()
                    refb = ref_arr[col]
                    q = float(qquals[rpos])  # per-base Phred

                    if base == refb and base in ("A", "C", "G", "T"):
                        # match
                        if mode == "quality_signed":
                            row[col] = +q
                        elif mode == "quality":
                            row[col] = q
                        elif mode == "binary":
                            row[col] = 1.0
                        elif mode == "ternary":
                            row[col] = 1.0
                        else:
                            raise ValueError(f"Unknown mode '{mode}'")
                    else:
                        # mismatch (including N or other IUPAC)
                        if mode == "quality_signed":
                            row[col] = -q
                        elif mode == "quality":
                            row[col] = 0.0
                        elif mode == "binary":
                            row[col] = -1.0
                        elif mode == "ternary":
                            row[col] = 0.0
                        else:
                            raise ValueError(f"Unknown mode '{mode}'")

            rows.append(row)
            read_names.append(read.query_name)
            count += 1
            if max_reads and count >= max_reads:
                break

    matrix = np.vstack(rows) if rows else np.empty((0, end - start), dtype=float)
    meta = {
        "read_names": read_names,
        "region": region,
        "mode": mode,
        "chrom": chrom,
        "start": start,
        "end": end,
    }
    return matrix, meta


def plot_alignment_heatmap(
    matrix: np.ndarray,
    meta: Dict,
    title: Optional[str] = None,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    show_colorbar: bool = True,
    show_clustered_heatmap: bool = True,
    show_read_names:  bool = False,
):
    """
    Plot the alignment heatmap with clean y-axis ticks at row centers.
    No changes to how rows are created—one row per record, just fixed y-axis.
    """
    if matrix.size == 0:
        raise ValueError("Empty matrix: nothing to plot.")

    nreads = matrix.shape[0]
    plt.figure(figsize=(10, max(10, nreads * 0.02)))
    im = plt.imshow(
        matrix,
        aspect="auto",
        interpolation="nearest",
        vmin=vmin,
        vmax=vmax,
    )

    if title is None:
        title = f"Alignment Heatmap [{meta.get('region','?')}] • mode={meta.get('mode','?')} \n Number of reads = {nreads}"
    plt.title(title)
    plt.xlabel("Reference position (offset in window)")
    plt.ylabel("Reads")

    # --- Fix: y-axis at row centers (0..N-1), no -0.5/1.5 edges ---
    ax = plt.gca()
    ax.set_yticks(np.arange(nreads))
    # If read names available and match row count, use them; otherwise index numbers
    read_names = meta.get("read_names")
    if show_read_names and isinstance(read_names, list) and len(read_names) == nreads:
        labels = [rn if len(rn) <= 30 else rn[:14] + "…" + rn[-14:] for rn in read_names]
        ax.set_yticklabels(labels, fontsize=8)
    else:
        ax.set_yticklabels([str(i) for i in range(nreads)], fontsize=8)
    # --------------------------------------------------------------
    if show_colorbar:
        plt.colorbar(im,cmap='plasma') # cmap='viridis') 
    plt.tight_layout()
    plt.show()

    # Create the clustermap and customize the colorbar
    if show_clustered_heatmap:
        if matrix.shape[0] >= 2:
            # sns.clustermap(G_revs_df.fillna(-1), col_cluster=False)
            g = sns.clustermap(pd.DataFrame(matrix).fillna(100), col_cluster=False, cmap='plasma')  # , cbar_kws={'label': "Y labet") # , 'orientation': 'vertical'
            g.ax_cbar.set_title(f"Clustered {title} \n Fill NA with +100 \n Number of reads = {nreads}", pad=16)
            # Display the plot
            plt.show()
        else:
            print("Not enough reads for Clustered Heatmap with " + str(matrix.shape[0])+ " reads.")


def plot_reads_quality_heatmap(
    bam_path: str | Path,
    region: str,
    ref_genome_path: str | Path,
    heatmap_title: str = "Read Quality Heatmap",
    mode: str = "quality_signed",
    max_reads: Optional[int] = None,
    primary_only: bool = True,
    min_mapq: Optional[int] = None,
    save_matrix: bool = False,
):
    """
    Convenience function to build and plot a quality heatmap for reads in a BAM over a region.

    # Example usage:
    region_str = "chr1:206586100-206586220"  # keep windows reasonably small
    original_bam_path = "/home/michalula/data/ont/t2t               
    bam = original_bam_path #  unedit_bam_path # removed_reads_bam_path # "/path/to/your.bam"
    ref_genome_path = Path('/home/michalula/data/ref_genomes/t2t_v2_0/up_chm13v2.0.fasta') # ref = "/path/to/reference.fa"  # must have .fai index
    plot_reads_quality_heatmap(bam, region_str, ref_genome_path, mode="quality_signed", max_reads=5000, primary_only=True, min_mapq=20, save_matrix=False)
    """
    matrix, meta = build_alignment_heatmap(
        bam_path=bam_path,
        region=region,
        reference_fasta=ref_genome_path,  # Not needed for quality heatmap
        mode=mode,  # Use "quality_signed" for signed or 'quality' mode for per-base Phred scores
        max_reads=max_reads,
        primary_only=primary_only,
        min_mapq=min_mapq,
        include_unmapped=False,
    )
    plot_alignment_heatmap(matrix, meta, heatmap_title, vmin=None, vmax=None)
    if save_matrix:
        np.save("alignment_heatmap_matrix.npy", matrix)
        


# import pysam
# import pandas as pd
# from pathlib import Path
def count_indels_and_mismatches(bam_path, ref_genome_path, region=None, output_csv=None):
    """
    Calculate the number of indels (insertions/deletions/soft clips), mismatches, and None values for each read in a BAM file.
    Saves the results into a table (CSV if output_csv is provided).

    Parameters
    ----------
    bam_path : str or Path
        Path to the BAM file.
    ref_genome_path : str or Path
        Path to the reference FASTA file (indexed).
    region : str, optional
        Genomic region 'chr:start-end' to restrict analysis.
    output_csv : str or Path, optional
        Path to save the output table as CSV.

    Output
    ------
    DataFrame with columns: ['read_name_str', 'num_indels', 'num_mismatches', 'num_nones']

    Goal
    ----
    For each read, count the number of indels, mismatches, and None values, and save the results.
    """
    bam_path = Path(bam_path)
    fa = pysam.FastaFile(str(ref_genome_path))

    read_lengths = []
    mapping_qualities = []
    avg_base_qualities = []

    mode_base_qualities = []
    median_base_qualities = []
    num_zeros_base_qualities = []
    less_10_base_qualities = []

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        if region:
            chrom, rng = region.split(":")
            start, end = map(int, rng.replace(",", "").split("-"))
            reads = bam.fetch(chrom, start, end)
            print(f"In function count_indels_and_mismatches: Processing region is {region}")
            print(f"In function count_indels_and_mismatches: Region length is {end - start}")
        else:
            reads = bam.fetch()


        results = []
        for read in reads:
            if read.is_unmapped or read.query_sequence is None:
                continue

            # Skip secondary/supplementary if you only want primary alignments
            if read.is_secondary or read.is_supplementary:
                continue

            read_lengths.append(read.query_length)
            mapping_qualities.append(read.mapping_quality)
            if read.query_qualities is not None:
                avg_base_qualities.append(np.mean(read.query_qualities))
                mode_result = stats.mode(read.query_qualities, keepdims=True) # keepdims=True returns array-like output; mode: The array of modal values. count: The number of times each corresponding modal value appears. # print(f"The mode is: {mode_result.mode[0]}") # print(f"The count is: {mode_result.count[0]}")
                mode_base_qualities.append(mode_result.mode[0])
                median_base_qualities.append(np.median(read.query_qualities))
                num_zeros_base_qualities.append(np.sum(np.array(read.query_qualities) == 0))
                less_10_base_qualities.append(np.sum(np.array(read.query_qualities) < 10))

            else:
                avg_base_qualities.append(np.nan)

            num_overlap_aligned_bases = read.get_overlap(start, end) # self, uint32_t start, uint32_t end)
                # return number of aligned bases of read overlapping the interval start and end on the reference sequence.
                # Return None if cigar alignment is not available.

            num_nones = 0
            # num_indels = 0
            num_inserts = 0
            num_mismatches = 0
            num_ambiguous = 0
            qseq = read.query_sequence
            # Get aligned pairs with sequence
            for rpos, refpos in read.get_aligned_pairs(matches_only=False, with_seq=False):
                if rpos is None and refpos is None:
                    raise ValueError("Both rpos and refpos are None, which should not happen.")

                elif rpos is None:
                    # rpos is None, it's a deletion in the read (relative to the reference).
                    # deletion relative to reference; keep as NaN (or 0, but NaN is clearer) 
                    # num_indels += 1
                    num_nones += 1

                elif refpos is None:
                    # If refpos is None, it's an insertion or soft clip in the read (relative to the reference).
                    # insertion relative to reference (or soft clip) → no reference column
                    num_inserts += 1
                     
                else:
                    # Both rpos and refpos are not None, the base is aligned (match or mismatch)
                    base = qseq[rpos].upper()
                    ref_base = fa.fetch(bam.get_reference_name(read.reference_id), refpos, refpos+1).upper()
                    if base != ref_base and base in "ACGT" and ref_base in "ACGT":
                        num_mismatches += 1
                    elif base not in "ACGT" or ref_base not in "ACGT":
                        num_ambiguous += 1 # count ambiguous bases as None
                
            results.append({
                "read_name_str": read.query_name,
                "read_lengths": read.query_length,
                "mapping_qualities": read.mapping_quality,
                "avg_base_qualities": avg_base_qualities[-1] if len(avg_base_qualities) > 0 else np.nan,
                "mode_base_qualities": mode_base_qualities[-1] if len(mode_base_qualities) > 0 else np.nan,
                "median_base_qualities": median_base_qualities[-1] if len(median_base_qualities) > 0 else np.nan,
                "num_overlap_aligned_bases": num_overlap_aligned_bases,
                "fraction_overlap_aligned": num_overlap_aligned_bases / (end - start) if region else None,
                # "num_indels": num_indels,
                "num_nones": num_nones,
                "fraction_nones": num_nones / len(read.query_sequence) if len(read.query_sequence) > 0 else None,
                "num_inserts": num_inserts,
                "fraction_inserts": num_inserts / len(read.query_sequence) if len(read.query_sequence) > 0 else None,
                "num_mismatches": num_mismatches,
                "fraction_mismatches": num_mismatches / len(read.query_sequence) if len(read.query_sequence) > 0 else None,
                "num_ambiguous": num_ambiguous,

            })
    fa.close()

    read_indel_mismatch_counts_df = pd.DataFrame(results)
    if output_csv:
        read_indel_mismatch_counts_df.to_csv(output_csv, index=False)
    return read_indel_mismatch_counts_df
# Example usage:
# bam_path = "/path/to/your.bam"
# ref_genome_path = "/path/to/reference.fasta"
# region = "chr1:100000-200000"  # Optional region
# output_csv = "/path/to/output.csv"  # Optional output CSV path
# df = count_indels_and_mismatches(bam_path, ref_genome_path, region, output_csv)
# print(df.head())          
