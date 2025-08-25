import pysam
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, Tuple, Dict, List
import seaborn as sns
import matplotlib.pyplot as plt


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
        title = f"Alignment Heatmap [{meta.get('region','?')}] • mode={meta.get('mode','?')}"
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
        # sns.clustermap(G_revs_df.fillna(-1), col_cluster=False)

        g = sns.clustermap(pd.DataFrame(matrix).fillna(-1), col_cluster=False, cmap='plasma')  # , cbar_kws={'label': "Y labet") # , 'orientation': 'vertical'

        g.ax_cbar.set_title("Clustered" + title, pad=16)

        # Display the plot
        plt.show()