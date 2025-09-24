import platform
import sys
import pysam
from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import os

from datetime import datetime
from dimelo import parse_bam, plot_reads, load_processed, plot_read_browser
import h5py


# Import the function from funcs_check_quality_bams.py
sys.path.append("/home/michalula/code/epiCausality/epiCode/utils/")
from funcs_check_quality_bams import count_indels_and_mismatches, plot_reads_quality_heatmap

def system_info():
    """Print system information."""
    print('System:', platform.system())
    print('Release:', platform.release())
    print('Version:', platform.version())
    print('Processor:', platform.processor())
    print('Python version:', sys.version)

def get_reference_sequence(ref_genome_file,  region_chr, region_start, region_end):
    """Fetch reference sequence from genome file."""
    if not ref_genome_file.exists():
        print(f"Reference genome file not found: {ref_genome_file}")
        sys.exit(1)
    try:
        ref_seq = pysam.FastaFile(ref_genome_file).fetch(region_chr, region_start-1, region_end-1) # removing 1 from the coord indexes s the reference genome starts with 1 and the pysam.FastaFile is 0-based and the coordinates are the index-1
        # ref_seq = pysam.FastaFile(ref_genome_file).fetch(region_chr, region_start, region_end) # removing 1 from the coord indexes s the reference genome starts with 1 and the pysam.FastaFile is 0-based and the coordinates are the index-1
        ref_seq_list = list(ref_seq)
        print(ref_seq)
        print(len(ref_seq))
        return ref_seq_list
    except Exception as e:
        print("Error fetching reference sequence:", e)
        return None

def create_output_directory(path):
    """Create output directory if it doesn't exist."""
    try:
        output_dir = Path(path)
        output_dir.mkdir(exist_ok=True)
        return output_dir
    except Exception as e:
        print("Error creating output directory:", e)
        return None

def extract_from_bam(experiment_name, bam_path, ref_genome_file, output_dir, window_size=None, threshold_mC=0.99, 
                    num_cores=32, region='chr1:206586162-206586192', motifs=['CG,0'], 
                    output_name='extract_output', save_fig=True):
    """Processes a BAM file using parse_bam.extract and plots the extracted reads."""
    try:
        # Parse region to calculate region length
        region_chr, region_coords = region.split(':')
        region_start, region_end = map(int, region_coords.split('-'))
        region_length = region_end - region_start
        print(f"Region length: {region_length}")

        extract_file, extract_region = parse_bam.extract(
            input_file=bam_path,
            output_name=output_name,
            ref_genome=ref_genome_file,
            output_directory=output_dir,
            regions=region,
            motifs=motifs,
            thresh=threshold_mC,
            window_size=window_size,
        )

        if threshold_mC == None: 
            print("Sort by strand reads post shuffling")
            fig_plot_browser = plot_read_browser.plot_read_browser(
                mod_file_name=extract_file,# mod_file_name: str | Path,
                region=region, #region: str,
                motifs=motifs, #motifs: list[str],
                thresh=threshold_mC, # thresh: int | float | None = None,
                single_strand = False, # : bool = False,
                sort_by=['shuffle', 'strand'], #: str | list[str] = "shuffle",
                hover = True, #: bool = True,
                )
            fig_plot_browser.update_layout(  
                title=f"{experiment_name}<br>Extracted Reads for {region}",
            )
            print("Do not suffle the reads!")
            fig_plot_browser = plot_read_browser.plot_read_browser(
                mod_file_name=extract_file,# mod_file_name: str | Path,
                region=region, #region: str,
                motifs=motifs, #motifs: list[str],
                thresh=threshold_mC, # thresh: int | float | None = None,
                single_strand = False, # : bool = False,
                # sort_by=['shuffle', 'strand'], #: str | list[str] = "shuffle",
                hover = True, #: bool = True,
                )
            fig_plot_browser.update_layout(  
                title=f"{experiment_name}<br>Extracted Reads for {region}",
            )
            # fig.show()
            if save_fig:
                output_html_path = Path(output_dir) / f"plot_browser_{region_length}bps_{experiment_name}_extract_reads_{region}.html"
                fig_plot_browser.write_html(str(output_html_path))
                print(f"Plot browser html figure saved to {output_html_path}")
            return extract_file, extract_region, fig_plot_browser
        else: 
            plot_reads.plot_reads(
                extract_file,
                region,
                motifs=motifs,
                window_size=window_size,
                sort_by=['shuffle', 'strand'],
                s=1
            )
            plt.xlabel(f'bp relative to {region}')
            plt.title(f"{experiment_name}\n Extracted Reads for {region} \n Sort reads by strand post shuffling")
            plt.show()

            plot_reads.plot_reads(
                extract_file,
                region,
                motifs=motifs,
                window_size=window_size,
                # sort_by=['shuffle', 'strand'],
                s=1
            )
            plt.xlabel(f'bp relative to {region}')
            plt.title(f"{experiment_name}\n Extracted Reads for {region} \n Do not suffle the reads")
            plt.show()
        return extract_file, extract_region
    
    except Exception as e:
        print("Error in BAM extraction:", e)
        return None, None

def process_extracted_reads_no_fully_unmethylated(extract_file, region, motifs, ref_seq_list):
    """
    Process extracted reads into a DataFrame.

    Warning: make sure that the ref_seq_list was created using the same region and reference genome using the function: 
        motifs=['CG,0']
        ref_seq_list = get_reference_sequence(ref_genome_v1_1_file, region_chr, region_start, region_end)  
    """
    try:
        reads, read_names, mods, region_dict = load_processed.readwise_binary_modification_arrays(
            file=extract_file,
            regions=region,
            motifs=motifs
        )
        mCG_reads_df = pd.DataFrame({
            'read_name': read_names,
            'mod': mods,
            'pos': reads
        }).explode('pos')

        # mCG_reads_df['pos_shifted'] = mCG_reads_df['pos'] + 15
        region_length = len(ref_seq_list)
        mCG_reads_df['pos_shifted'] = mCG_reads_df['pos'] + (region_length // 2)
        return mCG_reads_df, region_dict
    except Exception as e:
        print("Error processing extracted reads:", e)
        return None, None


# def process_extracted_reads_add_fully_unmethylated(extract_file, region, motifs, ref_seq_list):
#     """
#     Process extracted reads into a DataFrame, ensuring all reads (methylated and unmethylated) are included.
#     """
#     try:
#         # Extract methylation-modified positions
#         mod_coords, read_ids, mods, region_dict = load_processed.readwise_binary_modification_arrays(
#             file=extract_file,
#             regions=region,
#             motifs=motifs
#         )

#         # Get all read names (both methylated and unmethylated)
#         with h5py.File(extract_file, "r") as h5:
#             all_read_names = np.array(h5["read_name"], dtype=str)  # Extract all read names

#         # Convert read IDs to strings to avoid type mismatches
#         read_ids = np.array(read_ids, dtype=str)

#         # Create a DataFrame for methylated reads
#         mCG_reads_df = pd.DataFrame({
#             'read_name': read_ids,
#             'mod': mods,
#             'pos': mod_coords
#         })

#         # Ensure 'pos' is converted to numeric type
#         mCG_reads_df['pos'] = pd.to_numeric(mCG_reads_df['pos'], errors='coerce')

#         # Identify unmethylated reads (present in BAM but missing from mCG_reads_df)
#         methylated_reads = set(mCG_reads_df['read_name'])
#         unmethylated_reads = [read for read in all_read_names if read not in methylated_reads]

#         # Create a DataFrame for unmethylated reads (no positions)
#         unmethylated_df = pd.DataFrame({
#             'read_name': unmethylated_reads,
#             'mod': None,
#             'pos': np.nan  # No methylation site
#         })

#         # Combine both DataFrames
#         mCG_reads_df = pd.concat([mCG_reads_df, unmethylated_df], ignore_index=True)

#         # Ensure 'pos' is numeric for all rows
#         mCG_reads_df['pos'] = pd.to_numeric(mCG_reads_df['pos'], errors='coerce')

#         # Compute shifted positions
#         region_length = len(ref_seq_list)
#         mCG_reads_df['pos_shifted'] = mCG_reads_df['pos'].apply(
#             lambda x: int(x + (region_length // 2)) if not np.isnan(x) else -1
#         )

#         return mCG_reads_df, region_dict

#     except Exception as e:
#         print("Error processing extracted reads:", e)
#         return None, None

# def has_no_NONE_alignment_scores(row):
#     """
#     Check if all alignment scores for a read are not None.

#     Parameters
#     ----------
#     row : pandas.Series or dict
#         A row from a DataFrame representing a read, expected to contain an 'alignment_scores' field
#         which is a list or array of scores for each base in the region.

#     Output
#     ------
#     bool
#         Returns True if all alignment scores are not None, otherwise False.

#     Goal
#     ----
#     To filter and identify reads that have complete alignment information (no missing scores)
#     across the entire region of interest.
#     """
#     scores = row.get('alignment_scores')
#     if scores is None:
#         return False
#     # Check all scores are not None
#     return all(score is not None for score in scores)


# def count_alignment_scores(row):
#     """
#     Count the total number of alignment_scores that are None in the row.

#     Parameters
#     ----------
#     row : pandas.Series or dict
#         A row from a DataFrame representing a read, expected to contain an 'alignment_scores' field
#         which is a list or array of scores for each base in the region.

#     Output
#     ------
#     int
#         Returns the count of alignment scores that are None.

#     Goal
#     ----
#     To quantify the number of missing alignment scores for a read across the region of interest.
#     """
#     print(row)
#     scores = row.get('alignment_scores')
#     print('scores:', scores)
#     if scores is None:
#         return 0
#     return sum(score is None for score in scores)


def process_extracted_reads(extract_file, original_bam_path, region, motifs, 
                            ref_seq_list, ref_genome_path,
                            experiment_name, output_dir, 
                            keep_full_coverage_reads_only=True,
                            save_indels_mismatches_count_csv_path="indels_mismatches_count.csv",
                            threshold_fraction_overlap_aligned=0.5, threshold_fraction_mismatches=0.5, 
                            threshold_mapping_qualities=60, threshold_avg_base_qualities=20,
                            max_reads_plot=3000):
                            # indel_fraction_threshold=0.5, non_fraction_threshold=0.5):
    """
    Process extracted reads into a DataFrame, ensuring only reads that cover the full start to end DNA coordinates are included.
    This function filters out reads that don't span the complete region of interest.
    
    Parameters:
    -----------
    extract_file : str or Path
        Path to the HDF5 file containing extracted read data
    original_bam_path : str or Path
        Path to the corresponding original BAM file for counting indels and mismatches
    region : str
        Genomic region in format 'chr:start-end' (e.g., 'chr1:206586162-206586192')
    motifs : list
        List of motifs to search for (e.g., ['CG,0'])
    ref_seq_list : list
        Reference sequence as a list of nucleotides
    ref_genome_path : str or Path
        Path to the reference genome FASTA file
    experiment_name:    str
        Name of the experiment for labeling outputs
    output_dir : str or Path    
        Directory to save output files
    keep_full_coverage_reads_only : bool, optional
        If True, only keeps reads that cover the full start to end DNA coordinates. If False, keeps all reads.
        Default is True.
    save_indels_mismatches_count_csv_path : str, optional
        Path to save the CSV file with indel and mismatch counts per read. Default is "indels_mismatches_count.csv".
    
    threshold_fraction_overlap_aligned : float, optional
        Threshold for the fraction of overlap aligned bases allowed per read. Reads with a lower fraction will be removed.
        Default is 0.5 (50%).
    threshold_fraction_mismatches : float, optional
        Threshold for the fraction of mismatches allowed per read. Reads with a higher fraction will be removed.    
        Default is 0.5 (50%).

    threshold_mapping_qualities :  float, optional
        Threshold for the mapping quality allowed per read. Reads with a lower mapping quality will be removed.
        Default is 60.

    threshold_avg_base_qualities : float, optional
        Threshold for the average base quality allowed per read. Reads with a lower average base quality will be removed.
        Default is 20.

    # indel_fraction_threshold : float, optional  
    #     Threshold for the fraction of indels (deletions/insertions/soft clips) allowed per read. Reads with a higher fraction will be removed.
    #     Default is 0.5 (50%).
    # non_fraction_threshold : float, optional    
    #     Threshold for the fraction of Nones allowed per read. Reads with a higher fraction will be removed.     
    #     Default is 0.5 (50%).
    max_reads_plot :  int, optional
        Maximum number of reads to plot in the heatmap. Default is 3000.
        
    Returns:
    --------
    tuple : (DataFrame, dict)
        - DataFrame with filtered reads containing columns: read_name, mod, pos, pos_shifted
        - Dictionary containing region information
        
    Example:
    --------
    # Extract reads that cover the full region
    mCG_reads_df, region_dict = process_extracted_reads(
        extract_file='path/to/extracted_reads.h5',
        region='chr1:206586162-206586192',
        motifs=['CG,0'],
        ref_seq_list=['A', 'T', 'G', 'C', ...],
        keep_full_coverage_reads_only=True
    )
    """
    # try:
    # Extract methylation-modified positions
    mod_coords, read_ids, mods, region_dict = load_processed.readwise_binary_modification_arrays(
        file=extract_file,
        regions=region,
        motifs=motifs
    )

    # Get all read names and their coordinate information
    read_starts = None
    read_ends = None
    read_coords = None
    
    with h5py.File(extract_file, "r") as h5:
        all_read_names = np.array(h5["read_name"], dtype=str)  # Extract all read names
        
        # Print the keys in the .h5 file
        # print("Available keys in HDF5 file:", list(h5.keys()))
        
        # Try to get read coordinate information if available
        read_coords_available = False
        try:
            # Check if read coordinates are available in the HDF5 file
            if "read_start" in h5.keys() and "read_end" in h5.keys():
                read_starts = np.array(h5["read_start"])
                read_ends = np.array(h5["read_end"])
                read_coords_available = True
                print(f"Found read coordinates: {len(read_starts)} reads")
            elif "read_coords" in h5.keys():
                read_coords = np.array(h5["read_coords"])
                read_coords_available = True
                print(f"Found read_coords: {len(read_coords)} reads")
            else:
                print("Warning: Read coordinate information not found in HDF5 file.")
        except Exception as e:
            print(f"Warning: Could not access read coordinate information: {e}")

    # Parse the region to get start and end coordinates
    region_chr, region_coords = region.split(':')
    region_start, region_end = map(int, region_coords.split('-'))
    region_length = region_end - region_start

    # Convert read IDs to actual read names using the all_read_names array
    # read_ids are indices into the all_read_names array
    methylated_read_names = []
    for read_id in read_ids:
        try:
            read_index = int(read_id)
            if read_index < len(all_read_names):
                methylated_read_names.append(all_read_names[read_index])
            else:
                print(f"Warning: Read ID {read_id} is out of bounds for all_read_names array")
                methylated_read_names.append(f"unknown_{read_id}")
        except (ValueError, TypeError):
            print(f"Warning: Could not convert read_id {read_id} to integer")
            methylated_read_names.append(f"unknown_{read_id}")

    # Convert read IDs to strings to avoid type mismatches
    read_ids = np.array(read_ids, dtype=str)

    # Create a DataFrame for methylated reads
    mCG_reads_df = pd.DataFrame({
        'read_name_str': methylated_read_names,
        'read_name': read_ids,
        'read_id_number': read_ids,
        'mod': mods,
        'pos': mod_coords
    })
    
    print(f"Unique read names with methylation: {len(mCG_reads_df['read_name'].unique())}")

    read_indel_mismatch_counts_df = count_indels_and_mismatches(original_bam_path, ref_genome_path, region, output_csv=save_indels_mismatches_count_csv_path)

    # Filter reads based on coverage criteria
    if keep_full_coverage_reads_only:
        if read_coords_available:
            # Method 1: Use explicit read coordinates if available
            if read_starts is not None and read_ends is not None:
                # Create a mapping of read names to their coordinates
                read_coord_map = {}
                for i, read_name in enumerate(all_read_names):
                    read_coord_map[read_name] = (read_starts[i], read_ends[i])
                
                # Filter reads that cover the full region
                full_coverage_reads = []
                for read_name in all_read_names:
                    if read_name in read_coord_map:
                        read_start, read_end = read_coord_map[read_name]
                        # Check if read covers the full region (read starts before or at region start and ends after or at region end)
                        if read_start <= region_start and read_end >= region_end:
                            full_coverage_reads.append(read_name)
                
                # Check overlap between full coverage reads and reads with methylation data
                reads_with_methylation = set(mCG_reads_df['read_name_str'].unique())
                full_coverage_set = set(full_coverage_reads)
                overlap = reads_with_methylation.intersection(full_coverage_set)
                
                print(f"Found {len(full_coverage_reads)} reads with full coverage")
                print(f"Reads with methylation data: {len(reads_with_methylation)}")
                print(f"Overlap between full coverage and methylation: {len(overlap)}")
                
                # Filter the mCG_reads_df to only include full coverage reads
                mCG_reads_df = mCG_reads_df[mCG_reads_df['read_name_str'].isin(full_coverage_reads)]
                print(f"After full coverage filtering: {len(np.unique(mCG_reads_df['read_name_str']))} reads with methylation data")
               
                # NEED TO REMOVE READS THAT ARE ONLY PARTLY ALIGNED TO THE REGION.. 
                # Removing reads that have bases with None alignment scores 
                # count_alignment_scores_df = mCG_reads_df.apply(count_alignment_scores, axis=1)
                # print(f"Count of None alignment scores per read:\n{count_alignment_scores_df.value_counts()}")
                # print(f"Reads with less than 100 None alignment scores: {(count_alignment_scores_df < 100).sum()}")
                # print("Reads with 100 or more None alignment scores")
                # mCG_reads_df = mCG_reads_df[mCG_reads_df.apply(count_alignment_scores, axis=1) < 100]
                # # mCG_reads_df = mCG_reads_df[mCG_reads_df.apply(has_no_NONE_alignment_scores, axis=1)]
                # print(f"After removing reads that have bases with None alignment scores and full coverage filtering: {len(np.unique(mCG_reads_df['read_name_str']))} reads with methylation data")
                # Remove reads with more than 50% deletions in the region
                # Merge mCG_reads_df with read_indel_mismatch_counts_df on read_name/read_name_str
                reads_with_overlap_indel_mismatch_counts_df = mCG_reads_df.merge(
                    read_indel_mismatch_counts_df,
                    left_on='read_name_str',
                    right_on='read_name_str',
                    how='left'
                )
                # print(reads_with_indel_mismatch_counts_df)
                # print(f"After full coverage filtering: {len(np.unique(reads_with_indel_mismatch_counts_df['read_name_str']))} reads with methylation data")

                # Calculate fraction of aligned bases per read
                # reads_with_overlap_indel_mismatch_counts_df['overlap_aligned_fraction'] = reads_with_overlap_indel_mismatch_counts_df['num_overlap_aligned_bases'] / region_length
                # Keep only reads with indel_fraction <= indel_fraction_threshold
                # filtered_reads_with_overlap_indel_mismatch_counts_df = reads_with_overlap_indel_mismatch_counts_df.copy()


                # if save_filtered_reads_bam:
                pre_filtered_reads_bam_name = "pre_filtered_ROI_reads_" + experiment_name + ".bam"
                pre_filtered_output_bam_path=Path(output_dir, pre_filtered_reads_bam_name)
                # filtered_output_bam_path
                subset_BAM_by_read_IDs(original_bam_path, reads_with_overlap_indel_mismatch_counts_df, output_bam_path=pre_filtered_output_bam_path, index_output=True)
                print(f"Pre Filtered BAM with reads that have some covarage of the ROI written to \n {pre_filtered_reads_bam_name}")
                plot_reads_quality_heatmap(pre_filtered_output_bam_path, 
                                            region, 
                                            ref_genome_path,
                                            heatmap_title = "Pre filtering within ROI: Reads Quality Heatmap " + experiment_name,
                                            mode= "quality_signed",
                                            max_reads =  3000,
                                            primary_only = True,
                                            min_mapq = None,  
                                            save_matrix= False)


                #TODO: Plot histogram of the the reads overlap_aligned kengths 
                # plot_histogram(data, title, num_bins=16,  xaxis_title="mod_vector values")
                plot_bam_quality_metrics(pre_filtered_output_bam_path)

                filtered_reads_with_overlap_indel_mismatch_counts_df = reads_with_overlap_indel_mismatch_counts_df[reads_with_overlap_indel_mismatch_counts_df['fraction_overlap_aligned'] >= threshold_fraction_overlap_aligned].copy()
                print(f"After removing reads with <{threshold_fraction_overlap_aligned*100}% threshold_fraction_overlap_aligned: {len(np.unique(filtered_reads_with_overlap_indel_mismatch_counts_df['read_name_str']))} reads with methylation data")

                
                filtered_reads_with_overlap_indel_mismatch_counts_df = filtered_reads_with_overlap_indel_mismatch_counts_df[filtered_reads_with_overlap_indel_mismatch_counts_df['fraction_mismatches'] < threshold_fraction_mismatches].copy()
                print(f"After removing reads with >{threshold_fraction_mismatches*100}% threshold_fraction_mismatches: {len(np.unique(filtered_reads_with_overlap_indel_mismatch_counts_df['read_name_str']))} reads with methylation data")

                # Remove reads with mapping_qualities >= 60
                filtered_reads_with_overlap_indel_mismatch_counts_df = filtered_reads_with_overlap_indel_mismatch_counts_df[filtered_reads_with_overlap_indel_mismatch_counts_df['mapping_qualities'] >= threshold_mapping_qualities].copy()
                print(f"After removing reads with >{threshold_mapping_qualities} threshold_mapping_qualities: {len(np.unique(filtered_reads_with_overlap_indel_mismatch_counts_df['read_name_str']))} reads with methylation data")

                # Remove reads with avg_base_qualities > 20
                filtered_reads_with_overlap_indel_mismatch_counts_df = filtered_reads_with_overlap_indel_mismatch_counts_df[filtered_reads_with_overlap_indel_mismatch_counts_df['avg_base_qualities'] >= threshold_avg_base_qualities].copy()
                print(f"After removing reads with >{threshold_avg_base_qualities} threshold_avg_base_qualities: {len(np.unique(filtered_reads_with_overlap_indel_mismatch_counts_df['read_name_str']))} reads with methylation data")

                # # Calculate fraction of Nones per read
                # filtered_reads_with_indel_mismatch_counts_df['nones_fraction'] = filtered_reads_with_indel_mismatch_counts_df['num_nones'] / region_length
                # # Keep only reads with num_none <= indel_fraction_threshold
                # filtered_reads_with_none_indel_mismatch_counts_df = filtered_reads_with_indel_mismatch_counts_df.copy()
                # # filtered_reads_with_none_indel_mismatch_counts_df = filtered_reads_with_indel_mismatch_counts_df[filtered_reads_with_indel_mismatch_counts_df['nones_fraction'] <= non_fraction_threshold].copy()
                # print(f"After removing reads with >{non_fraction_threshold*100}% Nones: {len(np.unique(filtered_reads_with_none_indel_mismatch_counts_df['read_name_str']))} reads with methylation data")

                # Update mCG_reads_df to the filtered DataFrame
                mCG_reads_df = filtered_reads_with_overlap_indel_mismatch_counts_df

            elif read_coords is not None:
                # Method 2: Use read_coords if available (assuming it contains start/end pairs)
                full_coverage_reads = []
                for i, read_name in enumerate(all_read_names):
                    if i < len(read_coords):
                        read_start, read_end = read_coords[i]
                        if read_start <= region_start and read_end >= region_end:
                            full_coverage_reads.append(read_name)
                
                mCG_reads_df = mCG_reads_df[mCG_reads_df['read_name_str'].isin(full_coverage_reads)]
                print(f"Found {len(full_coverage_reads)} reads with full coverage")
                print(f"After full coverage filtering: {len(np.unique(mCG_reads_df['read_name_str']))} reads with methylation data")
                
                # # Removing reads that have bases with None alignment scores 
                # mCG_reads_df = mCG_reads_df[mCG_reads_df.apply(has_no_NONE_alignment_scores, axis=1)]
                # print(f"After removing reads that have bases with None alignment scores and full coverage filtering: {len(np.unique(mCG_reads_df['read_name_str']))} reads with methylation data")
                

        else:
            print("Warning: No coordinate information available for full coverage filtering. Keeping all reads with methylation information.")
    
    # # Debug: Show the mapping between read IDs and actual read names
    # print(f"Converted {len(read_ids)} read IDs to actual read names")
    # print(f"Sample mapping: read_id 0 -> {read_names[0] if len(read_names) > 0 else 'N/A'}")
    # Ensure 'pos' is numeric for all rows
    mCG_reads_df['pos'] = pd.to_numeric(mCG_reads_df['pos'], errors='coerce')

    # Compute shifted positions
    region_length = len(ref_seq_list)
    mCG_reads_df['pos_shifted'] = mCG_reads_df['pos'].apply(
        lambda x: int(x + (region_length // 2)) if not np.isnan(x) else -1
    )

    print(f"Final result: {len(mCG_reads_df)} reads with methylation information out of {len(all_read_names)} total reads")
    
    # if save_filtered_reads_bam:
    filtered_reads_bam_name = "filtered_reads_overlap_MORE_than_" + str(threshold_fraction_overlap_aligned) + "_" + experiment_name + ".bam"
    filtered_output_bam_path=Path(output_dir, filtered_reads_bam_name)
    # filtered_output_bam_path
    subset_BAM_by_read_IDs(original_bam_path, mCG_reads_df, output_bam_path=filtered_output_bam_path, index_output=True)
    print(f"Filtered BAM with reads that have mC and base overlap > {threshold_fraction_overlap_aligned} written to \n {filtered_output_bam_path}")
    plot_reads_quality_heatmap(original_bam_path, 
                                region, 
                                ref_genome_path,
                                heatmap_title = "Read Quality Heatmap for all original full BAM " + experiment_name,
                                mode= "quality_signed",
                                max_reads =  max_reads_plot,
                                primary_only = True,
                                min_mapq = None,  
                                save_matrix= False)

    plot_reads_quality_heatmap(filtered_output_bam_path, 
                                region, 
                                ref_genome_path,
                                heatmap_title = "Read Quality Heatmap for filtered reads " + experiment_name,
                                mode= "quality_signed",
                                max_reads = max_reads_plot,
                                primary_only = True,
                                min_mapq = None,  
                                save_matrix= False)

    plot_bam_quality_metrics(original_bam_path)
    plot_bam_quality_metrics(filtered_output_bam_path)
    
    return mCG_reads_df, region_dict

    # except Exception as e:
    #     print("Error processing extracted reads with full coverage filtering:", e)
    #     return None, None


def remove_low_methylated_reads(mCG_reads_df, threshold_percent=50):
    """
    Remove reads that have less than a specified percentage of the maximum number 
    of methylated CGs per read in the dataset.
    
    Parameters:
    -----------
    mCG_reads_df : pandas.DataFrame
        DataFrame containing read methylation data with columns: read_name, mod, pos, pos_shifted
    region_dict : dict
        Dictionary containing region information
    threshold_percent : float, optional
        Percentage threshold relative to the mean number of methylated CGs per read.
        Reads with fewer methylated CGs than this percentage will be removed.
        Default is 10 (10% of maximum).
        
    Returns:
    --------
    tuple : (DataFrame, dict)
        - Filtered DataFrame with low methylated reads removed
        
    Example:
    --------
    # Remove reads with less than 10% of max methylation
    filtered_mCG_reads_df = remove_low_methylated_reads(
        mCG_reads_df=mCG_reads_df,
        threshold_percent=10
    )
    """
    try:
        # Count methylated CGs per read
        # Convert mod column to numeric to enable mathematical operations
        mCG_reads_df['num_CG_methylated'] = pd.to_numeric(mCG_reads_df['mod'], errors='coerce').fillna(1)
        # mCG_reads_df['num_CG_methylated'] = 1
        
        # Group by read_name and count the number of methylation events (mod=1)
        methylation_counts = mCG_reads_df.groupby('read_name')['num_CG_methylated'].sum().reset_index()
        methylation_counts.columns = ['read_name', 'methylation_count']
        
        # # Find the maximum number of methylated CGs across all reads
        max_methylation = methylation_counts['methylation_count'].max()
        # # Find the median number of methylated CGs across all reads
        # max_methylation = methylation_counts['methylation_count'].median()
        # Find the mean number of methylated CGs across all reads
        mean_methylation = methylation_counts['methylation_count'].median()
        
        # Calculate the threshold (10% of maximum by default)
        threshold = mean_methylation * (threshold_percent / 100)
        
        print(f"Mean methylated CGs per read: {mean_methylation}")
        print(f"Maximum methylated CGs per read: {max_methylation}")
        print(f"Threshold ({threshold_percent}% of max): {threshold:.2f}")
        
        # Get reads that meet the threshold
        reads_to_keep = methylation_counts[methylation_counts['methylation_count'] >= threshold]['read_name'].tolist()
        reads_to_remove = methylation_counts[methylation_counts['methylation_count'] < threshold]['read_name'].tolist()
        
        # Filter the original DataFrame
        filtered_mCG_reads_df = mCG_reads_df[mCG_reads_df['read_name'].isin(reads_to_keep)].copy()
        remove_mCG_reads_df = mCG_reads_df[mCG_reads_df['read_name'].isin(reads_to_remove)].copy()
        
        print(f"Original number of reads: {len(mCG_reads_df['read_name'].unique())}")
        print(f"Number of reads after filtering: {len(filtered_mCG_reads_df['read_name'].unique())}")
        print(f"Removed {len(mCG_reads_df['read_name'].unique()) - len(filtered_mCG_reads_df['read_name'].unique())} reads")
        
        return filtered_mCG_reads_df, methylation_counts, remove_mCG_reads_df
        
    except Exception as e:
        print("Error in remove_low_methylated_reads:", e)
        return mCG_reads_df

import pysam
from pathlib import Path

def subset_BAM_by_read_IDs(bam_path, remove_mCG_reads_df, output_bam_path=None, index_output=True):
    """
    Create a BAM that contains only reads whose names are listed in remove_mCG_reads_df['read_name'].

    Parameters
    ----------
    bam_path : str | Path
        Path to the source BAM.
    remove_mCG_reads_df : pandas.DataFrame
        DataFrame with a column 'read_name' (the qname in the BAM). Duplicates are OK.
    output_bam_path : str | Path | None
        Where to write the subset BAM. If None, will write alongside the input as
        '<input_basename>.subset_by_remove_reads.bam'.
    index_output : bool
        If True, create a .bai index for the output BAM.

    Returns
    -------
    Path
        Path to the written subset BAM.
    """
    bam_path = Path(bam_path)

    if output_bam_path is None:
        output_bam_path = bam_path.with_suffix("")  # strip .bam
        output_bam_path = Path(str(output_bam_path) + ".subset_by_remove_reads.bam")
    else:
        output_bam_path = Path(output_bam_path)

    # Determine which column has the read names:
    # Prefer 'read_name' (typical). If not present, try a common alternative used upstream.
    if 'read_name_str' in remove_mCG_reads_df.columns:
        name_series = remove_mCG_reads_df['read_name_str']
    elif 'read_name' in remove_mCG_reads_df.columns:
        name_series = remove_mCG_reads_df['read_name']
    else:
        raise ValueError(
            "remove_mCG_reads_df must contain a 'read_name' or 'read_name_str' column with BAM qnames."
        )

    # Build a set for fast lookups; ensure strings
    target_names = set(map(str, name_series.dropna().unique()))
    if len(target_names) == 0:
        raise ValueError("remove_mCG_reads_df has no read names to subset.")

    # Open input and create output with the same header
    with pysam.AlignmentFile(bam_path, "rb") as in_bam:
        with pysam.AlignmentFile(output_bam_path, "wb", header=in_bam.header) as out_bam:
            # Iterate over all records (including unmapped) safely
            for aln in in_bam.fetch(until_eof=True):
                # query_name is the read's QNAME in BAM
                if aln.query_name in target_names:
                    out_bam.write(aln)

    # Optionally index the output BAM
    if index_output:
        # Remove existing index if present to avoid pysam errors
        bai = output_bam_path.with_suffix(output_bam_path.suffix + ".bai")
        if bai.exists():
            bai.unlink()
        pysam.index(str(output_bam_path))

    print(f"Subset BAM written to: {output_bam_path}")
    if index_output:
        print(f"Index written to: {output_bam_path}.bai")

    return output_bam_path



def bam_to_sam(bam_path, sam_path=None):
    # # Example usage:
    # bam_to_sam("input.bam", "output.sam")
    bam_path = Path(bam_path)
    if sam_path is None:
        sam_path = bam_path.with_suffix(".sam")
    else:
        sam_path = Path(sam_path)

    with pysam.AlignmentFile(bam_path, "rb") as bam_file:
        with pysam.AlignmentFile(sam_path, "w", header=bam_file.header) as sam_file:
            for read in bam_file.fetch(until_eof=True):
                sam_file.write(read)

    print(f"Converted BAM → SAM: {sam_path}")
    return sam_path


import pysam
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path

def plot_bam_quality_metrics(bam_path):
    # # Example usage
    # plot_bam_quality_metrics("input.bam")

    bam_path = Path(bam_path)

    read_lengths = []
    mapping_qualities = []
    avg_base_qualities = []

    with pysam.AlignmentFile(bam_path, "rb") as bam_file:
        for read in bam_file.fetch(until_eof=True):
            # Skip secondary/supplementary if you only want primary alignments
            if read.is_secondary or read.is_supplementary:
                continue

            read_lengths.append(read.query_length)
            mapping_qualities.append(read.mapping_quality)
            if read.query_qualities is not None:
                avg_base_qualities.append(np.mean(read.query_qualities))
            else:
                avg_base_qualities.append(np.nan)

    # Set plot style
    # sns.set(style="whitegrid", font_scale=1.2)
    sns.set(font_scale=1.2)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    

    # Read length distribution
    sns.histplot(read_lengths, kde=False, ax=axes[0], color="steelblue") # , bins=100
    axes[0].set_title("Read Length Distribution")
    axes[0].set_xlabel("Read Length (bp)")
    axes[0].set_ylabel("Count")
    #TODO: make the limit setting an input variable? or plot more distributions of the data
    axes[0].set_xlim(min(read_lengths), min(max(read_lengths), 50000))

    # Mapping quality distribution
    sns.histplot(mapping_qualities, kde=False, ax=axes[1], color="orange") #  bins=60,
    axes[1].set_title("Mapping Quality Distribution")
    axes[1].set_xlabel("Mapping Quality")
    axes[1].set_ylabel("Count")

    # Average base quality distribution
    sns.histplot(avg_base_qualities, kde=False, ax=axes[2], color="green") #  bins=50,
    axes[2].set_title("Average Base Quality per Read")
    axes[2].set_xlabel("Average Phred Quality Score")
    axes[2].set_ylabel("Count")

    # plt.title(f"Run Quality for bam file \n {os.path.basename(bam_path)}\n Total reads processed: {len(read_lengths)}")
    plt.suptitle(f"Run Quality for bam file \n {os.path.basename(bam_path)}\n Total reads processed: {len(read_lengths)}", fontsize=12) 

    fig.subplots_adjust(top=0.88) # Adjust this value as needed

    plt.tight_layout()
    plt.show()

    print(f"Total reads processed: {len(read_lengths)}")



# def visualize_data_old(mCG_reads_df):
#     """Generate visualizations for the data."""
#     try:
#         mCG_reads_df['read_name'].plot(kind='hist', bins=1600, title='#mC of individual Reads Distribution')
#         plt.gca().spines[['top', 'right']].set_visible(False)
#         plt.show()

#         sns.scatterplot(
#             data=mCG_reads_df,
#             x="pos",
#             y="read_name",
#             hue="mod",
#             s=0.5,
#             marker="s",
#             linewidth=0
#         )
#         # plt.xticks(ticks=np.arange(len(ref_seq_list)), labels=ref_seq_list, size=font_size) # 'small') #, rotation=90)
#         plt.xlabel('Position')
#         plt.ylabel('Read Name')
#         plt.show()
#     except Exception as e:
#         print("Error in visualization:", e)

def visualize_data(mCG_reads_df):
    """Generate visualizations for the data."""
    try:
        # Ensure 'pos' is numeric and drop NaNs
        mCG_reads_df = mCG_reads_df.copy()  # Avoid modifying the original DataFrame
        mCG_reads_df['pos'] = pd.to_numeric(mCG_reads_df['pos'], errors='coerce')
        mCG_reads_df = mCG_reads_df.dropna(subset=['pos'])  # Remove unmethylated reads for plotting

        # Histogram of read distribution
        mCG_reads_df['read_name'].value_counts().plot(kind='bar', title='#mC of individual Reads Distribution')
        plt.gca().spines[['top', 'right']].set_visible(False)
        plt.show()

        # Scatter plot of modifications
        sns.scatterplot(
            data=mCG_reads_df,
            x="pos",
            y="read_name",
            hue="mod",
            s=0.5,
            marker="s",
            linewidth=0
        )

        plt.xlabel('Position')
        plt.ylabel('Read Name')
        plt.show()
    except Exception as e:
        print("Error in visualization:", e)

def create_padded_reads_no_fully_unmethylated(mCG_reads_df, region_dict, region_length):
    """Generate padded reads matrix."""
    try:
        read_names_unique = np.unique(mCG_reads_df['read_name'])
        num_reads = len(read_names_unique)
        reads_dict = {name: i for i, name in enumerate(read_names_unique)}
        padded_reads = np.full((num_reads, region_length), np.nan)

        for i in range(len(mCG_reads_df['read_name'])):
            padded_reads[reads_dict[mCG_reads_df['read_name'][i]], mCG_reads_df['pos_shifted'][i]] = 1

        return padded_reads
    except Exception as e:
        print("Error creating padded reads matrix:", e)
        return None

def create_padded_reads(mCG_reads_df, region_dict, region_length):
    """Generate padded reads matrix, including reads with no methylation."""
    try:
        # Ensure 'pos_shifted' is numeric
        mCG_reads_df['pos_shifted'] = pd.to_numeric(mCG_reads_df['pos_shifted'], errors='coerce')

        read_names_unique = np.unique(mCG_reads_df['read_name'])
        num_reads = len(read_names_unique)
        reads_dict = {name: i for i, name in enumerate(read_names_unique)}
        padded_reads = np.full((num_reads, region_length), np.nan)

        for _, row in mCG_reads_df.iterrows():
            if row['pos_shifted'] >= 0:  # Ignore unmethylated reads (set as -1)
                padded_reads[reads_dict[row['read_name']], int(row['pos_shifted'])] = 1

        return padded_reads

    except Exception as e:
        print("Error creating padded reads matrix:", e)
        return None


def plot_padded_reads(padded_reads, ref_seq_list):
    """Plot padded reads matrix using matshow with x-ticks as reference sequence."""
    try:
        plt.figure(figsize=(10, 150))
        plt.matshow(padded_reads, fignum=1)
        plt.colorbar()
        plt.title("Padded Reads Matrix")

        # Scale font size: decreases as seq_length increases, but within reasonable bounds
        font_size = max(2, min(8, 500 / len(ref_seq_list)))  # Now it stays between 2 and 8

        plt.xticks(ticks=np.arange(len(ref_seq_list)), labels=ref_seq_list, size=font_size) # 'small') #, rotation=90)
        # plt.xlabel("Reference Sequence")

        plt.show()
    except Exception as e:
        print("Error plotting padded reads matrix:", e)


def save_padded_reads(padded_reads, output_dir, file_name):
    """Save padded reads as a NumPy array."""
    try:
        np.save(Path(output_dir, file_name), padded_reads)
        print(f"Padded reads saved to {file_name}")
    except Exception as e:
        print("Error saving padded reads:", e)


# =======================
# For unthresholded data plotting:

def plot_histogram(data, title, num_bins=16,  xaxis_title="mod_vector values"):
    # Compute the histogram
    hist, bin_edges = np.histogram(data, bins=num_bins)

    # Create the bar plot
    fig = go.Figure(data=[
        go.Bar(
            x=bin_edges[:-1],  # Start of each bin
            y=hist,            # Frequency in each bin
            width=np.diff(bin_edges),  # Width of each bin
            marker=dict(color='blue', opacity=0.7)
        )
    ])

    # Add labels and title
    fig.update_layout(
        title=title,
        xaxis_title= xaxis_title,
        yaxis_title="Frequency",
        bargap=0.1
    )

    return fig
    

# Define the number of bins
# num_bins = 50
# title=f"Distribution of Non-Zero mod_vector Values<br>Experiment: {experiment_name}<br>Region length: {region_length} [{region_str}]"
# fig_hist = plot_histogram(data=filtered_mod_vector_no0, num_bins=num_bins, title=title)
# # Show the plot
# fig_hist.show()

def plot_mov_values_percentages(filtered_mod_vector_no0, title, num_bins=16):
    # Compute the histogram
    hist, bin_edges = np.histogram(filtered_mod_vector_no0, bins=num_bins)

    # Normalize the histogram to percentages
    percentages = (hist / len(filtered_mod_vector_no0)) * 100

    # Create the bar plot
    fig = go.Figure(data=[
        go.Bar(
            x=bin_edges[:-1],  # Start of each bin
            y=percentages,     # Percentage in each bin
            width=np.diff(bin_edges),  # Width of each bin
            marker=dict(color='blue', opacity=0.7)
        )
    ])

    # Add labels and title
    fig.update_layout(
        title=title,
        xaxis_title="mod_vector values",
        yaxis_title="Percentage",
        bargap=0.1
    )

    # Show the plot
    fig.show()

    # Print the percentage values in each bin
    for i in range(len(percentages)):
        print(f"Bin {i + 1}: Range [{bin_edges[i]:.4f}, {bin_edges[i + 1]:.4f}) - Percentage: {percentages[i]:.2f}%")

# plot_mov_values_percentages(filtered_mod_vector_no0)

def parse_region(region_str):
    # Split the region string into chromosome and range
    region_chr, region_range = region_str.split(':')
    # Split the range into start and end
    region_start, region_end = map(int, region_range.split('-'))
    # Calculate the region length
    region_length = region_end - region_start
    return region_chr, region_start, region_end, region_length

# # Example usage
# region_chr, region_start, region_end, region_length = parse_region(region_str)
# print("region_chr:", region_chr)
# print("region_start:", region_start)
# print("region_end:", region_end)
# print("region_length:", region_length)


def plot_violin_mod_vector(mod_vector, title="Violin Plot of mod_vector Values"):
    """Plot a violin plot of the mod_vector values using seaborn."""
    plt.figure(figsize=(8, 4))
    sns.violinplot(y=mod_vector, inner="box", color="skyblue")
    plt.title(title)
    plt.ylabel("mod_vector value")
    plt.xlabel("")
    plt.tight_layout()
    plt.show()


def mod_vectors_noThreshold_analyze(
    experiment_name,
    extract_file, # bam_path,
    region_str,
    motifs,
    num_bins=16,
    # ref_genome_path,
    # output_dir,
    # threshold_mC=None,
    ):
    # extract_file, extract_region, fig_plot_browser = extract_from_bam(
    #     experiment_name=experiment_name,
    #     bam_path=bam_path,
    #     ref_genome_file=ref_genome_path,
    #     output_dir=output_dir,
    #     region=region_str,
    #     motifs=motifs,
    #     output_name='extracted_reads',
    #     threshold_mC=threshold_mC,
    # )
    sorted_read_tuples, readwise_datasets, region_dict = load_processed.read_vectors_from_hdf5(
                                            extract_file,  # extract_file, #     file: str | Path,
                                            motifs=motifs,  #     motifs: list[str],
                                            regions=region_str,  #     regions: str | Path | list[str | Path] | None = None,
                                            #     window_size: int | None = None,
                                            #     single_strand: bool = False,
                                            #     sort_by: str | list[str] = ["chromosome", "region_start", "read_start"],
                                            #     calculate_mod_fractions: bool = True,
                                        ) # -> tuple[list[tuple], list[str], dict | None]


    # Aggregate the second elements (mod_vector) from all tuples
    aggregated_mod_vector = [read[1] for read in sorted_read_tuples] # np.sum(, axis=0
    # print('aggregated_mod_vector', aggregated_mod_vector)                      

    flattened_mod_vector = np.concatenate(aggregated_mod_vector)
    # print("flattened_mod_vector", flattened_mod_vector)

    # Violin plot of the flattened_mod_vector
    plot_violin_mod_vector(flattened_mod_vector, title=f"Violin Plot of mod_vector Values\nExperiment: {experiment_name} [{region_str}]")

    # Filter out values equal to 0
    filtered_mod_vector_no0 = flattened_mod_vector[flattened_mod_vector != 0]

    # Define the number of bins
    # num_bins = 50
    region_chr, region_start, region_end, region_length = parse_region(region_str)
    title=f"Distribution of Non-Zero mod_vector Values<br>Experiment: {experiment_name}<br>Region length: {region_length} [{region_str}]"
    fig_hist = plot_histogram(data=filtered_mod_vector_no0, num_bins=num_bins, title=title)
    # Show the plot
    fig_hist.show()

    plot_mov_values_percentages(filtered_mod_vector_no0, 
                                title=f"Percentage Distribution of Non-Zero mod_vector Values<br>Experiment: {experiment_name}<br>Region length: {region_length} [{region_str}]",
                                num_bins=num_bins)

    return sorted_read_tuples, readwise_datasets, region_dict, aggregated_mod_vector, filtered_mod_vector_no0

# sorted_read_tuples, readwise_datasets, region_dict, aggregated_mod_vector, filtered_mod_vector_no0 = mod_vectors_analyze(
#         experiment_name,
#         extract_file, # bam_path,
#         region_str,
#         motifs
#     )


def main():
    """Main function to execute all tasks."""
    system_info()

    experiment_name = "unedited_T_primerES_nCATS"
    threshold_mC =  0.7 #  0.9 #0.99
    bam_path = "/home/michalula/data/cas9_nanopore/data/20241226_MR_nCATs_TcellsPrES_unedit_P2R9/passed_fast5/5mCG/to_t2t_v1_1/sort_align_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed.dna_r9.4.1_e8_sup@v3.3.5mCG.bam"

    date_today = datetime.today().strftime('%Y-%m-%d')

    ref_genome_v1_1_file = Path('/home/michalula/data/ref_genomes/to_t2t_v1_1/chm13.draft_v1.1.fasta')
    reg_genome_version = "t2t_v1_1"
    # t2t_v1_1_cd55_30bps = 'chr1:206586162-206586192'
    region_chr = 'chr1'

    # Expend window size
    expand_window_size = 16 # 50 # 50 #000 
    print("Expend window size by 2 * ", expand_window_size)
    region_start = 206586162 - expand_window_size
    region_end = 206586192 + expand_window_size + 1
    region_str = region_chr + ":" + str(region_start) + "-" + str(region_end) #'chr1:206586162-206586192'
    region_length = region_end - region_start
    print("region_length", region_length)


    save_padded_reads_name_np = f"padded_reads_{experiment_name}_mCthresh{threshold_mC}_{reg_genome_version}_{region_str}_{date_today}.npy"
    output_dir = create_output_directory("./dimelo_v2_output")

    motifs=['CG,0']
    ref_seq_list = get_reference_sequence(ref_genome_v1_1_file, region_chr, region_start, region_end)


    extract_file, extract_region = extract_from_bam(
        experiment_name=experiment_name,
        bam_path=bam_path,
        ref_genome_file=ref_genome_v1_1_file,
        output_dir=output_dir,
        region=region_str,
        motifs=motifs,
        output_name='extracted_reads',
        threshold_mC=threshold_mC,
    )

    keep_unmethylated_reads = False
    if extract_file:
        if keep_unmethylated_reads:
            mCG_reads_df, region_dict = process_extracted_reads(extract_file, region_str, motifs, ref_seq_list)
            visualize_data(mCG_reads_df)

            padded_reads = create_padded_reads(mCG_reads_df, region_dict, region_length)
        else:

            mCG_reads_df, region_dict = process_extracted_reads_no_fully_unmethylated(extract_file, region_str, motifs, ref_seq_list)
            visualize_data(mCG_reads_df)
            padded_reads = create_padded_reads_no_fully_unmethylated(mCG_reads_df, region_dict, region_length)
        
        if padded_reads is not None:
            plot_padded_reads(padded_reads, ref_seq_list)
            save_padded_reads(padded_reads, output_dir, save_padded_reads_name_np)
    
    
if __name__ == "__main__":
    main()
