import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from pathlib import Path
import sys
from datetime import datetime
import plotly.graph_objects as go

# Add custom function paths
sys.path.append("/home/michalula/code/epiCausality/epiCode/utils/")

from funcs_extract_mC_profiles_from_BAMs import (
    get_reference_sequence
)

def load_padded_reads(save_folder_path, save_name_np):
    """
    Load padded reads from a .npy file.

    Parameters:
        save_folder_path (str): Path to the save folder.
        save_name_np (str): Name of the .npy file.

    Returns:
        np.ndarray: Padded reads array.
    """
    try:
        file_path = Path(save_folder_path, save_name_np)
        padded_reads = np.load(file_path)
        print("Padded reads loaded successfully.")
        return padded_reads
    except Exception as e:
        print(f"Error loading padded reads: {e}")
        raise

def visualize_padded_reads(padded_reads):
    """
    Visualize the padded reads as a matrix plot.

    Parameters:
        padded_reads (np.ndarray): Array of padded reads.
    """
    try:
        plt.figure(figsize=(10, 50))
        plt.matshow(padded_reads, fignum=1)
        plt.title("Padded Reads Visualization")
        plt.show()
        print("Visualization completed.")
    except Exception as e:
        print(f"Error visualizing padded reads: {e}")
        raise

def process_region(ref_genome_file, region_chr, region_start, region_end):
    """
    Process a genomic region and retrieve reference sequences.

    Parameters:
        ref_genome_file (str): Path to the reference genome file.
        region_chr (str): Chromosome name.
        region_start (int): Start position of the region.
        region_end (int): End position of the region.

    Returns:
        list: Reference sequence for the specified region.
    """
    try:
        region_length = region_end - region_start
        print("Region length:", region_length)
        ref_seq_list = get_reference_sequence(ref_genome_file, region_chr, region_start, region_end)
        print("Reference sequence retrieved.")
        return ref_seq_list
    except Exception as e:
        print(f"Error processing region: {e}")
        raise

def generate_dataframe(padded_reads, ref_seq_list):
    """
    Generate a DataFrame from padded reads and reference sequences.

    Parameters:
        padded_reads (np.ndarray): Array of padded reads.
        ref_seq_list (list): Reference sequence list.

    Returns:
        pd.DataFrame: DataFrame with reference sequence as columns.
    """
    try:
        if padded_reads.shape[1] != len(ref_seq_list):
            raise ValueError(f"Length mismatch: padded_reads has {padded_reads.shape[1]} columns, but ref_seq_list has {len(ref_seq_list)} elements.")
        
        padded_reads_df = pd.DataFrame(padded_reads, columns=ref_seq_list)
        print("DataFrame created successfully.")
        return padded_reads_df
    except Exception as e:
        print(f"Error generating DataFrame: {e}")
        raise

def generate_cgs_all(padded_reads_df, ref_seq_list):
    """
    Generate the CGs_all DataFrame along with related DataFrames.

    Parameters:
        padded_reads_df (pd.DataFrame): DataFrame to analyze.
        ref_seq_list (list): Reference sequence list.

    Returns:
        tuple: CGs_all DataFrame, C_fwd_df, G_revs_df, CG_pair_idx.
    """
    try:
        seq_str = ''.join(ref_seq_list)
        CG_pair_idx = [i for i in range(len(seq_str) - 1) if seq_str[i] == 'C' and seq_str[i + 1] == 'G']
        print("CG Pair Indices:", CG_pair_idx)

        C_reads_df = padded_reads_df.iloc[:, CG_pair_idx]
        G_reads_df = padded_reads_df.iloc[:, [i + 1 for i in CG_pair_idx]]

        fwd_reads_bools = C_reads_df.sum(axis=1) != 0
        rvs_reads_bools = G_reads_df.sum(axis=1) != 0

        print("Forward reads:", sum(fwd_reads_bools))
        print("Reverse reads:", sum(rvs_reads_bools))

        C_fwd_df = C_reads_df[fwd_reads_bools]
        G_revs_df = G_reads_df[rvs_reads_bools]

        # Update column names
        C_fwd_df.columns = [f"C_{i+1}" for i in range(C_fwd_df.shape[1])]
        G_revs_df.columns = [f"G_{i+1}" for i in range(G_revs_df.shape[1])]

        CGs_all = pd.DataFrame(
            np.concatenate([np.array(C_fwd_df), np.array(G_revs_df)], axis=0),
            columns=[f"CG_{i+1}" for i in range(np.array(C_fwd_df).shape[1])]
        )
        return CGs_all, C_fwd_df, G_revs_df, CG_pair_idx, fwd_reads_bools, rvs_reads_bools

    except Exception as e:
        print(f"Error generating CGs_all: {e}")
        raise


def visualize_cgs_all(padded_reads_df, CGs_all, C_fwd_df, G_revs_df, CG_pair_idx, ref_seq_list, fwd_reads_bools, rvs_reads_bools):
    """
    Visualize the CGs_all DataFrame and related plots.

    Parameters:
        CGs_all (pd.DataFrame): CGs_all DataFrame.
        C_fwd_df (pd.DataFrame): Forward reads DataFrame.
        G_revs_df (pd.DataFrame): Reverse reads DataFrame.
        CG_pair_idx (list): CG pair indices.
        ref_seq_list (list): Reference sequence list.
        fwd_reads_bools (np.ndarray): Boolean array of forward reads.
        rvs_reads_bools (np.ndarray): Boolean array of reverse reads.
        padded_reads_df (pd.DataFrame): DataFrame of padded reads.
    """
    try:
        # Basic statistics
        print("DataFrame shape:", padded_reads_df.shape)
        print(padded_reads_df.describe())

        read_sums = np.nansum(padded_reads_df.values, axis=1)
        print("Read sums:", read_sums)
        print("Zero reads:", sum(read_sums == 0), ", Non-zero reads:", sum(read_sums != 0))

        mC_sums = np.nansum(padded_reads_df.values, axis=0)
        print("mC sums:", mC_sums)

        # Visualize mC sums as bar plot
        plt.figure(figsize=(10, 5))
        plt.bar(np.arange(len(mC_sums)), mC_sums)
        plt.xticks(range(len(ref_seq_list)), ref_seq_list, size='small')
        plt.title("mC Sums Bar Plot")
        plt.show()

        plt.figure(figsize=(10, 5))
        plt.scatter(np.arange(len(mC_sums)), mC_sums)        
        # plt.xticks(range(len(ref_seq_list)), ref_seq_list, size='small')
        plt.title("mC Sums Scatter Plot")
        plt.show()

        # Heatmap of CGs_all
        sns.heatmap(CGs_all.fillna(-1))
        plt.title(f"Concatinated CGs_all (Forward and Reverse) (bright = mC, dark = C)\n(Fwd: {sum(fwd_reads_bools)}, Rev: {sum(rvs_reads_bools)})")
        plt.show()

        # Clustered Heatmap of CGs_all
        sns.clustermap(CGs_all.fillna(-1), col_cluster=False)
        plt.title(f"Clustered CGs_all (Forward and Reverse) (bright = mC, dark = C) ClusterMap\n(Fwd: {sum(fwd_reads_bools)}, Rev: {sum(rvs_reads_bools)})")
        plt.show()

        # Compute sums and fractions
        CGs_all_sums = np.nansum(CGs_all.values, axis=0)
        CGs_all_on_fwd_C_sums = np.zeros(len(ref_seq_list))
        CGs_all_on_fwd_C_sums[CG_pair_idx] = CGs_all_sums

        plt.bar(np.arange(len(CGs_all_sums)), CGs_all_sums)
        plt.title("Sum of mCs (fwd + rvs)")
        plt.show()

        plt.bar(np.arange(len(CGs_all_on_fwd_C_sums)), CGs_all_on_fwd_C_sums)
        plt.xticks(range(len(ref_seq_list)), ref_seq_list, size='small')
        plt.title("Total sum of mCs (fwd + rvs) in the T cells Cas9 data")
        plt.show()

        mC_fracs = CGs_all_on_fwd_C_sums / len(CGs_all)
        plt.bar(np.arange(len(mC_fracs)), mC_fracs)
        plt.xticks(range(len(ref_seq_list)), ref_seq_list, size='small')
        plt.title("Fractions of mC [mC_sums / num_reads]")
        plt.show()

        # Cluster maps for filtered DataFrames
        sns.clustermap(C_fwd_df.fillna(-1), col_cluster=False)
        plt.title("Filtered Forward Reads ClusterMap")
        plt.show()

        sns.clustermap(G_revs_df.fillna(-1), col_cluster=False)
        plt.title("Filtered Reverse Reads ClusterMap")
        plt.show()

    except Exception as e:
        print(f"Error visualizing CGs_all: {e}")
        raise

def save_cgs_all(CGs_all, save_folder_path, save_cpg_name_np):
    """
    Save the CGs_all DataFrame to a file.

    Parameters:
        CGs_all (pd.DataFrame): CGs_all DataFrame to save.
        save_folder_path (str): Path to save processed data.
    """
    try: 
        np.save(Path(save_folder_path, save_cpg_name_np), CGs_all)
        print(f"CGs_all saved as {save_cpg_name_np} in {save_folder_path}")
    except Exception as e:
        print(f"Error saving CGs_all: {e}")
        raise


def analize_forward_reverse_CGs_pipeline(experiment_name, save_folder_path, save_padded_reads_name_np, ref_genome_file, region_chr, region_start, region_end):
    """
    Process the pipeline with the given constants.

    Parameters:
        experiment_name (str): Name of the experiment.
        save_folder_path (str): Path to save folder.
        save_padded_reads_name_np (str): Name of the padded reads file.
        ref_genome_file (str): Path to reference genome file.
        region_chr (str): Chromosome name.
        region_start (int): Start position of the region.
        region_end (int): End position of the region.

    Returns:
        tuple: CGs_all, C_fwd_df, G_revs_df, and padded_reads_df.
    """
    try:
        # Load padded reads
        padded_reads = load_padded_reads(save_folder_path, save_padded_reads_name_np)

        # Visualize padded reads
        visualize_padded_reads(padded_reads)

        # Process region
        ref_seq_list = process_region(ref_genome_file, region_chr, region_start, region_end)

        # Generate DataFrame
        padded_reads_df = generate_dataframe(padded_reads, ref_seq_list)

        # Generate CGs_all and related DataFrames
        CGs_all, C_fwd_df, G_revs_df, CG_pair_idx, fwd_reads_bools, rvs_reads_bools = generate_cgs_all(padded_reads_df, ref_seq_list)

        # Visualize CGs_all
        visualize_cgs_all(padded_reads_df, CGs_all, C_fwd_df, G_revs_df, CG_pair_idx, ref_seq_list, fwd_reads_bools, rvs_reads_bools)

        # File name generation
        fwd_count = sum(fwd_reads_bools)
        rvs_count = sum(rvs_reads_bools)
        date_today = datetime.today().strftime('%Y-%m-%d')
        save_cpg_name_np = f"CG_combined_{experiment_name}_numFWD{fwd_count}_numRVS{rvs_count}_{save_padded_reads_name_np}"#_{date_today}.npy"

        # Save CGs_all
        save_cgs_all(CGs_all, save_folder_path, save_cpg_name_np) 

        return CGs_all, C_fwd_df, G_revs_df, padded_reads_df

    except Exception as e:
        print(f"Error in process pipeline: {e}")
        raise


def main():
    """
    Main pipeline for processing and visualizing CpG units on forward and reverse strands 
    """
    try:
        # Define constants
        experiment_name = "unedited_T_primerES_nCATS"
        save_folder_path = "/home/michalula/code/epiCausality/epiCode/notebooks/dimelo_v2_output"
        save_padded_reads_name_np = "padded_reads.npy"
        ref_genome_file = "/home/michalula/data/ref_genomes/to_t2t_v1_1/chm13.draft_v1.1.fasta"
        region_chr = "chr1"
        region_start = 206586162
        region_end = 206586192

        # Process pipeline
        CGs_all, C_fwd_df, G_revs_df, padded_reads_df = analize_forward_reverse_CGs_pipeline(
            experiment_name, save_folder_path, save_padded_reads_name_np, 
            ref_genome_file, region_chr, region_start, region_end
        )

        print("Pipeline executed successfully (analize_forward_reverse_CGs_pipeline function)")
        return CGs_all, C_fwd_df, G_revs_df, padded_reads_df

    except Exception as e:
        print(f"Error in main pipeline (analize_forward_reverse_CGs_pipeline function): {e}")

if __name__ == "__main__":
    CGs_all, C_fwd_df, G_revs_df, padded_reads_df = main()
