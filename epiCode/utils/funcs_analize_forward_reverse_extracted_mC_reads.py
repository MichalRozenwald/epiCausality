import os
import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from pathlib import Path
from datetime import datetime
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
from IPython.display import display

# !  python3 -m pip install tensorflow
# !  python3 -m pip install keras
# ! python3 -m pip install 'scikit-learn'
# ! python3 -m pip install shap
# ! python3 -m  pip uninstall -y plotly kaleido
# ! python3 -m  pip install --upgrade plotly kaleido

pio.renderers.default = "vscode"


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

def generate_CGs_all(padded_reads_df, ref_seq_list, region_chr, region_start):
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
        # Calculate the genomic coordinates of the CGs
        CG_coordinates = [(region_start + idx) for idx in CG_pair_idx]
        # Print the genomic coordinates with CG order number
        for order, (idx, coord) in enumerate(zip(CG_pair_idx, CG_coordinates), start=1):
            print(f"CG_{order} at index {idx} has genomic coordinate: {region_chr}:{coord}")

        # Create a DataFrame with the CG index, position in the region, chromosome, and coordinate
        CG_info_df = pd.DataFrame({
            'Position_in_region': CG_pair_idx,
            'Chromosome': [region_chr] * len(CG_pair_idx),
            'Coordinate': CG_coordinates
        })
        CG_info_df['CG_number'] = CG_info_df.index + 1
        print('CG_info_df', CG_info_df)

        C_reads_df = padded_reads_df.iloc[:, CG_pair_idx]
        G_reads_df = padded_reads_df.iloc[:, [i - 1 for i in CG_pair_idx]]

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
        return CGs_all, C_fwd_df, G_revs_df, CG_pair_idx, CG_coordinates, CG_info_df, fwd_reads_bools, rvs_reads_bools

    except Exception as e:
        print(f"Error generating CGs_all: {e}")
        raise

    
def plot_mC_sums_bar_aggregated_per_base(mC_sums, ref_seq_list, title=f"mC Sums Bar Plot"):
    """
    Creates an interactive bar plot of mC_sums using Plotly.

    Parameters:
    - mC_sums: List or array of mC sums values.
    - ref_seq_list: List of reference sequence labels (same length as mC_sums).

    Returns:
    - A Plotly figure.

    # Example usage
    ref_seq_list = ['A', 'T', 'C', 'G', 'A', 'C', 'T', 'G']
    mC_sums = np.random.randint(0, 100, len(ref_seq_list))

    fig = plot_mC_sums_bar(mC_sums, ref_seq_list)
    fig.show()
    """
    if len(mC_sums) != len(ref_seq_list):
        raise ValueError("mC_sums and ref_seq_list must have the same length.")
    if mC_sums is None or ref_seq_list is None:
        raise ValueError("mC_sums or ref_seq_list is None!")
    if len(mC_sums) == 0 or len(ref_seq_list) == 0:
        raise ValueError("mC_sums or ref_seq_list is empty!")
    if len(mC_sums) != len(ref_seq_list):
        raise ValueError("Mismatch: mC_sums and ref_seq_list lengths are different!")

    # Scale font size: decreases as seq_length increases, but stays within reasonable bounds
    font_size = max(2, min(8, 300 / len(ref_seq_list)))

    # Create DataFrame for Plotly
    df = pd.DataFrame({
        "Reference Sequence": ref_seq_list,
        "mC Sums": mC_sums
    })
    if df is None or df.empty:
        raise ValueError("Error: DataFrame is empty or None, cannot generate plot")
    else:
        print("In plot_mC_sums_bar DataFrame columns:", df.columns)
        print("In plot_mC_sums_bar DataFrame shape:", df.shape)
        # print("In plot_mC_sums_bar DataFrame:", df)


    # Create interactive bar plot using plotly.express
    fig = px.bar(
        df, 
        x="Reference Sequence", 
        y="mC Sums",
        text_auto=True
    )
    # fig.update_xaxes(tickangle= -90)  # Rotate x-axis labels for better readability
    fig.update_xaxes(tickangle= 90, font=dict(size=4))  # Rotate x-axis labels for better readability

    if fig is None:
        raise ValueError("❌ Plotly failed to create a figure. Check input data!")
    # else:
    #     print("✅ Figure created successfully:", type(fig))

    # Update layout for better readability
    fig.update_layout(
        title=title,
        xaxis_title="Reference Sequence",
        yaxis_title="mC Sums",
        font=dict(size=5),
    )

    if fig is None:
        raise ValueError("❌  Post fig.update_layout Plotly failed to create a figure. Check input data!")
    else:
        print("✅ Figure created successfully  Post fig.update_layout() :", type(fig))

    # fig.show()  # Ensure this is executed to render the plot
    # pio.show(fig)

    # Show the figure using multiple methods
    try:
        fig.show()
        # pio.show(fig)
        # display(fig)
    except Exception as e:
        print("❌ Error displaying figure:", e)
        fig.write_html("mC_sums_plot.html")
        print("✅ Plot saved as 'mC_sums_plot.html'. Open it in a browser.")

    return fig

    
def plot_mC_sums_bar(mC_sums, ref_seq_list, title="mC Sums Bar Plot", yaxis_title="mC Sums"):
    """
    Creates an interactive bar plot of mC_sums using Plotly.

    Parameters:
    - mC_sums: List or array of mC sums values.
    - ref_seq_list: List of reference sequence labels (same length as mC_sums).

    Returns:
    - A Plotly figure.

    # Example usage
    ref_seq_list = ['A', 'T', 'C', 'G', 'A', 'C', 'T', 'G']
    mC_sums = np.random.randint(0, 100, len(ref_seq_list))

    fig = plot_mC_sums_bar(mC_sums, ref_seq_list)
    fig.show()
    """
    if len(mC_sums) != len(ref_seq_list):
        raise ValueError("mC_sums and ref_seq_list must have the same length.")
    if mC_sums is None or ref_seq_list is None:
        raise ValueError("mC_sums or ref_seq_list is None!")
    if len(mC_sums) == 0 or len(ref_seq_list) == 0:
        raise ValueError("mC_sums or ref_seq_list is empty!")
    if len(mC_sums) != len(ref_seq_list):
        raise ValueError("Mismatch: mC_sums and ref_seq_list lengths are different!")

    # Scale font size: decreases as seq_length increases, but stays within reasonable bounds
    font_size = max(2, min(8, 500 / len(ref_seq_list)))

    # Create DataFrame and add an index column
    df = pd.DataFrame({
        "Index": range(len(ref_seq_list)),  # Ensure sequence order
        "Reference Sequence": ref_seq_list,
        "mC Sums": mC_sums
    })

    if df is None or df.empty:
        raise ValueError("Error: DataFrame is empty or None, cannot generate plot")
    else:
        print("In plot_mC_sums_bar DataFrame columns:", df.columns)
        print("In plot_mC_sums_bar DataFrame shape:", df.shape)
        # print("In plot_mC_sums_bar DataFrame:", df)


    # Create bar plot with sequential x-axis
    fig = px.bar(df, 
                 x="Index",  # Use index instead of 'Reference Sequence' directly
                 y="mC Sums",
                 text_auto=True,
                #  labels={"Index": "Reference Sequence"}
                 )

    # # Update x-axis to show sequence labels but preserve order
    # fig.update_layout(
    #     title=title,
    #     xaxis=dict(
    #         tickmode="array",
    #         tickvals=df["Index"],
    #         ticktext=df["Reference Sequence"],  # Show original sequence
    #     ),
    #     yaxis_title=yaxis_title
    # )
    # fig.update_xaxes(tickangle= 90, font=dict(size=50))  # Rotate x-axis labels for better readability
        # Update x-axis to show sequence labels without grouping
    fig.update_layout(
        title=title,
        xaxis=dict(
            tickmode="array",
            tickvals=df["Index"],  # Set tick positions
            ticktext=df["Reference Sequence"],  # Show original sequence labels
        ),
        yaxis_title="mC Sums"
    )
    # ✅ Rotate x-axis labels and set font size
    fig.update_xaxes(tickangle=0, tickfont=dict(size=5))  # Adjust size as needed
    if fig is None:
        raise ValueError("❌ Plotly failed to create a figure. Check input data!")
    # else:
    #     print("✅ Figure created successfully:", type(fig))
    if fig is None:
        raise ValueError("❌  Post fig.update_layout Plotly failed to create a figure. Check input data!")
    else:
        print("✅ Figure created successfully  Post fig.update_layout() :", type(fig))
    # fig.show()  # Ensure this is executed to render the plot
    # pio.show(fig)

    # Show the figure using multiple methods
    try:
        fig.show()
        # pio.show(fig)
        # display(fig)
    except Exception as e:
        print("❌ Error displaying figure:", e)
        fig.write_html("mC_sums_plot.html")
        print("✅ Plot saved as 'mC_sums_plot.html'. Open it in a browser.")

    return fig


def plot_mCG_bars(CGs_all, CG_pair_idx, ref_seq_list, experiment_name):
    """
    Creates an interactive bar plot of mCG_sums using Plotly.

    Parameters:
    - mCG_sums: List or array of mCG sums values.
    - ref_seq_list: List of reference sequence labels (same length as mCG_sums).

    Returns:
    - A Plotly figure.
    """

    # Compute sums and fractions
    CGs_all_sums = np.nansum(CGs_all.values, axis=0)
    CGs_all_on_fwd_C_sums = np.zeros(len(ref_seq_list))
    CGs_all_on_fwd_C_sums[CG_pair_idx] = CGs_all_sums
    mC_fracs = CGs_all_sums / len(CGs_all)
    print("CGs_all_sums  =", CGs_all_sums)
    print("CGs_all_sums / len(CGs_all) =", CGs_all_sums / len(CGs_all))

    print("len(CGs_all) =", len(CGs_all))
    print("CGs_all.shape =", CGs_all.shape)
    print("len(CGs_all_on_fwd_C_sums) =", len(CGs_all_on_fwd_C_sums))
    print("CGs_all_on_fwd_C_sums.shape =", CGs_all_on_fwd_C_sums.shape)
    print("CGs_all_on_fwd_C_sums =", CGs_all_on_fwd_C_sums)
    print("CGs_all_on_fwd_C_sums / len(CGs_all) =", CGs_all_on_fwd_C_sums / len(CGs_all))

    # Scale font size: decreases as seq_length increases, but within reasonable bounds
    font_size = max(2, min(8, 500 / len(ref_seq_list)))  # Now it stays between 2 and 8

    plt.figure(figsize=(10, 5))
    plt.bar(np.arange(len(CGs_all_sums)), CGs_all_sums, snap=False)
    plt.title(f"{experiment_name}\nSum of mCs of CpG units (fwd + rvs)")
    plt.show()

    plt.figure(figsize=(10, 5))
    plt.bar(np.arange(len(mC_fracs)), mC_fracs, snap=False)
    if len(ref_seq_list) < 160:      
        plt.xticks(ticks=np.arange(CGs_all.shape[1]), labels=CGs_all.columns) #, size=font_size) # 'small') #, rotation=90)
    plt.title(f"{experiment_name}\nFractions of mC [mC_sums / num_reads], num_reads= {len(CGs_all)}  of CpG units  (fwd + rvs)")
    plt.show()

    plt.figure(figsize=(10, 5))
    plt.bar(np.arange(len(CGs_all_on_fwd_C_sums)), CGs_all_on_fwd_C_sums, snap=False) # , width=0.0001
    if len(ref_seq_list) < 160:       
        plt.xticks(ticks=np.arange(len(ref_seq_list)), labels=ref_seq_list, size=font_size) # 'small') #, rotation=90)
    plt.title(f"{experiment_name}\n Total sum of mCs (fwd + rvs) bap plot with reference seq")
    plt.show()

    mC_fracs = CGs_all_on_fwd_C_sums / len(CGs_all)
    plt.figure(figsize=(10, 5))
    plt.bar(np.arange(len(mC_fracs)), mC_fracs, snap=False) # , width=0.0001
    if len(ref_seq_list) < 160:               
        plt.xticks(ticks=np.arange(len(ref_seq_list)), labels=ref_seq_list, size=font_size) # 'small') #, rotation=90)
    plt.title(f"{experiment_name}\n Fractions of mC [mC_sums / num_reads]  with reference seq,  num_reads= {len(CGs_all)}")
    plt.show()

    if len(ref_seq_list) < 160:  
        # Visualize mC sums as bars from plotly
        plot_mC_sums_bar(CGs_all_sums, CGs_all.columns, title=f"{experiment_name}\nSum of mCs of CpG units (fwd + rvs)", yaxis_title="mC Sums (fwd + rvs)")
        # Visualize mC sums as bars from plotly
        plot_mC_sums_bar(CGs_all_sums / len(CGs_all), CGs_all.columns,
                    title=f"{experiment_name}<br> Fractions of mC [mC_sums / num_reads], num_reads= {len(CGs_all)}  of CpG units  (fwd + rvs)",
                    yaxis_title="mC Fractions (fwd + rvs)")
        
        # Visualize mC sums as bars from plotly
        plot_mC_sums_bar(CGs_all_on_fwd_C_sums, ref_seq_list, 
                    title=f"{experiment_name}<br> Total sum of mCs (fwd + rvs) bap plot with reference seq,  num_reads= {len(CGs_all)}", yaxis_title="mC Sums")
        # # Visualize mC sums as bars from plotly
        plot_mC_sums_bar(CGs_all_on_fwd_C_sums / len(CGs_all), ref_seq_list, 
                    title=f"{experiment_name}<br> Fractions of mC [mC_sums / num_reads] with reference seq, num_reads= {len(CGs_all)}", yaxis_title="mC Fractions")




def visualize_CGs_all(padded_reads_df, CGs_all, C_fwd_df, G_revs_df, CG_pair_idx, ref_seq_list, fwd_reads_bools, rvs_reads_bools, experiment_name):
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

    # Scale font size: decreases as seq_length increases, but within reasonable bounds
    font_size = max(2, min(8, 500 / len(ref_seq_list)))  # Now it stays between 2 and 8
    try:
        read_sums = np.nansum(padded_reads_df.values, axis=1)
        mC_sums = np.nansum(padded_reads_df.values, axis=0)
        # Print basic statistics
        print("DataFrame shape:", padded_reads_df.shape)
        print(padded_reads_df.describe())
        print("Read sums:", read_sums)
        print("Zero reads:", sum(read_sums == 0), ", Non-zero reads:", sum(read_sums != 0))
        print("mC sums = ", mC_sums)

        plt.figure(figsize=(10, 5))
        plt.bar(np.arange(len(mC_sums)), mC_sums, snap=False)
        if len(ref_seq_list) < 160:       
            # plt.xticks(range(len(ref_seq_list)), ref_seq_list, size='small')
            plt.xticks(ticks=np.arange(len(ref_seq_list)), labels=ref_seq_list, size=font_size) # 'small') #, rotation=90)
        plt.title(f"{experiment_name} mC Sums Bar Plot")
        plt.show()

        print("(mC_sums / len(CGs_all) =", mC_sums/ len(CGs_all))
        plt.figure(figsize=(10, 5))
        # mC_fracs = CGs_all_sums / len(CGs_all) 
        plt.scatter(np.arange(len(mC_sums)), mC_sums / len(CGs_all)) 
        if len(ref_seq_list) < 160:       
            plt.xticks(ticks=np.arange(len(ref_seq_list)), labels=ref_seq_list, size=font_size) # 'small') #, rotation=90)
        plt.title(f"{experiment_name}\n mC Fractions Scatter Plot  [mC_sums / num_reads], num_reads= {len(CGs_all)}")
        plt.show()

        if len(ref_seq_list) < 160: 
                # Visualize mC sums as bar plotly
                plot_mC_sums_bar(mC_sums, ref_seq_list, title=f"{experiment_name}<br> mC Sums Bar Plot (Plotly)", yaxis_title="mC Sums")
                # Visualize mC sums as bar plotly
                plot_mC_sums_bar(mC_sums/ len(CGs_all), ref_seq_list,
                            title=f"{experiment_name}<br> mC Fractions Scatter Plot  [mC_sums / num_reads], num_reads= {len(CGs_all)}", yaxis_title="mC Fractions")
                
        # Visualize mC sums as bars
        plot_mCG_bars(CGs_all, CG_pair_idx, ref_seq_list, experiment_name)


        # Cluster maps for filtered DataFrames
        sns.clustermap(C_fwd_df.fillna(-1), col_cluster=False)
        plt.title(f"{experiment_name}\nFiltered Forward Reads ClusterMap")
        plt.show()

        sns.clustermap(G_revs_df.fillna(-1), col_cluster=False)
        plt.title(f"{experiment_name}\nFiltered Reverse Reads ClusterMap")
        plt.show()

        # Heatmap of CGs_all
        sns.heatmap(CGs_all.fillna(-1))
        plt.title(f"{experiment_name}\nConcatinated CGs_all (Forward and Reverse) (bright = mC, dark = C)\n(Fwd: {sum(fwd_reads_bools)}, Rev: {sum(rvs_reads_bools)})")
        plt.show()

        # Clustered Heatmap of CGs_all
        sns.clustermap(CGs_all.fillna(-1), col_cluster=False)
        plt.title(f"{experiment_name}\nClustered CGs_all (Forward and Reverse) (bright = mC, dark = C) ClusterMap\n(Fwd: {sum(fwd_reads_bools)}, Rev: {sum(rvs_reads_bools)})")
        plt.show()

    except Exception as e:
        print(f"Error visualizing CGs_all: {e}")
        raise

def save_CGs_all(CGs_all, save_folder_path, save_cpg_name_np):
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


def save_CGs_coords_info(CG_info_df, save_folder_path, save_cpg_info_name_df):
    """
    Save the CGs_all DataFrame to a file.

    Parameters:
        CGs_all (pd.DataFrame): CGs_all DataFrame to save.
        save_folder_path (str): Path to save processed data.
    """
    try: 
        CG_info_df.to_csv(Path(save_folder_path, save_cpg_info_name_df))
        print(f"CGs_all saved as {save_cpg_info_name_df} in {save_folder_path}")
    except Exception as e:
        print(f"Error saving CGs_all: {e}")
        raise


def analize_forward_reverse_CGs_pipeline(experiment_name, save_folder_path, 
        save_padded_reads_name_np, ref_genome_file, region_chr, region_start, region_end, 
        do_save_CGs_coords_info=False, save_cpg_info_name="CG_info_df"):
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
        CGs_all, C_fwd_df, G_revs_df, CG_pair_idx, CG_coordinates, CG_info_df, \
            fwd_reads_bools, rvs_reads_bools = generate_CGs_all(padded_reads_df, ref_seq_list, region_chr, region_start)

        # Visualize CGs_all
        visualize_CGs_all(padded_reads_df, CGs_all, C_fwd_df, G_revs_df, CG_pair_idx, \
            ref_seq_list, fwd_reads_bools, rvs_reads_bools, experiment_name)

        # File name generation
        fwd_count = sum(fwd_reads_bools)
        rvs_count = sum(rvs_reads_bools)
        date_today = datetime.today().strftime('%Y-%m-%d')
        save_cpg_name_np = f"CG_{CGs_all.shape[1]}_{save_padded_reads_name_np[:-4]}_units_combined_numFWD{fwd_count}_numRVS{rvs_count}.npy" #_{date_today}.npy"

        # Save CGs_all
        save_CGs_all(CGs_all, save_folder_path, save_cpg_name_np) 
        save_cpg_info_name_df= save_cpg_info_name + f"_{experiment_name}_numFWD{fwd_count}_numRVS{rvs_count}_{date_today}.csv"
        if do_save_CGs_coords_info:
            save_CGs_coords_info(CG_info_df, save_folder_path, save_cpg_info_name_df)

        return CGs_all, C_fwd_df, G_revs_df, padded_reads_df, CG_pair_idx, CG_coordinates, CG_info_df

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
        threshold_mC =  0.7 #  0.9 #0.99
        bam_path = "/home/michalula/data/cas9_nanopore/data/20241226_MR_nCATs_TcellsPrES_unedit_P2R9/passed_fast5/5mCG/to_t2t_v1_1/sort_align_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed.dna_r9.4.1_e8_sup@v3.3.5mCG.bam"

        date_today = datetime.today().strftime('%Y-%m-%d')

        ref_genome_path = Path('/home/michalula/data/ref_genomes/to_t2t_v1_1/chm13.draft_v1.1.fasta')
        reg_genome_version = "t2t_v1_1"
        # t2t_v1_1_cd55_30bps = 'chr1:206586162-206586192'
        region_chr = 'chr1'

        # Expend window size
        expand_window_size = 16 # 0 # 50 # 50 #000
        expand_window_size
        print("Expend window size by 2 * ", expand_window_size)
        region_start = 206586162 - expand_window_size
        region_end = 206586192 + expand_window_size + 1
        region_str = region_chr + ":" + str(region_start) + "-" + str(region_end) #'chr1:206586162-206586192'
        region_length = region_end - region_start
        print("region_length", region_length)


        save_padded_reads_name_np = f"padded_reads_{experiment_name}_mCthresh{threshold_mC}_{reg_genome_version}_{region_str}_{date_today}.npy"
        output_dir = create_output_directory("./dimelo_v2_output")

        motifs=['CG,0']
        ref_seq_list = get_reference_sequence(ref_genome_path, region_chr, region_start, region_end)


        # Process pipeline
        CGs_all, C_fwd_df, G_revs_df, padded_reads_df, CG_pair_idx, CG_coordinates, CG_info_df = analize_forward_reverse_CGs_pipeline(
            experiment_name=experiment_name, save_folder_path=output_dir, 
            save_padded_reads_name_np=save_padded_reads_name_np, 
            ref_genome_file=ref_genome_path, region_chr=region_chr, region_start=region_start, region_end=region_end
        )

        print("Pipeline executed successfully (analize_forward_reverse_CGs_pipeline function)")
        return CGs_all, C_fwd_df, G_revs_df, padded_reads_df, CG_pair_idx, CG_coordinates, CG_info_df

    except Exception as e:
        print(f"Error in main pipeline (analize_forward_reverse_CGs_pipeline function): {e}") 
        
if __name__ == "__main__":
    CGs_all, C_fwd_df, G_revs_df, padded_reads_df = main()
