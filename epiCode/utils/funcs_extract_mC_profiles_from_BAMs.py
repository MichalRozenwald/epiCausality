import platform
import sys
import pysam
from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from datetime import datetime
from dimelo import parse_bam, plot_reads, load_processed
import h5py


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
        # ref_seq = pysam.FastaFile(ref_genome_file).fetch(region_chr, region_start-1, region_end-1) # removing 1 from the coord indexes s the reference genome starts with 1 and the pysam.FastaFile is 0-based and the coordinates are the index-1
        ref_seq = pysam.FastaFile(ref_genome_file).fetch(region_chr, region_start, region_end) # removing 1 from the coord indexes s the reference genome starts with 1 and the pysam.FastaFile is 0-based and the coordinates are the index-1
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

def extract_from_bam(bam_path, ref_genome_file, output_dir, window_size=None, threshold_mC=0.99, 
                    num_cores=32, regions='chr1:206586162-206586192', motifs=['CG,0'], 
                    output_name='extract_output'):
    """Processes a BAM file using parse_bam.extract and plots the extracted reads."""
    try:
        extract_file, extract_regions = parse_bam.extract(
            input_file=bam_path,
            output_name=output_name,
            ref_genome=ref_genome_file,
            output_directory=output_dir,
            regions=regions,
            motifs=motifs,
            thresh=threshold_mC,
            window_size=window_size,
        )

        plot_reads.plot_reads(
            extract_file,
            regions,
            motifs=motifs,
            window_size=window_size,
            sort_by=['shuffle', 'strand'],
            s=1
        )

        plt.xlabel(f'bp relative to {regions}')
        plt.show()

        return extract_file, extract_regions
    except Exception as e:
        print("Error in BAM extraction:", e)
        return None, None

def process_extracted_reads_no_fully_unmethylated(extract_file, regions, motifs, ref_seq_list):
    """
    Process extracted reads into a DataFrame.

    Warning: make sure that the ref_seq_list was created using the same region and reference genome using the function: 
        motifs=['CG,0']
        ref_seq_list = get_reference_sequence(ref_genome_v1_1_file, region_chr, region_start, region_end)  
    """
    try:
        reads, read_names, mods, regions_dict = load_processed.readwise_binary_modification_arrays(
            file=extract_file,
            regions=regions,
            motifs=motifs
        )
        reads_df = pd.DataFrame({
            'read_name': read_names,
            'mod': mods,
            'pos': reads
        }).explode('pos')

        # reads_df['pos_shifted'] = reads_df['pos'] + 15
        region_length = len(ref_seq_list)
        reads_df['pos_shifted'] = reads_df['pos'] + (region_length // 2)
        return reads_df, regions_dict
    except Exception as e:
        print("Error processing extracted reads:", e)
        return None, None

import pandas as pd
import numpy as np
import h5py

def process_extracted_reads(extract_file, regions, motifs, ref_seq_list):
    """
    Process extracted reads into a DataFrame, ensuring all reads (methylated and unmethylated) are included.
    """
    try:
        # Extract methylation-modified positions
        mod_coords, read_ids, mods, regions_dict = load_processed.readwise_binary_modification_arrays(
            file=extract_file,
            regions=regions,
            motifs=motifs
        )

        # Get all read names (both methylated and unmethylated)
        with h5py.File(extract_file, "r") as h5:
            all_read_names = np.array(h5["read_name"], dtype=str)  # Extract all read names

        # Convert read IDs to strings to avoid type mismatches
        read_ids = np.array(read_ids, dtype=str)

        # Create a DataFrame for methylated reads
        reads_df = pd.DataFrame({
            'read_name': read_ids,
            'mod': mods,
            'pos': mod_coords
        })

        # Ensure 'pos' is converted to numeric type
        reads_df['pos'] = pd.to_numeric(reads_df['pos'], errors='coerce')

        # Identify unmethylated reads (present in BAM but missing from reads_df)
        methylated_reads = set(reads_df['read_name'])
        unmethylated_reads = [read for read in all_read_names if read not in methylated_reads]

        # Create a DataFrame for unmethylated reads (no positions)
        unmethylated_df = pd.DataFrame({
            'read_name': unmethylated_reads,
            'mod': None,
            'pos': np.nan  # No methylation site
        })

        # Combine both DataFrames
        reads_df = pd.concat([reads_df, unmethylated_df], ignore_index=True)

        # Ensure 'pos' is numeric for all rows
        reads_df['pos'] = pd.to_numeric(reads_df['pos'], errors='coerce')

        # Compute shifted positions
        region_length = len(ref_seq_list)
        reads_df['pos_shifted'] = reads_df['pos'].apply(
            lambda x: int(x + (region_length // 2)) if not np.isnan(x) else -1
        )

        return reads_df, regions_dict

    except Exception as e:
        print("Error processing extracted reads:", e)
        return None, None



def visualize_data_old(reads_df):
    """Generate visualizations for the data."""
    try:
        reads_df['read_name'].plot(kind='hist', bins=1600, title='Read Names Distribution')
        plt.gca().spines[['top', 'right']].set_visible(False)
        plt.show()

        sns.scatterplot(
            data=reads_df,
            x="pos",
            y="read_name",
            hue="mod",
            s=0.5,
            marker="s",
            linewidth=0
        )
        # plt.xticks(ticks=np.arange(len(ref_seq_list)), labels=ref_seq_list, size=font_size) # 'small') #, rotation=90)
        plt.xlabel('Position')
        plt.ylabel('Read Name')
        plt.show()
    except Exception as e:
        print("Error in visualization:", e)

def visualize_data(reads_df):
    """Generate visualizations for the data."""
    try:
        # Ensure 'pos' is numeric and drop NaNs
        reads_df = reads_df.copy()  # Avoid modifying the original DataFrame
        reads_df['pos'] = pd.to_numeric(reads_df['pos'], errors='coerce')
        reads_df = reads_df.dropna(subset=['pos'])  # Remove unmethylated reads for plotting

        # Histogram of read distribution
        reads_df['read_name'].value_counts().plot(kind='bar', title='Read Names Distribution')
        plt.gca().spines[['top', 'right']].set_visible(False)
        plt.show()

        # Scatter plot of modifications
        sns.scatterplot(
            data=reads_df,
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

def create_padded_reads_no_fully_unmethylated(reads_df, regions_dict, region_length):
    """Generate padded reads matrix."""
    try:
        read_names_unique = np.unique(reads_df['read_name'])
        num_reads = len(read_names_unique)
        reads_dict = {name: i for i, name in enumerate(read_names_unique)}
        padded_reads = np.full((num_reads, region_length), np.nan)

        for i in range(len(reads_df['read_name'])):
            padded_reads[reads_dict[reads_df['read_name'][i]], reads_df['pos_shifted'][i]] = 1

        return padded_reads
    except Exception as e:
        print("Error creating padded reads matrix:", e)
        return None

def create_padded_reads(reads_df, regions_dict, region_length):
    """Generate padded reads matrix, including reads with no methylation."""
    try:
        # Ensure 'pos_shifted' is numeric
        reads_df['pos_shifted'] = pd.to_numeric(reads_df['pos_shifted'], errors='coerce')

        read_names_unique = np.unique(reads_df['read_name'])
        num_reads = len(read_names_unique)
        reads_dict = {name: i for i, name in enumerate(read_names_unique)}
        padded_reads = np.full((num_reads, region_length), np.nan)

        for _, row in reads_df.iterrows():
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
    ref_seq_list = get_reference_sequence(ref_genome_v1_1_file, region_chr, region_start, region_end)


    extract_file, extract_regions = extract_from_bam(
        bam_path=bam_path,
        ref_genome_file=ref_genome_v1_1_file,
        output_dir=output_dir,
        regions=region_str,
        motifs=motifs,
        output_name='extracted_reads',
        threshold_mC=threshold_mC,
    )

    keep_unmethylated_reads = False
    if extract_file:
        if keep_unmethylated_reads:
            reads_df, regions_dict = process_extracted_reads(extract_file, region_str, motifs, ref_seq_list)
            visualize_data(reads_df)

            padded_reads = create_padded_reads(reads_df, regions_dict, region_length)
        else:

            reads_df, regions_dict = process_extracted_reads_no_fully_unmethylated(extract_file, region_str, motifs, ref_seq_list)
            visualize_data(reads_df)
            padded_reads = create_padded_reads_no_fully_unmethylated(reads_df, regions_dict, region_length)
        
        if padded_reads is not None:
            plot_padded_reads(padded_reads, ref_seq_list)
            save_padded_reads(padded_reads, output_dir, save_padded_reads_name_np)
    
    
if __name__ == "__main__":
    main()
