import platform
import sys
import pysam
from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from dimelo import parse_bam, plot_reads, load_processed

def system_info():
    """Print system information."""
    print('System:', platform.system())
    print('Release:', platform.release())
    print('Version:', platform.version())
    print('Processor:', platform.processor())
    print('Python version:', sys.version)

def get_reference_sequence(ref_genome_file, region):
    """Fetch reference sequence from genome file."""
    try:
        ref_seq = pysam.FastaFile(ref_genome_file).fetch(*region)
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

def process_extracted_reads(extract_file, regions, motifs):
    """Process extracted reads into a DataFrame."""
    try:
        reads, read_names, mods, regions_dict = load_processed.readwise_binary_modification_arrays(
            file=extract_file,
            regions=regions,
            motifs=motifs
        )
        df = pd.DataFrame({
            'read_name': read_names,
            'mod': mods,
            'pos': reads
        }).explode('pos')

        df['pos_shifted'] = df['pos'] + 15
        return df, regions_dict
    except Exception as e:
        print("Error processing extracted reads:", e)
        return None, None

def visualize_data(df):
    """Generate visualizations for the data."""
    try:
        df['read_name'].plot(kind='hist', bins=1600, title='Read Names Distribution')
        plt.gca().spines[['top', 'right']].set_visible(False)
        plt.show()

        sns.scatterplot(
            data=df,
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

def create_padded_reads(df, regions_dict):
    """Generate padded reads matrix."""
    try:
        read_names_unique = np.unique(df['read_name'])
        num_reads = len(read_names_unique)
        reads_dict = {name: i for i, name in enumerate(read_names_unique)}
        padded_reads = np.full((num_reads, 30), np.nan)

        for i in range(len(df['read_name'])):
            padded_reads[reads_dict[df['read_name'][i]], df['pos_shifted'][i]] = 1

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
        plt.xticks(ticks=np.arange(len(ref_seq_list)), labels=ref_seq_list, size='small') #, rotation=90)
        # plt.xticks(range(len(seq_list)), seq_list, size='small')
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

    ref_genome_v1_1_file = Path('/home/michalula/data/ref_genomes/to_t2t_v1_1/chm13.draft_v1.1.fasta')
    region = ('chr1', 206586162, 206586192)
    # region_str = 'chr1:206586162-206586192'
    t2t_v1_1_cd55_30bps = 'chr1:206586162-206586192'
    motifs=['CG,0']

    ref_seq_list = get_reference_sequence(ref_genome_v1_1_file, region)

    output_dir = create_output_directory("./dimelo_v2_output")
    
    unedited_bam_path = "/home/michalula/data/cas9_nanopore/data/20241226_MR_nCATs_TcellsPrES_unedit_P2R9/passed_fast5/5mCG/to_t2t_v1_1/sort_align_trim_20241226_MR_nCATs_TcellsPrES_unedit_P2R9_passed.dna_r9.4.1_e8_sup@v3.3.5mCG.bam"

    extract_file, extract_regions = extract_from_bam(
        bam_path=unedited_bam_path,
        ref_genome_file=ref_genome_v1_1_file,
        output_dir=output_dir,
        regions=t2t_v1_1_cd55_30bps,
        motifs=motifs,
        output_name='extracted_reads'
    )

    if extract_file:
        df, regions_dict = process_extracted_reads(extract_file, t2t_v1_1_cd55_30bps, motifs)
        visualize_data(df)

        padded_reads = create_padded_reads(df, regions_dict)
        if padded_reads is not None:
            plot_padded_reads(padded_reads, ref_seq_list)
            save_padded_reads(padded_reads, output_dir, "padded_reads.npy")

if __name__ == "__main__":
    main()
