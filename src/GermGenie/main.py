import argparse
import os
import glob
import subprocess
from typing import Literal
from GermGenie.pipeline import EMU
from plotly import express as px
from plotly import graph_objects as go
import gzip

import pandas as pd
from Bio import SeqIO

def get_name(fastq: str) -> str:
    """Get the name of a file"""
    return os.path.basename(fastq).split(".")[0]

def create_output_dirs(output_dir: str, subsample: bool = False) -> str:
    """Create output directories

    Args:
        output_dir (str): Base output directory
        subsample (bool, optional): Create subsample directory. Defaults to False.

    Returns:
        str: Path to emu directory
    """
    emu_dir: str = os.path.join(output_dir, "emu")
    
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
        os.mkdir(emu_dir)
        if subsample:
            os.mkdir(os.path.join(output_dir, "subsample"))
    
    return emu_dir

def find_input_files(input_dir: str) -> list[str]:
    """Find input fastq.gz files in a directory

    Args:
        input_dir (str): Path to directory containing fastq.gz files

    Returns:
        list[str]: List of fastq.gz files
    """
    infiles: list[str] = glob.glob(os.path.join(input_dir, "*.gz"))    
    if len(infiles) < 1:
        print("No input files found...") # TODO: replace with logging
        os.abort()
    else:
        return infiles
    
def subsample_fastq(nreads: int, fastq: str, outdir: str) -> str|None:
    """Subsample a fastq file

    Args:
        nreads (int): Number of reads to subsample
        fastq (str): Path to fastq file
        outdir (str): Path to output directory

    Returns:
        str|None: Path to subsampled fastq file, or None if an error occurred
    """
    name: str = get_name(fastq)
    out: str = os.path.join(outdir, f"{name}.fastq.gz")
        
    try:
        subprocess.run(['gunzip', '-c', fastq, '|', 'seqtk', 'sample', '-', nreads, '|', 'gzip', '>', out], shell=True)
        return out
    except subprocess.CalledProcessError as e:
        print(f"{e}") # TODO: replace with logging
        return None    

def run_emu(fastq: str, db: str, threads: int, output_dir: str) -> str|None:
    """Run EMU on a fastq file

    Args:
        fastq (str): Path to fastq file
        db (str): Path to EMU compatible database
        threads (int): Number of threads to use
        output_dir (str): Path to output directory

    Returns:
        str|None: stdout from EMU, or None if an error occurred
    """
    name: str = get_name(fastq)
    
    try:
        stdout = subprocess.check_output(
        ["emu", "abundance", fastq, '--db', db, '--threads', threads, '--output-dir', output_dir, '--output-basename', name],
        shell=True, text=True, stderr=subprocess.PIPE
    )
    except subprocess.CalledProcessError as e:
        print(f'Error processing {name}') # TODO: replace with logging
        print(e)
        print(e.stderr)
        return None
    
    return stdout
        

class ReadAssignment:
    """Class for storing read assignments from emu stdout"""	
    def __init__(self):
        self.sample: list[str]
        self.assigned: list[int]
        self.unassigned: list[int]
    
    def add(self, filename, stdout: str) -> None:
        """Add a read assignment from emu stdout

        Args:
            sample (str): Name of the sample
            assigned (str): Number of assigned reads
            unassigned (str): Number of unassigned reads
        """
        self.sample.append(get_name(filename))
        self.assigned.append(int(stdout.split('\n')[1].split(' ')[-1]))
        self.unassigned.append(int(stdout.split('\n')[0].split(' ')[-1]))
    
    def get_dict(self) -> dict[str, list]:
        """Get data as a dictionary"""
        return {
            "sample": self.sample,
            "assigned": self.assigned,
            "unassigned": self.unassigned
        }
        
def count_reads(fastq: str) -> int:
    """Count the number of reads in a fastq file

    Args:
        fq (str): path to fastq file

    Returns:
        int: number of reads
    """
    with gzip.open(fastq, "rt") as fq:
        line_count: int = sum(1 for _ in fq)
    read_count: int = line_count // 4
    return read_count

def concatenate_results(output_dir: str) -> pd.DataFrame:
    """Concatenate abundance tsv files, add a column for the sample name

    Args:
        output_dir (str): Path to directory containing abundance.tsv files

    Returns:
        pd.DataFrame: Concatenated dataframe
    """
    dfs: list[pd.DataFrame] = []
    for file in glob.glob(os.path.join(output_dir, "*abundance.tsv")):
        df: pd.DataFrame = pd.read_csv(file, sep="\t")
        df["sample"] = get_name(file)
        dfs.append(df)
    
    return pd.concat(dfs)
        
def parse_abundances(df: pd.DataFrame, threshold: int, level: Literal['genus', 'species']) -> pd.DataFrame:
    """Parse abundance data to remove low abundance taxa and group them as 'other'

    Args:
        df (pd.DataFrame): dataframe containing multi-sample abundance data
        threshold (int): Minimum abundance threshold
        level (Literal['genus', 'species']): Taxonomic level to parse

    Returns:
        pd.DataFrame: Parsed dataframe
    """
    df = df[['sample', 'abundance', level]]
    df = df[df['abundance'] > threshold]
    df.loc[df['abundance'] < threshold, level] = f'Other {'genera' if level == 'genus' else level} < {threshold}%'
    return df

def plot(df: pd.DataFrame) -> go.Figure:
    """Plot relative abundances

    Args:
        df (pd.DataFrame): Dataframe containing abundance data

    Returns:
        go.Figure: Plotly figure
    """
    
    fig = go.Figure(px.bar(
        df,
        x="sample",
        y="abundance",
        color=list(df.columns)[1],
        color_discrete_sequence=px.colors.qualitative.Dark24,
        title=f"Relative Abundances of {(list(df.columns)[1]).capitalize()}",
        labels={"sample": "Sample Name", "abundance": "Relative Abundance (%)"},
    ))
    
    return fig

def main():
    parser = argparse.ArgumentParser(
        "GermGenie",
        description="EMU wrapper for analyzing and plotting relative abundance from 16S data",
        epilog="Developed by Daan Brackel & Sander Boden @ ATLS-Avans",
    )
    parser.add_argument(
        "fastq", help="Path to folder containing gzipped fastq files", type=str
    )
    parser.add_argument(
        "output",
        help="Path to directory to place results (created if not exists.)",
        type=str,
    )
    parser.add_argument(
        "db", help="Path to EMU database", 
        type=str
        )
    parser.add_argument(
        "--threads",
        "-t",
        help="Number of threads to use for EMU classification (defaults to 2)",
        type=int,
        default=2,
    )
    parser.add_argument(
        "--threshold",
        "-T",
        help="Percent abundance threshold. Abundances below threshold will be shown as 'other' (defaults to 1 percent)",
        type=int,
        default=1,
    )
    parser.add_argument(
        "--tsv",
        help="Write abundances to tsv file (abundances.tsv)",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--nreads", "--nr",
        action="store_true",
        default=False,
        help="Visualize number of reads per sample in barplot",
    )
    parser.add_argument(
        "--subsample",
        "-s",
        help="Subsample fastq files to a specific number of reads. defaults to None (use all data)",
        type=int,
        default=None,
    )

    # args = parser.parse_args()

    # analysis = EMU(args.fastq, args.output, args.db, args.threads, args.threshold, args.nreads, args.subsample)

    # if args.tsv:
    #     analysis.df.to_csv(
    #         os.path.join(args.output, "abundances.tsv"), sep="\t", index=False
    #     )


if __name__ == "__main__":
    main()
