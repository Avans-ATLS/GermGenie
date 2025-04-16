import argparse
import os
import glob
import subprocess
from typing import List, Dict, Tuple
from plotly import express as px
from plotly import graph_objects as go

import pandas as pd


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


def find_input_files(input_dir: str) -> List[str]:
    """Find input fastq.gz files in a directory

    Args:
        input_dir (str): Path to directory containing fastq.gz files

    Returns:
        List[str]: List of fastq.gz files
    """
    # Find all fastq.gz files
    infiles: List[str] = glob.glob(os.path.join(input_dir, "*.gz"))
    if len(infiles) < 1:
        # If no files are found, abort
        print("No input files found...")  # TODO: replace with logging
        os.abort()
    else:
        # Return list of files
        return infiles

def subsample_fastq(nreads: int, fastq: str, outdir: str) -> str:
    """Subsample a fastq file

    Args:
        nreads (int): Number of reads to subsample
        fastq (str): Path to fastq file
        outdir (str): Path to output directory

    Returns:
        str: Path to subsampled fastq file, or input fastq path if an error occurred
    """
    name: str = os.path.basename(fastq).split(".")[0] 
    out: str = os.path.join(outdir, f"{name}_subsampled.fastq.gz")
    pass


def run_emu(fastq: str, db: str, threads: int, output_dir: str) -> str:
    """Run EMU on a fastq file

    Args:
        fastq (str): Path to fastq file
        db (str): Path to EMU compatible database
        threads (int): Number of threads to use
        output_dir (str): Path to output directory

    Returns:
        str: stdout from emu, or empty string if an error occurred
    """
    name: str = get_name(fastq)

    try:
        cmd = [
            "emu",
            "abundance",
            "--db",
            db,
            "--threads",
            str(threads),
            "--output-dir",
            output_dir,
            "--output-basename",
            name,
            fastq,
        ]
        print(" ".join(cmd))
        emu = subprocess.run(
            cmd,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        print(f"Error processing {name}")  # TODO: replace with logging
        print(e.stderr)
        print(e.stdout)
        return ""

    return emu.stdout


class ReadMapping:
    """Class for storing read statistics from emu stdout"""

    def __init__(self):
        self.sample: List[str] = []
        self.mapped: List[int] = []
        self.unmapped: List[int] = []
        self.unclassified_mapped: List[int] = []

    def add(self, sample, stdout: str) -> None:
        """Add a read stats from emu stdout

        Args:
            sample (str): Name of the sample
        """
        self.sample.append(sample)
        self.unmapped.append(int(stdout.split("\n")[0].split(" ")[-1]))
        self.mapped.append(int(stdout.split("\n")[1].split(" ")[-1]))
        self.unclassified_mapped.append(int(stdout.split("\n")[2].split(" ")[-1]))

    def to_dict(self) -> Dict[str, List]:
        """Get data as a dictionary"""
        return {
            "sample": self.sample,
            "mapped": self.mapped,
            "unmapped": self.unmapped,
            "unclassified_mapped": self.unclassified_mapped,
        }

def concatenate_results(output_dir: str) -> pd.DataFrame:
    """Concatenate abundance tsv files, add a column for the sample name

    Args:
        output_dir (str): Path to directory containing abundance.tsv files

    Returns:
        pd.DataFrame: Concatenated dataframe
    """
    dfs: List[pd.DataFrame] = []
    # Loop through all abundance files
    for file in glob.glob(os.path.join(output_dir, "*abundance.tsv")):
        # Read file and add sample name
        df: pd.DataFrame = pd.read_csv(file, sep="\t", skipfooter=2, engine="python")
        df["sample"] = "_".join(get_name(file).split("_")[:-1])
        # Append to list
        dfs.append(df)
    # Concatenate all dataframes and return
    return pd.concat(dfs)


def parse_abundances(df: pd.DataFrame, threshold: int, level: str) -> pd.DataFrame:
    """Parse abundance data to remove low abundance taxa and group them as 'other'

    Args:
        df (pd.DataFrame): dataframe containing multi-sample abundance data
        threshold (int): Minimum abundance threshold
        level (str): Taxonomic level to parse ('genus' or 'species)

    Returns:
        pd.DataFrame: Parsed dataframe
    """
    # Grab relevant columns
    df = df[["sample", "abundance", level]]
    # Convert abundance from fraction to percentage
    df = df.copy()
    df["abundance"] *= 100
    # Rename low abundance taxa to other
    other = f"Other {'genera' if level == 'genus' else level} < {threshold}%"
    df.loc[df["abundance"] < threshold, level] = other
    # Group by sample and taxon and sum percentages
    df = df.groupby(["sample", level]).sum().reset_index()
    return df


def plot(df: pd.DataFrame) -> go.Figure:
    """Plot relative abundances

    Args:
        df (pd.DataFrame): Dataframe containing abundance data

    Returns:
        go.Figure: Plotly figure
    """
    fig = go.Figure(
        px.bar(
            df,
            x="sample",
            y="abundance",
            color="species",
            color_discrete_sequence=px.colors.qualitative.Dark24,
            title=f"Relative Abundances of {(list(df.columns)[1]).capitalize()}",
            labels={"sample": "Sample Name", "abundance": "Relative Abundance (%)"},
        )
    )

    return fig


def plot_reads_mapped(readstats: ReadMapping) -> Tuple[go.Figure, pd.DataFrame]:
    """Plot number of reads mapped, unmapped, and unclassified"""
    # Create dataframe from readstats object
    stats = pd.DataFrame(readstats.to_dict())
    # Melt dataframe into longform
    stats = stats.melt(id_vars="sample", var_name="status", value_name="reads")
    fig = px.bar(
        stats,
        x="sample",
        y="reads",
        color="status",
        title="Read Mapping Statistics",
        labels={"sample": "Sample Name", "reads": "Number of Reads"},
    )
    return fig, stats

def cli() -> argparse.Namespace:
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
        "db", 
        help="Path to EMU database", 
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
        "--nreads",
        "-nr",
        action="store_true",
        default=False,
        help="Visualize number of reads per sample in barplot",
    )
    parser.add_argument( ## BROKEN :(
        "--subsample",
        "-s",
        help="WARNING: DO NOT USE !!!", # Subsample fastq files to a specific number of reads. defaults to None (use all data)",
        type=int,
        default=None,
    )
    
    parser.add_argument(
        '--top_n', 
        '-tn', 
        type=int, 
        default=0, 
        help="Number of top taxa to plot. 0 for all taxa.")

    return parser.parse_args()

def main() -> None:
    args = cli()
    # instantiate assignments
    readstats = ReadMapping()

    # Create output directories
    emu_dir: str = create_output_dirs(args.output, True if args.subsample else False)

    # Find input files
    infiles: List[str] = find_input_files(args.fastq)
    print(f'Found {len(infiles)} input files')

    # Subsample fastq files
    if args.subsample:
        infiles = [
            subsample_fastq(args.subsample, f, os.path.join(args.output, "subsample"))
            for f in infiles
        ]

    # Run EMU on fastq files
    for f in infiles:
        emu_stdout = run_emu(f, args.db, args.threads, emu_dir)
        if emu_stdout == "":
            print(f"Error processing {f}")  # TODO: replace with logging
            continue
        else:
            readstats.add(get_name(f), emu_stdout)

    # Plot read mapping statistics
    if args.nreads:
        fig, readsdf = plot_reads_mapped(readstats)
        fig.write_html(os.path.join(args.output, "read_mapping.html"))
        # Write to tsv
        if args.tsv:
            readsdf.pivot(
            index="sample",
            columns=["status"],
            values="reads",
        ).to_csv(os.path.join(args.output, "read_mapping.tsv"), sep="\t")

    # Concatenate results
    df = concatenate_results(emu_dir)

    # Parse data
    df = parse_abundances(df, args.threshold, "species")
    if args.top_n > 0:
        # Filter to top n taxa, group the rest as "other"
        top_taxa = df.groupby("species")["abundance"].sum().nlargest(args.top_n).index
        df.loc[~df["species"].isin(top_taxa), "species"] = f"Other species"
        df = df.groupby(["sample", "species"]).sum().reset_index()
        
    # Plot data
    fig = plot(df)
    fig.write_html(os.path.join(args.output, "relative_abundances.html"))

    # Write to tsv
    if args.tsv:
        df.to_csv(os.path.join(args.output, "abundances.tsv"), sep="\t", index=False)
        
    print("Done!")
        
if __name__ == "__main__":
    main()