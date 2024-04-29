import subprocess
import os
import glob

import pandas as pd
import plotly.express as px


class EMU:
    def __init__(
        self,
        input_dir: str,
        output_dir: str,
        db: str,
        threads: int = 2,
        threshold: int = 1,
    ) -> None:
        """Run EMU on samples, write barplots

        Args:
            input_dir (str): Directory containing (gzipped) fastq files
            output_dir (str): Directory for writing results
            db (str): Path to EMU database
            threads (int, optional): Number of threads. Defaults to 2.
            threshold (int, optional): %-abundance threshold for visualizing taxa. Defaults to 1.
        """
        # assign std vars
        self.input: str = input_dir
        self.output: str = output_dir
        self.t: int = threads
        self.db: str = db
        self.threshold: int = threshold
        self.assignment: dict[str, list] = {
            "sample": [],
            "assigned": [],
            "unassigned": [],
        }
        # create path for emu results
        self.emu_dir: str = os.path.join(self.output, "emu")
        # make output dirs
        if not os.path.exists(self.output):
            os.mkdir(self.output)
            os.mkdir(self.emu_dir)

        # search for input files
        infiles = glob.glob(os.path.join(self.input, "*.gz"))
        if len(infiles) < 1:
            print("No input files found...")
            os.abort()

        # run emu on input files
        for fastq in sorted(infiles):
            self.run_emu(fastq)
        # merge emu results
        self.df: pd.DataFrame = self.merge_results()
        # plot species
        fig = self.plot(self.species_data())
        fig.write_html(os.path.join(self.output, "species_bar.html"))
        # plot genus
        fig = self.plot(self.genus_data())
        fig.write_html(os.path.join(self.output, "genus_bar.html"))

        # write table of assigned/unassigned reads
        pd.DataFrame(self.assignment).to_csv(
            os.path.join(self.output, "assignments.tsv"), sep="\t", index=False
        )

    def run_emu(self, fq: str) -> None:
        """Run emu on a folder of fastq.gz files

        Args:
            fq (str): path to fastq file
        """
        name = fq.split("/")[-1].split(".")[0]
        command = f"emu abundance {fq} --db {self.db} --output-basename {name} --output-dir {self.emu_dir} --threads {self.t}"
        run = subprocess.run(command, shell=True, text=True, stdout=subprocess.PIPE)
        stdout = run.stdout.split("\n")

        self.assignment["sample"].append(name)
        self.assignment["assigned"].append(stdout[1].split(" ")[-1])
        self.assignment["unassigned"].append(stdout[0].split(" ")[-1])

    def merge_results(self) -> pd.DataFrame:
        """Merge abundance tsv's from EMU

        Returns:
            pd.DataFrame: Table of genus and species abundances
        """

        def process_file(infile: str) -> pd.DataFrame:
            data = pd.read_csv(infile, sep="\t")
            # Select only the required columns
            processed = data[["abundance", "species", "genus"]]
            return processed

        # get list of tsv's to merge
        emu_results = glob.glob(os.path.join(self.emu_dir, "*abundance.tsv"))
        # init empty dataframe
        df = pd.DataFrame(columns=["sample", "species", "genus", "abundance"])
        for tsv in emu_results:
            # append data to dataframe
            name = tsv.split("/")[-1].split("_")[0]
            tmpdf = process_file(tsv)
            tmpdf["sample"] = name
            df = pd.concat([df, tmpdf])

        return df

    def species_data(self):
        # Species plot
        species_df = self.df.drop("genus", axis=1)  # Drop genus column
        species_df["abundance"] *= 100  # Multiply abundance by 100 to get percentage
        species_df = species_df.groupby(["sample", "species"]).sum().reset_index()
        # Aggregate species with <1% abundance into 'Other species <1%'
        species_df.loc[
            species_df["abundance"] < self.threshold, "species"
        ] = f"Other species <{self.threshold}%"
        species_df = species_df.groupby(["sample", "species"]).sum().reset_index()

        return species_df

    def genus_data(self):
        # Genus plot
        genus_df = self.df.drop("species", axis=1)  # Drop species column
        genus_df["abundance"] *= 100  # Multiply abundance by 100 to get percentage
        genus_df = genus_df.groupby(["sample", "genus"]).sum().reset_index()
        # Aggregate genera with <1% abundance into 'Other genera <1%'
        genus_df.loc[
            genus_df["abundance"] < self.threshold, "genus"
        ] = f"Other genera <{self.threshold}%"
        genus_df = genus_df.groupby(["sample", "genus"]).sum().reset_index()

        return genus_df

    def plot(self, df: pd.DataFrame):
        return px.bar(
            df,
            x="sample",
            y="abundance",
            color=list(df.columns)[1],
            title=f"Relative Abundances of {(list(df.columns)[1]).capitalize()}",
            labels={"sample": "Sample Name", "abundance": "Relative Abundance (%)"},
        )
