# GermGenie

GermGenie was specifically designed to analyse 16S data from clinical FFPE specimens, however it can be used to analyse any bacterial sample. GermGenie outputs stacked barplot showing the abundance of every species in your sample. By setting an abundance threshold, any species below the threshold will be added to an 'other' category (>1% by default).  
This tool was designed with Oxford Nanopore sequencing reads (ONT), and was not tested with any other sequencing data. The input should be a folder containing one or more samples in a fastq.gz format.

## Dependencies
The pipeline is based on [EMU](https://github.com/treangenlab/emu). Optional QC is performed with [chopper](https://github.com/wdecoster/chopper). Data is visualized using the [Plotly library](https://plotly.com).


# Installation

Follow EMU's installation instructions from the [repo](https://github.com/treangenlab/emu).
If you want to filter based on qualityscores or length, install chopper with the instructions from the [repo](https://github.com/wdecoster/chopper?tab=readme-ov-file#installation)
After installing EMU, install GermGenie in the same conda environment.
```bash
python -m pip install GermGenie
```

# Usage
```bash
usage: GermGenie [-h] [--version] [--threads THREADS]
                 [--threshold THRESHOLD] [--tsv] [--nreads]
                 [--subsample SUBSAMPLE] [--top_n TOP_N]
                 [--min-length MIN_LENGTH] [--max-length MAX_LENGTH]
                 [--min-quality MIN_QUALITY]
                 fastq output db

EMU wrapper for analyzing and plotting relative abundance from 16S
data

positional arguments:
  fastq                 Path to folder containing gzipped fastq
                        files
  output                Path to directory to place results (created
                        if not exists.)
  db                    Path to EMU database

options:
  -h, --help            show this help message and exit
  --version             Show program's version number and exit
  --threads THREADS, -t THREADS
                        Number of threads to use for EMU
                        classification (defaults to 2)
  --threshold THRESHOLD, -T THRESHOLD
                        Percent abundance threshold. Abundances
                        below threshold will be shown as 'other'
                        (defaults to 1 percent)
  --tsv                 Write abundances to tsv file
                        (abundances.tsv)
  --nreads, -nr         Visualize number of reads per sample in
                        barplot
  --subsample SUBSAMPLE, -s SUBSAMPLE
                        WARNING: DO NOT USE !!!
  --top_n TOP_N, -tn TOP_N
                        Number of top taxa to plot. 0 for all taxa.
  --min-length MIN_LENGTH, -mil MIN_LENGTH
                        Minimum length of reads to keep. Default is
                        to keep all reads.
  --max-length MAX_LENGTH, -mal MAX_LENGTH
                        Maximum length of reads to keep. Default is
                        to keep all reads.
  --min-quality MIN_QUALITY, -miq MIN_QUALITY
                        Minimum average Phred quality score of reads
                        to keep. Default is to keep all reads.

Developed by Daan Brackel, Birgit Rijvers & Sander Boden @ ATLS-
Avans
```
