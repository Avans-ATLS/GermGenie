# GermGenie

GermGenie was specifically designed to analyse 16S data from clinical FFPE specimens, however it can be used to analyse any bacterial sample. GermGenie outputs stacked barplot showing the abundance of every species/genus in your sample. By setting an abundance threshold, any species/genus below the threshold will be added to an 'other' category.  
This tool was designed with Oxford Nanopore sequencing reads (ONT), and was not tested with any other sequencing data. The input should be a folder containing one or more samples in a fastq.gz format.

## Dependencies
The pipeline is based on [EMU](https://github.com/treangenlab/emu). Data is visualized using the [Plotly library](https://plotly.com).


# Installation

Follow EMU's installation instructions from the [repo](https://github.com/treangenlab/emu).  
After installing EMU, install GermGenie in the same conda environment.
```bash
python -m pip install GermGenie
```



# Usage
```bash
usage: GermGenie [-h] [--threads THREADS] [--threshold THRESHOLD] [--tsv]
                 fastq output db

EMU wrapper for analyzing and plotting relative abundance from 16S data

positional arguments:
  fastq                 Path to folder containing gzipped fastq files
  output                Path to directory to place results (created if not
                        exists.)
  db                    Path to EMU database

optional arguments:
  -h, --help            show this help message and exit
  --threads THREADS, -t THREADS
                        Number of threads to use for EMU classification
                        (defaults to 2)
  --threshold THRESHOLD, -T THRESHOLD
                        Percent abundance threshold. Abundances below
                        threshold will be shown as 'other' (defaults to 1
                        percent)
  --tsv                 Write abundances to tsv file (abundances.tsv)

Developed by Daan Brackel & Sander Boden @ ATLS-Avans
```
