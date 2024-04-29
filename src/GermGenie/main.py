import argparse
import os

from GermGenie.pipeline import EMU


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
    parser.add_argument("db", help="Path to EMU database", type=str)
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
    )

    args = parser.parse_args()

    analysis = EMU(args.fastq, args.output, args.db, args.threads, args.threshold)

    if args.tsv:
        analysis.df.to_csv(
            os.path.join(args.output, "abundances.tsv"), sep="\t", index=False
        )


if __name__ == "__main__":
    main()
