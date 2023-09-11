import sys

sys.stderr = open(snakemake.log[0], "w")
from pathlib import Path
import pandas as pd
from math import ceil
from collections import Counter

DELIM = ","
SIGFIG = 4


def summary(counts: Counter) -> tuple[int, int, int, int, float, float, float]:
    fns = counts["FN"]
    tps = counts["TP"]
    fps = counts["FP"]
    tns = counts["TN"]

    precision = round(tps / (tps + fps), SIGFIG)
    sn = round(tps / (tps + fns), SIGFIG)
    f1 = round((2 * tps) / ((2 * tps) + fps + fns), SIGFIG)

    return fns, fps, tns, tps, sn, precision, f1


def main():
    data = {}
    for p in map(Path, snakemake.input.classifications):
        tool = ".".join(p.name.split(".")[1:3])
        df = pd.read_csv(p, sep="\t")

        genus_counts = Counter(df["genus_classification"])
        species_counts = Counter(df["species_classification"])
        mtb_counts = Counter(df["mtb_classification"])
        mtbc_counts = Counter(df["mtbc_classification"])
        counts = dict(
            genus=genus_counts, species=species_counts, mtb=mtb_counts, mtbc=mtbc_counts
        )
        data[tool] = counts

    # Genus summary
    with open(snakemake.output.genus_summary, "w") as fd_out:
        print(
            DELIM.join(
                [
                    "tool",
                    "Seconds",
                    "Max. Memory (MB)",
                    "FN",
                    "FP",
                    "TN",
                    "TP",
                    "Recall",
                    "Precision",
                    "F-score",
                ]
            ),
            file=fd_out,
        )
        for p in map(Path, snakemake.input.benchmarks):
            tool = ".".join(p.parts[-3:-1])

            df = pd.read_csv(p, sep="\t")
            secs = str(ceil(list(df["s"])[0]))
            rss = str(ceil(list(df["max_rss"])[0]))

            counts = data[tool].get("genus")
            res = list(map(str, summary(counts)))
            print(DELIM.join([tool, secs, rss, *res]), file=fd_out)

    # Species summary
    with open(snakemake.output.species_summary, "w") as fd_out:
        print(
            DELIM.join(
                [
                    "tool",
                    "Seconds",
                    "Max. Memory (MB)",
                    "FN",
                    "FP",
                    "TN",
                    "TP",
                    "Recall",
                    "Precision",
                    "F-score",
                ]
            ),
            file=fd_out,
        )
        for p in map(Path, snakemake.input.benchmarks):
            tool = ".".join(p.parts[-3:-1])

            df = pd.read_csv(p, sep="\t")
            secs = str(ceil(list(df["s"])[0]))
            rss = str(ceil(list(df["max_rss"])[0]))

            counts = data[tool].get("species")
            res = list(map(str, summary(counts)))
            print(DELIM.join([tool, secs, rss, *res]), file=fd_out)

    # MTB summary
    with open(snakemake.output.mtb_summary, "w") as fd_out:
        print(
            DELIM.join(
                [
                    "tool",
                    "Seconds",
                    "Max. Memory (MB)",
                    "FN",
                    "FP",
                    "TN",
                    "TP",
                    "Recall",
                    "Precision",
                    "F-score",
                ]
            ),
            file=fd_out,
        )
        for p in map(Path, snakemake.input.benchmarks):
            tool = ".".join(p.parts[-3:-1])

            df = pd.read_csv(p, sep="\t")
            secs = str(ceil(list(df["s"])[0]))
            rss = str(ceil(list(df["max_rss"])[0]))

            counts = data[tool].get("mtb")
            res = list(map(str, summary(counts)))
            print(DELIM.join([tool, secs, rss, *res]), file=fd_out)

    # MTBC summary
    with open(snakemake.output.mtbc_summary, "w") as fd_out:
        print(
            DELIM.join(
                [
                    "tool",
                    "Seconds",
                    "Max. Memory (MB)",
                    "FN",
                    "FP",
                    "TN",
                    "TP",
                    "Recall",
                    "Precision",
                    "F-score",
                ]
            ),
            file=fd_out,
        )
        for p in map(Path, snakemake.input.benchmarks):
            tool = ".".join(p.parts[-3:-1])

            df = pd.read_csv(p, sep="\t")
            secs = str(ceil(list(df["s"])[0]))
            rss = str(ceil(list(df["max_rss"])[0]))

            counts = data[tool].get("mtbc")
            res = list(map(str, summary(counts)))
            print(DELIM.join([tool, secs, rss, *res]), file=fd_out)


main()
