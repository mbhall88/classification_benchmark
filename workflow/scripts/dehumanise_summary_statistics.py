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
    n_reads = None
    for p in map(Path, snakemake.input.classifications):
        tool = p.name.split(".", maxsplit=1)[1].rsplit(".", maxsplit=2)[0]
        df = pd.read_csv(p, sep="\t")
        if n_reads is None:
            n_reads = len(df)
        else:
            assert len(df) == n_reads, f"different number of reads in {tool}"

        counts = Counter(df["classification"])
        data[tool] = counts

    with open(snakemake.output.summary, "w") as fd_out:
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
            if "kraken" in str(p):
                suf = p.parts[-2]
                tool = f"kraken.{suf}"
            else:
                tool = p.parts[-2]

            df = pd.read_csv(p, sep="\t")
            secs = str(ceil(list(df["s"])[0]))
            rss = str(ceil(list(df["max_rss"])[0]))

            counts = data[tool]
            res = list(map(str, summary(counts)))
            print(DELIM.join([tool, secs, rss, *res]), file=fd_out)


main()
