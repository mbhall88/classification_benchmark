import sys

sys.stderr = open(snakemake.log[0], "w")
from pathlib import Path
import pandas as pd
from collections import Counter

DELIM = ","


def summary(counts: Counter) -> tuple[int, int, int, int, float, float, float]:
    fns = counts["FN"]
    tps = counts["TP"]
    fps = counts["FP"]
    tns = counts["TN"]

    precision = round(tps / (tps + fps), 5)
    sn = round(tps / (tps + fns), 5)
    f1 = round((2 * tps)/((2*tps) + fps + fns), 5)

    return fns, fps, tns, tps, sn, precision, f1


def main():
    data = {}
    is_illumina = "illumina" in snakemake.output.summary
    n_reads = None
    for p in map(Path, snakemake.input.classifications):
        tool = p.name.split(".", maxsplit=1)[1].rsplit(".", maxsplit=1)[0]
        df = pd.read_csv(p, sep="\t")
        if n_reads is None:
            n_reads = len(df)
        else:
            if is_illumina and tool.startswith("kraken"):
                assert len(df) == n_reads / 2, f"different number of reads in {tool}"
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
            secs = str(round(list(df["s"])[0], 1))
            rss = str(list(df["max_rss"])[0])

            counts = data[tool]
            res = list(map(str, summary(counts)))
            print(DELIM.join([tool, secs, rss, *res]), file=fd_out)


main()
