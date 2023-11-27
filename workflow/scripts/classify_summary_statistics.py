import sys

sys.stderr = open(snakemake.log[0], "w")
from collections import Counter
from math import ceil, sqrt
from pathlib import Path

import pandas as pd
from scipy import stats

DELIM = ","
SIGFIG = 4
RSS_SIGFIG = 1


def confidence_interval(n_s: int, n_f: int, conf: float = 0.95) -> tuple[float, float]:
    """Calculate the Wilson score interval.
    Equation take from https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Wilson_score_interval
    n_s: Number of successes or, in the case of confusion matrix statistics, the numerator
    n_f: Number of failures or, in the case of confusion matrix statistics, the denominator minus the numerator
    conf: the confidence level. i.e. 0.95 is 95% confidence
    """
    n = n_f + n_s
    z = stats.norm.ppf(1 - (1 - conf) / 2)  # two-sided
    z2 = z**2
    nz2 = n + z2
    A = (n_s + (0.5 * z2)) / nz2
    B = z / nz2
    C = sqrt(((n_s * n_f) / n) + (z2 / 4))
    CI = B * C
    return A - CI, A + CI


def summary(counts: Counter) -> tuple[int, int, int, int, float, float, float]:
    fns = counts["FN"]
    tps = counts["TP"]
    fps = counts["FP"]
    tns = counts["TN"]

    precision = round(tps / (tps + fps), SIGFIG)
    precision_lwr_bound, precision_upr_bound = confidence_interval(n_s=tps, n_f=fps)
    precision_lwr_bound = round(precision_lwr_bound, SIGFIG)
    precision_upr_bound = round(precision_upr_bound, SIGFIG)
    sn = round(tps / (tps + fns), SIGFIG)
    sn_lwr_bound, sn_upr_bound = confidence_interval(n_s=tps, n_f=fns)
    sn_lwr_bound = round(sn_lwr_bound, SIGFIG)
    sn_upr_bound = round(sn_upr_bound, SIGFIG)
    f1 = round((2 * tps) / ((2 * tps) + fps + fns), SIGFIG)
    f1_lwr_bound, f1_upr_bound = confidence_interval(n_s=tps, n_f=fps + fns)
    f1_lwr_bound = round(f1_lwr_bound, SIGFIG)
    f1_upr_bound = round(f1_upr_bound, SIGFIG)

    return (
        fns,
        fps,
        tns,
        tps,
        (sn, sn_lwr_bound, sn_upr_bound),
        (precision, precision_lwr_bound, precision_upr_bound),
        (f1, f1_lwr_bound, f1_upr_bound),
    )


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
                    "Max. Memory (GB)",
                    "FN",
                    "FP",
                    "TN",
                    "TP",
                    "Recall",
                    "Recall Lower Bound",
                    "Recall Upper Bound",
                    "Precision",
                    "Precision Lower Bound",
                    "Precision Upper Bound",
                    "F-score",
                    "F-score Lower Bound",
                    "F-score Upper Bound",
                ]
            ),
            file=fd_out,
        )
        for p in map(Path, snakemake.input.benchmarks):
            tool = ".".join(p.parts[-3:-1])

            df = pd.read_csv(p, sep="\t")
            # average the seconds
            secs = str(ceil(df["s"].mean()))

            # convert rss (MB) to GB and use the maximum of the repeats
            rss = str(round(df["max_rss"].max() / 1024, RSS_SIGFIG))

            counts = data[tool].get("genus")
            res = []
            for i in summary(counts):
                if isinstance(i, tuple):
                    res.extend(map(str, i))
                else:
                    res.append(str(i))
            print(DELIM.join([tool, secs, rss, *res]), file=fd_out)

    # Species summary
    with open(snakemake.output.species_summary, "w") as fd_out:
        print(
            DELIM.join(
                [
                    "tool",
                    "Seconds",
                    "Max. Memory (GB)",
                    "FN",
                    "FP",
                    "TN",
                    "TP",
                    "Recall",
                    "Recall Lower Bound",
                    "Recall Upper Bound",
                    "Precision",
                    "Precision Lower Bound",
                    "Precision Upper Bound",
                    "F-score",
                    "F-score Lower Bound",
                    "F-score Upper Bound",
                ]
            ),
            file=fd_out,
        )
        for p in map(Path, snakemake.input.benchmarks):
            tool = ".".join(p.parts[-3:-1])

            df = pd.read_csv(p, sep="\t")
            # average the seconds
            secs = str(ceil(df["s"].mean()))
            # convert rss (MB) to GB and use the maximum of the repeats
            rss = str(round(df["max_rss"].max() / 1024, RSS_SIGFIG))

            counts = data[tool].get("species")
            res = []
            for i in summary(counts):
                if isinstance(i, tuple):
                    res.extend(map(str, i))
                else:
                    res.append(str(i))
            print(DELIM.join([tool, secs, rss, *res]), file=fd_out)

    # MTB summary
    with open(snakemake.output.mtb_summary, "w") as fd_out:
        print(
            DELIM.join(
                [
                    "tool",
                    "Seconds",
                    "Max. Memory (GB)",
                    "FN",
                    "FP",
                    "TN",
                    "TP",
                    "Recall",
                    "Recall Lower Bound",
                    "Recall Upper Bound",
                    "Precision",
                    "Precision Lower Bound",
                    "Precision Upper Bound",
                    "F-score",
                    "F-score Lower Bound",
                    "F-score Upper Bound",
                ]
            ),
            file=fd_out,
        )
        for p in map(Path, snakemake.input.benchmarks):
            tool = ".".join(p.parts[-3:-1])

            df = pd.read_csv(p, sep="\t")
            # average the seconds
            secs = str(ceil(df["s"].mean()))
            # convert rss (MB) to GB and use the maximum of the repeats
            rss = str(round(df["max_rss"].max() / 1024, RSS_SIGFIG))

            counts = data[tool].get("mtb")
            res = []
            for i in summary(counts):
                if isinstance(i, tuple):
                    res.extend(map(str, i))
                else:
                    res.append(str(i))
            print(DELIM.join([tool, secs, rss, *res]), file=fd_out)

    # MTBC summary
    with open(snakemake.output.mtbc_summary, "w") as fd_out:
        print(
            DELIM.join(
                [
                    "tool",
                    "Seconds",
                    "Max. Memory (GB)",
                    "FN",
                    "FP",
                    "TN",
                    "TP",
                    "Recall",
                    "Recall Lower Bound",
                    "Recall Upper Bound",
                    "Precision",
                    "Precision Lower Bound",
                    "Precision Upper Bound",
                    "F-score",
                    "F-score Lower Bound",
                    "F-score Upper Bound",
                ]
            ),
            file=fd_out,
        )
        for p in map(Path, snakemake.input.benchmarks):
            tool = ".".join(p.parts[-3:-1])

            df = pd.read_csv(p, sep="\t")
            # average the seconds
            secs = str(ceil(df["s"].mean()))
            # convert rss (MB) to GB and use the maximum of the repeats
            rss = str(round(df["max_rss"].max() / 1024, RSS_SIGFIG))

            counts = data[tool].get("mtbc")
            res = []
            for i in summary(counts):
                if isinstance(i, tuple):
                    res.extend(map(str, i))
                else:
                    res.append(str(i))
            print(DELIM.join([tool, secs, rss, *res]), file=fd_out)


main()
