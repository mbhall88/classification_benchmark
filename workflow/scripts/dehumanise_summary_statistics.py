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
TOOL_NAMES = {
    "sra": "HRRT",
    "hostile": "Hostile",
}


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


def f_score_interval(
    tps: int, fps: int, fns: int, conf: float = 0.95
) -> tuple[float, float]:
    """Calculate the Wilson score interval for the F-score.
    tps: Number of true positives
    fps: Number of false positives
    fns: Number of false negatives
    conf: the confidence level. i.e. 0.95 is 95% confidence
    """
    precision_lwr_bound, precision_upr_bound = confidence_interval(n_s=tps, n_f=fps)
    recall_lwr_bound, recall_upr_bound = confidence_interval(n_s=tps, n_f=fns)

    f1_lwr_bound = (
        2
        * (precision_lwr_bound * recall_lwr_bound)
        / (precision_lwr_bound + recall_lwr_bound)
    )
    f1_upr_bound = (
        2
        * (precision_upr_bound * recall_upr_bound)
        / (precision_upr_bound + recall_upr_bound)
    )

    return f1_lwr_bound, f1_upr_bound


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
    f1_lwr_bound, f1_upr_bound = f_score_interval(tps=tps, fps=fps, fns=fns)
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
    n_reads = None
    for p in map(Path, snakemake.input.classifications):
        tool = p.name.split(".", maxsplit=1)[1].rsplit(".", maxsplit=2)[0]
        if tool in TOOL_NAMES:
            tool = TOOL_NAMES[tool]
        else:
            tool = tool.replace(".", " ")
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
                    "Method",
                    "Rate (reads/sec)",
                    "Peak Memory (GB)",
                    "FN",
                    "FP",
                    "TN",
                    "TP",
                    "Recall (95% CI)",
                    "Precision (95% CI)",
                    "F-score (95% CI)",
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
            if tool in TOOL_NAMES:
                tool = TOOL_NAMES[tool]
            else:
                tool = tool.replace(".", " ")

            df = pd.read_csv(p, sep="\t")

            # average the seconds
            secs = df["s"].mean()
            rate = str(ceil(n_reads / secs))

            # convert rss (MB) to GB and use the maximum of the repeats
            rss = str(round(df["max_rss"].max() / 1024, RSS_SIGFIG))

            counts = data[tool]

            res = []
            for i in summary(counts):
                if isinstance(i, tuple):
                    res.append(f"{i[0]} ({i[1]}-{i[2]})")
                else:
                    res.append(str(i))
            print(DELIM.join([tool, rate, rss, *res]), file=fd_out)


main()
