import sys

sys.stderr = open(snakemake.log[0], "w")

from pathlib import Path

import pandas as pd

DELIM = ","


def main():
    predf = pd.read_table(snakemake.input.preclassification[0])
    pre_covg = predf["coverage"][0]
    pre_depth = predf["meandepth"][0]

    with open(snakemake.output[0], "w") as fp:
        fp.write(
            f"tool{DELIM}db{DELIM}tech{DELIM}pre_covg{DELIM}pre_depth{DELIM}post_covg{DELIM}post_depth{DELIM}rel_covg{DELIM}rel_depth\n"
        )

        for p in map(Path, snakemake.input.postclassifications):
            tool = p.parts[-2]
            db = p.name.split(".")[0]
            tech = snakemake.wildcards.tech
            postdf = pd.read_table(p)
            post_covg = postdf["coverage"][0]
            post_depth = postdf["meandepth"][0]

            rel_covg = post_covg / pre_covg
            rel_depth = post_depth / pre_depth

            fp.write(
                f"{tool}{DELIM}{db}{DELIM}{tech}{DELIM}{pre_covg}{DELIM}{pre_depth}{DELIM}{post_covg}{DELIM}{post_depth}{DELIM}{rel_covg}{DELIM}{rel_depth}\n"
            )


if __name__ == "__main__":
    main()
