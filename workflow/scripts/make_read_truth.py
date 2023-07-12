import csv
import sys

sys.stderr = open(snakemake.log[0], "w")

import pysam

DELIM = "\t"
UNMAPPED = {"junk_seq", "random_seq"}
COLUMNS = [
    "strain_id",
    "strain",
    "species_id",
    "species",
    "genus_id",
    "genus",
]


def main():
    truth = dict()
    with open(snakemake.input.acc2tax, newline="") as fd:
        reader = csv.DictReader(fd, delimiter="\t")
        for row in reader:
            truth[row["accession"]] = {c: row[c] for c in COLUMNS}

    with open(snakemake.output.metadata, "w") as fd_out:
        print(
            DELIM.join(["read_id", *COLUMNS]),
            file=fd_out,
        )

        with pysam.FastxFile(snakemake.input.reads) as fh:
            for read in fh:
                read_acc = read.comment.split()[0].split(",")[0]
                if read_acc in UNMAPPED:
                    row = [read.name, "", "", "", "", "", ""]
                else:
                    tax = truth.get(read_acc)
                    if tax is None:
                        raise KeyError(
                            f"Couldn't find accession {read_acc} for read {read.name} in accession truth"
                        )
                    row = [read.name, *[tax[c] for c in COLUMNS]]
                print(DELIM.join(row), file=fd_out)


main()
