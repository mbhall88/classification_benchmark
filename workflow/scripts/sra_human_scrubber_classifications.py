import csv
import sys

sys.stderr = open(snakemake.log[0], "w")

import pysam

DELIM = "\t"
FN = "FN"
TN = "TN"
FP = "FP"
TP = "TP"
HUMAN_SPECIES_ID = "9606"

COLUMNS = [
    "strain_id",
    "strain",
    "species_id",
    "species",
    "genus_id",
    "genus",
]


def classify_kept(read_id, truth) -> str:
    species_id = truth.get(read_id, dict()).get("species_id")
    if species_id is None:
        raise KeyError(f"{read_id} is not in truth")

    if species_id == HUMAN_SPECIES_ID:
        return FN
    else:
        return TN


def classify_scrubbed(read_id, truth) -> str:
    species_id = truth.get(read_id, dict()).get("species_id")
    if species_id is None:
        raise KeyError(f"{read_id} is not in truth")

    if species_id == HUMAN_SPECIES_ID:
        return TP
    else:
        return FP


def main():
    truth = dict()
    with open(snakemake.input.truth, newline="") as fd:
        reader = csv.DictReader(fd, delimiter="\t")
        for row in reader:
            truth[row["read_id"]] = {c: row[c] for c in COLUMNS}

    with open(snakemake.output.classification, "w") as fd_out:
        print(f"read_id{DELIM}classification", file=fd_out)
        kept = set()
        with pysam.FastxFile(snakemake.input.reads) as fh:
            for entry in fh:
                if entry.name in kept:
                    raise ValueError(f"Seen {entry.name} multiple times in kept reads")
                kept.add(entry.name)
                clf = classify_kept(entry.name, truth)
                print(f"{entry.name}{DELIM}{clf}", file=fd_out)

        scrubbed = set()
        with pysam.FastxFile(snakemake.input.removed) as fh:
            for entry in fh:
                if entry.name in scrubbed:
                    raise ValueError(
                        f"Seen {entry.name} multiple times in scrubbed reads"
                    )
                if entry.name in kept:
                    raise ValueError(f"Seen {entry.name} in scrubbed and kept reads")
                scrubbed.add(entry.name)
                clf = classify_scrubbed(entry.name, truth)
                print(f"{entry.name}{DELIM}{clf}", file=fd_out)

    read_ids = set(truth.keys())
    seen_read_ids = kept.union(scrubbed)

    if seen_read_ids != read_ids:
        diff = seen_read_ids.symmetric_difference(read_ids)
        raise ValueError(
            f"{len(read_ids)} reads in truth and {len(seen_read_ids)} in output. These read IDs were not found in both the input and output\n{diff}"
        )


main()
