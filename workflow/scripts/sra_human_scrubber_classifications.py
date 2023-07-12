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


def classify_kept(acc, truth) -> str:
    species_id = truth.get(acc, dict()).get("species_id")
    if species_id is None:
        raise KeyError(f"{acc} is not in truth")

    if species_id == HUMAN_SPECIES_ID:
        return FN
    else:
        return TN


def classify_scrubbed(acc, truth) -> str:
    species_id = truth.get(acc, dict()).get("species_id")
    if species_id is None:
        raise KeyError(f"{acc} is not in truth")

    if species_id == HUMAN_SPECIES_ID:
        return TP
    else:
        return FP


def main():
    with pysam.FastxFile(snakemake.input.all_reads) as fh:
        read_ids = {entry.name for entry in fh}

    truth = dict()
    with open(snakemake.input.truth, newline="") as fd:
        reader = csv.DictReader(fd, delimiter="\t")
        for row in reader:
            truth[row["accession"]] = {
                "species_id": row["species_id"],
                "genus_id": row["genus_id"],
            }

    with open(snakemake.output.classification, "w") as fd_out:
        print(f"read_id{DELIM}classification", file=fd_out)
        kept = set()
        with pysam.FastxFile(snakemake.input.reads) as fh:
            for entry in fh:
                if entry.name in kept:
                    raise ValueError(f"Seen {entry.name} multiple times in kept reads")
                kept.add(entry.name)
                read_acc = entry.comment.split()[0].split(",")[0]
                clf = classify_kept(read_acc, truth)
                print(f"{entry.name}{DELIM}{clf}", file=fd_out)

            scrubbed = set()
            with pysam.FastxFile(snakemake.input.removed) as fh:
                for entry in fh:
                    if entry.name in scrubbed:
                        raise ValueError(
                            f"Seen {entry.name} multiple times in scrubbed reads"
                        )
                    if entry.name in kept:
                        raise ValueError(
                            f"Seen {entry.name} in scrubbed and kept reads"
                        )
                    scrubbed.add(entry.name)
                    read_acc = entry.comment.split()[0].split(",")[0]
                    clf = classify_scrubbed(read_acc, truth)
                    print(f"{entry.name}{DELIM}{clf}", file=fd_out)

    seen_read_ids = kept.intersection(scrubbed)

    if seen_read_ids != read_ids:
        diff = seen_read_ids.symmetric_difference(read_ids)
        raise ValueError(
            f"These read IDs were not found in both the input and output\n{diff}"
        )


main()
