import sys

sys.stderr = open(snakemake.log[0], "w")

import pysam

DELIM = "\t"
FN = "FN"
TN = "TN"
FP = "FP"
TP = "TP"
NA = "NA"
HUMAN_SPECIES_ID = "human"


def classify_kept(read_id, truth) -> str:
    species_id = truth.get(read_id)
    if species_id is None:
        raise KeyError(f"{read_id} is not in truth")

    if species_id == HUMAN_SPECIES_ID:
        return FN
    else:
        return TN


def classify_scrubbed(read_id, truth) -> str:
    species_id = truth.get(read_id)
    if species_id is None:
        raise KeyError(f"{read_id} is not in truth")

    if species_id == HUMAN_SPECIES_ID:
        return TP
    else:
        return FP


def main():
    truth = dict()
    with open(snakemake.input.truth) as fd:
        for row in map(str.rstrip, fd):
            read_id, species_id = row.split(DELIM)
            truth[read_id] = species_id

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

        if "illumina" in snakemake.output.classification:
            with pysam.FastxFile(snakemake.input.reads2) as fh:
                for entry in fh:
                    if entry.name in kept:
                        raise ValueError(
                            f"Seen {entry.name} multiple times in kept reads"
                        )
                    kept.add(entry.name)
                    clf = classify_kept(entry.name, truth)
                    print(f"{entry.name}{DELIM}{clf}", file=fd_out)

        truth_read_ids = set(truth.keys())
        not_seen_read_ids = truth_read_ids - kept
        for read_id in not_seen_read_ids:
            clf = classify_scrubbed(read_id, truth)
            print(f"{read_id}{DELIM}{clf}", file=fd_out)


main()
