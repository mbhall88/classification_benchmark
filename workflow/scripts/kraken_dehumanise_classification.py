import sys

sys.stderr = open(snakemake.log[0], "w")

import csv

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


def main():
    truth = dict()
    is_illumina = "illumina" in snakemake.output.classification
    with open(snakemake.input.truth, newline="") as fd:
        reader = csv.DictReader(fd, delimiter="\t")
        for row in reader:
            truth[row["read_id"]] = {c: row[c] for c in COLUMNS}

    with open(snakemake.output.classification, "w") as fd_out:
        print(f"read_id{DELIM}classification", file=fd_out)
        seen = set()
        with open(snakemake.input.classification) as fd_in:
            for line in map(str.rstrip, fd_in):
                fields = line.split("\t")
                read_id = fields[1]
                if read_id in seen:
                    raise ValueError(f"Seen {read_id} multiple times")
                else:
                    seen.add(read_id)

                if is_illumina:
                    read_tax = truth.get(f"{read_id}/1")
                    read_tax2 = truth.get(f"{read_id}/2")
                    assert read_tax == read_tax2, read_id
                else:
                    read_tax = truth.get(read_id)
                if read_tax is None:
                    raise KeyError(f"{read_id} not in truth")
                read_is_human = read_tax["species_id"] == HUMAN_SPECIES_ID
                taxid = fields[2]

                if taxid == HUMAN_SPECIES_ID:
                    if read_is_human:
                        clf = TP
                    else:
                        clf = FP
                else:
                    if read_is_human:
                        clf = FN
                    else:
                        clf = TN

                if is_illumina:
                    print(f"{read_id}/1{DELIM}{clf}", file=fd_out)
                    print(f"{read_id}/2{DELIM}{clf}", file=fd_out)
                else:
                    print(f"{read_id}{DELIM}{clf}", file=fd_out)

    truth_read_ids = set(truth.keys())
    if is_illumina:
        new_seen = set()
        for readid in seen:
            for i in [1, 2]:
                new_seen.add(f"{readid}/{i}")
        seen = new_seen
    diff = truth_read_ids.symmetric_difference(seen)
    if len(diff) > 0:
        raise ValueError(f"Got read differences between truth and PAF\n{diff}")


main()
