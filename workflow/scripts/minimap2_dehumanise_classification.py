import sys

sys.stderr = open(snakemake.log[0], "w")

from pafpy import PafFile
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
    with open(snakemake.input.truth, newline="") as fd:
        reader = csv.DictReader(fd, delimiter="\t")
        for row in reader:
            truth[row["read_id"]] = {c: row[c] for c in COLUMNS}

    with open(snakemake.output.classification, "w") as fd_out:
        print(f"read_id{DELIM}classification", file=fd_out)
        seen = set()
        with PafFile(snakemake.input.alignment) as paf:
            for record in paf:
                read_id = record.qname
                if read_id in seen:
                    continue
                else:
                    seen.add(read_id)
                read_tax = truth.get(read_id)
                if read_tax is None:
                    raise KeyError(f"{read_id} not in truth")
                read_is_human = read_tax["species_id"] == HUMAN_SPECIES_ID
                if record.is_unmapped():
                    if read_is_human:
                        clf = FN
                    else:
                        clf = TN
                else:
                    if read_is_human:
                        clf = TP
                    else:
                        clf = FP

                print(f"{read_id}{DELIM}{clf}", file=fd_out)

    truth_read_ids = set(truth.keys())
    diff = truth_read_ids.symmetric_difference(seen)
    if len(diff) > 0:
        raise ValueError(f"Got read differences between truth and PAF\n{diff}")


main()
