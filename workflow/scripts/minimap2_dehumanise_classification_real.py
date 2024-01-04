import sys

sys.stderr = open(snakemake.log[0], "w")

from pafpy import PafFile

DELIM = "\t"
FN = "FN"
TN = "TN"
FP = "FP"
TP = "TP"
NA = "NA"
HUMAN_SPECIES_ID = "human"


def main():
    truth = dict()
    with open(snakemake.input.truth) as fd:
        for row in map(str.rstrip, fd):
            read_id, species_id = row.split(DELIM)
            truth[read_id] = species_id

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
                read_is_human = read_tax == HUMAN_SPECIES_ID

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
