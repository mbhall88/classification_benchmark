import sys
from dataclasses import dataclass
from typing import Optional

sys.stderr = open(snakemake.log[0], "w")

import csv
from taxonomy import Taxonomy, TaxonomyError

DELIM = "\t"
FN = "FN"
TN = "TN"
FP = "FP"
TP = "TP"
NA = "NA"

COLUMNS = [
    "strain_id",
    "strain",
    "species_id",
    "species",
    "genus_id",
    "genus",
]

MTB_TAXID = "1773"
MTBC_TAXID = "77643"
MYCO_TAXIDS = {
    "1763",  # Mycobacterium
    "670516",  # Mycobacteroides
    "2126281",  # Mycolicibacillus
    "1073531",  # Mycolicibacter
    "1866885",  # Mycolicibacterium
}


def is_mtbc(taxid: str, taxtree) -> bool:
    lineage = taxtree.lineage(taxid)
    for node in lineage:
        if node.id == MTBC_TAXID:
            return True

    return False


def is_mtb(taxid: str, taxtree) -> bool:
    lineage = taxtree.lineage(taxid)
    for node in lineage:
        if node.id == MTB_TAXID:
            return True

    return False


def is_mycobacterium(taxid: str, taxtree: Taxonomy) -> bool:
    lineage = taxtree.lineage(taxid)
    for node in lineage:
        if node.id in MYCO_TAXIDS:
            return True

    return False


def is_genus_correct(true_taxid: str, called_taxid: str, taxtree: Taxonomy) -> bool:
    true_genus = taxtree.parent(true_taxid, at_rank="genus")
    called_genus = taxtree.parent(called_taxid, at_rank="genus")
    if true_genus is None:
        raise ValueError(f"True taxid {true_taxid} is above genus level")

    if called_genus is None:
        return False

    return called_genus == true_genus


def is_species_correct(true_taxid: str, called_taxid: str, taxtree: Taxonomy) -> bool:
    if is_mtbc(true_taxid, taxtree) and is_mtbc(called_taxid, taxtree):
        return True

    true_species = taxtree.parent(true_taxid, at_rank="species")
    called_species = taxtree.parent(called_taxid, at_rank="species")

    if true_species is None:
        raise ValueError(f"True taxid {true_taxid} is above species level")

    if called_species is None:
        return False

    if is_mtbc(true_taxid, taxtree) and is_mtbc(called_taxid, taxtree):
        return True

    return called_species == true_species


@dataclass
class TaxonomyNode:
    id: str
    name: str
    parent: Optional["TaxonomyNode"]
    rank: str


def main():
    taxdir = str(Path(snakemake.input.nodes).parent)
    taxtree = Taxonomy.from_ncbi(taxdir)

    truth = dict()
    is_illumina = snakemake.wildcards.tech == "illumina"
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

                truth_taxid = read_tax["species_id"]
                called_taxid = fields[2]

                # todo check genus and species level classifications

                if not truth_taxid and snakemake.params.ignore_unmapped:
                    clf = NA
                elif not is_mycobacterium(truth_taxid, taxtree):
                    if
                elif taxid == HUMAN_SPECIES_ID:
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
