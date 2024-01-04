import sys
from dataclasses import dataclass
from functools import cache
from typing import Optional, List

sys.stderr = open(snakemake.log[0], "w")

from pathlib import Path
import csv
from taxonomy import Taxonomy, TaxonomyError
from Bio import Entrez


Entrez.email = "michael.hall2@unimelb.edu.au"
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

MTB = "mtb"
ZYMO = "zymo"
HUMAN = "human"
MTB_TAXID = "1773"
UNCLASSIFIED = "0"
MTBC_TAXID = "77643"
MYCO_TAXIDS = {
    "1763",  # Mycobacterium
    "670516",  # Mycobacteroides
    "2126281",  # Mycolicibacillus
    "1073531",  # Mycolicibacter
    "1866885",  # Mycolicibacterium
}


@dataclass
class TaxonomyNode:
    id: str
    name: str
    parent: Optional["TaxonomyNode"]
    rank: str


@cache
def fetch_taxonomy(taxid: str) -> List[TaxonomyNode]:
    handle = Entrez.efetch(db="taxonomy", id=taxid)
    record = Entrez.read(handle)[0]
    lineage = []
    prev_node = None
    for lineage_info in record["LineageEx"]:
        _id = lineage_info["TaxId"]
        name = lineage_info["ScientificName"]
        rank = lineage_info["Rank"]
        parent = prev_node
        node = TaxonomyNode(id=_id, name=name, rank=rank, parent=parent)
        lineage.append(node)
        prev_node = node

    return lineage


def is_mtbc(taxid: str, taxtree) -> bool:
    if taxid == UNCLASSIFIED:
        return False

    if taxid == MTB:
        return True

    if taxid in (ZYMO, HUMAN):
        return False

    try:
        lineage = taxtree.lineage(taxid)
    except TaxonomyError:
        lineage = fetch_taxonomy(taxid)

    for node in lineage:
        if node.id == MTBC_TAXID:
            return True

    return False


def is_mtb(taxid: str, taxtree) -> bool:
    if taxid == UNCLASSIFIED:
        return False

    if taxid == MTB:
        return True

    if taxid in (ZYMO, HUMAN):
        return False

    try:
        lineage = taxtree.lineage(taxid)
    except TaxonomyError:
        lineage = fetch_taxonomy(taxid)

    for node in lineage:
        if node.id == MTB_TAXID:
            return True

    return False


def is_mycobacterium(taxid: str, taxtree: Taxonomy) -> bool:
    if taxid == UNCLASSIFIED:
        return False

    if taxid == MTB:
        return True

    if taxid in (ZYMO, HUMAN):
        return False

    try:
        lineage = taxtree.lineage(taxid)
    except TaxonomyError:
        lineage = fetch_taxonomy(taxid)

    for node in lineage:
        if node.id in MYCO_TAXIDS:
            return True

    return False


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
    is_illumina = "illumina" in snakemake.output.classification
    with open(snakemake.input.truth) as fd:
        for row in map(str.rstrip, fd):
            read_id, species_id = row.split(DELIM)
            truth[read_id] = species_id

    with open(snakemake.output.classification, "w") as fd_out:
        print(
            DELIM.join(["read_id", "mtb_classification"]),
            file=fd_out,
        )
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

                truth_taxid = read_tax
                called_taxid = fields[2]

                if not is_mycobacterium(truth_taxid, taxtree):
                    if is_mycobacterium(called_taxid, taxtree):
                        if is_mtb(called_taxid, taxtree):
                            mtb_clf = FP
                        else:
                            mtb_clf = TN
                    else:
                        mtb_clf = TN
                else:  # truth is mycobacterium
                    if not is_mycobacterium(called_taxid, taxtree):
                        if is_mtb(truth_taxid, taxtree):
                            mtb_clf = FN
                        else:
                            mtb_clf = NA
                    else:  # called is mycobacterium
                        if is_mtb(truth_taxid, taxtree):
                            if is_mtb(called_taxid, taxtree):
                                mtb_clf = TP
                            else:
                                mtb_clf = FN

                row = [mtb_clf]
                if is_illumina:
                    print(
                        DELIM.join([f"{read_id}/1", *row]),
                        file=fd_out,
                    )
                    print(
                        DELIM.join([f"{read_id}/2", *row]),
                        file=fd_out,
                    )
                else:
                    print(
                        DELIM.join([read_id, *row]),
                        file=fd_out,
                    )


main()
