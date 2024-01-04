import sys
from dataclasses import dataclass
from functools import cache
from typing import List, Optional

sys.stderr = open(snakemake.log[0], "w")

import csv
from pathlib import Path

from Bio import Entrez
from pafpy import PafFile
from taxonomy import Taxonomy, TaxonomyError

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
HUMAN_SPECIES_ID = "9606"
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

    acc2taxid = dict()
    with open(snakemake.input.acc_truth) as fd:
        reader = csv.DictReader(fd, delimiter="\t")
        for row in reader:
            acc2taxid[row["accession"]] = {c: row[c] for c in COLUMNS}

    truth = dict()
    is_illumina = "illumina" in snakemake.output.classification
    with open(snakemake.input.truth) as fd:
        for row in map(str.rstrip, fd):
            read_id, species_id = row.split(DELIM)
            truth[read_id] = species_id

    with open(snakemake.output.classification, "w") as fd_out:
        print(
            DELIM.join(
                [
                    "read_id",
                    "mtb_classification",
                ]
            ),
            file=fd_out,
        )

        clfs = dict()

        with PafFile(snakemake.input.aln) as paf:
            for record in paf:
                read_id = record.qname
                if is_illumina:
                    read_id = read_id[:-2]

                if is_illumina:
                    read_tax = truth.get(f"{read_id}/1")
                    read_tax2 = truth.get(f"{read_id}/2")
                    assert read_tax == read_tax2, read_id
                else:
                    read_tax = truth.get(read_id)
                if read_tax is None:
                    raise KeyError(f"{read_id} not in truth")

                truth_taxid = read_tax

                if not truth_taxid:  # is probably a virus
                    print(
                        f"ERROR: {read_id} has no species ID",
                        file=sys.stderr,
                    )
                    sys.exit(1)

                called_taxid = acc2taxid.get(record.tname, {}).get("species_id", "0")

                if not truth_taxid:
                    mtb_clf = NA
                elif not is_mycobacterium(truth_taxid, taxtree):
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
                    else:
                        if is_mtb(truth_taxid, taxtree):
                            if is_mtb(called_taxid, taxtree):
                                mtb_clf = TP
                            else:
                                mtb_clf = FN

                current_clfs = clfs.get(read_id)
                new_clfs = [mtb_clf]
                if current_clfs is None:
                    clfs[read_id] = new_clfs
                else:
                    updated_clfs = []
                    for current_c, new_c in zip(current_clfs, new_clfs):
                        if current_c == new_c:
                            updated_clfs.append(current_c)
                        elif TP in (current_c, new_c):
                            updated_clfs.append(TP)
                        elif FP in (current_c, new_c):
                            updated_clfs.append(FP)
                        else:
                            raise NotImplementedError(
                                f"Got an unexpected classification comparison. The current is {current_c} and the new is {new_c}. QNAME {record.qname} TNAME {record.tname}"
                            )
                    assert len(updated_clfs) == len(new_clfs), record
                    clfs[read_id] = updated_clfs

        non_human_read_ids = set()
        for rid, read_tax in truth.items():
            taxid = read_tax
            if taxid != HUMAN:
                non_human_read_ids.add(rid)

        # make sure all read ids in truth have been classified, if not, print all missing read ids to stderr
        missing_read_ids = non_human_read_ids - set(clfs.keys())
        if missing_read_ids:
            print(
                f"WARNING: {len(missing_read_ids)} read ids were not classified. They are:",
                file=sys.stderr,
            )
            for read_id in missing_read_ids:
                print(read_id, file=sys.stderr)

        for read_id, entires in clfs.items():
            if is_illumina:
                print(
                    DELIM.join([f"{read_id}/1", *entires]),
                    file=fd_out,
                )
                print(
                    DELIM.join([f"{read_id}/2", *entires]),
                    file=fd_out,
                )
            else:
                print(
                    DELIM.join([read_id, *entires]),
                    file=fd_out,
                )


main()
