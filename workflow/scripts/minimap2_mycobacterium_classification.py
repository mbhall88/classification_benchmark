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

    try:
        lineage = taxtree.lineage(taxid)
    except TaxonomyError:
        lineage = fetch_taxonomy(taxid)

    for node in lineage:
        if node.id in MYCO_TAXIDS:
            return True

    return False


def is_genus_correct(true_taxid: str, called_taxid: str, taxtree: Taxonomy) -> bool:
    if called_taxid == UNCLASSIFIED:
        return False

    true_genus = taxtree.parent(true_taxid, at_rank="genus")
    called_genus = taxtree.parent(called_taxid, at_rank="genus")
    if true_genus is None:
        raise ValueError(f"True taxid {true_taxid} is above genus level")

    if called_genus is None:
        return False

    return called_genus == true_genus


def is_species_correct(true_taxid: str, called_taxid: str, taxtree: Taxonomy) -> bool:
    if called_taxid == UNCLASSIFIED:
        return False

    # if is_mtbc(true_taxid, taxtree) and is_mtbc(called_taxid, taxtree):
    #     return True

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

    acc2taxid = dict()
    with open(snakemake.input.acc_truth) as fd:
        reader = csv.DictReader(fd, delimiter="\t")
        for row in reader:
            acc2taxid[row["accession"]] = {c: row[c] for c in COLUMNS}

    truth = dict()
    is_illumina = snakemake.wildcards.tech == "illumina"
    with open(snakemake.input.truth, newline="") as fd:
        reader = csv.DictReader(fd, delimiter="\t")
        for row in reader:
            truth[row["read_id"]] = {c: row[c] for c in COLUMNS}

    with open(snakemake.output.classification, "w") as fd_out:
        print(
            DELIM.join(
                [
                    "read_id",
                    "genus_classification",
                    "species_classification",
                    "mtb_classification",
                    "mtbc_classification",
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

                truth_taxid = read_tax["species_id"]

                if not truth_taxid:  # is probably a virus
                    genus_name = read_tax["genus"]
                    print(
                        f"WARN: {read_id} has no species ID, but has genus {genus_name}",
                        file=sys.stderr,
                    )
                    truth_taxid = read_tax["genus_id"]

                called_taxid = acc2taxid.get(record.tname, {}).get("species_id", "0")

                if not truth_taxid:
                    genus_clf = NA
                    species_clf = NA
                    mtb_clf = NA
                    mtbc_clf = NA
                elif not is_mycobacterium(truth_taxid, taxtree):
                    if is_mycobacterium(called_taxid, taxtree):
                        genus_clf = FP
                        species_clf = FP
                        if is_mtb(called_taxid, taxtree):
                            mtb_clf = FP
                        else:
                            mtb_clf = TN
                        if is_mtbc(called_taxid, taxtree):
                            mtbc_clf = FP
                        else:
                            mtbc_clf = TN
                    else:
                        genus_clf = TN
                        species_clf = TN
                        mtb_clf = TN
                        mtbc_clf = TN
                else:  # truth is mycobacterium
                    if not is_mycobacterium(called_taxid, taxtree):
                        genus_clf = FN
                        species_clf = FN
                        if is_mtb(truth_taxid, taxtree):
                            mtb_clf = FN
                        else:
                            mtb_clf = NA
                        if is_mtbc(truth_taxid, taxtree):
                            mtbc_clf = FN
                        else:
                            mtbc_clf = NA
                    else:
                        if is_genus_correct(truth_taxid, called_taxid, taxtree):
                            genus_clf = TP
                        else:
                            genus_clf = FN
                        if is_species_correct(truth_taxid, called_taxid, taxtree):
                            species_clf = TP
                        else:
                            species_clf = FN
                        if is_mtb(truth_taxid, taxtree):
                            if is_mtb(called_taxid, taxtree):
                                mtb_clf = TP
                            else:
                                mtb_clf = FN
                        if is_mtbc(truth_taxid, taxtree):
                            if is_mtbc(called_taxid, taxtree):
                                mtbc_clf = TP
                            else:
                                mtbc_clf = FN

                current_clfs = clfs.get(read_id)
                new_clfs = [genus_clf, species_clf, mtb_clf, mtbc_clf]
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
            taxid = read_tax["species_id"]
            if taxid != HUMAN_SPECIES_ID:
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
