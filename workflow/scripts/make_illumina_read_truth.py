import csv
import sys
from dataclasses import dataclass
from functools import cache
from typing import Optional, List

from Bio import Entrez

sys.stderr = open(snakemake.log[0], "w")

import pysam
from taxonomy import Taxonomy, TaxonomyError
from pathlib import Path

DELIM = "\t"
UNMAPPED = {"random"}
COLUMNS = [
    "strain_id",
    "strain",
    "species_id",
    "species",
    "genus_id",
    "genus",
]


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


def main():
    taxdir = str(Path(snakemake.input.nodes).parent)
    taxtree = Taxonomy.from_ncbi(taxdir)
    truth = dict()
    with open(snakemake.input.acc2tax, newline="") as fd:
        reader = csv.DictReader(fd, delimiter="\t")
        for row in reader:
            truth[row["accession"]] = {c: row[c] for c in COLUMNS}

    with open(snakemake.output.metadata, "w") as fd_out:
        print(
            DELIM.join(["read_id", *COLUMNS]),
            file=fd_out,
        )

        with pysam.FastxFile(snakemake.input.r1) as fh:
            for read in fh:
                read_acc = read.name.rsplit("-", maxsplit=1)[0]
                row = [read.name, "", "", "", "", "", ""]
                if read_acc not in UNMAPPED:
                    tax = truth.get(read_acc)
                    if tax is None:
                        if read_acc.startswith("kraken"):
                            taxid = read_acc.split("|")[1]

                            try:
                                lineage = taxtree.lineage(taxid)
                            except TaxonomyError:
                                lineage = fetch_taxonomy(taxid)

                            for l in lineage:
                                if l.rank == "strain":
                                    row[1] = l.id
                                    row[2] = l.name
                                elif l.rank == "species":
                                    row[3] = l.id
                                    row[4] = l.name
                                elif l.rank == "genus":
                                    row[5] = l.id
                                    row[6] = l.name
                        else:
                            raise KeyError(
                                f"Couldn't find accession {read_acc} for read {read.name} in accession truth"
                            )
                    else:
                        for i, c in enumerate(COLUMNS):
                            row[i + 1] = tax[c]

                print(DELIM.join(row), file=fd_out)

        with pysam.FastxFile(snakemake.input.r2) as fh:
            for read in fh:
                read_acc = read.name.rsplit("-", maxsplit=1)[0]
                row = [read.name, "", "", "", "", "", ""]
                if read_acc not in UNMAPPED:
                    tax = truth.get(read_acc)
                    if tax is None:
                        if read_acc.startswith("kraken"):
                            taxid = read_acc.split("|")[1]

                            try:
                                lineage = taxtree.lineage(taxid)
                            except TaxonomyError:
                                lineage = fetch_taxonomy(taxid)

                            for l in lineage:
                                if l.rank == "strain":
                                    row[1] = l.id
                                    row[2] = l.name
                                elif l.rank == "species":
                                    row[3] = l.id
                                    row[4] = l.name
                                elif l.rank == "genus":
                                    row[5] = l.id
                                    row[6] = l.name
                        else:
                            raise KeyError(
                                f"Couldn't find accession {read_acc} for read {read.name} in accession truth"
                            )
                    else:
                        for i, c in enumerate(COLUMNS):
                            row[i + 1] = tax[c]

                print(DELIM.join(row), file=fd_out)


main()
