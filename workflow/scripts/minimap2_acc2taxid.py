import gzip
import json
import re
import sys
from dataclasses import dataclass
from functools import cache
from typing import Optional, List

sys.stderr = open(snakemake.log[0], "w")

from pathlib import Path
from taxonomy import Taxonomy, TaxonomyError
from Bio import Entrez


Entrez.email = "michael.hall2@unimelb.edu.au"
DELIM = "\t"
MTB_TAXID = "1773"
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
def accession2taxid(acc: str) -> str:
    handle = Entrez.esearch(db="nucleotide", term=acc)
    record = Entrez.read(handle)
    gi = record["IdList"][0]
    handle = Entrez.esummary(db="nucleotide", id=gi, retmode="json")
    result = json.load(handle)["result"]
    taxid = result[gi]["taxid"]
    return str(taxid)


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

    with open(snakemake.input.metadata) as fd_in, open(
        snakemake.output.truth, "w"
    ) as fd_out:
        print(
            DELIM.join(
                [
                    "accession",
                    "strain_id",
                    "strain",
                    "species_id",
                    "species",
                    "genus_id",
                    "genus",
                ]
            ),
            file=fd_out,
        )

        for row in map(str.rstrip, fd_in):
            category, _, identifier = row.split("\t")
            if category == "Bacteria":
                continue
            if "|" in identifier:
                acc = identifier.split("|")[-1]
            else:
                acc = identifier

            try:
                taxid = accession2taxid(acc)
            except Exception as err:
                print(f"Failed to get taxid for {identifier}", file=sys.stderr)
                raise err

            row = [identifier, "", "", "", "", "", ""]

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

            print(DELIM.join(row), file=fd_out)

        with gzip.open(snakemake.input.fasta, mode="rt") as fasta_fd:
            for line in fasta_fd:
                if not line[0] == ">":
                    continue
                identifier = line[1:].split()[0]
                m = re.match(r"^N\d{4}$", identifier)
                if m:
                    taxid = MTB_TAXID
                else:
                    try:
                        taxid = accession2taxid(acc)
                    except Exception as err:
                        print(f"Failed to get taxid for {identifier}", file=sys.stderr)
                        raise err

                row = [identifier, "", "", "", "", "", ""]

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

                print(DELIM.join(row), file=fd_out)


main()
