import gzip
import json
import sys
from dataclasses import dataclass
from typing import Optional, List

sys.stderr = open(snakemake.log[0], "w")

from taxonomy import Taxonomy, TaxonomyError
from Bio import Entrez
from pathlib import Path
from functools import cache
import re

Entrez.email = "michael.hall2@unimelb.edu.au"

DELIM = "\t"
MTB_TAXID = "1773"
HUMAN_TAXID = "9606"
TAXID_REGEX = re.compile(r"taxid=(?P<id>\d+)")


@cache
def accession2taxid(acc: str) -> str:
    handle = Entrez.esearch(db="nucleotide", term=acc)
    record = Entrez.read(handle)
    gi = record["IdList"][0]
    handle = Entrez.esummary(db="nucleotide", id=gi, retmode="json")
    result = json.load(handle)["result"]
    taxid = result[gi]["taxid"]
    return str(taxid)


@dataclass
class TaxonomyNode:
    id: str
    name: str
    parent: Optional["TaxonomyNode"]
    rank: str


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


def is_header(s: str) -> bool:
    return s[0] == ">"


def main():
    taxdir = str(Path(snakemake.input.nodes).parent)
    taxtree = Taxonomy.from_ncbi(taxdir)

    with open(snakemake.output.metadata, "w") as fd_out:
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
        for file in snakemake.input.refs:
            organism = Path(file).name.split(".")[0]
            if organism == "library":
                organism = Path(file).parts[-2].capitalize()
                if organism == "Viral":
                    organism = "Virus"
            elif organism == "hg02886":
                organism = "Human"


            if file.endswith("gz"):
                openf = gzip.open
            else:
                openf = open

            with openf(file, mode="rt") as fd_in:
                for header in filter(is_header, fd_in):
                    seqid = header[1:].split()[0]
                    taxid_match = TAXID_REGEX.search(header)
                    if taxid_match:
                        taxid = taxid_match.group("id")
                    elif organism == "Human":
                        taxid = HUMAN_TAXID
                    else:
                        taxid = seqid.split("|")[1]

                    row = [seqid, "", "", "", "", "", ""]

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
