import gzip
import sys
from dataclasses import dataclass
from typing import Optional, List

sys.stderr = open(snakemake.log[0], "w")

from taxonomy import Taxonomy, TaxonomyError
from Bio import Entrez, SeqIO
from pathlib import Path
from functools import cache
from urllib.error import HTTPError

Entrez.email = "michael.hall2@unimelb.edu.au"

DELIM = "\t"
MTB_TAXID = "1773"
HUMAN_TAXID = "9606"
MTB_ACCESSIONS = {
    "MTB_an",
    "N000",
    "N003",
    "N005",
    "N005",
    "N007",
    "N009",
    "N013",
    "N014",
    "N015",
    "N015",
    "N015",
    "N117",
    "N117",
    "N120",
    "N121",
    "N127",
    "N128",
}


@cache
def accession2taxid(acc: str) -> str:
    handle = Entrez.efetch(db="nucleotide", id=acc, retmode="text", rettype="gb")
    record = SeqIO.read(handle, "gb")
    handle = Entrez.esearch(db="taxonomy", term=record.annotations["organism"])
    record = handle.read()
    if len(record["IdList"]) > 1:
        raise ValueError(
            f"Got more than one IdList for {acc}\n{record['IdList']}\n{record}"
        )
    taxid = record["IdList"][0]
    return taxid


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

            if file.endswith("gz"):
                openf = gzip.open
            else:
                openf = open

            with openf(file, mode="rt") as fd_in:
                for header in filter(is_header, fd_in):
                    seqid = header[1:].split()[0]
                    if seqid in MTB_ACCESSIONS:
                        taxid = MTB_TAXID
                    elif organism == "Human":
                        taxid = HUMAN_TAXID
                    elif organism == "NTM":
                        if "|" in seqid:
                            seqid = seqid.split("|")[-1]
                        try:
                            taxid = accession2taxid(seqid)
                        except Exception as err:
                            print(
                                f"Failed to fetch taxid ({taxid}) for {seqid} from {organism} {header}",
                                file=sys.stderr,
                            )
                            raise err
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
