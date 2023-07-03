import sys
sys.stderr = open(snakemake.log[0], "w")

from taxonomy import Taxonomy
from Bio import Entrez, SeqIO
from pathlib import Path
from functools import cache

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
    record = Entrez.read(handle)
    if len(record["IdList"]) > 1:
        raise ValueError(
            f"Got more than one IdList for {acc}\n{record['IdList']}\n{record}"
        )
    taxid = record["IdList"][0]
    return taxid


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
            with open(file) as fd_in:
                for header in filter(is_header, fd_in):
                    seqid = header[1:].split()[0]
                    if seqid in MTB_ACCESSIONS:
                        taxid = MTB_TAXID
                    elif organism == "Human":
                        taxid = HUMAN_TAXID
                    elif organism == "NTM":
                        if "|" in seqid:
                            seqid = seqid.split("|")[-1]
                        taxid = accession2taxid(seqid)
                    else:
                        taxid = accession2taxid(seqid.split("|")[1])

                    row = [seqid, "", "", "", "", "", ""]

                    lineage = taxtree.lineage(taxid)
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
