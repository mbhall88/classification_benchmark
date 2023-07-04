import sys

sys.stderr = open(snakemake.log[0], "w")

import gzip
from pathlib import Path
from taxonomy import Taxonomy
import pysam

MTBC_TAXID = "77643"


def is_mtbc(lineage) -> bool:
    return any(l.id == MTBC_TAXID for l in lineage)


def eprint(msg):
    print(msg, file=sys.stderr)


def main():
    eprint("Loading taxonomy...")

    dump_dir = str(Path(snakemake.input.nodes).parent)
    taxtree = Taxonomy.from_ncbi(dump_dir)

    eprint("Loading gramtools assembly info...")

    gramtools = dict()
    with open(snakemake.input.lineage_info) as fd:
        _ = next(fd)
        for row in map(str.rstrip, fd):
            strain, lineage, genbank, taxid = row.split(",")
            gramtools[strain] = (lineage, genbank, taxid)

    eprint("Extracting TB sequences...")

    mtbc_file = snakemake.output.mtbc_fasta
    if mtbc_file.endswith(".gz"):
        mtbc_fd = gzip.open(mtbc_file, "wt")
    else:
        mtbc_fd = open(mtbc_file, "w")

    with pysam.FastxFile(snakemake.input.tb_asm) as fh:
        for entry in fh:
            if entry.name in gramtools:
                info = gramtools[entry.name]
                lineage = info[0]
                genbank = info[1]
                taxid = info[2]
            elif entry.name == "NC_000962.3":
                taxid = "83332"
                genbank = "AL123456"
                lineage = "4"
            else:
                continue

            header = f">{entry.name} lineage={lineage} genbank={genbank} taxid={taxid} {entry.comment}"
            print(header, file=mtbc_fd)
            print(entry.sequence.rstrip(), file=mtbc_fd)

    ntm_file = snakemake.output.ntm_fasta
    if ntm_file.endswith(".gz"):
        ntm_fd = gzip.open(ntm_file, "wt")
    else:
        ntm_fd = open(ntm_file, "w")

    eprint("Splitting mycobacteria file...")

    with pysam.FastxFile(snakemake.input.fasta) as fh:
        for entry in fh:
            taxid = entry.name.split("|")[1]
            if is_mtbc(taxtree.lineage(taxid)):
                print(str(entry).rstrip(), file=mtbc_fd)
            else:
                print(str(entry).rstrip(), file=ntm_fd)

    eprint("Done")


main()
