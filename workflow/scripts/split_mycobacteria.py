import sys

sys.stderr = open(snakemake.log[0], "w")

import gzip
from pathlib import Path
from taxonomy import Taxonomy
import pysam

MTB_TAXID = "1773"


def is_mtbc(lineage) -> bool:
    return any(l.id == MTBC_TAXID for l in lineage)


def eprint(msg):
    print(msg, file=sys.stderr)


def main():
    eprint("Loading summary...")

    asm_dir = Path(snakemake.input.genomes)
    summary_file = asm_dir / "assembly_summary.txt"

    if not summary_file.exists():
        raise FileNotFoundError(f"{summary_file} does not exist")

    acc2tax = dict()
    with open(summary_file) as fp:
        for line in fp:
            fields = line.split("\t")
            acc = fields[0]
            taxid = fields[6]
            acc2tax[acc] = taxid

    asm_files = list(asm_dir.rglob("*.fna*"))

    if len(acc2tax) != len(asm_files):
        raise ValueError(f"Different number of assembly files {asm_files} and accessions {acc2tax}")

    eprint("Extracting sequences...")

    mtbc_file = snakemake.output.mtbc_fasta
    if mtbc_file.endswith(".gz"):
        mtbc_fd = gzip.open(mtbc_file, "wt")
    else:
        mtbc_fd = open(mtbc_file, "w")

    ntm_file = snakemake.output.ntm_fasta
    if ntm_file.endswith(".gz"):
        ntm_fd = gzip.open(ntm_file, "wt")
    else:
        ntm_fd = open(ntm_file, "w")

    n_written = 0
    for p in asm_files:
        for acc, taxid in acc2tax.items():
            if acc in p.name:
                if taxid == MTB_TAXID:
                    out_fp = mtbc_fd
                else:
                    out_fp = ntm_fd

                with pysam.FastxFile(str(p)) as fh:
                    for entry in fh:
                        header = f">{entry.name} taxid={taxid} {entry.comment}"
                        print(header, file=out_fp)
                        print(entry.sequence.rstrip(), file=out_fp)
                n_written += 1
                break

    if n_written != len(acc2tax):
        raise ValueError(f"Wrote different amount of asms ({n_written}) to expected ({len(acc2tax)})")

    eprint("Done")


main()
