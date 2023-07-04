import sys

sys.stderr = open(snakemake.log[0], "w")
from pathlib import Path
from pysam import FastaFile
import gzip
from collections import Counter, defaultdict


def eprint(msg):
    print(msg, file=sys.stderr)


def invert_dict(d: dict) -> dict:
    d_inv = defaultdict(list)
    for k, v in d.items():
        d_inv[v].append(k)

    return d_inv


def load_metadata(path: str) -> dict[str, list[str]]:
    """Maps organism->list of reference names"""
    metadata: dict[str, str] = {}
    organisms = []
    with open(path) as fd:
        for line in map(str.rstrip, fd):
            organism, _, seqname = line.split("\t")
            if seqname in metadata:
                raise KeyError(f"{seqname} is in reference multiple times")
            metadata[seqname] = organism
            organisms.append(organism)

    eprint(f"Got the following organism group counts:\n{Counter(organism)}")
    eprint(f"Loaded {len(metadata)} sequences from reference")
    return invert_dict(metadata)


def main():
    fa_reader = FastaFile(snakemake.input.fasta)

    eprint("Loading metadata...")
    metadata = load_metadata(snakemake.input.metadata)

    for outpath in snakemake.output.refs:
        if outpath.endswith(".gz"):
            openf = gzip.open
        else:
            openf = open

        organism = Path(outpath).name.split(".")[0]
        eprint(f"Extracting {organism} sequences...")

        with openf(outpath, mode="wt") as fd_out:
            contig_names = metadata[organism]
            for contig in contig_names:
                seq = fa_reader.fetch(contig)
                header = f">{contig}"
                if organism in ("Bacteria", "TB", "NTM"):
                    header += " circular=true"
                print(header, file=fd_out)
                print(seq, file=fd_out)

    eprint("Finished")


main()
