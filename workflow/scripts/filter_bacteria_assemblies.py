import sys

sys.stderr = open(snakemake.log[0], "w")

import pysam
from collections import defaultdict
import random
from pathlib import Path


def main():
    random.seed(88)

    exclude_genera = set(snakemake.params.exclude)
    genera_map = defaultdict(list)

    with pysam.FastxFile(snakemake.input.fasta) as fd:
        for entry in fd:
            genus = entry.comment.split()[0]
            if (
                genus in exclude_genera
                or len(entry.sequence) < snakemake.params.min_length
            ):
                continue
            genera_map[genus].append(f"{entry.name} {entry.comment}")

    for genus in genera_map.copy():
        ids = genera_map[genus]
        if len(ids) < snakemake.params.min_asm:
            del genera_map[genus]
        elif len(ids) > snakemake.params.max_asm:
            sub = random.sample(ids, k=snakemake.params.max_asm)
            genera_map[genus] = sub

    print(
        f"{len(genera_map)} genera passed filtering. Extracting sequences...",
        file=sys.stderr,
    )

    outdir = Path(snakemake.output.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    fa = pysam.FastaFile(snakemake.input.fasta)

    for genus, headers in genera_map.items():
        with open(outdir / f"{genus}.fa", "w") as fd:
            for header in headers:
                seqid = header.split()[0]
                seq = fa.fetch(seqid)
                print(f">{header} circular=true", file=fd)
                print(seq, file=fd)


main()
