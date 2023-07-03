import sys

sys.stderr = open(snakemake.log[0], "w")

import pysam
from collections import defaultdict
import random
from pathlib import Path
import re


def only_alphabetic(s: str) -> str:
    return re.sub(r"[^a-zA-Z]", "", s)


def extract_genus_from_comment(comment: str) -> str:
    fields = comment.split()
    genus = only_alphabetic(fields[0])
    if genus.isupper():
        genus = only_alphabetic(fields[1])

    return genus


def main():
    random.seed(88)

    include_genera = set(snakemake.params.include)
    genera_map = defaultdict(list)

    with pysam.FastxFile(snakemake.input.fasta) as fd:
        for entry in fd:
            genus = extract_genus_from_comment(entry.comment)
            if (
                genus not in include_genera
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
        genus_dir = outdir / genus
        genus_dir.mkdir()

        for header in headers:
            seqid = header.split()[0]
            p = genus_dir / f"{seqid}.fa"
            with open(p, mode="w") as fd:
                seq = fa.fetch(seqid)
                print(f">{header} circular=true", file=fd)
                print(seq, file=fd)


main()
