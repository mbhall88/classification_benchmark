import sys

sys.stderr = open(snakemake.log[0], "w")

import gzip
from pathlib import Path


def main():
    metadata = dict()
    groups = set()
    with open(snakemake.input.metadata) as fp:
        for row in map(str.rstrip, fp):
            if not row:
                continue
            group, is_contam, contig = row.split("\t")
            metadata[contig] = (group, int(is_contam))
            groups.add(group)

    print(f"{len(metadata)} contigs in metadata", file=sys.stderr)

    groups = sorted(groups)
    # move TB group to be the first group
    groups.insert(0, groups.pop(groups.index("TB")))

    outdir = Path(snakemake.output.fasta_dir).resolve()
    outdir.mkdir(exist_ok=True, parents=True)

    counter = 0

    with gzip.open(snakemake.input.fasta, "rt") as fp:
        seq = ""
        name = ""
        for line in fp:
            if line.startswith(">"):
                if name:
                    assert seq
                    p = outdir / f"{name}.fa"
                    with open(p, "w") as fa:
                        print(f">{name}", file=fa)
                        fa.write(seq)
                        counter += 1
                    seq = ""

                name = line.split()[0][1:]
                if name not in metadata:
                    raise KeyError(f"Missing contig {name} from metadata file")
            else:
                seq += line

        assert seq
        assert name
        p = outdir / f"{name}.fa"
        with open(p, "w") as fa:
            print(f">{name}", file=fa)
            fa.write(seq)
            counter += 1

    print(f"Wrote {counter} fasta files", file=sys.stderr)

    # write file list
    file_list = []
    for contig, (group, _is_contam) in metadata.items():
        p = outdir / f"{contig}.fa"
        assert p.exists(), p
        i = groups.index(group) + 1
        file_list.append((i, p, group))

    with open(snakemake.output.file_list, "w") as fd:
        for i, p, group in sorted(file_list):
            p = outdir / f"{contig}.fa"
            assert p.exists(), p
            i = groups.index(group) + 1
            print(f"{str(p)} {i} {group}", file=fd)


main()
