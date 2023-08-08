import sys
import gzip
import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("indir")
    parser.add_argument("-s", "--asm-summary")
    parser.add_argument("-T", "--taxon", help="Taxon ID to apply to all assemblies")
    parser.add_argument(
        "-x", "--exclude", help="Comma-separated list of assembly accessions to exclude"
    )
    parser.add_argument(
        "-r",
        "--recursive",
        help="Search indir recursively for assemblies",
        action="store_true",
    )
    parser.add_argument("-o", "--output", default="-")

    args = parser.parse_args()

    acc2file = dict()

    if args.recursive:
        files = Path(args.indir).rglob("*.fna*")
    else:
        files = Path(args.indir).glob("*.fna*")

    exclude = set()
    if args.exclude:
        exclude = {s.strip() for s in args.exclude.split(",")}

    for p in files:
        acc = "_".join(p.name.split("_", maxsplit=2)[:2])
        assert acc not in acc2file, p
        acc2file[acc] = p.resolve()

    acc2tax = dict()
    asm_summary = Path(sys.argv[2])
    if not asm_summary.exists():
        all_tax = sys.argv[2]
    else:
        with open(sys.argv[2]) as fp:
            for line in fp:
                fields = line.split("\t")
                acc = fields[0]
                assert acc in acc2file, acc
                assert acc not in acc2tax
                taxid = fields[5]
                acc2tax[acc] = taxid

    n_processed = 0

    if args.output == "-":
        out_fp = sys.stdout
    else:
        p = Path(args.output)
        if p.suffix == ".gz":
            out_fp = gzip.open(p, mode="wt")
        else:
            out_fp = open(p, mode="w")

    for acc, p in acc2file.items():
        n_processed += 1
        if n_processed % 50 == 0:
            print(f"Processed {n_processed}", file=sys.stderr)

        if acc in exclude:
            print(f"Exclusing {acc}", file=sys.stderr)
            continue
        if acc2tax:
            taxid = acc2tax[acc]
        else:
            taxid = all_tax

        if p.suffix == ".gz":
            fp = gzip.open(p, mode="rt")
        else:
            fp = open(p)

        for line in map(str.rstrip, fp):
            if line[0] == ">":
                contig = line.split()[0][1:].split("_")[-1]
                contig = f"kraken:taxid|{taxid}|{contig}"
                comment = " ".join(line.split()[1:])
                header = f">{contig} {comment}"
                print(header, file=out_fp)
            else:
                print(line, file=out_fp)

        fp.close()

    out_fp.close()


if __name__ == "__main__":
    main()
