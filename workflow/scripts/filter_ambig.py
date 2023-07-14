import pysam
import sys


def ambig_prop(seq: str) -> float:
    n_acgt = sum(1 for c in seq if c.upper() not in "ACGT")
    return n_acgt / len(seq)


def main():
    if len(sys.argv) < 3:
        raise ValueError("USAGE: filter_ambig.py <fastq> <max. ambig %>")
    infile = sys.argv[1]
    max_ambig = float(sys.argv[2])

    with pysam.FastxFile(infile) as fh:
        for read in fh:
            if ambig_prop(read.sequence) > max_ambig:
                continue
            sys.stdout.write(str(read))


if __name__ == "__main__":
    main()
