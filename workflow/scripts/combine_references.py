import sys

sys.stderr = open(snakemake.log[0], "w")

import gzip
import shutil


def main():
    out_fasta = snakemake.output.fasta
    out_metadata = snakemake.output.metadata
    shutil.copy2(snakemake.input.decontam_metadata, out_metadata)

    add_metadata = []

    # add ancestral genome
    is_contam = "0"
    group = "TB"
    with gzip.open(snakemake.input.ancestral, "rt") as f_in:
        for line in map(str.rstrip, f_in):
            if line.startswith(">"):
                name = line.split()[0][1:]
                add_metadata.append([group, is_contam, name])

    # do the same thing for the gramtools assemblies
    for file in snakemake.input.gramtools_asms:
        with gzip.open(file, "rt") as f_in:
            for line in map(str.rstrip, f_in):
                if line.startswith(">"):
                    name = line.split()[0][1:]
                    add_metadata.append([group, is_contam, name])

    with open(out_metadata, "a") as metadata_fp:
        for row in add_metadata:
            print("\t".join(row), file=metadata_fp)

    # combine all of the fastas into one
    with open(out_fasta, "wb") as wfd:
        for f in [
            snakemake.input.decontam_fa,
            snakemake.input.ancestral,
            *snakemake.input.gramtools_asms,
        ]:
            if f.endswith(".gz"):
                openf = gzip.open
            else:
                openf = open
            with openf(f, "rb") as fd:
                shutil.copyfileobj(fd, wfd)


main()
