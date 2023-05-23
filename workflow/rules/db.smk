from humanfriendly import parse_size


def infer_kraken_memory(wildcards, attempt):
    if wildcards.size == "full":
        mb = 40 * GB
    else:
        bytes = parse_size(wildcards.size)
        mb = bytes / (10**6) + (4 * GB)  # add 4gb for wriggle room
    return int(mb) * attempt


def infer_max_db_size_opt(wildcards):
    if wildcards.size == "full":
        return ""
    else:
        return f"--max-db-size {parse_size(wildcards.size)}"


def infer_minimizer_spaces(wildcards):
    l = int(wildcards.l)
    s = int(l / 4)
    return f"--minimizer-spaces {s}"


rule build_kraken_database:
    output:
        db=directory(RESULTS / "kraken/db/k{k}/l{l}/{size}"),
    log:
        LOGS / "build_kraken_database/k{k}/l{l}/{size}.log",
    threads: 16
    resources:
        mem_mb=infer_kraken_memory,
        runtime="1d",
    benchmark:
        BENCH / "kraken/build/k{k}/l{l}/{size}.tsv"
    params:
        opts="--standard --kmer-len {k} --minimizer-len {l}",
        max_db_size=infer_max_db_size_opt,
        spaces=infer_minimizer_spaces,
    container:
        CONTAINERS["kraken"]
    shell:
        "kraken2-build {params.opts} {params.max_db_size} {params.spaces} --threads {threads} --db {output.db} &> {log}"



rule build_decontamination_db:
    output:
        fasta=RESOURCES / "decontamination/remove_contam.fa.gz",
        metadata=RESOURCES / "decontamination/remove_contam.tsv",
    params:
        script=SCRIPTS / "download_tb_reference_files.pl",
        outdir=lambda wildcards, output: Path(output.fasta).parent,
    container:
        CONTAINERS["clockwork"]
    log:
        LOGS / "build_decontamination_db.log",
    shadow:
        "shallow"
    shell:
        "perl {params.script} {params.outdir} &> {log}"