organisms = ["Human", "Virus", "Bacteria", "NTM", "TB"]


rule separate_db_by_organism:
    input:
        metadata=rules.combine_references.output.metadata,
        fasta=rules.combine_references.output.fasta,
        faidx=rules.faidx_db.output.faidx,
    output:
        refs=[RESULTS / f"simulate/references/{org}.fa.gz" for org in organisms],
    log:
        LOGS / "separate_db_by_organism.log",
    resources:
        runtime="1h",
        mem_mb=int(4 * GB),
    container:
        CONTAINERS["pysam"]
    script:
        SCRIPTS / "separate_db_by_organism.py"


def simulate_options(wildcards):
    opts = []
    if wildcards.read_type == "Unmapped":
        junk = 25
        random_reads = 75
        chimera = 1.0
    else:
        junk = 0
        random_reads = 0
        chimera = 0.5
    opts.append(f"--junk_reads {junk}")
    opts.append(f"--random_reads {random_reads}")
    opts.append(f"--chimeras {chimera}")

    if wildcards.read_type in ("NTM", "TB"):
        length = "4000,3000"
        opts.append(f"--length {length}")

    return " ".join(opts)


def calculate_total_bases(wildcards):
    quantity = config["simulate"]["quantity"]
    n_bases = parse_size(quantity)

    proportions = config["simulate"]["proportions"]
    prop = proportions.get(wildcards.read_type)
    if prop is None:
        raise ValueError(f"No proportion for {wildcards.read_type}")
    return n_bases * prop


def infer_simulate_input(wildcards):
    if wildcards.read_type == "Unmapped":
        # just pick smallest fasta as we wont actually take sequences from it
        return str(RESULTS / "simulate/references/TB.fa.gz")
    else:
        return str(RESULTS / f"simulate/references/{wildcards.read_type}.fa.gz")


rule simulate_nanopore_reads:
    input:
        reference=infer_simulate_input,
    output:
        fastq=RESULTS / "simulate/reads/{read_type}.fq.gz",
    log:
        LOGS / "simulate_nanopore_reads/{read_type}.log",
    threads: 8
    resources:
        mem_mb=int(8 * GB),
        runtime="2d",
    container:
        CONTAINERS["badread"]
    params:
        opts=simulate_options,
        total_bases=calculate_total_bases,
    shell:
        """
        exec 2> {log}
        bases=$(({params.total_bases} / {threads}))
        for p in $(seq 1 {threads}); do
            reads=temp_"$p".fq
            badread simulate --reference {input.reference} --quantity $bases {opts} > $reads &
        done
        wait
        cat temp_*.fq | gzip > {output.fastq}
        rm temp_*.fq
        """


rule combine_simulated_reads:
    input:
        fastqs=[
            RESULTS / f"simulate/reads/{read_type}.fq.gz"
            for read_type in config["simulate"]["proportions"]
        ],
    output:
        reads=RESULTS / "simulate/reads/metagenome.fq.gz",
    log:
        LOGS / "combine_simulated_reads.log",
    resources:
        runtime="5m",
    shell:
        "cat {input.fastqs} > {output.reads} 2> {log}"