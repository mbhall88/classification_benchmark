rule index_kraken_bacteria_library:
    input:
        fasta=rules.download_kraken_bacteria_db.output.fasta,
    output:
        idx=RESULTS / "kraken/db/library/bacteria/library.fna.fai",
    log:
        LOGS / "index_kraken_bacteria_library.log",
    resources:
        runtime="30m",
    container:
        CONTAINERS["samtools"]
    shell:
        "samtools faidx {input.fasta} 2> {log}"


rule filter_bacteria_assemblies:
    input:
        faidx=rules.index_kraken_bacteria_library.output.idx,
        fasta=rules.index_kraken_bacteria_library.input.fasta,
    output:
        outdir=temp(directory(RESULTS / "simulate/references/genera")),
    log:
        LOGS / "filter_bacteria_assemblies.log",
    resources:
        runtime="1h",
        mem_mb=int(4 * GB),
    container:
        CONTAINERS["pysam"]
    params:
        min_length=50_000,
        min_asm=10,
        max_asm=1_000,
        exclude=[
            "Mycobacterium",
            "Mycobacteroides",
            "Mycolicibacter",
            "Mycolicibacterium",
            "Mycolicibacillus",
        ],
    script:
        SCRIPTS / "filter_bacteria_assemblies.py"

# todo: cite https://doi.org/10.1093/gigascience/giac022
rule download_koref:
    output:
        fasta=RESULTS / "simulate/references/Human.fa.gz",
    log:
        LOGS / "download_koref.log",
    resources:
        runtime="30m",
    conda:
        ENVS / "download.yaml"
    params:
        url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/020/497/085/GCA_020497085.1_KOREF_S1v2.1/GCA_020497085.1_KOREF_S1v2.1_genomic.fna.gz",
        min_length=10_000
    shell:
        "(wget -nv {params.url} -O - | seqtk seq -L {params.min_length} - | gzip) > {output.fasta} 2> {log}"



rule download_held_out_mycobacterium:
    output:
        genomes=temp(directory(RESULTS / "simulate/reference/mycobacteria/2023-08-09/")),
    log:
        LOGS / "download_held_out_mycobacterium.log",
    resources:
        runtime="20m",
        mem_mb=int(2 * GB),
    threads: 2
    conda:
        ENVS / "genome_updater.yaml"
    params:
        opts='-m -a -f "genomic.fna.gz" -g "bacteria" -d "refseq,genbank"',
        filter=rules.prepare_mycobacterium_for_kraken.params.exclude,
        outdir=lambda wildcards, output: Path(output.genomes).parent,
        version=lambda wildcards, output: Path(output.genomes).parts[-1],
    shell:
        """
        genome_updater.sh {params.opts} -o {params.outdir} -b {params.version} \
            -t {threads} -F "1:{params.filter}" &> {log}
        """


rule split_mycobacteria:
    input:
        genomes=rules.download_held_out_mycobacterium.output.genomes,
    output:
        ntm_fasta=RESULTS / "simulate/references/NTM.fa.gz",
        mtbc_fasta=RESULTS / "simulate/references/MTBC.fa.gz",
    log:
        LOGS / "split_mycobacteria.log",
    resources:
        runtime="5m",
        mem_mb=int(2 * GB),
    conda:
        ENVS / "split_mycobacteria.yaml"
    script:
        SCRIPTS / "split_mycobacteria.py"


rule reduce_bacteria_assemblies:
    input:
        indir=rules.filter_bacteria_assemblies.output.outdir,
        script=SCRIPTS / "dereplicator.py",
    output:
        outdir=temp(directory(RESULTS / "simulate/references/refined_genera")),
    log:
        LOGS / "reduce_bacteria_assemblies.log",
    threads: 16
    resources:
        runtime="6h",
        mem_mb=int(16 * GB),
    conda:
        ENVS / "derep.yaml"
    params:
        opts="-f 0.1",
    shell:
        """
        exec &> {log}
        for dir in {input.indir}/*/
        do
            dir=${{dir%/}}
            genus=${{dir##*/}}
            echo "Dereplicating $genus"
            outdir={output.outdir}/$genus
            python {input.script} {params.opts} --threads {threads} $dir $outdir
        done
        """


rule combine_bacteria_assemblies:
    input:
        asmdir=rules.reduce_bacteria_assemblies.output.outdir,
    output:
        fasta=RESULTS / "simulate/references/Bacteria.fa.gz",
    log:
        LOGS / "combine_bacteria_assemblies.log",
    resources:
        mem_mb=GB,
        runtime="2h",
    container:
        CONTAINERS["rs_utils"]
    shell:
        "fd -e fa -X gzip -c \; . {input.asmdir} 2> {log} > {output.fasta}"


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

    if wildcards.read_type in ("NTM", "MTBC", "Virus"):
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


def infer_simulate_input(read_type):
    if read_type == "Unmapped":
        # just pick smallest fasta as we wont actually take sequences from it
        return str(RESULTS / "simulate/references/MTBC.fa.gz")
    elif read_type in ("NTM", "MTBC", "Bacteria", "Human"):
        return str(RESULTS / f"simulate/references/{read_type}.fa.gz")
    elif read_type == "Virus":
        return str(RESULTS / f"kraken/db/library/viral/library.fna")
    else:
        raise NotImplementedError(f"Don't know what reference for {read_type}")


simulate_mem = {
    "Bacteria": 128 * GB,
    "Human": 128 * GB,
    "Virus": 16 * GB,
    "MTBC": 8 * GB,
    "NTM": 16 * GB,
    "Unmapped": 8 * GB,
}


rule simulate_nanopore_reads:
    input:
        reference=lambda wildcards: infer_simulate_input(wildcards.read_type),
    output:
        fastq=RESULTS / "simulate/reads/{read_type}.ont.fq.gz",
    log:
        LOGS / "simulate_nanopore_reads/{read_type}.log",
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * simulate_mem[wildcards.read_type],
        runtime="4h",
    container:
        CONTAINERS["badread"]
    params:
        opts=simulate_options,
        total_bases=calculate_total_bases,
    shadow:
        "shallow"
    shell:
        """
        exec 2> {log}
        bases=$(python -c "print(int({params.total_bases} / {threads}))")
        tmpd=$(mktemp -d)
        for p in $(seq 1 {threads}); do
            reads="$tmpd/$p.fq"
            badread simulate --reference {input.reference} --quantity $bases {params.opts} > $reads &
        done
        wait
        cat $tmpd/*.fq | gzip > {output.fastq}
        rm -rf $tmpd
        """


rule combine_simulated_reads:
    input:
        fastqs=[
            RESULTS / f"simulate/reads/{read_type}.ont.fq.gz"
            for read_type in config["simulate"]["proportions"]
            if read_type != "Unmapped"
        ],
        script=SCRIPTS / "filter_ambig.py",
    output:
        reads=RESULTS / "simulate/reads/metagenome.ont.fq.gz",
    log:
        LOGS / "combine_simulated_reads.log",
    resources:
        runtime="30m",
    conda:
        ENVS / "combine_simulated_reads.yaml"
    params:
        opts="-L 500",
        max_ambig=0,
    shell:
        "(zcat {input.fastqs} | seqtk seq {params.opts} - | python {input.script} - {params.max_ambig} | gzip) > {output.reads} 2> {log}"


ILLUMINA_READ_LENGTH = 150


def count_contigs(path):
    if str(path).endswith("gz"):
        import gzip

        fopen = gzip.open
    else:
        fopen = open

    with fopen(path, mode="rt") as fp:
        return sum(1 for l in fp if l[0] == ">")


def calculate_number_reads(wildcards, input):
    quantity = config["simulate"]["quantity"]
    n_bases = parse_size(quantity) / 5  # use a fifth of the number of bases in nanopore

    proportions = config["simulate"]["proportions"]
    prop = proportions.get(wildcards.read_type)
    if prop is None:
        raise ValueError(f"No proportion for {wildcards.read_type}")
    required_bases = n_bases * prop
    total_reads = required_bases / ILLUMINA_READ_LENGTH / 2  # divide by 2 as paired
    return total_reads


rule simulate_illumina_reads:
    input:
        reference=lambda wildcards: infer_simulate_input(wildcards.read_type),
    output:
        r1=RESULTS / "simulate/reads/{read_type}_R1.illumina.fq.gz",
        r2=RESULTS / "simulate/reads/{read_type}_R2.illumina.fq.gz",
    log:
        LOGS / "simulate_illumina_reads/{read_type}.log",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(8 * GB),
        runtime="1h",
    conda:
        ENVS / "art.yaml"
    params:
        opts=f"-na -ss MSv3 -l {ILLUMINA_READ_LENGTH} -s 10 -m {int(ILLUMINA_READ_LENGTH + 100)} -p",
        n_reads=calculate_number_reads,
    shadow:
        "shallow"
    shell:
        """
        exec 2> {log}
        if [ {wildcards.read_type} == "Unmapped" ]; then
            ref=$(mktemp -u)
            fastaq make_random_contigs 1 100000000 $ref
            sed -i "1s/.*/>random/" $ref
        elif file {input.reference} | grep -q compressed; then
            ref=$(mktemp -u)
            gzip -dc {input.reference} > $ref
        else
            ref={input.reference}
        fi
        n_contigs=$(grep -c '^>' $ref)
        n_reads=$(python -c "print(max(round(int({params.n_reads}) / int($n_contigs)), 1))")
        outprefix=$(mktemp -u)
            
        art_illumina {params.opts} -i $ref -o $outprefix -c $n_reads
        gzip -c "${{outprefix}}1.fq" > {output.r1}
        gzip -c "${{outprefix}}2.fq" > {output.r2}
        """


rule combine_illumina_simulated_reads:
    input:
        r1s=sorted(
            [
                RESULTS / f"simulate/reads/{read_type}_R1.illumina.fq.gz"
                for read_type in config["simulate"]["proportions"]
                if read_type != "Unmapped"
            ]
        ),
        r2s=sorted(
            [
                RESULTS / f"simulate/reads/{read_type}_R2.illumina.fq.gz"
                for read_type in config["simulate"]["proportions"]
                if read_type != "Unmapped"
            ]
        ),
        script=SCRIPTS / "filter_ambig.py",
    output:
        r1=RESULTS / "simulate/reads/metagenome_R1.illumina.fq.gz",
        r2=RESULTS / "simulate/reads/metagenome_R2.illumina.fq.gz",
    log:
        LOGS / "combine_illumina_simulated_reads.log",
    resources:
        runtime="30m",
    conda:
        ENVS / "combine_illumina_simulated_reads.yaml"
    params:
        max_ambig=0,
    shell:
        """
        exec 2> {log}
        tmp1=$(mktemp -u --suffix=_1.fq)
        tmp2=$(mktemp -u --suffix=_2.fq)
        (zcat {input.r1s} | python {input.script} - {params.max_ambig}) > $tmp1
        (zcat {input.r2s} | python {input.script} - {params.max_ambig}) > $tmp2
        td=$(mktemp -d)
        seqkit pair --id-regexp '^(\S+)\/[12]' -1 $tmp1 -2 $tmp2 -O $td
        gzip -c $td/$(basename $tmp1) > {output.r1}
        gzip -c $td/$(basename $tmp2) > {output.r2}
        """
