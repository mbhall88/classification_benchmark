use rule combine_simulated_reads as extract_dehumanised_ont_reads with:
    input:
        fastqs=[
            RESULTS / f"simulate/reads/{read_type}.ont.fq.gz"
            for read_type in config["simulate"]["proportions"]
            if read_type not in ("Unmapped", "Human")
        ],
    output:
        reads=RESULTS / "dehumanise/metagenome.dehumanised.ont.fq.gz",
    log:
        LOGS / "extract_dehumanised_ont_reads.log",


use rule combine_illumina_simulated_reads as extract_dehumanised_illumina_reads with:
    input:
        r1s=sorted(
            [
                RESULTS / f"simulate/reads/{read_type}_R1.illumina.fq.gz"
                for read_type in config["simulate"]["proportions"]
                if read_type not in ("Unmapped", "Human")
            ]
        ),
        r2s=sorted(
            [
                RESULTS / f"simulate/reads/{read_type}_R2.illumina.fq.gz"
                for read_type in config["simulate"]["proportions"]
                if read_type not in ("Unmapped", "Human")
            ]
        ),
    output:
        r1=RESULTS / "dehumanise/metagenome_R1.dehumanised.illumina.fq.gz",
        r2=RESULTS / "dehumanise/metagenome_R2.dehumanised.illumina.fq.gz",
    log:
        LOGS / "extract_dehumanised_illumina_reads.log"


def infer_classify_reads(wildcards):
    if wildcards.tech == "ont":
        return RESULTS / "dehumanise/metagenome.dehumanised.ont.fq.gz"
    elif wildcards.tech == "illumina":
        return [
            RESULTS / f"dehumanise/metagenome_R{i}.dehumanised.illumina.fq.gz"
            for i in [1, 2]
        ]
    else:
        raise ValueError(f"Don't recognise tech {wildcards.tech}")


rule remove_held_out_assemblies_clockwork:
    input:
        db=rules.faidx_db.output.fasta,
        ntm=rules.split_mycobacteria.output.ntm_fasta,
        mtb=rules.split_mycobacteria.output.mtbc_fasta,
    output:
        db=RESULTS / "classify/minimap2/db/db.fa.gz",
    log:
        LOGS / "remove_held_out_assemblies_clockwork.log",
    resources:
        mem_mb=int(4 * GB),
        runtime="10m",
    container:
        CONTAINERS["seqkit"]
    shell:
        """
        exec 2> {log}
        rmids=$(mktemp -u)
        seqkit seq -n -i {input.mtb} > $rmids
        seqkit seq -n -i {input.ntm} >> $rmids
        seqkit grep -v -f $rmids -o {output.db} {input.db}
        """


rule index_heldout_clockwork_with_minimap2:
    input:
        fasta=rules.remove_held_out_assemblies_clockwork.output.db,
    output:
        index=RESULTS / "classify/minimap2/db/db.{preset}.mmi",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(32 * GB),
        runtime="5m",
    log:
        LOGS / "index_heldout_clockwork_with_minimap2/{preset}.log",
    threads: 4
    container:
        CONTAINERS["minimap2"]
    params:
        opts="-x {preset} -I 12G",
    shell:
        "minimap2 {params.opts} -t {threads} -d {output.index} {input.fasta} 2> {log}"


rule remove_held_out_assemblies_mycobacterium:
    input:
        db=rules.create_minimap2_mycobacterium_db.output.db,
        ntm=rules.split_mycobacteria.output.ntm_fasta,
        mtb=rules.split_mycobacteria.output.mtbc_fasta,
    output:
        db=RESULTS / "classify/minimap2/db/Mycobacterium.rep.fa.gz",
    log:
        LOGS / "remove_held_out_assemblies_mycobacterium.log",
    resources:
        mem_mb=int(4 * GB),
        runtime="10m",
    container:
        CONTAINERS["seqkit"]
    shell:
        """
        exec 2> {log}
        rmids=$(mktemp -u)
        seqkit seq -n -i {input.mtb} > $rmids
        seqkit seq -n -i {input.ntm} >> $rmids
        seqkit grep -v -f $rmids -o {output.db} {input.db}
        """


rule index_heldout_mycobacterium_db_minimap2:
    input:
        fasta=rules.remove_held_out_assemblies_mycobacterium.output.db,
    output:
        index=RESULTS / "classify/minimap2/db/Mycobacterium.rep.{preset}.mmi",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(16 * GB),
        runtime="5m",
    log:
        LOGS / "index_heldout_mycobacterium_db_minimap2/{preset}.log",
    threads: 4
    container:
        CONTAINERS["minimap2"]
    params:
        opts="-x {preset} -I 12G",
    shell:
        "minimap2 {params.opts} -t {threads} -d {output.index} {input.fasta} 2> {log}"


def infer_minimap2_db(wildcards):
    preset = PRESETS[wildcards.tech]
    if wildcards.db == "clockwork":
        return RESULTS / f"classify/minimap2/db/db.{preset}.mmi"
    elif wildcards.db == "mtbc":
        return RESULTS / f"db/GTDB_genus_Mycobacterium/MTB.{preset}.mmi"
    elif wildcards.db == "mycobacterium":
        return RESULTS / f"classify/minimap2/db/Mycobacterium.rep.{preset}.mmi"
    else:
        raise ValueError(f"Don't recognise db {wildcards.db}")


PRESETS = {"ont": "map-ont", "illumina": "sr"}


rule minimap2_classify:
    input:
        reads=infer_classify_reads,
        db=infer_minimap2_db,
    output:
        aln=RESULTS / "classify/minimap2/aln.{db}.{tech}.paf",
    log:
        LOGS / "minimap2_classify/{db}/{tech}.log",
    resources:
        runtime="2h",
        mem_mb=int(20 * GB),
    threads: 4
    container:
        CONTAINERS["minimap2"]
    benchmark:
        repeat(BENCH / "classify/minimap2/{db}/{tech}.tsv", config["benchmark_repeats"])
    params:
        opts="--secondary=no -c",
        preset=lambda wildcards: PRESETS[wildcards.tech],
    shell:
        """
        minimap2 {params.opts} -x {params.preset} -o {output.aln} -t {threads} \
            {input.db} {input.reads} 2> {log}
        """


CLASSIFY_DBS = {
    "standard": RESULTS / "db/kraken/standard",
    "standard-8": RESULTS / "db/kraken/standard-8",
    "mycobacterium": RESULTS / "db/GTDB_genus_Mycobacterium/kraken/db/hash.k2d",
}


rule kraken_classify:
    input:
        reads=infer_classify_reads,
        db=lambda wildcards: CLASSIFY_DBS[wildcards.db],
    output:
        report=RESULTS / "classify/kraken/{db}.{tech}.k2report",
        out=RESULTS / "classify/kraken/{db}.{tech}.k2",
    log:
        LOGS / "kraken_classify/{db}/{tech}.log",
    threads: 4
    resources:
        runtime="20m",
        mem_mb=lambda wildcards: int(80 * GB)
        if wildcards.db == "standard"
        else int(12 * GB),
    container:
        CONTAINERS["kraken"]
    benchmark:
        repeat(BENCH / "classify/kraken/{db}/{tech}.tsv", config["benchmark_repeats"])
    params:
        db=lambda wildcards, input: Path(input.db).parent
        if wildcards.db == "mycobacterium"
        else input.db,
        opts="--minimum-hit-groups 3 --report-minimizer-data",
        tech_opts=lambda wildcards: "--paired" if wildcards.tech == "illumina" else "",
    shell:
        """
        kraken2 {params.opts} {params.tech_opts} --threads {threads} --db {params.db} \
            --report {output.report} --output {output.out} {input.reads} 2> {log}
        """
