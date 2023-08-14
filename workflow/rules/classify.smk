BEST_KLIB = "HPRC"


rule extract_dehumanised_ont_reads:
    input:
        reads=rules.combine_simulated_reads.output.reads,
        classification=RESULTS
        / f"dehumanise/kraken/classify/{BEST_KLIB}/metagenome.ont.k2",
    output:
        reads=RESULTS / "dehumanise/metagenome.dehumanised.ont.fq.gz",
    log:
        LOGS / "extract_dehumanised_ont_reads.log",
    resources:
        runtime="30m",
    container:
        CONTAINERS["seqtk"]
    shell:
        """
        (awk -F'\t' '$1=="U" {{print $2}}' {input.classification} \
            | seqtk subseq {input.reads} - \
            | gzip) > {output.reads} 2> {log}
        """


rule extract_dehumanised_illumina_reads:
    input:
        reads=rules.combine_illumina_simulated_reads.output.r1,
        reads2=rules.combine_illumina_simulated_reads.output.r2,
        classification=RESULTS
        / f"dehumanise/kraken/classify/{BEST_KLIB}/metagenome.illumina.k2",
    output:
        reads=RESULTS / "dehumanise/metagenome_R1.dehumanised.illumina.fq.gz",
        reads2=RESULTS / "dehumanise/metagenome_R2.dehumanised.illumina.fq.gz",
    log:
        LOGS / "extract_dehumanised_illumina_reads.log",
    resources:
        runtime="30m",
    container:
        CONTAINERS["seqtk"]
    shell:
        """
        (awk -F'\t' '$1=="U" {{print $2"/1"}}' {input.classification} \
            | seqtk subseq {input.reads} - \
            | gzip) > {output.reads} 2> {log}
        (awk -F'\t' '$1=="U" {{print $2"/2"}}' {input.classification} \
            | seqtk subseq {input.reads2} - \
            | gzip) > {output.reads2} 2>> {log}
        """


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


def infer_minimap2_db(wildcards):
    preset = PRESETS[wildcards.tech]
    if wildcards.db == "clockwork":
        return RESULTS / f"db/minimap2/db.{preset}.mmi"
    elif wildcards.db == "mtbc":
        return RESULTS / f"db/GTDB_genus_Mycobacterium/MTB.{preset}.mmi"
    elif wildcards.db == "mycobacterium":
        return RESULTS / f"db/GTDB_genus_Mycobacterium/Mycobacterium.rep.{preset}.mmi"
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
        mem_mb=int(16 * GB),
    threads: 4
    container:
        CONTAINERS["minimap2"]
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
    params:
        opts="--minimum-hit-groups 3 --report-minimizer-data",
        tech_opts=lambda wildcards: "--paired" if wildcards.tech == "illumina" else "",
        db=lambda wildcards, input: Path(input.db).parent
        if wildcards.db == "mycobacterium"
        else input.db,
    shell:
        """
        kraken2 {params.opts} {params.tech_opts} --threads {threads} --db {params.db} \
            --report {output.report} --output {output.out} {input.reads} 2> {log}
        """
