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
    if wildcards.db == "clockwork":
        return RESULTS / "db/db.fa.gz"
    elif wildcards.db == "mtbc":
        return RESULTS / "db/GTDB_genus_Mycobacterium/MTB.fna.gz"
    elif wildcards.db == "mycobacterium":
        return RESULTS / "db/GTDB_genus_Mycobacterium/Mycobacterium.rep.fna.gz"
    else:
        raise ValueError(f"Don't recognise db {wildcards.db}")


PRESETS = {"ont": "map-ont", "illumina": "sr"}


rule minimap2_classify_clockwork:
    input:
        reads=infer_classify_reads,
        db=rules.faidx_db.output.fasta,
    output:
        aln=RESULTS / "classify/minimap2/aln.{db}.{tech}.paf",
    log:
        LOGS / "minimap2_classify_clockwork/{db}/{tech}.log",
    resources:
        runtime="2h",
        mem_mb=int(16 * GB),
    container:
        CONTAINERS["minimap2"]
    params:
        opts="--secondary=no -c",
        preset=lambda wildcards: PRESETS[wildcards.tech],
    shell:
        "minimap2 {params.opts} -x {params.preset} -o {output.aln} {input.db} {input.reads} 2> {log}"


# todo: kraken - default db
# todo: kraken - default small db
# todo: kraken - Myco db
