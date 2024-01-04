REPEAT = config["benchmark_repeats"]

use rule combine_simulated_reads as extract_dehumanised_ont_reads with:
    input:
        fastqs=[
            RESULTS / f"simulate/reads/{read_type}.ont.fq.gz"
            for read_type in config["simulate"]["proportions"]
            if read_type not in ("Unmapped", "Human")
        ],
        script=SCRIPTS / "filter_ambig.py",
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
        script=SCRIPTS / "filter_ambig.py",
    output:
        r1=RESULTS / "dehumanise/metagenome_R1.dehumanised.illumina.fq.gz",
        r2=RESULTS / "dehumanise/metagenome_R2.dehumanised.illumina.fq.gz",
    log:
        LOGS / "extract_dehumanised_illumina_reads.log"


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


rule minimap2_classify:
    input:
        reads=infer_classify_reads,
        db=infer_minimap2_db,
    output:
        aln=RESULTS / "classify/minimap2/aln.{db}.{tech}.paf",
    log:
        LOGS / "minimap2_classify/{db}/{tech}.log",
    resources:
        runtime=f"{2 * REPEAT}h",
        mem_mb=int(30 * GB),
    threads: 4
    container:
        CONTAINERS["minimap2"]
    benchmark:
        repeat(BENCH / "classify/minimap2/{db}/{tech}.tsv", REPEAT)
    params:
        opts="--secondary=no -c --paf-no-hit",
        preset=lambda wildcards: PRESETS[wildcards.tech],
    shell:
        """
        minimap2 {params.opts} -x {params.preset} -o {output.aln} -t {threads} \
            {input.db} {input.reads} 2> {log}
        """


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
        runtime=f"{20 * REPEAT}m",
        mem_mb=lambda wildcards: int(80 * GB)
        if wildcards.db == "standard"
        else int(12 * GB),
    container:
        CONTAINERS["kraken"]
    benchmark:
        repeat(BENCH / "classify/kraken/{db}/{tech}.tsv", REPEAT)
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
