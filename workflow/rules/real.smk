"""This snakemake file contains the rules for the analysis of the pseudo-real data set."""


def infer_accession(wildcards, tech):
    return config["real_data"][wildcards.sample][tech]


rule download_real_data_ont:
    output:
        DATA / "ont/{sample}.fq.gz",
    log:
        LOGS / "download_real_data_ont/{sample}.log",
    wildcard_constraints:
        sample="|".join([s for s in config["real_data"].keys() if s != "human"]),
    resources:
        mem_mb=int(GB * 32),
        runtime="6h",
    conda:
        ENVS / "download_real_data.yaml"
    params:
        accession=lambda wildcards: infer_accession(wildcards, "ont"),
        max_bases=config["max_real_bases"],
        min_len=500,
        seed=23,
    shadow:
        "shallow"
    shell:
        """
        exec 2> {log}
        tdir=$(mktemp -d)
        fastq-dl -a {params.accession} -o "$tdir"
        tmpfq="$tdir/"tmp.fq
        seqkit seq -m {params.min_len} -o $tmpfq "$tdir/"{params.accession}.fastq.gz
        rasusa -i $tmpfq -o {output[0]} -b {params.max_bases} -s {params.seed}
        rm -rf "$tdir"
        """


rule download_real_data_illumina:
    output:
        [DATA / f"illumina/{{sample}}_{i}.fq.gz" for i in [1, 2]],
    log:
        LOGS / "download_real_data_illumina/{sample}.log",
    wildcard_constraints:
        sample="|".join([s for s in config["real_data"].keys() if s != "human"]),
    resources:
        mem_mb=int(GB * 32),
        runtime="6h",
    conda:
        ENVS / "download_real_data.yaml"
    params:
        accession=lambda wildcards: infer_accession(wildcards, "illumina"),
        max_bases=config["max_real_bases"],
        seed=23,
    shadow:
        "shallow"
    shell:
        """
        exec 2> {log}
        tdir=$(mktemp -d)
        fastq-dl -a {params.accession} -o "$tdir"
        tmp1="$tdir/"{params.accession}_1.fastq.gz
        tmp2="$tdir/"{params.accession}_2.fastq.gz
        # td=$(mktemp -d)
        # seqkit pair --id-regexp '^(\S+(\s\d+)?)[\/\.][12]' -1 $tmp1 -2 $tmp2 -O $td
        # tmp1="$td/"$(basename $tmp1)
        # tmp2="$td/"$(basename $tmp2)
        rasusa -i $tmp1 -i $tmp2 -o {output[0]} -o {output[1]} -b {params.max_bases} -s {params.seed}
        # rm -rf "$td"
        rm -rf "$tdir"
        """


rule create_ont_truthsheet_and_combine:
    input:
        mtb_reads=DATA / "ont/mtb.fq.gz",
        zymo_reads=DATA / "ont/zymo.fq.gz",
        human_reads=DATA / "ont/human.fq.gz",
    output:
        truth=RESULTS / "real/truth/read2truth.ont.tsv",
        all_reads=RESULTS / "real/reads/all.ont.fq.gz",
        non_human_reads=RESULTS / "real/reads/nonhuman.ont.fq.gz",
    log:
        LOGS / "create_ont_truthsheet.log",
    resources:
        mem_mb=int(GB * 4),
        runtime="1h",
    container:
        CONTAINERS["seqkit"]
    shell:
        """
        exec 2> {log}
        seqkit seq -ni {input.mtb_reads} | sed 's/$/\tmtb/' > {output.truth}
        seqkit seq -ni {input.zymo_reads} | sed 's/$/\tzymo/' >> {output.truth}
        seqkit seq -ni {input.human_reads} | sed 's/$/\thuman/' >> {output.truth}
        seqkit seq -o {output.all_reads} {input.mtb_reads} {input.zymo_reads} {input.human_reads}
        seqkit seq -o {output.non_human_reads} {input.mtb_reads} {input.zymo_reads}
        """


rule create_illumina_truthsheet_and_combine:
    input:
        mtb_reads1=DATA / "illumina/mtb_1.fq.gz",
        mtb_reads2=DATA / "illumina/mtb_2.fq.gz",
        zymo_reads1=DATA / "illumina/zymo_1.fq.gz",
        zymo_reads2=DATA / "illumina/zymo_2.fq.gz",
        human_reads1=DATA / "illumina/human_1.fq.gz",
        human_reads2=DATA / "illumina/human_2.fq.gz",
    output:
        truth=RESULTS / "real/truth/read2truth.illumina.tsv",
        all_reads1=RESULTS / "real/reads/all_1.illumina.fq.gz",
        all_reads2=RESULTS / "real/reads/all_2.illumina.fq.gz",
        non_human_reads1=RESULTS / "real/reads/nonhuman_1.illumina.fq.gz",
        non_human_reads2=RESULTS / "real/reads/nonhuman_2.illumina.fq.gz",
    log:
        LOGS / "create_illumina_truthsheet.log",
    resources:
        mem_mb=int(GB * 4),
        runtime="1h",
    container:
        CONTAINERS["seqkit"]
    shell:
        """
        exec 2> {log}
        # rename reads
        mtb1="$(mktemp -u)_1.fq.gz"
        mtb2="$(mktemp -u)_2.fq.gz"
        seqkit replace -p .+ -r "mtb_{{nr}}/1" -o $mtb1 {input.mtb_reads1}
        seqkit replace -p .+ -r "mtb_{{nr}}/2" -o $mtb2 {input.mtb_reads2}
        zymo1="$(mktemp -u)_1.fq.gz"
        zymo2="$(mktemp -u)_2.fq.gz"
        seqkit replace -p .+ -r "zymo_{{nr}}/1" -o $zymo1 {input.zymo_reads1}
        seqkit replace -p .+ -r "zymo_{{nr}}/2" -o $zymo2 {input.zymo_reads2}
        human1="$(mktemp -u)_1.fq.gz"
        human2="$(mktemp -u)_2.fq.gz"
        seqkit replace -p .+ -r "human_{{nr}}/1" -o $human1 {input.human_reads1}
        seqkit replace -p .+ -r "human_{{nr}}/2" -o $human2 {input.human_reads2}
        # create truth sheet
        seqkit seq -ni $mtb1 $mtb2 | sed 's/$/\tmtb/' > {output.truth}
        seqkit seq -ni $zymo1 $zymo2 | sed 's/$/\tzymo/' >> {output.truth}
        seqkit seq -ni $human1 $human2 | sed 's/$/\thuman/' >> {output.truth}
        # combine reads
        seqkit seq -o {output.all_reads1} $mtb1 $zymo1 $human1
        seqkit seq -o {output.all_reads2} $mtb2 $zymo2 $human2
        seqkit seq -o {output.non_human_reads1} $mtb1 $zymo1
        seqkit seq -o {output.non_human_reads2} $mtb2 $zymo2
        rm $mtb1 $mtb2 $zymo1 $zymo2 $human1 $human2
        """


def infer_dehumanise_real_reads(wildcards):
    if wildcards.tech == "illumina":
        return [RESULTS / f"real/reads/all_{i}.{wildcards.tech}.fq.gz" for i in [1, 2]]
    else:
        return RESULTS / f"real/reads/all.{wildcards.tech}.fq.gz"


rule kraken_dehumanise_real:
    input:
        reads=infer_dehumanise_real_reads,
        db=infer_kraken_db,
    output:
        report=RESULTS / "real/dehumanise/kraken/{lib}/{tech}.k2report",
        out=RESULTS / "real/dehumanise/kraken/{lib}/{tech}.k2",
    log:
        LOGS / "kraken_dehumanise_real/{lib}/{tech}.log",
    threads: 4
    resources:
        runtime=f"{20 * REPEAT}m",
        mem_mb=int(10 * GB),
    container:
        CONTAINERS["kraken"]
    benchmark:
        repeat(BENCH / "real/dehumanise/kraken/{lib}/{tech}.tsv", REPEAT)
    params:
        opts="--minimum-hit-groups 3 --report-minimizer-data",
        tech_opts=lambda wildcards: "--paired" if wildcards.tech == "illumina" else "",
    shell:
        """
        kraken2 {params.opts} {params.tech_opts} --threads {threads} --db {input.db} \
            --report {output.report} --output {output.out} {input.reads} 2> {log}
        """


use rule sra_human_scrubber as sra_human_scrubber_real with:
    input:
        reads=RESULTS / "real/reads/all.ont.fq.gz",
    output:
        reads=RESULTS / "real/dehumanise/sra/metagenome.scrubbed.ont.fq.gz",
        removed=RESULTS / "real/dehumanise/sra/metagenome.removed.ont.fq.gz",
    log:
        LOGS / "sra_human_scrubber_real.log",
    benchmark:
        repeat(BENCH / "real/dehumanise/sra/ont.tsv", REPEAT)


use rule sra_human_scrubber_illumina as sra_human_scrubber_illumina_real with:
    input:
        r1=RESULTS / f"real/reads/all_1.illumina.fq.gz",
        r2=RESULTS / f"real/reads/all_2.illumina.fq.gz",
    output:
        reads1=RESULTS / "real/dehumanise/sra/metagenome_R1.scrubbed.illumina.fq.gz",
        reads2=RESULTS / "real/dehumanise/sra/metagenome_R2.scrubbed.illumina.fq.gz",
        removed1=RESULTS / "real/dehumanise/sra/metagenome_R1.removed.illumina.fq.gz",
        removed2=RESULTS / "real/dehumanise/sra/metagenome_R2.removed.illumina.fq.gz",
    log:
        LOGS / "sra_human_scrubber_illumina_real.log",
    benchmark:
        repeat(BENCH / "real/dehumanise/sra/illumina.tsv", REPEAT)


use rule minimap2_human_scrubber as minimap2_human_scrubber_real with:
    input:
        reads=rules.sra_human_scrubber_real.input.reads,
        index=RESULTS / "db/chm13.map-ont.mm2",
    output:
        aln=RESULTS / "real/dehumanise/minimap2/aln.ont.paf",
    log:
        LOGS / "minimap2_human_scrubber_real.log",
    benchmark:
        repeat(BENCH / "real/dehumanise/minimap2/ont.tsv", REPEAT)


use rule minimap2_human_scrubber_illumina as minimap2_human_scrubber_illumina_real with:
    input:
        r1=rules.sra_human_scrubber_illumina_real.input.r1,
        r2=rules.sra_human_scrubber_illumina_real.input.r2,
        index=RESULTS / "db/chm13.sr.mm2",
    output:
        aln=RESULTS / "real/dehumanise/minimap2/aln.illumina.paf",
    log:
        LOGS / "minimap2_human_scrubber_illumina_real.log",
    benchmark:
        repeat(BENCH / "real/dehumanise/minimap2/illumina.tsv", REPEAT)


use rule miniwinnow_human_scrubber as miniwinnow_human_scrubber_real with:
    input:
        reads=rules.sra_human_scrubber_real.input.reads,
        ref=rules.download_chm13.output.fasta,
        aln=rules.minimap2_human_scrubber_real.output.aln,
    output:
        aln=RESULTS / "real/dehumanise/miniwinnow/aln.ont.paf",
    log:
        LOGS / "miniwinnow_human_scrubber_real.log",
    benchmark:
        repeat(BENCH / "real/dehumanise/miniwinnow/ont.tsv", REPEAT)


use rule hostile_human_scrubber as hostile_human_scrubber_real with:
    input:
        reads=rules.sra_human_scrubber_real.input.reads,
    output:
        reads=RESULTS / "real/dehumanise/hostile/all.ont.clean.fastq.gz",
    log:
        LOGS / "hostile_human_scrubber_real.log",
    benchmark:
        repeat(BENCH / "real/dehumanise/hostile/ont.tsv", REPEAT)


use rule hostile_human_scrubber_illumina as hostile_human_scrubber_illumina_real with:
    input:
        r1=rules.sra_human_scrubber_illumina_real.input.r1,
        r2=rules.sra_human_scrubber_illumina_real.input.r2,
    output:
        reads1=RESULTS / "real/dehumanise/hostile/all_1.illumina.clean_1.fastq.gz",
        reads2=RESULTS / "real/dehumanise/hostile/all_2.illumina.clean_2.fastq.gz",
    log:
        LOGS / "hostile_human_scrubber_illumina_real.log",
    benchmark:
        repeat(BENCH / "real/dehumanise/hostile/illumina.tsv", REPEAT)


rule sra_human_scrubber_classifications_real:
    input:
        reads=rules.sra_human_scrubber_real.output.reads,
        removed=rules.sra_human_scrubber_real.output.removed,
        truth=rules.create_ont_truthsheet_and_combine.output.truth,
    output:
        classification=RESULTS / "real/dehumanise/classifications.sra.ont.tsv",
    log:
        LOGS / "sra_human_scrubber_classifications_real.log",
    resources:
        runtime="15m",
        mem_mb=int(8 * GB),
    container:
        CONTAINERS["pysam"]
    script:
        SCRIPTS / "sra_human_scrubber_classifications_real.py"


classification_base_mem = int(32 * GB)


rule sra_human_scrubber_classifications_illumina_real:
    input:
        reads=rules.sra_human_scrubber_illumina_real.output.reads1,
        reads2=rules.sra_human_scrubber_illumina_real.output.reads2,
        removed=rules.sra_human_scrubber_illumina_real.output.removed1,
        removed2=rules.sra_human_scrubber_illumina_real.output.removed2,
        truth=rules.create_illumina_truthsheet_and_combine.output.truth,
    output:
        classification=RESULTS / "real/dehumanise/classifications.sra.illumina.tsv",
    log:
        LOGS / "sra_human_scrubber_classifications_illumina_real.log",
    resources:
        runtime="15m",
        mem_mb=lambda wildcards, attempt: classification_base_mem * attempt,
    container:
        CONTAINERS["pysam"]
    params:
        ignore_unmapped=True,
    script:
        SCRIPTS / "sra_human_scrubber_classifications_real.py"


rule minimap2_dehumanise_classification_real:
    input:
        truth=RESULTS / "real/truth/read2truth.{tech}.tsv",
        alignment=RESULTS / "real/dehumanise/minimap2/aln.{tech}.paf",
    output:
        classification=RESULTS / "real/dehumanise/classifications.minimap2.{tech}.tsv",
    log:
        LOGS / "minimap2_dehumanise_classification_real/{tech}.log",
    resources:
        runtime="15m",
        mem_mb=lambda wildcards, attempt: classification_base_mem * attempt,
    conda:
        ENVS / "minimap2_dehumanise_classification.yaml"
    script:
        SCRIPTS / "minimap2_dehumanise_classification_real.py"


rule miniwinnow_dehumanise_classification_real:
    input:
        truth=rules.create_ont_truthsheet_and_combine.output.truth,
        alignment=rules.miniwinnow_human_scrubber_real.output.aln,
    output:
        classification=RESULTS / "real/dehumanise/classifications.miniwinnow.ont.tsv",
    log:
        LOGS / "miniwinnow_dehumanise_classification_real.log",
    resources:
        runtime="10m",
        mem_mb=lambda wildcards, attempt: classification_base_mem * attempt,
    conda:
        ENVS / "minimap2_dehumanise_classification.yaml"
    script:
        SCRIPTS / "miniwinnow_dehumanise_classification_real.py"


rule kraken_dehumanise_classification_real:
    input:
        truth=RESULTS / "real/truth/read2truth.{tech}.tsv",
        classification=rules.kraken_dehumanise_real.output.out,
    output:
        classification=RESULTS
        / "real/dehumanise/classifications.kraken.{lib}.{tech}.tsv",
    log:
        LOGS / "kraken_dehumanise_classification_real/{lib}/{tech}.log",
    resources:
        runtime="15m",
        mem_mb=lambda wildcards, attempt: classification_base_mem * attempt,
    container:
        CONTAINERS["python"]
    script:
        SCRIPTS / "kraken_dehumanise_classification_real.py"


rule hostile_human_scrubber_classifications_real:
    input:
        reads=rules.hostile_human_scrubber_real.output.reads,
        truth=rules.create_ont_truthsheet_and_combine.output.truth,
    output:
        classification=RESULTS / "real/dehumanise/classifications.hostile.ont.tsv",
    log:
        LOGS / "hostile_human_scrubber_classifications_real.log",
    resources:
        runtime="15m",
        mem_mb=lambda wildcards, attempt: classification_base_mem * attempt,
    container:
        CONTAINERS["pysam"]
    script:
        SCRIPTS / "hostile_human_scrubber_classification_real.py"


rule hostile_human_scrubber_classifications_illumina_real:
    input:
        reads=rules.hostile_human_scrubber_illumina_real.output.reads1,
        reads2=rules.hostile_human_scrubber_illumina_real.output.reads2,
        truth=rules.create_illumina_truthsheet_and_combine.output.truth,
    output:
        classification=RESULTS / "real/dehumanise/classifications.hostile.illumina.tsv",
    log:
        LOGS / "hostile_human_scrubber_classifications_illumina_real.log",
    resources:
        runtime="15m",
        mem_mb=lambda wildcards, attempt: classification_base_mem * attempt,
    container:
        CONTAINERS["pysam"]
    script:
        SCRIPTS / "hostile_human_scrubber_classification_real.py"


tools = [
    "sra",
    "hostile",
    "minimap2",
    "miniwinnow",
    "kraken.default",
    "kraken.HPRC",
]


def infer_real_benchmarks(wildcards):
    benchmarks = []
    for tool in tools:
        if tool.startswith("kraken"):
            tool = tool.replace(".", "/")
        if wildcards.tech == "illumina" and tool == "miniwinnow":
            continue
        benchmarks.append(BENCH / f"real/dehumanise/{tool}/{wildcards.tech}.tsv")

    return benchmarks


def infer_real_classifications(wildcards):
    clfs = []
    for tool in tools:
        if wildcards.tech == "illumina" and tool == "miniwinnow":
            continue
        clfs.append(
            RESULTS / f"real/dehumanise/classifications.{tool}.{wildcards.tech}.tsv"
        )

    return clfs


use rule dehumanise_summary_statistics as dehumanise_summary_statistics_real with:
    input:
        classifications=infer_real_classifications,
        benchmarks=infer_real_benchmarks,
    output:
        summary=RESULTS / "real/dehumanise/summary.{tech}.csv",
    log:
        LOGS / "dehumanise_summary_statistics_real/{tech}.log",


# ===============================================
# ================ MTB CLASSIFICATION ===========


use rule minimap2_classify as minimap2_classify_real with:
    input:
        reads=infer_classify_real_reads,
        db=infer_minimap2_db_real,
    output:
        aln=RESULTS / "real/classify/minimap2/aln.{db}.{tech}.paf",
    log:
        LOGS / "minimap2_classify_real/{db}/{tech}.log",
    benchmark:
        repeat(BENCH / "real/classify/minimap2/{db}/{tech}.tsv", REPEAT)


use rule kraken_classify as kraken_classify_real with:
    input:
        reads=infer_classify_real_reads,
        db=lambda wildcards: CLASSIFY_DBS[wildcards.db],
    output:
        report=RESULTS / "real/classify/kraken/{db}.{tech}.k2report",
        out=RESULTS / "real/classify/kraken/{db}.{tech}.k2",
    log:
        LOGS / "kraken_classify_real/{db}/{tech}.log",
    benchmark:
        repeat(BENCH / "real/classify/kraken/{db}/{tech}.tsv", REPEAT)


def infer_truthsheet_real(wildcards):
    if wildcards.tech == "illumina":
        return RESULTS / "real/truth/read2truth.illumina.tsv"
    else:
        return RESULTS / "real/truth/read2truth.ont.tsv"


rule kraken_mycobacterium_classification_real:
    input:
        truth=infer_truthsheet_real,
        classification=rules.kraken_classify_real.output.out,
        names=rules.download_kraken_taxonomy.output.names,
        nodes=rules.download_kraken_taxonomy.output.nodes,
    output:
        classification=RESULTS / "real/classify/classifications.kraken.{db}.{tech}.tsv",
    log:
        LOGS / "kraken_mycobacterium_classification_real/{db}/{tech}.log",
    resources:
        runtime="15m",
        mem_mb=int(16 * GB),
    conda:
        ENVS / "accession2taxonomy.yaml"
    script:
        SCRIPTS / "kraken_mycobacterium_classification_real.py"


rule minimap2_mycobacterium_classification_real:
    input:
        truth=infer_truthsheet_real,
        aln=rules.minimap2_classify_real.output.aln,
        names=rules.download_kraken_taxonomy.output.names,
        nodes=rules.download_kraken_taxonomy.output.nodes,
        acc_truth=rules.minimap2_acc2taxid.output.truth,
    output:
        classification=RESULTS / "real/classify/classifications.minimap2.{db}.{tech}.tsv",
    log:
        LOGS / "minimap2_mycobacterium_classification_real/{db}/{tech}.log",
    resources:
        runtime="1h",
        mem_mb=int(16 * GB),
    conda:
        ENVS / "minimap2_mycobacterium_classification.yaml"
    script:
        SCRIPTS / "minimap2_mycobacterium_classification_real.py"


def infer_classify_benchmarks_real(wildcards):
    benchmarks = []
    tech = wildcards.tech
    for db in ["standard", "standard-8", "mycobacterium"]:
        benchmarks.append(BENCH / f"real/classify/kraken/{db}/{tech}.tsv")
    for db in ["clockwork", "mtbc", "mycobacterium"]:
        benchmarks.append(BENCH / f"real/classify/minimap2/{db}/{tech}.tsv")

    return benchmarks


def infer_classify_classifications_real(wildcards):
    clfs = []
    tech = wildcards.tech
    for db in ["standard", "standard-8", "mycobacterium"]:
        clfs.append(RESULTS / f"real/classify/classifications.kraken.{db}.{tech}.tsv")
    for db in ["clockwork", "mtbc", "mycobacterium"]:
        clfs.append(RESULTS / f"real/classify/classifications.minimap2.{db}.{tech}.tsv")

    return clfs


rule classify_summary_statistics_real:
    input:
        classifications=infer_classify_classifications_real,
        benchmarks=infer_classify_benchmarks_real,
    output:
        mtb_summary=RESULTS / "real/classify/summary.mtb.{tech}.csv",
    log:
        LOGS / "classify_summary_statistics_real/{tech}.log",
    resources:
        runtime="5m",
    conda:
        ENVS / "datasci.yaml"
    script:
        SCRIPTS / "classify_summary_statistics_real.py"