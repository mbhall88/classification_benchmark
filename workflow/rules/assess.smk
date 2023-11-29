rule accession2taxonomy:
    input:
        names=rules.download_kraken_taxonomy.output.names,
        nodes=rules.download_kraken_taxonomy.output.nodes,
        refs=[
            infer_simulate_input(read_type)
            for read_type in config["simulate"]["proportions"]
            if read_type != "Unmapped"
        ],
    output:
        metadata=RESULTS / "assess/acc2tax.tsv",
    log:
        LOGS / "accession2taxonomy.log",
    resources:
        runtime="2h",
        mem_mb=int(4 * GB),
    conda:
        ENVS / "accession2taxonomy.yaml"
    script:
        SCRIPTS / "accession2taxid.py"


rule make_read_truth:
    input:
        acc2tax=rules.accession2taxonomy.output.metadata,
        nodes=rules.download_kraken_taxonomy.output.nodes,
        reads=rules.combine_simulated_reads.output.reads,
    output:
        metadata=RESULTS / "assess/read2taxonomy.ont.tsv",
    log:
        LOGS / "make_read_truth.log",
    resources:
        runtime="2h",
    conda:
        ENVS / "make_read_truth.yaml"
    script:
        SCRIPTS / "make_read_truth.py"


rule make_illumina_read_truth:
    input:
        acc2tax=rules.accession2taxonomy.output.metadata,
        nodes=rules.download_kraken_taxonomy.output.nodes,
        r1=rules.combine_illumina_simulated_reads.output.r1,
        r2=rules.combine_illumina_simulated_reads.output.r2,
    output:
        metadata=RESULTS / "assess/read2taxonomy.illumina.tsv",
    log:
        LOGS / "make_illumina_read_truth.log",
    resources:
        runtime="2h",
    conda:
        ENVS / "make_read_truth.yaml"
    script:
        SCRIPTS / "make_illumina_read_truth.py"


rule sra_human_scrubber_classifications:
    input:
        reads=rules.sra_human_scrubber.output.reads,
        removed=rules.sra_human_scrubber.output.removed,
        truth=rules.make_read_truth.output.metadata,
    output:
        classification=RESULTS / "dehumanise/classifications.sra.ont.tsv",
    log:
        LOGS / "sra_human_scrubber_classifications.log",
    resources:
        runtime="15m",
    container:
        CONTAINERS["pysam"]
    params:
        ignore_unmapped=True,
    script:
        SCRIPTS / "sra_human_scrubber_classifications.py"


rule sra_human_scrubber_classifications_illumina:
    input:
        reads=rules.sra_human_scrubber_illumina.output.reads1,
        reads2=rules.sra_human_scrubber_illumina.output.reads2,
        removed=rules.sra_human_scrubber_illumina.output.removed1,
        removed2=rules.sra_human_scrubber_illumina.output.removed2,
        truth=rules.make_illumina_read_truth.output.metadata,
    output:
        classification=RESULTS / "dehumanise/classifications.sra.illumina.tsv",
    log:
        LOGS / "sra_human_scrubber_classifications_illumina.log",
    resources:
        runtime="15m",
        mem_mb=int(4 * GB),
    container:
        CONTAINERS["pysam"]
    params:
        ignore_unmapped=True,
    script:
        SCRIPTS / "sra_human_scrubber_classifications.py"


rule minimap2_dehumanise_classification:
    input:
        truth=RESULTS / "assess/read2taxonomy.{tech}.tsv",
        alignment=RESULTS / "dehumanise/minimap2/metagenome.aln.{tech}.paf",
    output:
        classification=RESULTS / "dehumanise/classifications.minimap2.{tech}.tsv",
    log:
        LOGS / "minimap2_dehumanise_classification/{tech}.log",
    resources:
        runtime="15m",
        mem_mb=int(4 * GB),
    conda:
        ENVS / "minimap2_dehumanise_classification.yaml"
    params:
        ignore_unmapped=True,
    script:
        SCRIPTS / "minimap2_dehumanise_classification.py"


rule kraken_dehumanise_classification:
    input:
        truth=RESULTS / "assess/read2taxonomy.{tech}.tsv",
        classification=RESULTS / "dehumanise/kraken/classify/{lib}/metagenome.{tech}.k2",
    output:
        classification=RESULTS / "dehumanise/classifications.kraken.{lib}.{tech}.tsv",
    log:
        LOGS / "kraken_dehumanise_classification/{lib}/{tech}.log",
    resources:
        runtime="15m",
        mem_mb=int(4 * GB),
    container:
        CONTAINERS["python"]
    params:
        ignore_unmapped=True,
    script:
        SCRIPTS / "kraken_dehumanise_classification.py"


rule miniwinnow_dehumanise_classification:
    input:
        truth=rules.make_read_truth.output.metadata,
        alignment=rules.miniwinnow_human_scrubber.output.aln,
    output:
        classification=RESULTS / "dehumanise/classifications.miniwinnow.ont.tsv",
    log:
        LOGS / "miniwinnow_dehumanise_classification.log",
    resources:
        runtime="10m",
    conda:
        ENVS / "minimap2_dehumanise_classification.yaml"
    params:
        ignore_unmapped=True,
    script:
        SCRIPTS / "miniwinnow_dehumanise_classification.py"


rule hostile_human_scrubber_classifications:
    input:
        reads=rules.hostile_human_scrubber.output.reads,
        truth=rules.make_read_truth.output.metadata,
    output:
        classification=RESULTS / "dehumanise/classifications.hostile.ont.tsv",
    log:
        LOGS / "hostile_human_scrubber_classifications.log",
    resources:
        runtime="15m",
        mem_mb=int(4 * GB),
    container:
        CONTAINERS["pysam"]
    params:
        ignore_unmapped=True,
    script:
        SCRIPTS / "hostile_human_scrubber_classification.py"


rule hostile_human_scrubber_classifications_illumina:
    input:
        reads=rules.hostile_human_scrubber_illumina.output.reads1,
        reads2=rules.hostile_human_scrubber_illumina.output.reads2,
        truth=rules.make_illumina_read_truth.output.metadata,
    output:
        classification=RESULTS / "dehumanise/classifications.hostile.illumina.tsv",
    log:
        LOGS / "hostile_human_scrubber_classifications_illumina.log",
    resources:
        runtime="15m",
        mem_mb=int(4 * GB),
    container:
        CONTAINERS["pysam"]
    params:
        ignore_unmapped=True,
    script:
        SCRIPTS / "hostile_human_scrubber_classification.py"


tools = [
    "sra",
    "hostile",
    "minimap2",
    "miniwinnow",
    "kraken.default",
    "kraken.HPRC",
]


def infer_benchmarks(wildcards):
    benchmarks = []
    for tool in tools:
        if tool.startswith("kraken"):
            tool = tool.replace(".", "/")
        if wildcards.tech == "illumina" and tool == "miniwinnow":
            continue
        benchmarks.append(BENCH / f"dehumanise/{tool}/{wildcards.tech}.tsv")

    return benchmarks


def infer_classifications(wildcards):
    clfs = []
    for tool in tools:
        if wildcards.tech == "illumina" and tool == "miniwinnow":
            continue
        clfs.append(RESULTS / f"dehumanise/classifications.{tool}.{wildcards.tech}.tsv")

    return clfs


rule dehumanise_summary_statistics:
    input:
        classifications=infer_classifications,
        benchmarks=infer_benchmarks,
    output:
        summary=RESULTS / "dehumanise/summary.{tech}.csv",
    log:
        LOGS / "dehumanise_summary_statistics/{tech}.log",
    resources:
        runtime="5m",
    conda:
        ENVS / "datasci.yaml"
    script:
        SCRIPTS / "dehumanise_summary_statistics.py"


rule kraken_mycobacterium_classification:
    input:
        truth=RESULTS / "assess/read2taxonomy.{tech}.tsv",
        classification=rules.kraken_classify.output.out,
        names=rules.download_kraken_taxonomy.output.names,
        nodes=rules.download_kraken_taxonomy.output.nodes,
    output:
        classification=RESULTS / "classify/classifications.kraken.{db}.{tech}.tsv",
    log:
        LOGS / "kraken_mycobacterium_classification/{db}/{tech}.log",
    resources:
        runtime="15m",
        mem_mb=int(2 * GB),
    conda:
        ENVS / "accession2taxonomy.yaml"
    params:
        ignore_unmapped=True,
    script:
        SCRIPTS / "kraken_mycobacterium_classification.py"


rule minimap2_acc2taxid:
    input:
        metadata=rules.combine_references.output.metadata,
        names=rules.download_kraken_taxonomy.output.names,
        nodes=rules.download_kraken_taxonomy.output.nodes,
        fasta=rules.create_minimap2_mycobacterium_db.output.db,
        acc2taxid=CONFIG / "acc2taxid.csv",
    output:
        truth=RESULTS / "assess/minimap2.db.acc2tax.tsv",
    log:
        LOGS / "minimap2_acc2taxid.log",
    resources:
        runtime="1w",
    conda:
        ENVS / "minimap2_acc2taxid.yaml"
    script:
        SCRIPTS / "minimap2_acc2taxid.py"


rule minimap2_mycobacterium_classification:
    input:
        truth=RESULTS / "assess/read2taxonomy.{tech}.tsv",
        aln=rules.minimap2_classify.output.aln,
        names=rules.download_kraken_taxonomy.output.names,
        nodes=rules.download_kraken_taxonomy.output.nodes,
        acc_truth=rules.minimap2_acc2taxid.output.truth,
    output:
        classification=RESULTS / "classify/classifications.minimap2.{db}.{tech}.tsv",
    log:
        LOGS / "minimap2_mycobacterium_classification/{db}/{tech}.log",
    resources:
        runtime="15m",
        mem_mb=int(2 * GB),
    conda:
        ENVS / "minimap2_mycobacterium_classification.yaml"
    params:
        ignore_unmapped=False,
    script:
        SCRIPTS / "minimap2_mycobacterium_classification.py"


def infer_classify_benchmarks(wildcards):
    benchmarks = []
    tech = wildcards.tech
    for db in ["standard", "standard-8", "mycobacterium"]:
        benchmarks.append(BENCH / f"classify/kraken/{db}/{tech}.tsv")
    for db in ["clockwork", "mtbc", "mycobacterium"]:
        benchmarks.append(BENCH / f"classify/minimap2/{db}/{tech}.tsv")

    return benchmarks


def infer_classify_classifications(wildcards):
    clfs = []
    tech = wildcards.tech
    for db in ["standard", "standard-8", "mycobacterium"]:
        clfs.append(RESULTS / f"classify/classifications.kraken.{db}.{tech}.tsv")
    for db in ["clockwork", "mtbc", "mycobacterium"]:
        clfs.append(RESULTS / f"classify/classifications.minimap2.{db}.{tech}.tsv")

    return clfs


rule classify_summary_statistics:
    input:
        classifications=infer_classify_classifications,
        benchmarks=infer_classify_benchmarks,
    output:
        genus_summary=RESULTS / "classify/summary.genus.{tech}.csv",
        species_summary=RESULTS / "classify/summary.species.{tech}.csv",
        mtb_summary=RESULTS / "classify/summary.mtb.{tech}.csv",
        mtbc_summary=RESULTS / "classify/summary.mtbc.{tech}.csv",
    log:
        LOGS / "classify_summary_statistics/{tech}.log",
    resources:
        runtime="5m",
    conda:
        ENVS / "datasci.yaml"
    script:
        SCRIPTS / "classify_summary_statistics.py"
