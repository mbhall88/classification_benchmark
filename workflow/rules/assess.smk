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
        runtime="20m",
    container:
        CONTAINERS["pysam"]
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
        runtime="20m",
    container:
        CONTAINERS["pysam"]
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
        runtime="20m",
    conda:
        ENVS / "minimap2_dehumanise_classification.yaml"
    script:
        SCRIPTS / "minimap2_dehumanise_classification.py"


rule kraken_dehumanise_classification:
    input:
        truth=RESULTS / "assess/read2taxonomy.{tech}.tsv",
        classification=RESULTS
        / "dehumanise/kraken/classify/k{k}l{l}/metagenome.{tech}.k2",
    output:
        classification=RESULTS / "dehumanise/classifications.kraken.k{k}l{l}.{tech}.tsv",
    log:
        LOGS / "kraken_dehumanise_classification/k{k}l{l}/{tech}.log",
    resources:
        runtime="20m",
    container:
        CONTAINERS["python"]
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
        runtime="10m",
    container:
        CONTAINERS["pysam"]
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
        runtime="10m",
    container:
        CONTAINERS["pysam"]
    script:
        SCRIPTS / "hostile_human_scrubber_classification.py"


tools = [
    "sra",
    "hostile",
    "minimap2",
    "miniwinnow",
    "kraken.k21l14",
    "kraken.k27l18",
    "kraken.k35l31",
]

benchmarks = []
for tool in tools:
    if tool.startswith("kraken"):
        tool = tool.replace(".", "/")
    benchmarks.append(BENCH / f"dehumanise/{tool}/ont.tsv")

tech_tool_combos = []
for tech in ["illumina", "ont"]:
    for tool in tools:
        if tech == "illumina" and tool == "miniwinnow":
            continue
        tech_tool_combos.append((tech, tool))


rule dehumanise_summary_statistics:
    input:
        classifications=[
            RESULTS / f"dehumanise/classifications.{tool}.{tech}.tsv"
            for tech, tool in tech_tool_combos
        ],
        benchmarks=benchmarks,
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
