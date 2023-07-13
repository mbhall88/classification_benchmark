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
        metadata=RESULTS / "assess/read2taxonomy.tsv",
    log:
        LOGS / "make_read_truth.log",
    resources:
        runtime="2h",
    conda:
        ENVS / "make_read_truth.yaml"
    script:
        SCRIPTS / "make_read_truth.py"


rule sra_human_scrubber_classifications:
    input:
        reads=rules.sra_human_scrubber.output.reads,
        removed=rules.sra_human_scrubber.output.removed,
        truth=rules.make_read_truth.output.metadata,
    output:
        classification=RESULTS / "dehumanise/classifications.sra.tsv",
    log:
        LOGS / "sra_human_scrubber_classifications.log",
    resources:
        runtime="20m",
    container:
        CONTAINERS["pysam"]
    script:
        SCRIPTS / "sra_human_scrubber_classifications.py"


rule minimap2_dehumanise_classification:
    input:
        truth=rules.make_read_truth.output.metadata,
        alignment=rules.minimap2_human_scrubber.output.aln,
    output:
        classification=RESULTS / "dehumanise/classifications.minimap2.tsv",
    log:
        LOGS / "minimap2_dehumanise_classification.log",
    resources:
        runtime="20m",
    conda:
        ENVS / "minimap2_dehumanise_classification.yaml"
    script:
        SCRIPTS / "minimap2_dehumanise_classification.py"


rule kraken_dehumanise_classification:
    input:
        truth=rules.make_read_truth.output.metadata,
        classification=rules.kraken_human_classify.output.out
    output:
        classification=RESULTS / "dehumanise/classifications.kraken.k{k}l{l}.tsv",
    log:
        LOGS / "kraken_dehumanise_classification/k{k}l{l}.log",
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
        classification=RESULTS / "dehumanise/classifications.miniwinnow.tsv",
    log:
        LOGS / "miniwinnow_dehumanise_classification.log",
    resources:
        runtime="10m",
    conda:
        ENVS / "minimap2_dehumanise_classification.yaml"
    script:
        SCRIPTS / "miniwinnow_dehumanise_classification.py"
