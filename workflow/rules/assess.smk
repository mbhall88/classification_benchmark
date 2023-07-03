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
