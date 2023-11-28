"""This snakemake file contains the rules for the analysis of the pseudo-real data set."""

rule download_and_subset_pod5:
    output:
        pod5_dir=directory(RESULTS / "pod5/{sample}/pod5/"),
    log:
        LOGS / "download_and_subset_pod5/{sample}.log",
    threads: 4
    resources:
        mem_mb=int(10 * GB),
        runtime="1w"
    params:
        url=lambda wildcards: config["real_data"]["human"][wildcards.sample]["ont"],
        n_reads=50_000,
    shadow:
        "shallow"
    conda:
        ENVS / "download_and_subset_pod5.yaml"
    script:
        SCRIPTS / "download_and_subset_pod5.sh"