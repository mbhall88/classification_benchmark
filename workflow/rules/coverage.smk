rule download_h37rv:
    output:
        RESULTS / "db/h37rv.fa.gz",
    log:
        LOGS / "download_h37rv.log",
    resources:
        mem_mb=500,
        runtime="2m",
    container:
        CONTAINERS["base"]
    shell:
        "wget -q -O {output} 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz' 2> {log}"


rule mtb_coverage_simulations_preclassification:
    input:
        reference=rules.download_h37rv.output,
        reads=infer_classify_reads,
    output:
        RESULTS / "coverage/preclassification/{tech}.tsv",
    log:
        LOGS / "mtb_coverage_simulations_preclassification/{tech}.log",
    resources:
        mem_mb=4 * GB,
        runtime="20m",
    threads: 4
    conda:
        ENVS / "coverage.yaml"
    shell:
        """
        (minimap2 -ax map-ont -t {threads} {input.reference} {input.reads} |
            samtools sort -@ {threads} | 
            samtools coverage -o {output} - ) 2> {log}
        """


use rule mtb_coverage_simulations_preclassification as mtb_coverage_real_preclassification with:
    input:
        reference=rules.download_h37rv.output,
        reads=infer_classify_real_reads,
    output:
        RESULTS / "real/coverage/preclassification/{tech}.tsv",
    log:
        LOGS / "mtb_coverage_real_preclassification/{tech}.log",


rule mtb_coverage_simulations_postclassification:
    input:
        reference=rules.download_h37rv.output,
        reads=infer_classify_reads,
        classifications=RESULTS / "classify/classifications.{tool}.{db}.{tech}.tsv",
    output:
        RESULTS / "coverage/postclassification/{tech}/{tool}/{db}.tsv",
    log:
        LOGS / "mtb_coverage_simulations_postclassification/{tech}/{tool}/{db}.log",
    resources:
        mem_mb=4 * GB,
        runtime="20m",
    threads: 4
    conda:
        ENVS / "coverage.yaml"
    shell:
        """
        (csvtk grep -tf mtb_classification -p TP {input.classifications} |
            csvtk cut -tf read_id |
            csvtk del-header |
            seqkit grep -f - {input.reads} |
            minimap2 -ax map-ont -t {threads} {input.reference} - |
            samtools sort -@ {threads} |
            samtools coverage -o {output} -) 2> {log}
        """


use rule mtb_coverage_simulations_postclassification as mtb_coverage_real_postclassification with:
    input:
        reference=rules.download_h37rv.output,
        reads=infer_classify_real_reads,
        classifications=RESULTS / "real/classify/classifications.{tool}.{db}.{tech}.tsv",
    output:
        RESULTS / "real/coverage/postclassification/{tech}/{tool}/{db}.tsv",
    log:
        LOGS / "mtb_coverage_real_postclassification/{tech}/{tool}/{db}.log",


rule relative_coverage_change_simulations:
    input:
        preclassification=rules.mtb_coverage_simulations_preclassification.output,
        postclassifications=infer_coverage_postclassifications_simulations,
    output:
        RESULTS / "coverage/relative_change/{tech}.csv",
    log:
        LOGS / "relative_coverage_change_simulations/{tech}.log",
    resources:
        runtime="2m",
        mem_mb=500,
    script:
        SCRIPTS / "relative_coverage.py"


use rule relative_coverage_change_simulations as relative_coverage_change_real with:
    input:
        preclassification=rules.mtb_coverage_real_preclassification.output,
        postclassifications=infer_coverage_postclassifications_real,
    output:
        RESULTS / "real/coverage/relative_change/{tech}.csv",
    log:
        LOGS / "relative_coverage_change_real/{tech}.log",
