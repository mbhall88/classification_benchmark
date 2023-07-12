rule sra_human_scrubber:
    input:
        reads=rules.combine_simulated_reads.output.reads,
    output:
        reads=RESULTS / "dehumanise/sra/metagenome.scrubbed.ont.fq.gz",
        removed=RESULTS / "dehumanise/sra/metagenome.removed.ont.fq.gz",
    log:
        LOGS / "sra_human_scrubber.log",
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(16 * GB),
        runtime="1h",
    benchmark:
        BENCH / "dehumanise/sra/ont.tsv"
    container:
        CONTAINERS["sra_human_scrubber"]
    params:
        opts="-x -r",
    shadow:
        "shallow"
    shell:
        """
        exec 2> {log}
        tmpout=$(mktemp -u --suffix='.fq')
        tmpin=$(mktemp -u --suffix='.fq')
        gzip -dc {input.reads} > "$tmpin"
        scrub.sh {params.opts} -p {threads} -o "$tmpout" "$tmpin" >> {log}
        gzip -c "$tmpout" > {output.reads}
        gzip -c "$tmpin".removed_spots > {output.removed}
        """


rule minimap2_human_scrubber:
    input:
        reads=rules.sra_human_scrubber.input.reads,
        ref=rules.download_chm13.output.fasta,
    output:
        aln=RESULTS / "dehumanise/minimap2/metagenome.aln.ont.paf",
    log:
        LOGS / "minimap2_human_scrubber.log",
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(12 * GB),
        runtime="30m",
    benchmark:
        BENCH / "dehumanise/minimap2/ont.tsv"
    container:
        CONTAINERS["minimap2"]
    params:
        opts="-x map-ont --paf-no-hit",
    shell:
        "minimap2 {params.opts} -t {threads} -o {output.aln} {input.ref} {input.reads} 2> {log}"


rule kraken_human_classify:
    input:
        db=rules.build_kraken_human_database.output.db,
        reads=rules.sra_human_scrubber.input.reads,
    output:
        report=RESULTS / "dehumanise/kraken/classify/k{k}l{l}/metagenome.ont.k2report",
        out=RESULTS / "dehumanise/kraken/classify/k{k}l{l}/metagenome.ont.k2",
    log:
        LOGS / "kraken_human_classify/k{k}l{l}/ont.log",
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(8 * GB),
        runtime="1h",
    benchmark:
        BENCH / "dehumanise/kraken/k{k}l{l}/ont.tsv"
    container:
        CONTAINERS["kraken"]
    params:
        opts="--minimum-hit-groups 3 --report-minimizer-data",
    shell:
        """
        kraken2 {params.opts} --threads {threads} --db {input.db} --report {output.report} \
            --output {output.out} {input.reads} 2> {log}
        """
