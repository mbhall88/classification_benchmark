REPEAT = 1


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
        repeat(BENCH / "dehumanise/sra/ont.tsv", REPEAT)
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


rule sra_human_scrubber_illumina:
    input:
        r1=rules.combine_illumina_simulated_reads.output.r1,
        r2=rules.combine_illumina_simulated_reads.output.r2,
    output:
        reads1=RESULTS / "dehumanise/sra/metagenome_R1.scrubbed.illumina.fq.gz",
        reads2=RESULTS / "dehumanise/sra/metagenome_R2.scrubbed.illumina.fq.gz",
        removed1=RESULTS / "dehumanise/sra/metagenome_R1.removed.illumina.fq.gz",
        removed2=RESULTS / "dehumanise/sra/metagenome_R2.removed.illumina.fq.gz",
    log:
        LOGS / "sra_human_scrubber_illumina.log",
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(16 * GB),
        runtime="1h",
    benchmark:
        repeat(BENCH / "dehumanise/sra/illumina.tsv", REPEAT)
    container:
        CONTAINERS["sra_human_scrubber"]
    params:
        opts="-x -r",
    shadow:
        "shallow"
    shell:
        """
        exec 2> {log}
        tmpout1=$(mktemp -u --suffix='.fq')
        tmpout2=$(mktemp -u --suffix='.fq')
        tmpin1=$(mktemp -u --suffix='.fq')
        tmpin2=$(mktemp -u --suffix='.fq')
        gzip -dc {input.r1} > "$tmpin1"
        gzip -dc {input.r2} > "$tmpin2"
        scrub.sh {params.opts} -p {threads} -o "$tmpout1" "$tmpin1" >> {log}
        scrub.sh {params.opts} -p {threads} -o "$tmpout2" "$tmpin2" >> {log}
        gzip -c "$tmpout1" > {output.reads1}
        gzip -c "$tmpout2" > {output.reads2}
        gzip -c "$tmpin1".removed_spots > {output.removed1}
        gzip -c "$tmpin2".removed_spots > {output.removed2}
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
        repeat(BENCH / "dehumanise/minimap2/ont.tsv", REPEAT)
    container:
        CONTAINERS["minimap2"]
    params:
        opts="-x map-ont --paf-no-hit",
    shell:
        "minimap2 {params.opts} -t {threads} -o {output.aln} {input.ref} {input.reads} 2> {log}"


rule minimap2_human_scrubber_illumina:
    input:
        r1=rules.sra_human_scrubber_illumina.input.r1,
        r2=rules.sra_human_scrubber_illumina.input.r2,
        ref=rules.download_chm13.output.fasta,
    output:
        aln=RESULTS / "dehumanise/minimap2/metagenome.aln.illumina.paf",
    log:
        LOGS / "minimap2_human_scrubber_illumina.log",
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(12 * GB),
        runtime="30m",
    benchmark:
        repeat(BENCH / "dehumanise/minimap2/illumina.tsv", REPEAT)
    container:
        CONTAINERS["minimap2"]
    params:
        opts="-x sr --paf-no-hit",
    shell:
        "minimap2 {params.opts} -t {threads} -o {output.aln} {input.ref} {input.r1} {input.r2} 2> {log}"


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
        mem_mb=lambda wildcards, attempt: attempt * int(6 * GB),
        runtime=lambda wildcards, attempt: f"{30 * attempt}m",
    benchmark:
        repeat(BENCH / "dehumanise/kraken/k{k}l{l}/ont.tsv", REPEAT)
    container:
        CONTAINERS["kraken"]
    params:
        opts="--minimum-hit-groups 3 --report-minimizer-data",
    shell:
        """
        kraken2 {params.opts} --threads {threads} --db {input.db} --report {output.report} \
            --output {output.out} {input.reads} 2> {log}
        """


rule miniwinnow_human_scrubber:
    input:
        reads=rules.sra_human_scrubber.input.reads,
        ref=rules.download_chm13.output.fasta,
        aln=rules.minimap2_human_scrubber.output.aln,
    output:
        aln=RESULTS / "dehumanise/miniwinnow/metagenome.aln.ont.paf",
    log:
        LOGS / "miniwinnow_human_scrubber.log",
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(16 * GB),
        runtime="30m",
    benchmark:
        repeat(BENCH / "dehumanise/miniwinnow/ont.tsv", REPEAT)
    conda:
        ENVS / "miniwinnow.yaml"
    shadow:
        "shallow"
    params:
        opts="-x map-ont --paf-no-hit",
    shell:
        """
        exec 2> {log}
        meryl count k=15 output merylDB {input.ref}
        meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt
        awk -F'\t' '$5=="*"' {input.aln} | cut -f1 | seqkit grep -f - {input.reads} | \
        winnowmap {params.opts} -W repetitive_k15.txt -t {threads} -o {output.aln} {input.ref} -
        """


rule hostile_human_scrubber:
    input:
        reads=rules.sra_human_scrubber.input.reads,
    output:
        reads=RESULTS / "dehumanise/hostile/metagenome.ont.clean.fastq.gz",
    log:
        LOGS / "hostile_human_scrubber.log",
    resources:
        mem_mb=int(14 * GB),
        runtime="1h",
    threads: 4
    shadow:
        "shallow"
    benchmark:
        repeat(BENCH / "dehumanise/hostile/ont.tsv", REPEAT)
    container:
        CONTAINERS["hostile"]
    params:
        opts="--aligner minimap2 --force",
        outdir=lambda wildcards, output: Path(output.reads).parent,
    shell:
        "hostile clean {params.opts} --out-dir {params.outdir} --threads {threads} --fastq1 {input.reads} &> {log}"
