from humanfriendly import parse_size


def infer_kraken_memory(wildcards, attempt):
    if wildcards.size == "full":
        mb = 40 * GB
    else:
        bytes = parse_size(wildcards.size)
        mb = bytes / (10**6) + (4 * GB)  # add 4gb for wriggle room
    return int(mb) * attempt


def infer_max_db_size_opt(wildcards):
    if wildcards.size == "full":
        return ""
    else:
        return f"--max-db-size {parse_size(wildcards.size)}"


def infer_minimizer_spaces(wildcards):
    l = int(wildcards.l)
    s = int(l / 4)
    return f"--minimizer-spaces {s}"


rule build_kraken_database:
    output:
        db=directory(RESULTS / "kraken/db/k{k}/l{l}/{size}"),
    log:
        LOGS / "build_kraken_database/k{k}/l{l}/{size}.log",
    threads: 16
    resources:
        mem_mb=infer_kraken_memory,
        runtime="1d",
    benchmark:
        BENCH / "kraken/build/k{k}/l{l}/{size}.tsv"
    params:
        opts="--standard --kmer-len {k} --minimizer-len {l} --use-ftp",
        max_db_size=infer_max_db_size_opt,
        spaces=infer_minimizer_spaces,
    container:
        CONTAINERS["kraken"]
    shell:
        "kraken2-build {params.opts} {params.max_db_size} {params.spaces} --threads {threads} --db {output.db} &> {log}"


rule build_decontamination_db:
    output:
        fasta=temp(RESULTS / "db/remove_contam.fa.gz"),
        metadata=temp(RESULTS / "db/remove_contam.tsv"),
    resources:
        runtime="90m",
    params:
        script=SCRIPTS / "download_tb_reference_files.pl",
        outdir=lambda wildcards, output: Path(output.fasta).parent,
    container:
        CONTAINERS["tb_decontam"]
    log:
        LOGS / "build_decontamination_db.log",
    shadow:
        "shallow"
    shell:
        "perl {params.script} {params.outdir} &> {log}"


rule download_ancestral_genome:
    output:
        asm=RESULTS / "db/mtb_ancestor.fa.gz",
    log:
        LOGS / "download_ancestral_genome.log",
    resources:
        mem_mb=int(0.5 * GB),
        runtime="3m",
    params:
        url=config["ancestral_genome"]["url"],
        md5=config["ancestral_genome"]["md5"],
    container:
        CONTAINERS["base"]
    shadow:
        "shallow"
    shell:
        """
        wget -O asm.fa {params.url} 2> {log}
        hash=$(md5sum asm.fa | cut -d' ' -f1)
        if [ "$hash" != {params.md5} ]; then
          echo "MD5 does not match" >> {log}
          exit 1
        fi
        gzip -c asm.fa > {output.asm} 2>> {log} 
        """


rule download_mtb_lineage_refs:
    output:
        asms=expand(str(RESULTS / "db/{sample}.fasta.gz"), sample=GRAMTOOLS_SAMPLES),
    log:
        LOGS / "download_mtb_lineage_refs.log",
    resources:
        mem_mb=int(0.5 * GB),
        runtime="5m",
    params:
        url=config["mtb_gramtools"]["url"],
        md5=config["mtb_gramtools"]["md5"],
        parent=lambda wildcards, output: Path(output.asms[0]).parent,
    container:
        CONTAINERS["base"]
    shadow:
        "shallow"
    shell:
        """
        exec 2> {log}
        wget -O mtb_asms.tar {params.url}
        hash=$(md5sum mtb_asms.tar | cut -d' ' -f1)
        if [ "$hash" != {params.md5} ]; then
          echo "MD5 does not match" >> {log}
          exit 1
        fi
        tar xf mtb_asms.tar
        cd analysis/input_data/mtuberculosis/pacb_ilmn/pacb_assemblies || exit 1
        for f in *.fasta.gz;
        do
            dst={params.parent}/$f
            mv $f $dst
        done
        """


rule combine_references:
    input:
        decontam_fa=rules.build_decontamination_db.output.fasta,
        decontam_metadata=rules.build_decontamination_db.output.metadata,
        ancestral=rules.download_ancestral_genome.output.asm,
        gramtools_asms=rules.download_mtb_lineage_refs.output.asms,
    output:
        fasta=RESULTS / "db/db.fa.gz",
        metadata=RESULTS / "db/db.tsv",
    log:
        LOGS / "combine_references.log",
    resources:
        runtime="10m",
    container:
        CONTAINERS["python"]
    script:
        SCRIPTS / "combine_references.py"


rule index_db_with_minimap2:
    input:
        fasta=rules.combine_references.output.fasta,
    output:
        index=RESULTS / "db/db.fa.gz.map-ont.mmi",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(32 * GB),
        runtime="5m",
    log:
        LOGS / "index_decontam_db_with_minimap2.log",
    threads: 4
    container:
        CONTAINERS["minimap2"]
    params:
        opts="-x map-ont -I 12G --idx-no-seq",
    shell:
        "minimap2 {params.opts} -t {threads} -d {output.index} {input.fasta} 2> {log}"


rule prepare_spumoni_index:
    input:
        fasta=rules.combine_references.output.fasta,
        metadata=rules.combine_references.output.metadata,
    output:
        file_list=RESULTS / "spumoni/db/file_list.txt",
        fasta_dir=temp(directory(RESULTS / "spumoni/db/individual_fastas")),
    log:
        LOGS / "prepare_spumoni_index.log",
    resources:
        runtime="15m",
    container:
        CONTAINERS["python"]
    script:
        SCRIPTS / "prepare_spumoni_index.py"


SPUMONI_EXTS = [".ms", ".slp", ".msnulldb", ".spumoni", ".pmlnulldb"]


rule build_spumoni_index:
    input:
        file_list=rules.prepare_spumoni_index.output.file_list,
        fasta_dir=rules.prepare_spumoni_index.output.fasta_dir,
    output:
        multiext(str(RESULTS / "spumoni/db/db"), *SPUMONI_EXTS),
    log:
        LOGS / "build_spumoni_index.log",
    resources:
        mem_mb=lambda wildcards, attempt: int(12 * GB) * attempt,
        runtime="2d",
    container:
        CONTAINERS["spumoni"]
    params:
        opts="-M -P -d -m",
        prefix=lambda wildcards, output: Path(output[0]).with_suffix(""),
    shell:
        "spumoni build {params.opts} -i {input.file_list} -o {params.prefix} &> {log}"
