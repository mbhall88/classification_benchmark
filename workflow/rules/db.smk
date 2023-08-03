

rule download_kraken_taxonomy:
    output:
        names=RESULTS / "kraken/db/taxonomy/names.dmp",
        nodes=RESULTS / "kraken/db/taxonomy/nodes.dmp",
    resources:
        mem_mb=1000,
        runtime="12h",
    container:
        CONTAINERS["kraken"]
    params:
        db=lambda wildcards, output: Path(output.names).parent.parent,
    log:
        LOGS / "download_kraken_taxonomy.log",
    shell:
        """
        k2 download-taxonomy --db {params.db} 2> {log}
        """


rule download_kraken_bacteria_db:
    input:
        rules.download_kraken_taxonomy.output.names,
    output:
        fasta=RESULTS / "kraken/db/library/bacteria/library.fna",
    resources:
        mem_mb=int(8 * GB),
        runtime="1w",
    container:
        CONTAINERS["kraken"]
    params:
        db=lambda wildcards, input: Path(input[0]).parent.parent,
    log:
        LOGS / "download_kraken_bacteria_db.log",
    shell:
        """
        k2 download-library --db {params.db} --library bacteria 2> {log}
        """


rule download_kraken_archaea_db:
    input:
        rules.download_kraken_taxonomy.output.names,
    output:
        fasta=RESULTS / "kraken/db/library/archaea/library.fna",
    resources:
        mem_mb=int(8 * GB),
        runtime="1w",
    container:
        CONTAINERS["kraken"]
    params:
        db=lambda wildcards, input: Path(input[0]).parent.parent,
    log:
        LOGS / "download_kraken_archaea_db.log",
    shell:
        """
        k2 download-library --db {params.db} --library archaea 2> {log}
        """


rule download_kraken_viral_db:
    input:
        rules.download_kraken_taxonomy.output.names,
    output:
        fasta=RESULTS / "kraken/db/library/viral/library.fna",
    resources:
        mem_mb=int(8 * GB),
        runtime="1w",
    container:
        CONTAINERS["kraken"]
    params:
        db=lambda wildcards, input: Path(input[0]).parent.parent,
    log:
        LOGS / "download_kraken_viral_db.log",
    shell:
        """
        k2 download-library --db {params.db} --library viral 2> {log}
        """


rule download_kraken_human_db:
    input:
        rules.download_kraken_taxonomy.output.names,
    output:
        fasta=RESULTS / "kraken/db/library/human/library.fna",
    resources:
        mem_mb=int(8 * GB),
        runtime="1w",
    container:
        CONTAINERS["kraken"]
    params:
        db=lambda wildcards, input: Path(input[0]).parent.parent,
    log:
        LOGS / "download_kraken_human_db.log",
    shell:
        """
        k2 download-library --db {params.db} --library human 2> {log}
        """


rule download_kraken_plasmid_db:
    input:
        rules.download_kraken_taxonomy.output.names,
    output:
        fasta=RESULTS / "kraken/db/library/plasmid/library.fna",
    resources:
        mem_mb=int(8 * GB),
        runtime="1w",
    container:
        CONTAINERS["kraken"]
    params:
        db=lambda wildcards, input: Path(input[0]).parent.parent,
    log:
        LOGS / "download_kraken_plasmid_db.log",
    shell:
        """
        k2 download-library --db {params.db} --library plasmid 2> {log}
        """


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
    input:
        libs=[RESULTS / f"kraken/db/library/{lib}/library.fna" for lib in kraken_libs],
    output:
        db=directory(RESULTS / "kraken/k{k}/l{l}/{size}/db"),
    log:
        LOGS / "build_kraken_database/k{k}/l{l}/{size}.log",
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt: int(80 * GB) * attempt,
        runtime=lambda wildcards, attempt: f"{2*attempt}d",
    benchmark:
        BENCH / "kraken/build/k{k}/l{l}/{size}.tsv"
    params:
        opts="--kmer-len {k} --minimizer-len {l}",
        max_db_size=infer_max_db_size_opt,
        spaces=infer_minimizer_spaces,
        original_db=lambda wildcards, input: Path(input.libs[0]).parent.parent.parent,
    container:
        CONTAINERS["kraken"]
    shell:
        """
        echo "Copying original db" > {log}
        cp -r {params.original_db} {output.db} 2>> {log}
        echo "Finished copying db" >> {log}
        k2 build {params.opts} {params.max_db_size} {params.spaces} --threads {threads} --db {output.db} &>> {log}
        """


rule build_kraken_human_database:
    input:
        libs=rules.download_kraken_human_db.output.fasta,
        taxonomy=rules.download_kraken_taxonomy.output.names,
    output:
        db=directory(RESULTS / "dehumanise/kraken/db/k{k}/l{l}/db"),
    log:
        LOGS / "build_kraken_human_database/k{k}/l{l}.log",
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: int(16 * GB) * attempt,
        runtime=lambda wildcards, attempt: f"{attempt * 4}h",
    params:
        opts="--kmer-len {k} --minimizer-len {l}",
        spaces=infer_minimizer_spaces,
    container:
        CONTAINERS["kraken"]
    shell:
        """
        echo "Copying original db" > {log}
        mkdir -p {output.db}/library 2>> {log}
        cp -r $(dirname {input.taxonomy}) {output.db}/taxonomy 2>> {log}
        cp -r $(dirname {input.libs}) {output.db}/library/ 2>> {log}
        echo "Finished copying db" >> {log}
        k2 build {params.opts} {params.spaces} --threads {threads} --db {output.db} &>> {log}
        """


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
        fasta=RESULTS / "db/db.fa",
        metadata=RESULTS / "db/db.tsv",
    log:
        LOGS / "combine_references.log",
    resources:
        runtime="10m",
    container:
        CONTAINERS["python"]
    script:
        SCRIPTS / "combine_references.py"


rule faidx_db:
    input:
        fasta=rules.combine_references.output.fasta,
    output:
        fasta=RESULTS / "db/db.fa.gz",
        faidx=RESULTS / "db/db.fa.gz.gzi",
    log:
        LOGS / "faidx_db.log",
    resources:
        runtime="30m",
    container:
        CONTAINERS["samtools"]
    shell:
        """
        bgzip --index -f {input} 2> {log}
        """


rule index_db_with_minimap2:
    input:
        fasta=rules.faidx_db.output.fasta,
    output:
        index=RESULTS / "db/minimap2/db.fa.gz.map-ont.mmi",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(32 * GB),
        runtime="5m",
    benchmark:
        BENCH / "minimap2/index.tsv"
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
        fasta=rules.faidx_db.output.fasta,
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


rule download_chm13:
    output:
        fasta=RESULTS / "db/chm13.fa.gz",
    log:
        LOGS / "download_chm13.log",
    resources:
        runtime="10m",
    container:
        CONTAINERS["base"]
    params:
        # chm13 plus HLA taken from https://github.com/bede/hostile
        url="https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/human-t2t-hla.fa.gz",
    shell:
        "wget {params.url} -O {output.fasta} 2> {log}"


rule download_hg02886:
    output:
        fasta=RESULTS / "db/hg02886.fa.gz",
    log:
        LOGS / "download_hg02886.log",
    resources:
        runtime="10m",
    container:
        CONTAINERS["base"]
    params:
        # from https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_018470455.1/
        # cite https://www.nature.com/articles/s41586-023-05896-x
        url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/018/470/455/GCA_018470455.1_HG02886.pri.mat.f1_v2/GCA_018470455.1_HG02886.pri.mat.f1_v2_genomic.fna.gz",
    shell:
        "wget -nv {params.url} -O {output.fasta} 2> {log}"


rule download_GTDB_mycobacterium:
    output:
        genomes=directory(RESULTS / "db/GTDB_genus_Mycobacterium/2023-08-03/files"),
        asm_summary=RESULTS
        / "db/GTDB_genus_Mycobacterium/2023-08-03/assembly_summary.txt",
        taxonomy=RESULTS
        / "db/GTDB_genus_Mycobacterium/2023-08-03/bac120_taxonomy.tsv.gz",
    log:
        LOGS / "download_GTDB_mycobacterium.log",
    resources:
        runtime="1d",
        mem_mb=int(4 * GB),
    threads: 8
    conda:
        ENVS / "genome_updater.yaml"
    params:
        opts='-A 0 -m -a -M "gtdb" -f "genomic.fna.gz" -g "bacteria" -d "refseq,genbank"',
        taxon='-T "g__Mycobacterium"',
        outdir=lambda wildcards, output: Path(output.genomes).parent.parent,
        version=lambda wildcards, output: Path(output.genomes).parts[-2],
    shell:
        """
        genome_updater.sh {params.opts} {params.taxon} -o {params.outdir} -b {params.version} -t {threads} &> {log}
        """


rule mycobacterium_full_tree:
    input:
        genomes=rules.download_GTDB_mycobacterium.output.genomes,
    output:
        tree=RESULTS / "db/GTDB_genus_Mycobacterium/tree.full.dnd",
        matrix=RESULTS / "db/GTDB_genus_Mycobacterium/distmatrix.full.tsv",
    log:
        LOGS / "mycobacterium_full_tree.log",
    resources:
        runtime="6h",
        mem_mb=int(8 * GB),
    threads: 6
    container:
        CONTAINERS["mashtree"]
    shell:
        """
        fofn=$(mktemp -u)
        find $(realpath {input.genomes}) -type f > $fofn
        mashtree --file-of-files $fofn --numcpus {threads} --outtree {output.tree} --outmatrix {output.matrix} 2> {log}
        """


rule dereplicate_mycobacterium_assemblies:
    input:
        genomes=rules.download_GTDB_mycobacterium.output.genomes,
        script=SCRIPTS / "dereplicator.py",
    output:
        outdir=directory(RESULTS / "db/GTDB_genus_Mycobacterium/dereplicate"),
    log:
        LOGS / "dereplicate_mycobacterium_assemblies.log",
    threads: 16
    resources:
        runtime="6h",
        mem_mb=int(16 * GB),
    conda:
        ENVS / "derep.yaml"
    params:
        opts="-d 0.01 --verbose",
    shell:
        """
        find {input.genomes} -type f | wc l > {log}
        python {input.script} {params.opts} --threads {threads} {input.genomes} {output.outdir} 2>> {log}
        find {output.outdir} -type f | wc l >> {log}
        """

rule mycobacterium_rep_tree:
    input:
        genomes=rules.dereplicate_mycobacterium_assemblies.output.outdir,
    output:
        tree=RESULTS / "db/GTDB_genus_Mycobacterium/tree.rep.dnd",
        matrix=RESULTS / "db/GTDB_genus_Mycobacterium/distmatrix.rep.tsv",
    log:
        LOGS / "mycobacterium_rep_tree.log",
    resources:
        runtime="1h",
        mem_mb=int(8 * GB),
    threads: 6
    container:
        CONTAINERS["mashtree"]
    shell:
        """
        fofn=$(mktemp -u)
        find $(realpath {input.genomes}) -type f > $fofn
        mashtree --file-of-files $fofn --numcpus {threads} --outtree {output.tree} --outmatrix {output.matrix} 2> {log}
        """


