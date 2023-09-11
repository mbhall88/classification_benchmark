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



def infer_max_db_size_opt(wildcards):
    if wildcards.size == "full":
        return ""
    else:
        return f"--max-db-size {parse_size(wildcards.size)}"


def infer_minimizer_spaces(wildcards):
    l = int(wildcards.l)
    s = int(l / 4)
    return f"--minimizer-spaces {s}"


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


KRAKEN_URLS = {
    "standard": "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20230605.tar.gz",
    "standard-8": "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20230605.tar.gz"
}

rule download_kraken_db:
    output:
        db=directory(RESULTS / "db/kraken/{db}")
    log:
        LOGS / "download_kraken_db/{db}.log"
    resources:
        runtime="1h"
    container:
        CONTAINERS["base"]
    params:
        url=lambda wildcards: KRAKEN_URLS[wildcards.db]
    shell:
        """
        exec &> {log}
        tmpf=$(mktemp -u)
        trap "rm -f $tmpf" EXIT
        wget {params.url} -O $tmpf
        mkdir -p {output.db}
        tar xzf $tmpf -C {output.db}
        rm {output.db}/database*
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
        index=RESULTS / "db/minimap2/db.{preset}.mmi",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(32 * GB),
        runtime="5m",
    log:
        LOGS / "index_decontam_db_with_minimap2/{preset}.log",
    threads: 4
    container:
        CONTAINERS["minimap2"]
    params:
        opts="-x {preset} -I 12G",
    shell:
        "minimap2 {params.opts} -t {threads} -d {output.index} {input.fasta} 2> {log}"


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

rule index_chm13_minimap2:
    input:
        fasta=rules.download_chm13.output.fasta
    output:
        index=RESULTS / "db/chm13.{preset}.mm2"
    log:
        LOGS / "index_chm13_minimap2/{preset}.log"
    resources:
        mem_mb=int(14 * GB),
        runtime="5m"
    container:
        CONTAINERS["minimap2"]
    shell:
        "minimap2 -d {output.index} -x {wildcards.preset} {input.fasta} &> {log}"

rule download_human_pangenome_assemblies:
    input:
        summary=CONFIG / "hprc_assembly_summary.txt"
    output:
        genomes=temp(
            directory(RESULTS / "db/HPRC/2023-08-08/files")
        ),
    log:
        LOGS / "download_human_pangenome_assemblies.log",
    resources:
        runtime="2h",
        mem_mb=int(4 * GB),
    threads: 8
    conda:
        ENVS / "genome_updater.yaml"
    params:
        opts='-m -a -f "genomic.fna.gz"',
        outdir=lambda wildcards, output: Path(output.genomes).parent.parent,
        version=lambda wildcards, output: Path(output.genomes).parts[-2],
    shell:
        """
        genome_updater.sh {params.opts} -o {params.outdir} -b {params.version} \
            -t {threads} -e "{input.summary}" &> {log}
        """



rule prepare_human_pangenome_for_kraken:
    input:
        genomes=rules.download_human_pangenome_assemblies.output.genomes,
        script=SCRIPTS / "prepare_kraken_fasta.py",
    output:
        fasta=temp(RESULTS / "db/HPRC/genomes.fna"),
    log:
        LOGS / "prepare_human_pangenome_for_kraken.log",
    resources:
        runtime="2h",
    container:
        CONTAINERS["python"]
    params:
        opts="-r -T 9606",
    shell:
        "python {input.script} {params.opts} -o {output.fasta} {input.genomes} 2> {log}"


rule build_human_pangenome_kraken_db:
    input:
        fasta=rules.prepare_human_pangenome_for_kraken.output.fasta,
    output:
        db=directory(RESULTS / "db/HPRC/kraken/db")
    log:
        LOGS / "build_human_pangenome_kraken_db.log",
    resources:
        runtime="1d",
        mem_mb=int(16 * GB),
    threads: 16
    container:
        CONTAINERS["kraken"]
    shell:
        """
        exec 2> {log}
        # remove annoying perl warnings
        export LANGUAGE=en_US.UTF-8
        export LC_ALL=en_US.UTF-8
        export LANG=en_US.UTF-8
        export LC_CTYPE=en_US.UTF-8
        >&2 echo "Downloading taxonomy..."
        kraken2-build --download-taxonomy --db {output.db}
        >&2 echo "Adding to library..."
        kraken2-build --add-to-library {input.fasta} --db {output.db}
        >&2 echo "Building..."
        kraken2-build --build --db {output.db} --threads {threads}
        #>&2 echo "Cleaning..."
        #k2 clean --db {output.db}
        """


rule download_GTDB_mycobacterium:
    output:
        genomes=temp(
            directory(RESULTS / "db/GTDB_genus_Mycobacterium/genomes/files")
        ),
        asm_summary=RESULTS
        / "db/GTDB_genus_Mycobacterium/genomes/assembly_summary.txt",
        taxonomy=RESULTS
        / "db/GTDB_genus_Mycobacterium/genomes/bac120_taxonomy.tsv.gz",
    log:
        LOGS / "download_GTDB_mycobacterium.log",
    resources:
        runtime="2h",
        mem_mb=int(4 * GB),
    threads: 32
    conda:
        ENVS / "genome_updater.yaml"
    params:
        opts='-A "species:1" -m -a -M "gtdb" -f "genomic.fna.gz" -g "bacteria" -d "refseq"',
        taxon='-T "f__Mycobacteriaceae,g__Klebsiella,g__Escherichia,g__Enterobacter,g__Salmonella,g__Streptococcus,g__Staphylococcus,g__Pseudomonas,g__Xanthomonas,g__Bifidobacterium"',
        outdir=lambda wildcards, output: Path(output.genomes).parent.parent,
        version=lambda wildcards, output: Path(output.genomes).parts[-2],
    shell:
        """
        genome_updater.sh {params.opts} {params.taxon} -o {params.outdir} -b {params.version} -t {threads} &> {log}
        """


rule prepare_mycobacterium_for_kraken:
    input:
        genomes=rules.download_GTDB_mycobacterium.output.genomes,
        script=SCRIPTS / "prepare_kraken_fasta.py",
        summary=rules.download_GTDB_mycobacterium.output.asm_summary,
    output:
        fasta=temp(RESULTS / "db/GTDB_genus_Mycobacterium/genomes.fna"),
    log:
        LOGS / "prepare_mycobacterium_for_kraken.log",
    resources:
        runtime="2h",
    container:
        CONTAINERS["python"]
    params:
        opts="-r",
        exclude=",".join([
            "GCF_932530395.1",  # MTB
            "GCF_017190695.1",  # M. abscessus
            "GCF_020735285.1",  # M. avium
            "GCA_014701265.1",  # M. kansasii
            "GCF_000013925.1",  # M. ulcerans
            "GCF_016756075.1",  # M. intracellulare
            "GCF_010727125.1",  # M. terrae
            "GCF_001307545.1",  # M. fortuitum
        ]),
    shell:
        "python {input.script} -x {params.exclude} {params.opts} -o {output.fasta} -s {input.summary} {input.genomes} 2> {log}"


rule download_and_add_mycobacterial_taxonomy:
    input:
        names=rules.download_kraken_taxonomy.output.names,
        genomes=rules.prepare_mycobacterium_for_kraken.output.fasta
    output:
        library=directory(RESULTS / "db/GTDB_genus_Mycobacterium/kraken/db/library"),
        taxonomy=directory(RESULTS / "db/GTDB_genus_Mycobacterium/kraken/db/taxonomy"),
    resources:
        mem_mb=int(8 * GB),
        runtime="4h",
    container:
        CONTAINERS["kraken"]
    log:
        LOGS / "download_and_add_mycobacterial_taxonomy.log",
    shell:
        """
        exec &> {log}
        mkdir -p {output.taxonomy}
        cp $(dirname {input.names})/*.dmp {output.taxonomy}
        kraken2-build --add-to-library {input.genomes} --db $(dirname {output.library})
        """


rule build_mycobacterium_kraken_db:
    input:
        lib=rules.download_and_add_mycobacterial_taxonomy.output.library
    output:
        hash=RESULTS / "db/GTDB_genus_Mycobacterium/kraken/db/hash.k2d",
    log:
        LOGS / "build_mycobacterium_kraken_db.log",
    resources:
        runtime="1d",
        mem_mb=int(32 * GB),
    threads: 32
    container:
        CONTAINERS["kraken"]
    shell:
        """
        exec &> {log}
        kraken2-build --build --db $(dirname {output.hash}) --threads {threads}
        #>&2 echo "Cleaning..."
        #k2 clean --db $(dirname {output.hash})
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
        runtime=lambda wildcards, attempt: f"{6 * attempt}h",
        mem_mb=lambda wildcards, attempt: attempt * int(16 * GB),
    threads: 6
    container:
        CONTAINERS["mashtree"]
    shell:
        """
        fofn=$(mktemp -u)
        find $(realpath {input.genomes}) -type f > $fofn
        mashtree --file-of-files $fofn --numcpus {threads} --outtree {output.tree} --outmatrix {output.matrix} 2> {log}
        """

rule create_minimap2_mycobacterium_db:
    input:
        extras=rules.download_mtb_lineage_refs.output.asms
    output:
        db=RESULTS / "db/GTDB_genus_Mycobacterium/Mycobacterium.rep.fna.gz"
    log:
        LOGS / "create_minimap2_mycobacterium_db.log"
    resources:
        runtime="30m"
    threads: 4
    conda:
        ENVS / "genome_updater.yaml"
    params:
        exclude=rules.prepare_mycobacterium_for_kraken.params.exclude,
        opts='-d "refseq" -g "bacteria" -T "g__Mycobacterium" -f "genomic.fna.gz" -M "gtdb" -A 1 -m'
    shell:
        """
        exec 2> {log}
        cat {input.extras} > {output.db}
        tmpd=$(mktemp -d)
        version="version"
        genome_updater.sh {params.opts} -t {threads} -o $tmpd -b "$version"
        (find $tmpd -type f -name '*.fna.gz' -print0 | xargs -0 cat) >> {output.db}
        """


rule index_mycobacterium_db_minimap2:
    input:
        fasta=rules.create_minimap2_mycobacterium_db.output.db,
    output:
        index=RESULTS / "db/GTDB_genus_Mycobacterium/Mycobacterium.rep.{preset}.mmi",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(16 * GB),
        runtime="5m",
    log:
        LOGS / "index_mycobacterium_db_minimap2/{preset}.log",
    threads: 4
    container:
        CONTAINERS["minimap2"]
    params:
        opts="-x {preset} -I 12G",
    shell:
        "minimap2 {params.opts} -t {threads} -d {output.index} {input.fasta} 2> {log}"


rule create_minimap2_mtb_db:
    input:
        extras=rules.download_mtb_lineage_refs.output.asms
    output:
        db=RESULTS / "db/GTDB_genus_Mycobacterium/MTB.fna.gz"
    log:
        LOGS / "create_minimap2_mtb_db.log"
    resources:
        runtime="5m"
    params:
        url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/195/955/GCA_000195955.2_ASM19595v2/GCA_000195955.2_ASM19595v2_genomic.fna.gz"
    shell:
        """
        cat {input.extras} > {output.db} 2> {log}
        wget {params.url} -O - >> {output.db} 2>> {log}
        """


rule index_mtb_db_minimap2:
    input:
        fasta=rules.create_minimap2_mtb_db.output.db,
    output:
        index=RESULTS / "db/GTDB_genus_Mycobacterium/MTB.{preset}.mmi",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(2 * GB),
        runtime="5m",
    log:
        LOGS / "index_mtb_db_minimap2/{preset}.log",
    threads: 4
    container:
        CONTAINERS["minimap2"]
    params:
        opts="-x {preset}",
    shell:
        "minimap2 {params.opts} -t {threads} -d {output.index} {input.fasta} 2> {log}"

