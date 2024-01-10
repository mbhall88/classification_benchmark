# Pangenome databases provide superior host removal and mycobacteria classification from clinical metagenomic data

> Hall, M, Coin, L., Pangenome databases provide superior host removal and mycobacteria classification from clinical metagenomic data. bioRxiv 2023. doi: [10.1101/2023.09.18.558339][doi]

Benchmarking different ways of doing read (taxonomic) classification, with a focus on
removal of contamination and classification of _M. tuberculosis_ reads.

This repository contains the code and snakemake pipeline to build/download the
databases, obtain all results from [the paper][doi], along with accompanying configuration
files.

Custom databases have all been uploaded to Zenodo, along with the simulated reads:

- Nanopore simulated metagenomic reads - <https://doi.org/10.5281/zenodo.8339788>
- Illumina simulated metagenomic reads - <https://doi.org/10.5281/zenodo.8339790>
- Nanopore and Illumina artificial real reads - <https://doi.org/10.5281/zenodo.10472796>
- Kraken2 database built from the Human Pangenome Reference Consortium
  genomes - <https://doi.org/10.5281/zenodo.8339731>
- Kraken2 database built from the kraken2 Human
  library - <https://doi.org/10.5281/zenodo.8339699>
- Kraken2 database built from a *Mycobacterium* representative set of
  genomes - <https://doi.org/10.5281/zenodo.8339821>
- A (fasta) database of representative genomes from the *Mycobacterium*
  genus - <https://doi.org/10.5281/zenodo.8339940>
- A (fasta) database of *M. tuberculosis* genomes from a variety of
  lineages - <https://doi.org/10.5281/zenodo.8339947>
- The fasta file built from the [Clockwork](https://github.com/iqbal-lab-org/clockwork)
  decontamination pipeline - <https://doi.org/10.5281/zenodo.8339802>

## Example usage

We provide some usage examples showing how to download the databases and then use them
on your reads.

### Human read removal

The method we found to give the best balance of runtime, memory usage, and precision and
recall was kraken2 with a database built from the Human Pangenome Reference Consortium
genomes.

This example has been wrapped into a standalone tool called [`nohuman`](https://github.com/mbhall88/nohuman/) which takes a fastq as input and returns a fastq with human reads removed.

#### Download human database

```
mkdir HPRC_db/
cd HPRC_db
URL="https://zenodo.org/record/8339732/files/k2_HPRC_20230810.tar.gz"
wget "$URL"
tar -xzf k2_HPRC_20230810.tar.gz
rm k2_HPRC_20230810.tar.gz
```

#### Run kraken2 with HPRC database

You'll need [kraken2](https://github.com/DerrickWood/kraken2) installed for this step.

```
kraken2 --threads 4 --db HPRC_db/ --output classifications.tsv reads.fq
```

If you are using Illumina reads, a slight adjustment is needed

```
kraken2 --paired --threads 4 --db HPRC_db/ --output classifications.tsv reads_1.fq reads_2.fq
```

#### Extract non-human reads

You'll need [seqkit](https://github.com/shenwei356/seqkit) installed for this step

For Nanopore data

```
awk -F'\t' '$1=="U" {print $2}' classifications.tsv | \
  seqkit grep -f - -o reads.depleted.fq reads.fq
```

For Illumina data

```
awk -F'\t' '$1=="U" {print $2}' classifications.tsv > ids.txt
seqkit grep --id-regexp '^(\S+)/[12]' -f ids.txt -o reads_1.depleted.fq reads_1.fq
seqkit grep --id-regexp '^(\S+)/[12]' -f ids.txt -o reads_2.depleted.fq reads_2.fq
```

### *M. tuberculosis* classification/enrichment

For this step we recommend either [minimap2](https://github.com/lh3/minimap2) or kraken2
with a *Mycobacterium* genus database. We leave it to the user to decide which approach
they prefer based on the results in our manuscript.

#### Download databases

```
mkdir Mycobacterium_db
cd Mycobacterium_db
# download database for use with minimap2
URL="https://zenodo.org/record/8339941/files/Mycobacterium.rep.fna.gz"
wget "$URL"
IDS_URL="https://zenodo.org/record/8343322/files/mtb.ids"
wget "$IDS_URL"
# download kraken database
URL="https://zenodo.org/record/8339822/files/k2_Mycobacterium_20230817.tar.gz"
wget "$URL"
tar -xzf k2_Mycobacterium_20230817.tar.gz
rm k2_Mycobacterium_20230817.tar.gz
```

#### Classify reads

**minimap2**

```
# nanopore
minimap2 --secondary=no -c -t 4 -x map-ont -o reads.aln.paf Mycobacterium_db/Mycobacterium.rep.fna.gz reads.depleted.fq
# illumina
minimap2 --secondary=no -c -t 4 -x sr -o reads.aln.paf Mycobacterium_db/Mycobacterium.rep.fna.gz reads_1.depleted.fq reads_2.depleted.fq
```

**kraken2**

```
# nanopore
kraken2 --db Mycobacterium_db --threads 4 --report myco.kreport --output classifications.myco.tsv reads.depleted.fq
# illumina
kraken2 --db Mycobacterium_db --paired --threads 4 --report myco.kreport --output classifications.myco.tsv reads_1.depleted.fq reads_2.depleted.fq
```

#### Extract *M. tuberculosis* reads

**minimap2**

```
# nanopore
grep -Ff Mycobacterium_db/mtb.ids reads.aln.paf | cut -f1 | \
  seqkit grep -f - -o reads.enriched.fq reads.depleted.fq
# illumina
grep -Ff Mycobacterium_db/mtb.ids reads.aln.paf | cut -f1 > keep.ids
seqkit grep -f keep.ids -o reads_1.enriched.fq reads_1.depleted.fq
seqkit grep -f keep.ids -o reads_2.enriched.fq reads_2.depleted.fq
```

**kraken2**

We'll use
the [`extract_kraken_reads.py` script](https://github.com/jenniferlu717/KrakenTools#extract_kraken_readspy)
for this

```
# nanopore
python extract_kraken_reads.py -k classifications.myco.tsv -1 reads.depleted.fq -o reads.enriched.fq -t 1773 -r myco.kreport --include-children
# illumina
python extract_kraken_reads.py -k classifications.myco.tsv -1 reads_1.depleted.fq -2 reads_2.depleted.fq -o reads_1.enriched.fq -o2 reads_2.enriched.fq -t 1773 -r myco.kreport --include-children
```

[doi]: https://doi.org/10.1101/2023.09.18.558339 
