# Contamination removal and classification of *Mycobacterium tuberculosis* metagenomic samples
Benchmarking different ways of doing read (taxonomic) classification, with a focus on removal of contamination and classification of _M. tuberculosis_ reads.

This repository contains the code and snakemake pipeline to build/download the databases, obtain all results from the paper, along with accompanying configuration files.

Custom databases have all been uploaded to Zenodo, along with the simulated reads:

- Nanopore simulated metagenomic reads - https://doi.org/10.5281/zenodo.8339788
- Illumina simulated metagenomic reads - https://doi.org/10.5281/zenodo.8339790
- Kraken2 database built from the Human Pangenome Reference Consortium genomes - https://doi.org/10.5281/zenodo.8339731
- Kraken2 database built from the kraken2 Human library - https://doi.org/10.5281/zenodo.8339699
- Kraken2 database built from a *Mycobacterium* representative set of genomes - https://doi.org/10.5281/zenodo.8339821
- A (fasta) database of representative genomes from the *Mycobacterium* genus - https://doi.org/10.5281/zenodo.8339940
- A (fasta) database of *M. tuberculosis* genomes from a variety of lineages - https://doi.org/10.5281/zenodo.8339947
- The fasta file built from the [Clockwork](https://github.com/iqbal-lab-org/clockwork) decontamination pipeline - https://doi.org/10.5281/zenodo.8339802

