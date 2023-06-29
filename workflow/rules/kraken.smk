"""From https://www.nature.com/articles/s41596-022-00738-y#Sec26
Following removal of host DNA, the remaining reads undergo classification using Kraken2Uniq.
The --report-minimizer-data flag forces Kraken 2 to provide unique k-mer counts per classification.
Run Kraken2 using the commands below. We use eight threads for faster run times, and
--minimum-hit-groups 3 for increased classification precision (i.e., fewer false positives).

$ kraken2 --db k2protocol_db --threads 8 --minimum-hit-groups 3 \
--report-minimizer-data --report SRR12486971.k2report \
--paired SRR12486971_1.fastq SRR12486971_2.fastq > SRR12486971.kraken2

Use https://github.com/jenniferlu717/KrakenTools#extract_kraken_readspy to extraxt reads
classified as a user-specified taxID

"""