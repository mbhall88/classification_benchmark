#!/usr/bin/env bash
set -euxo pipefail

# kingfisher get -r ERR3241786 ERR3239713 ERR3243073 \
#     -m ena-ascp ena-ftp prefetch --download-threads 8 --check-md5sums

for acc in ERR3241786 ERR3239713 ERR3243073; do
    seqkit pair -1 ${acc}_1.fastq.gz -2 ${acc}_2.fastq.gz --id-regexp '^(\S+(\s\d+)?)[\/\.][12]'
    rm ${acc}_1.fastq.gz ${acc}_2.fastq.gz
    rasusa -i ${acc}_1.paired.fastq.gz ${acc}_2.paired.fastq.gz -o ${acc}_1.fq.gz ${acc}_2.fq.gz -b 1g -s 222
    rm ${acc}_1.paired.fastq.gz ${acc}_2.paired.fastq.gz
done

seqkit seq ERR3241786_1.fq.gz ERR3239713_1.fq.gz ERR3243073_1.fq.gz -o human_1.fq.gz
rm ERR3241786_1.fq.gz ERR3239713_1.fq.gz ERR3243073_1.fq.gz
seqkit seq ERR3241786_2.fq.gz ERR3239713_2.fq.gz ERR3243073_2.fq.gz -o human_2.fq.gz
rm ERR3241786_2.fq.gz ERR3239713_2.fq.gz ERR3243073_2.fq.gz
