#usr/bin/bash

# Section __: Data retrieval
mkdir FASTQ
fastq-dump -O FASTQ --split-files --gzip SRR606187 SRR606188 SRR606189 SRR606190 
