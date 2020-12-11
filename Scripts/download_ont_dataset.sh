#usr/bin/bash

# Section __: Data retrieval
mkdir FASTQ

# Parse tasic_fastq_curated.txt and download FASTQ files
while IFS=$'\t' read -r ACC LINK NAME;do
 	wget $LINK -P FASTQ
done < ont_fastq_curated.txt
