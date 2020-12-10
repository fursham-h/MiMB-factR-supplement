#usr/bin/bash

# Section __: Data retrieval
mkdir FASTQ

# Parse tasic_fastq_curated.txt and download FASTQ files
while IFS=$'\t' read -r NAME CLASS LINK CLASSNUM;do
 	wget $LINK -P FASTQ
 	tar -xvf "FASTQ/$NAME.fastq.tar" -C FASTQ
 	rm FASTQ/$NAME.fastq.tar
done < tasic_fastq_curated.txt
