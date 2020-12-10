#usr/bin/bash

# Section __: Data retrieval
mkdir FASTQ

# Parse tasic_fastq_curated.txt and download FASTQ files
while IFS=$'\t' read -r NAME CLASS LINK CLASSNUM;do
 	wget $LINK -P FASTQ
 	tar -xvf "FASTQ/$NAME.fastq.tar" -C FASTQ
 	rm FASTQ/$NAME.fastq.tar
done < tasic_fastq_curated.txt

# Merging

# Create arrays to store filenames and groupnames
declare -A fastqArray
declare -A fastqNames

while IFS=$'\t' read -r NAME CLASS LINK CLASSNUM;do
	fastqArray["$CLASSNUM.1"]="${fastqArray["$CLASSNUM.1"]} FASTQ/${NAME}_R1.fastq.gz"
	fastqArray["$CLASSNUM.2"]="${fastqArray["$CLASSNUM.2"]} FASTQ/${NAME}_R2.fastq.gz"
	fastqNames["$CLASSNUM.1"]=${CLASS}_R1
	fastqNames["$CLASSNUM.2"]=${CLASS}_R2
done < tasic_fastq_curated.txt

# Loop groups and concatenate FASTQ by direction
for (( n=1; n<=4; n++ ));do
	cat ${fastqArray["$n.1"]} > FASTQ/${fastqNames["$n.1"]}.fastq.gz
	cat ${fastqArray["$n.2"]} > FASTQ/${fastqNames["$n.2"]}.fastq.gz
done