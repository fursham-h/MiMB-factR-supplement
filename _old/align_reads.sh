#usr/bin/bash


# Section __: Read alignment

# Download Hisat2 indices for mm10 if absent
if [ ! -d "mm10" ]; then
	wget https://genome-idx.s3.amazonaws.com/hisat/mm10_genome.tar.gz
	tar -xf "mm10_genome.tar.gz"
	rm mm10_genome.tar.gz
fi


# pre-step: detect cpu cores
cores=`nproc`

# Align reads to mm10 genome
mkdir Hisat2_SAMs
for fastq in "$@";do
	if [ -f "FASTQ/${fastq}_R1.fastq.gz" ]; then
		hisat2 -p $cores --dta -x mm10/genome -1 FASTQ/"$fastq"_R1.fastq.gz -2 FASTQ/"$fastq"_R2.fastq.gz -S Hisat2_SAMs/$fastq.sam
	else
		hisat2 -p $cores --dta -x mm10/genome -1 FASTQ/"$fastq"_1.fastq.gz -2 FASTQ/"$fastq"_2.fastq.gz -S Hisat2_SAMs/$fastq.sam
	fi
done

# Convert SAMs to BAMs
mkdir Hisat2_sorted_BAMs
for filename in Hisat2_SAMs/*.sam; do
	file=${filename##*/}
	file=${file%.sam}
	samtools view -@ $cores -Su $filename | samtools sort -@ $cores -o Hisat2_sorted_BAMs/$file.bam
done