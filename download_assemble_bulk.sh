#usr/bin/bash

# Section _____: Retrieving data
mkdir FASTQ
fastq-dump -O FASTQ --split-files --gzip SRR606187 SRR606188 SRR606189 SRR606190 


# Section _____: Aligning reads

# Download mm10 index if absent
if [ ! -d "mm10" ]; then
	wget https://genome-idx.s3.amazonaws.com/hisat/mm10_genome.tar.gz
	tar -xf "mm10_genome.tar.gz"
	rm mm10_genome.tar.gz
fi

# Align reads to mm10 genome
cores=`nproc`
mkdir Hisat2_SAMs
for fastq in SRR606187 SRR606188 SRR606189 SRR606190;do
	hisat2 -p $cores --dta -x mm10/genome -1 FASTQ/"$fastq"_1.fastq.gz -2 FASTQ/"$fastq"_2.fastq.gz -S Hisat2_SAMs/$fastq.sam
done

# Convert SAMs to BAMs
mkdir Hisat2_sorted_BAMs
for file in SRR606187 SRR606188 SRR606189 SRR606190; do
	samtools view -@ $cores -Su Hisat2_SAMs/$file.sam | samtools sort -@ $cores -o Hisat2_sorted_BAMs/$file.bam
done


# Section ____: Assembling transcriptome

# Download reference transcriptome, if absent
if [ ! -f "gencode.vM25.annotation.gtf" ]; then
	wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
	gzip -d gencode.vM25.annotation.gtf.gz
fi

# Assemble transcriptome using StringTie
mkdir Stringtie_gtf
for file in SRR606187 SRR606188 SRR606189 SRR606190; do
	stringtie Hisat2_sorted_BAMs/$file.bam -p $cores -o Stringtie_gtf/$file.gtf -G gencode.vM25.annotation.gtf
done

# Merge transcriptome
stringtie --merge -G gencode.vM25.annotation.gtf -o bulk_merged.gtf Stringtie_gtf/SRR606187.gtf Stringtie_gtf/SRR606188.gtf Stringtie_gtf/SRR606189.gtf Stringtie_gtf/SRR606190.gtf
