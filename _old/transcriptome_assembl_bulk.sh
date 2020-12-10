#usr/bin/bash

# create subdirectories for experiment and detect cores
mkdir FASTQ
mkdir Hisat2_mm10_index
mkdir Hisat2_SAMs
mkdir Hisat2_sorted_BAMs
mkdir Stringtie_gtf
cores=`nproc`


# download fastq files from SRA
cd FASTQ
fastq-dump SRR606187 SRR606188 SRR606189 SRR606190 # Add fastq-dump to list of programs


# download Hisat2 indices for mm10 and run alignment
cd ..
wget https://genome-idx.s3.amazonaws.com/hisat/mm10_genome.tar.gz -P Hisat2_mm10_index
tar -xf "Hisat2_mm10_index/mm10_genome.tar.gz" -C Hisat2_mm10_index
rm Hisat2_mm10_index/mm10_genome.tar.gz

for fastq in SRR606187 SRR606188 SRR606189 SRR606190;do
	hisat2 -p $cores --dta -x Hisat2_mm10_index/mm10/genome -1 FASTQ/"$fastq"_1.fastq -2 FASTQ/"$fastq"_2.fastq -S Hisat2_SAMs/$fastq.sam
done


# convert SAMs to BAMs
for filename in Hisat2_SAMs/*.sam; do
	file=${filename##*/}
	file=${file%.sam}
	samtools view -@ $cores -Su $filename | samtools sort -@ $cores -o Hisat2_sorted_BAMs/$file.bam
done
rm -r Hisat2_SAMs


# download reference transcriptome and assemble transcripts using stringtie
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
gzip -d gencode.vM25.annotation.gtf.gz

for filename in Hisat2_sorted_BAMs/*.bam; do
	file=${filename##*/}
	file=${file%.bam}
	stringtie $filename -p $cores -o Stringtie_gtf/$file.gtf -G gencode.vM25.annotation.gtf
done


# merge assembled transcriptome
declare -a array_exp
for filename in Stringtie_gtf/*.gtf; do
	array_exp+=($filename)
done
stringtie --merge -G gencode.vM25.annotation.gtf -o burge_merged.gtf ${array_exp[*]}
