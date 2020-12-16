#usr/bin/bash

# Section __: Data retrieval
mkdir FASTQ

# Parse tasic_fastq_curated.txt and download FASTQ files
while IFS=$'\t' read -r ACC LINK NAME;do
 	wget $LINK -P FASTQ
done < ont_fastq_curated.txt


if [ ! -f "gencode.vM25.annotation.gtf" ]; then
	wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz
	gzip -d GRCm38.primary_assembly.genome.fa.gz
fi


cores=`nproc`
mkdir minimap2_SAMs
for fastq in ERR2680377 ERR3363658_1 ERR3363660_1; do
	minimap2 -t $cores -ax splice -o minimap2_SAMs/$fastq.sam GRCm38.primary_assembly.genome.fa FASTQ/$fastq.fastq.gz
done


minimap2 -t $cores -ax map-ont -o minimap2_SAMs/$fastq.sam GRCm38.primary_assembly.genome.fa FASTQ/$fastq.fastq.gz





# Convert SAMs to BAMs
mkdir minimap2_sorted_BAMs
for file in ERR2680377 ERR3363658_1 ERR3363660_1; do
	samtools view -@ $cores -Su minimap2_SAMs/$file.sam | samtools sort -@ $cores -o minimap2_sorted_BAMs/$file.bam
done


# Section ____: Assembling transcriptome

# Download reference transcriptome, if absent
if [ ! -f "gencode.vM25.annotation.gtf" ]; then
	wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
	gzip -d gencode.vM25.annotation.gtf.gz
fi

# Assemble transcriptome using StringTie
mkdir Stringtie_gtf
for file in ERR2680377 ERR3363658_1 ERR3363660_1; do
	stringtie minimap2_sorted_BAMs/$file.bam -p $cores -o Stringtie_gtf/$file.gtf -G gencode.vM25.annotation.gtf
done

# Merge transcriptome
stringtie --merge -G gencode.vM25.annotation.gtf -o ont_merged.gtf Stringtie_gtf/ERR2680377.gtf Stringtie_gtf/ERR3363658_1.gtf Stringtie_gtf/ERR3363660_1.gtf 

