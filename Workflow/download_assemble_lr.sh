#usr/bin/bash

# Section __: Retrieve data
mkdir FASTQ_lr

# Parse tasic_fastq_curated.txt and download FASTQ files
while IFS=$'\t' read -r ACC LINK NAME;do
 	wget $LINK -P FASTQ_lr
done < lr_fastq.txt


if [ ! -f "gencode.vM25.annotation.gtf" ]; then
	wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz
	gzip -d GRCm38.primary_assembly.genome.fa.gz
fi

# Section __: Align reads
cores=`nproc`
mkdir minimap2_SAMs_lr
for fastq in ERR2680377 ERR3363658_1 ERR3363660_1; do
	minimap2 -t $cores -ax splice -o minimap2_SAMs_lr/$fastq.sam GRCm38.primary_assembly.genome.fa FASTQ_lr/$fastq.fastq.gz
done




# Convert SAMs to BAMs
mkdir minimap2_sorted_BAMs_lr
for file in ERR2680377 ERR3363658_1 ERR3363660_1; do
	samtools view -@ $cores -Su minimap2_SAMs_lr/$file.sam | samtools sort -@ $cores -o minimap2_sorted_BAMs_lr/$file.bam
done


# Section ____: Assemble transcriptome

# Download reference transcriptome, if absent
if [ ! -f "gencode.vM25.annotation.gtf" ]; then
	wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
	gzip -d gencode.vM25.annotation.gtf.gz
fi

# Assemble transcriptome using StringTie
mkdir Stringtie_gtf_lr
for file in ERR2680377 ERR3363658_1 ERR3363660_1; do
	stringtie minimap2_sorted_BAMs_lr/$file.bam -p $cores -o Stringtie_gtf_lr/$file.gtf -G gencode.vM25.annotation.gtf
done

# Merge transcriptome
stringtie --merge -G gencode.vM25.annotation.gtf -o lr_merged.gtf Stringtie_gtf_lr/ERR2680377.gtf Stringtie_gtf_lr/ERR3363658_1.gtf Stringtie_gtf_lr/ERR3363660_1.gtf 

