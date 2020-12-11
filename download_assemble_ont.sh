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
	minimap2 -t $cores -ax map-ont -o minimap_SAMs/$fastq.sam GRCm38.primary_assembly.genome.fa FASTQ/$fastq.fastq.gz
done


# Convert SAMs to BAMs
mkdir minimap2_sorted_BAMs
for file in ERR2680377 ERR3363658_1 ERR3363660_1; do
	samtools view -@ $cores -Su minimap2_SAMs/$file.sam | samtools sort -@ $cores -o minimap2_sorted_BAMs/$file.bam
done


