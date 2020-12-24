#usr/bin/bash

# Section 3.3.3: Retrieve data
mkdir FASTQ_lr
while IFS=$'\t' read -r ACC LINK NAME;do
 	wget $LINK -P FASTQ_lr
done < lr_fastq.txt

# Section 3.3.4: Download reference genome
if [ ! -f "GRCm38.primary_assembly.genome.fa" ]; then
	wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz
	gzip -d GRCm38.primary_assembly.genome.fa.gz
fi

# Section 3.3.5: Align reads
cores=`nproc`
mkdir minimap2_SAMs_lr
for fastq in ERR2680377 ERR3363658_1 ERR3363660_1; do
	minimap2 -t $cores -ax splice -o minimap2_SAMs_lr/$fastq.sam GRCm38.primary_assembly.genome.fa FASTQ_lr/$fastq.fastq.gz
done

# Section 3.3.6: Convert SAM to BAM
mkdir minimap2_sorted_BAMs_lr
for file in ERR2680377 ERR3363658_1 ERR3363660_1; do
	samtools view -@ $cores -Su minimap2_SAMs_lr/$file.sam | samtools sort -@ $cores -o minimap2_sorted_BAMs_lr/$file.bam
done

# Section 3.3.7: Download reference annotation
if [ ! -f "gencode.vM25.annotation.gtf" ]; then
	wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
	gzip -d gencode.vM25.annotation.gtf.gz
fi

# Section 3.3.8: Assemble transcriptome
mkdir Stringtie_gtf_lr
for file in ERR2680377 ERR3363658_1 ERR3363660_1; do
	stringtie minimap2_sorted_BAMs_lr/$file.bam -p $cores -o Stringtie_gtf_lr/$file.gtf -G gencode.vM25.annotation.gtf
done

# Section 3.3.9: Merge transcriptome and compress output
stringtie --merge -G gencode.vM25.annotation.gtf Stringtie_gtf_lr/ERR2680377.gtf Stringtie_gtf_lr/ERR3363658_1.gtf Stringtie_gtf_lr/ERR3363660_1.gtf | gzip > lr_merged.gtf.gz 
