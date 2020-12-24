#usr/bin/bash

# Section 3.1.2: Retrieve data
mkdir FASTQ_bulk
fastq-dump -O FASTQ_bulk --split-files --gzip SRR606187 SRR606188 SRR606189 SRR606190 

# Section 3.1.3: Download index
if [ ! -d "mm10" ]; then
	wget https://genome-idx.s3.amazonaws.com/hisat/mm10_genome.tar.gz
	tar -xf "mm10_genome.tar.gz"
	rm mm10_genome.tar.gz
fi

# Section 3.1.4: Align reads
cores=`nproc`
mkdir Hisat2_SAMs_bulk
for fastq in SRR606187 SRR606188 SRR606189 SRR606190;do
	hisat2 -p $cores --dta -x mm10/genome -1 FASTQ_bulk/"$fastq"_1.fastq.gz -2 FASTQ_bulk/"$fastq"_2.fastq.gz -S Hisat2_SAMs_bulk/$fastq.sam
done

# Section 3.1.5: Convert SAM to BAM
mkdir Hisat2_sorted_BAMs_bulk
for file in SRR606187 SRR606188 SRR606189 SRR606190; do
	samtools view -@ $cores -Su Hisat2_SAMs_bulk/$file.sam | samtools sort -@ $cores -o Hisat2_sorted_BAMs_bulk/$file.bam
done

# Section 3.1.6: Download reference annotation
if [ ! -f "gencode.vM25.annotation.gtf" ]; then
	wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
	gzip -d gencode.vM25.annotation.gtf.gz
fi

# Section 3.1.7: Assemble transcriptome
mkdir Stringtie_gtf_bulk
for file in SRR606187 SRR606188 SRR606189 SRR606190; do
	stringtie Hisat2_sorted_BAMs_bulk/$file.bam -p $cores -o Stringtie_gtf_bulk/$file.gtf -G gencode.vM25.annotation.gtf
done

# Section 3.1.8: Merge transcriptome and compress
stringtie --merge -G gencode.vM25.annotation.gtf Stringtie_gtf_bulk/SRR606187.gtf Stringtie_gtf_bulk/SRR606188.gtf Stringtie_gtf_bulk/SRR606189.gtf Stringtie_gtf_bulk/SRR606190.gtf | gzip > bulk_merged.gtf.gz
