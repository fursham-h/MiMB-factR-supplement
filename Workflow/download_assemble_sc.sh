#usr/bin/bash

# Section _____: Retrieving data


## Parse tasic_fastq_curated.txt and download FASTQ files
mkdir FASTQ_sc
while IFS=$'\t' read -r NAME CLASS LINK CLASSNUM;do
 	wget $LINK -P FASTQ_sc
 	tar -xvf "FASTQ_sc/$NAME.fastq.tar" -C FASTQ_sc
 	rm FASTQ_sc/$NAME.fastq.tar
done < sc_fastq_curated.txt

## Merge cell from same clusters

### Create arrays to store filenames and groupnames
declare -A fastqArray
declare -A fastqNames

while IFS=$'\t' read -r NAME CLASS LINK CLASSNUM;do
	fastqArray["$CLASSNUM.1"]="${fastqArray["$CLASSNUM.1"]} FASTQ_sc/${NAME}_R1.fastq.gz"
	fastqArray["$CLASSNUM.2"]="${fastqArray["$CLASSNUM.2"]} FASTQ_sc/${NAME}_R2.fastq.gz"
	fastqNames["$CLASSNUM.1"]=${CLASS}_R1
	fastqNames["$CLASSNUM.2"]=${CLASS}_R2
done < sc_fastq_curated.txt

### Loop groups and concatenate FASTQ by direction
for (( n=1; n<=4; n++ ));do
	cat ${fastqArray["$n.1"]} > FASTQ_sc/${fastqNames["$n.1"]}.fastq.gz
	cat ${fastqArray["$n.2"]} > FASTQ_sc/${fastqNames["$n.2"]}.fastq.gz
done



# Section _____: Aligning reads

# Download mm10 index if absent
if [ ! -d "mm10" ]; then
	wget https://genome-idx.s3.amazonaws.com/hisat/mm10_genome.tar.gz
	tar -xf "mm10_genome.tar.gz"
	rm mm10_genome.tar.gz
fi

# Align reads to mm10 genome
cores=`nproc`
mkdir Hisat2_SAMs_sc
for fastq in Glutamatergic GABAergic Endothelial Astrocyte;do
	hisat2 -p $cores --dta -x mm10/genome -1 FASTQ_sc/"$fastq"_R1.fastq.gz -2 FASTQ_sc/"$fastq"_R2.fastq.gz -S Hisat2_SAMs_sc/$fastq.sam
done

# Convert SAMs to BAMs
mkdir Hisat2_sorted_BAMs_sc
for file in Glutamatergic GABAergic Endothelial Astrocyte; do
	samtools view -@ $cores -Su Hisat2_SAMs_sc/$file.sam | samtools sort -@ $cores -o Hisat2_sorted_BAMs_sc/$file.bam
done


# Section ____: Assembling transcriptome

# Download reference transcriptome, if absent
if [ ! -f "gencode.vM25.annotation.gtf" ]; then
	wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
	gzip -d gencode.vM25.annotation.gtf.gz
fi

# Assemble transcriptome using StringTie
mkdir Stringtie_gtf_sc
for file in Glutamatergic GABAergic Endothelial Astrocyte; do
	stringtie Hisat2_sorted_BAMs_sc/$file.bam -p $cores -o Stringtie_gtf_sc/$file.gtf -G gencode.vM25.annotation.gtf
done

# Merge transcriptome
stringtie --merge -G gencode.vM25.annotation.gtf -o sc_merged.gtf Stringtie_gtf_sc/Glutamatergic.gtf Stringtie_gtf_sc/GABAergic.gtf Stringtie_gtf_sc/Endothelial.gtf Stringtie_gtf_sc/Astrocyte.gtf

