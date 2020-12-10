#usr/bin/bash

# pre-step: detect cpu cores
cores=`nproc`


# Section __: Transcriptome assembly

# Download reference transcriptome, if absent
if [ ! -f "gencode.vM25.annotation.gtf" ]; then
	wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
	gzip -d gencode.vM25.annotation.gtf.gz
fi


# Assemble transcriptome using StringTie
mkdir Stringtie_gtf
for file in "$@"; do
	stringtie Hisat2_sorted_BAMs/$file.bam -p $cores -o Stringtie_gtf/$file.gtf -G gencode.vM25.annotation.gtf
done

# Section __: Merge assembled transcriptome
declare -a array_exp
for filename in "$@"; do
	array_exp+="Stringtie_gtf/${filename}.gtf "
done
stringtie --merge -G gencode.vM25.annotation.gtf ${array_exp[*]}





