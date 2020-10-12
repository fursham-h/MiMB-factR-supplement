#usr/bin/bash

mkdir .tmp
cores=`nproc`

# parse tasic_bam_curated.txt and download BAM files to processed_BAMs directory
 while IFS=$'\t' read -r NAME CLASS LINK CLASSNUM;do
 	wget $LINK -P .tmp
 done < tasic_bam_curated.txt

# untar all files 
cd .tmp
for files in *.bam.tar;do 
	tar -xvf "$files" 
done

# merge BAM files from samples of the same class
cd ..
mkdir merged_BAMs

declare -A bamArray
declare -A bamNames
while IFS=$'\t' read -r NAME CLASS LINK CLASSNUM;do
	bamArray[$CLASSNUM]="${bamArray[$CLASSNUM]} .tmp/$NAME.bam"
	bamNames[$CLASSNUM]=$CLASS
done < tasic_bam_curated.txt

for (( n=1; n<=4; n++ ));do
	samtools merge -@ $cores merged_BAMs/${bamNames[$n]}.bam ${bamArray[$n]}
	#samtools index -@ 30 merged_BAMs/${bamNames[$n]}.bam merged_BAMs/${fastqNames[$n]}.bam.bai
done

rm -r .tmp

