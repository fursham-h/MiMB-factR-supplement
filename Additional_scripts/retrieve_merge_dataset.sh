#usr/bin/bash

# Create directories and detect CPU cores
mkdir .tmp
mkdir .tmpbams
mkdir merged_BAMs
cores=`nproc`

# Create Arrays
declare -A bamArray
declare -A bamNames

# Parse tasic_bam_curated.txt and download BAM files
while IFS=$'\t' read -r NAME CLASS LINK CLASSNUM;do
 	wget $LINK -P .tmp
 	tar -xvf ".tmp/$NAME.bam.tar" -C .tmp
 	rm .tmp/$NAME.bam.tar
 	bamArray[$CLASSNUM]="${bamArray[$CLASSNUM]} .tmp/$NAME.bam"
	bamNames[$CLASSNUM]=$CLASS
done < tasic_bam_curated.txt

# Merge bam files and clean up tmp folder
for (( n=1; n<=4; n++ ));do
	samtools merge -@ $cores merged_BAMs/${bamNames[$n]}.bam ${bamArray[$n]}
done
rm -r .tmp

# # Correct header and clean up tmpbams folder
# for bam in .tmpbams/*.bam;do
# 	file=${bam##*/}
# 	samtools reheader header_corrected.sam $bam > merged_BAMs/$file
# done
# rm -r .tmpbams
