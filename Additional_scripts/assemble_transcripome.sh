#usr/bin/bash


# Assemble transcripts using StringTie
mkdir Assembled_transcriptome
for (( n=1; n<=4; n++ ));do
	stringtie merged_BAMs/${bamNames[$n]}.bam \
	 -p $cores -o Assembled_transcriptome/${bamNames[$n]}.gtf -G rsem_GRCm38.p3.gtf
done
