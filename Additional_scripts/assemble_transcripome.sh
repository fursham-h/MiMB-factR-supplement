#usr/bin/bash

# Create directories and detect CPU cores
mkdir Assembled_transcriptome
cores=`nproc`

# retrieve reference annotation
wget http://celltypes.brain-map.org/api/v2/well_known_file_download/502999254
unzip 502999254
rm 502999254

# Assemble transcripts using StringTie
for bam in merged_BAMs/*.bam;do
	file=${bam##*/}
	file=${file%.bam}
	stringtie $bam -p $cores -o Assembled_transcriptome/$file.gtf -G rsem_GRCm38.p3.gtf
done

# Merge each transcriptome
stringtie --merge -o assembled_merged_tasic.gtf Assembled_transcriptome/*.gtf -G rsem_GRCm38.p3.gtf

# Special script to rename seqlevels
