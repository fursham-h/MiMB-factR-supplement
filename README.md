# Publication resources
This repository holds key supplementary materials to our paper **"Functional annotation of custom transcriptomes"**. In this publication, we described the methods to assemble custom transcriptomes from different RNA-sequencing (RNA-seq) experiments (bulk, single-cell and long-read) and annotate the functional output of novel isoforms using our custom R package, *factR*.

### Contents
This site contains the following sub-folders:

1. *"Workflow"* containing shell scripts and text files to build custom transcriptomes from bulk, single-cell or long-read RNA-seq experiments. 
2. *"Custom transcriptome"* containing pre-built custom transcriptomes from bulk, single-cell or long-read RNA-seq experiments.

### Instructions

#### Download
Contents of this repository can be cloned to your local directory using Git or downloaded as a ZIP file by clicking the green button "Code" found at the top of this page.

#### Dependencies
Below is a list of dependencies required to execute the workflow.

1. [*factR*](https://fursham-h.github.io/factR)
2. [SRA toolkit](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit)
3. [HISAT2](http://daehwankimlab.github.io/hisat2/download/)
4. [Minimap2](https://github.com/lh3/minimap2)
5. [Samtools](www.htslib.org)
6. [StringTie2](https://ccb.jhu.edu/software/stringtie/#install)

#### Usage
The shell scripts in folder *"Workflow"* contain code to download and assemble transcriptomes from three RNA-seq datasets (bulk: download_assemble_bulk.sh; single-cell: download_assemble_sc.sh; long-read: download_assemble_lr.sh). Dependencies [2-6] listed [above](#dependencies) are to be installed on your computer and added to PATH. To execute any one these scripts, open terminal and type the following :

```bash
cd Path/to/Workflow
bash ./download_assemble_bulk.sh 
#or
bash ./download_assemble_sc.sh 
#or
bash ./download_assemble_lr.sh
```

The above scripts will output custom transcriptomes as compressed GTF files, identical to those found in folder "Custom transcriptome". These GTFs can be further processed to annotate functions of newly-discovered transcripts using tools from our custom R package *factR*. To perform this, refer to "Methods" section of our publication or to *factR*'s [vignette walkthrough](https://fursham-h.github.io/factR/articles/factR.html)

### Citing us
If you use the materials from this repository, please cite _____

### Contact
You may email Fursham Hamid (fursham.hamid@kcl.ac.uk) or Eugene Makeyev (eugene.makeyev@kcl.ac.uk) for any enquiries with regards to the publication.
