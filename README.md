# Publication resources
This repository holds key supplementary materials to our paper **"Functional annotation of custom transcriptomes"**. In this publication, we described the methods to assemble custom transcriptomes from different generations of RNA-sequencing (RNA-seq) experiements and annotate the functional output of novel isoforms using our custom R package, *factR*.

### Contents
This site contains the following sub-folders:

1. *"Workflow"* containing shell scripts and text files to build custom transcriptomes from bulk, single-cell or long-read RNA-seq experiments. 
2. *"Custom transcriptome"* containing pre-built custom transcriptomes from bulk, single-cell or long-read RNA-seq experiments.
3. *"Vignette"* containing a walkthrough of functional annotation analyses of custom transcriptomes assembled from bulk, single-cell or long-read RNA-seq experiments. (in preparation)

### Instructions

#### Download
Contents of this repository can be cloned to your local directory using Git or downloaded as a ZIP file by clicking the green button "Code" found at the top of this page.

#### Dependencies
Below is a list of dependencies required to execute the workflow.

1. [*factR*](https://github.com/fursham-h/factR)
2. [SRA toolkit](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit)
3. [HISAT2](http://daehwankimlab.github.io/hisat2/download/)
4. [Minimap2](https://github.com/lh3/minimap2)
5. [Samtools](www.htslib.org)
6. [StringTie2](https://ccb.jhu.edu/software/stringtie/#install)

#### Usage
The shell scripts in folder *"Workflow"* are to be executed in command-line terminal. We advise installing required dependencies and adding them to PATH first. To run the script named "download_assemble_bulk.sh" (for example), type the following in terminal:

```bash
cd Path/to/Workflow
bash download_assemble_bulk.sh
```

The compressed GTF files in folder "Custom transcriptome" serve as starting materials to demonstrate the usability of our custom R package *factR*. Please refer to the Methods section of our publication for a general workflow, or the PDFs in the *"Vignette"* folder for detailed workflow.

### Citation
If you use the materials from this repository, please cite _____

### Contact
You may email Fursham Hamid (fursham.hamid@kcl.ac.uk) or Eugene Makeyev (eugene.makeyev@kcl.ac.uk) for any enquiries with regards to the publication.
