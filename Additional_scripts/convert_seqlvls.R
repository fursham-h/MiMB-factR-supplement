# Check dependencies
if(!"dplyr" %in% installed.packages()[,"Package"]){
  install.packages("dplyr")
}
library(dplyr, quietly = T)
# Read files
input.gtf <- read.delim("assembled_merged_tasic.gtf", header = F)
ncbi_accession <- read.delim("RefSeq_genomic_accessions.txt")

# Convert seqlevels
suppressMessages(input.gtf %>% as.data.frame() %>% 
  dplyr::left_join(ncbi_accession) %>% 
  dplyr::select(Vnew, V2:V9) %>% 
  dplyr::filter(!is.na(Vnew)) %>% 
  write.table("assembled_merged_tasic.gtf", sep = "\t", quote = F, row.names = F, col.names = F))
  