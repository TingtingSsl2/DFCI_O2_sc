# title: making geneID and geneSymbol table from Mouse gtf file 
# author: Tingting Zhao
# data: 09/02/2022
# usage: Rscript makeGeneAnnoTable.R 

args <- commandArgs(trailingOnly = TRUE)
refdir = args[1]
geneName_file = args[2]

# read in gene name file
geneAnno <- read.delim(paste0(refdir, geneName_file), header = T, col.names=c("geneID", "geneSymbol"))

# replace gene ID to genes without gene symbols
# for (i in 1:length(geneAnno$geneSymbol)){
#   if (geneAnno$geneSymbol[i]==" gene_source FlyBase"){
#     geneAnno$geneSymbol[i] <- geneAnno$geneID[i]
#   }
# }

# save gene names as a table
write.csv(geneAnno, paste0(refdir, "geneAnnotationTable.csv"))
