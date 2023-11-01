# tximport code is from https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
# deseq code is from from https://bioc.ism.ac.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/beginner.pdf


library(tximport)
library(tidyverse)


# run tximport

mRNAdata <- read.csv("mRNA-Early.txt",header=TRUE,sep=" ")

tx2gene <- cbind(mRNAdata$transcript_id,mRNAdata$gene_id)
colnames(tx2gene) <- c("TXNAME","GENEID")
tx2gene <- as_tibble(tx2gene)


files <- c("PERI_656.txt","POST_656.txt","PERI_1208.txt","POST_1208.txt","PERI_1174.txt","POST_1174.txt",
           "PERI_2453.txt","POST_2453.txt","PERI_3159.txt", "POST_3159.txt", "PERI_3408.txt", "POST_3408.txt",
           "PERI_4180.txt", "POST_4180.txt")

names(files) <- c("PERI_656","POST_656","PERI_1208","POST_1208","PERI_1174","POST_1174",
                  "PERI_2453","POST_2453","PERI_3159", "POST_3159", "PERI_3408", "POST_3408", "PERI_4180", "POST_4180")

txi <- tximport(files, type = "none",  tx2gene=tx2gene, txIdCol = "id", abundanceCol="abundance", countsCol="counts", lengthCol="length", importer=read.csv)


library(DESeq2)

# run DERSeq2
# code is from 

# set up DESeq2 dataset
phenodata <- read.csv("PhenoData.txt",sep="\t",header=TRUE) # sample IDs
groups <- as.character(phenodata$type)  # time point, PERI vs POST menopause
exp_factor <- as.character(phenodata$sample)  # person column 
group_levels <- levels(as.factor(groups))
design <- data.frame(condition=as.factor(groups), exp_factor=exp_factor)
rownames(design) <- colnames(txi$counts)
dds <- DESeqDataSetFromTximport(txi, design, ~ exp_factor + condition)

# DESeq calculate 
dds <- DESeq(dds)  

# display results 

res <- results(dds)


