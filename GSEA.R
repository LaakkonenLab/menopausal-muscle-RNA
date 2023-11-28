# This code is adapted from https://stephenturner.github.io/deseq-to-fgsea/
# and from https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html

library(tidyverse)
library(DESeq2)
library(org.Hs.eg.db)
library(msigdbr)
library(fgsea) 
library(data.table)

Annotation <- read.csv("annotations.txt",header=TRUE,sep="\t")

Counts <- read.csv("mRNA-transcript-counts-Early.txt",header=TRUE,sep=" ")


# read protein coding transcript list

PClist <- read.csv("ProteinCoding_Early.txt",header=TRUE,sep="\t")


# Choose only Protein Coding transcripts

PCcounts <- subset(Counts, rownames(Counts)%in%PClist$Transcript_ID)


# set up DESeq2 dataset

phenodata <- read.csv("PhenoDataEarly.txt",sep="\t",header=TRUE) # sample IDs
groups <- as.character(phenodata$type)  # perimenopause vs postmenopause
exp_factor <- as.character(phenodata$sample)  # person ID column 
group_levels <- levels(as.factor(groups))
design <- data.frame(condition=as.factor(groups), exp_factor=exp_factor)
rownames(design) <- colnames(countdata)
dds <- DESeqDataSetFromMatrix(countData=PCcounts, colData=design, design = ~ exp_factor + condition)


# DESeq calculate 
dds <- DESeq(dds)  


# display results 

res <- results(dds)
res2 <- results(dds, tidy=TRUE)
resOrdered <- res2[order(res2$padj),]  # the data has column names "row", "stat", and baseMean that we will use

annotated_list <- inner_join(resOrdered, Annotation, by=c("row"="transcript_id")) 


# make ranked list for GSEA

res3 <- annotated_list %>% 
  dplyr::select(gene_id, stat, baseMean) %>% 
  na.omit() %>% 
  group_by(gene_id) %>%
  arrange( -baseMean) %>%   ## choose transcript with maximum basemean
  slice(1) %>%
  summarize(stat=stat)


ranks <- deframe(res3)

########

# FGSEA ANALYSIS

########

# select database

# GO BP
My_sets = msigdbr(species = "human", category = "C5", subcategory = "GO:BP")

# REACTOME
My_sets = msigdbr(species = "human", category = "C2", subcategory = "CP:REACTOME")

# KEGG
My_sets = msigdbr(species = "human", category = "C2", subcategory = "CP:KEGG")


# Use the gene sets data frame for fgsea

msigdbr_list = split(x = My_sets$ensembl_gene, f = My_sets$gs_name)

fgseaRes <- fgsea(pathways = msigdbr_list, stats = ranks, minSize=15, maxSize=500)


# annotate gene IDs to gene symbols

fgseaRes[, leadingEdge := mapIdsList(
  x=org.Hs.eg.db, 
  keys=leadingEdge,
  keytype="ENSEMBL", 
  column="SYMBOL")]


# save results

fwrite(fgseaResMain, file="EarlyMT-GSEA-BP.txt", sep="\t", sep2=c("", " ", ""))



