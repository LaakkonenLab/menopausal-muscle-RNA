# This code is taken from the mixOmics vignette

library(mixOmics)
library(tidyverse)

browseVignettes("mixOmics")


# load data

mir <- read.csv("miRNA.txt",row.names = 1, header=TRUE,sep="\t" )  
mRNA <- read.csv("mRNA.txt",row.names = 1, header=TRUE,sep="\t" )  
lnc <- read.csv("lncRNA.txt",row.names = 1, header=TRUE,sep="\t" )  

mirDEG <- read.csv("miR2-DEG.txt", header=FALSE)
mRNADEG <- read.csv("DEG-mRNA-peri2.txt", header=TRUE,sep="\t")
lncDEG <- read.csv("lnc2-DEG.txt", header=FALSE)



X <- list(mRNA = mRNA,
          mir = mir, 
          lncRNA = lnc)

Y <- c(rep("Early",7),rep("Late",17))


# calculate block PLS-DA

MyResult.diablo <- block.plsda(X, Y)
plotVar(MyResult.diablo)


# calculate sparse block PLS-DA
list.keepX <- list(mRNA = c(10, 10), mir = c(10,10), lncRNA = c(10, 10))

MyResult.diablo <- block.splsda(X, Y, keepX=list.keepX)


plotIndiv(MyResult.diablo)
plotVar(MyResult.diablo)


#relevance network, which is built on the similarity matrix

network(MyResult.diablo, blocks = c(1,2,3),
        color.node = c('orange', 'yellow', 'blue'), 
        cutoff = 0.7, save = 'jpeg', name.save = 'Network1')




