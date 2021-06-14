#source("http://bioconductor.org/biocLite.R")
#biocLite(c("XML", "CNTools", "TCGA2STAT"))
#biocLite("TCGA2STAT")
library("TCGA2STAT")

### Create the main working directory
md = 'tmp'
dir.create(md, showWarnings=TRUE)

### Choose cancer, create necessary folders, and switch directory
# cancer
can <- "THCA"

# working directory
wd <- file.path(md, can)
dir.create(wd, showWarnings=TRUE)
setwd(wd)

# additional folders
dir.create('intermediate_files', showWarnings=TRUE)
dir.create('results', showWarnings=TRUE)

############################
### Raw data preparation ###
############################

### Read in the data
# note: always remove blood normal samples (if present)

rnaseq <- getTCGA(disease=can, data.type="RNASeq2", type="count", clinical=TRUE)
if (dim(rnaseq$dat.clean <- rnaseq$dat[, -grep("-10A-", colnames(rnaseq$dat))])[2] == 0) { rnaseq$dat.clean <- rnaseq$dat }

meth <- getTCGA(disease=can, data.type="Methylation", type="450K", clinical=TRUE)
if (dim(meth$dat.clean <- meth$dat[, -grep("-10A-", colnames(meth$dat))])[2] == 0) { meth$dat.clean <- meth$dat } else { meth$dat.clean <- meth$dat }

cnvsnp <- getTCGA(disease=can, data.type="CNV_SNP", clinical=TRUE , filter=c("X", "Y"))
if (dim(cnvsnp$dat.clean <- cnvsnp$dat[, -grep("-10A-", colnames(cnvsnp$dat))])[2] == 0) { cnvsnp$dat.clean <- cnvsnp$dat }


##########################
### Data preprocessing ###
##########################

### --- Methylation --- ###

# remove the sex chromosomes
mt.g <- as.data.frame(meth$cpgs, stringsAsFactors=FALSE)
mt <- as.data.frame(meth$dat.clean, stringsAsFactors=FALSE)
meth.g <- merge(mt.g,mt,by="row.names")
meth.g <- na.omit(meth.g)
meth.g <- meth.g[-grep("X|Y", meth.g$Chromosome),]
rownames(meth.g) <- meth.g$Row.names

# median the beta values per gene
meth.g <- meth.g[,c(2,5:length(colnames(meth.g)))]
m.unique <- aggregate(meth.g[2:length(colnames(meth.g))], list(meth.g$Gene_Symbol), median)
m.unique$Group.1 <- NULL
meth.samples <- colnames(m.unique)
m.unique$gene <- unique(meth.g$Gene_Symbol)

# recalculate beta into m-values
beta <- m.unique[,meth.samples]

M <- log2(beta/(1-beta))
meth.g = cbind(m.unique[,'gene'], M)
colnames(meth.g)[1] <- 'gene'
rownames(meth.g) <- meth.g$gene

meth$dat.clean <- meth.g

##########################
### Sample preparation ###
##########################

### --- Split by sample type (tumor, normal) --- ###
rnaseq.bySample <- SampleSplit(rnaseq$dat.clean) # 501 patients in $primary.tumor
meth.bySample <- SampleSplit(meth$dat.clean)     # 503 patients in $primary.tumor
cnvsnp.bySample <- SampleSplit(cnvsnp$dat.clean) # 499 patients in $primary.tumor

### --- Combine multiple OMICs profiles: rnaseq, meth, cvnsnp --- ###
# for primary tumor

d2.T <- OMICSBind(dat1 = rnaseq.bySample$primary.tumor, dat2 = meth.bySample$primary.tumor)
d3.T <- OMICSBind(dat1 = d2.T$merged.data, dat2 = cnvsnp.bySample$primary.tumor)

## -- Preprocess each data type separately

## rnaseq
rnaseqRows <- rownames(d3.T$merged.data)[grep("^d1.d1.", rownames(d3.T$merged.data))] # rnaseq
rnaseqT <- d3.T$merged.data[rnaseqRows,]
rnaseqRows <- gsub("d1.d1.", "", rnaseqRows)
rownames(rnaseqT) <- rnaseqRows

## meth
methRows <- rownames(d3.T$merged.data)[grep("^d1.d2.", rownames(d3.T$merged.data))] # meth
methT <- d3.T$merged.data[methRows,]
methRows <- gsub("d1.d2.", "", methRows)
rownames(methT) <- methRows

## cnvsnp
cnvsnpRows <- rownames(d3.T$merged.data)[grep("^d2.", rownames(d3.T$merged.data))] # cnvsnp
cnvsnpT <- d3.T$merged.data[cnvsnpRows,]
cnvsnpRows <- gsub("d2.", "", cnvsnpRows)
rownames(cnvsnpT) <- cnvsnpRows

## clinical data
clinicalData <- rnaseq$clinical # we focus on common only


## -- save data
dataTypes <- paste(c('rnaseq', 'meth', 'cnvsnp'),collapse='.')
save(rnaseqT, methT, cnvsnpT, clinicalData, file=paste0('intermediate_files/', dataTypes, ".T.comPatients.", can, ".rda"))

