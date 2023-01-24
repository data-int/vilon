##########################
## -- load libraries -- ##

## parallel execution
library("foreach")
library("doMC")

## set the main working directory
md = 'tmp'
can <- "vilon.online"
wd <- file.path(md, can)
dir.create(wd,recursive=T)
setwd(wd)

cfgTable <- as.data.frame(t(read.table("config.txt",sep="\t",row.names=1,
                                       col.names=c("NULL","value"))));

ncores <- 8 # choose the number of available threads
registerDoMC(ncores) # change to your number of CPU cores


#######################
## -- main script -- ##

## collect info
datDir <- 'intermediate_files/'

files <- dir(datDir,ignore.case=T,pattern="*.csv");
fileTypes <- sub(".*[.]","",sub("[.]csv","",files,ignore.case=T));
isClinical <- grepl("clinical",fileTypes,ignore.case=T);
## which data should be Voom-TMM normalized
dataTypes <- fileTypes[!isClinical];
nData <- 1:length(dataTypes)
dataTypesVoom <- ifelse(grepl("rna-?seq",dataTypes,ignore.case=T),1,0);
names(dataTypesVoom) <- dataTypes;

Data<-list();
for (i in seq(length(fileTypes))) {
    if (isClinical[i]) {
        cat("Reading clinical survival data from",files[i],"\n");
        Data[["clinicalData"]] <- try(read.csv(file.path(datDir,files[i]),
                                               head=T));
    } else {
        cat("Reading",fileTypes[i],"data from",files[i],
            ifelse(dataTypesVoom[fileTypes[i]]>0,
                   "for Voom/TMM normalization\n","\n"));
        Data[[paste0(fileTypes[i],"T")]] <- try(read.csv(file.path(datDir,
                                                                   files[i]),
                                                         head=T));
    }
}
if (any(sapply(Data,function(le){inherits(le, "try-error")}))) {
    stop("Error reading CSV input files - aborting.");
}


cachefile <- paste0(datDir,
                    paste(dataTypes, collapse='.'),
                    ".T.comPatients.",can,".rda");
metafile <- paste0(datDir,"meta",".T.comPatients.",can,".rda");
dataenv <- as.environment(Data);
cat("Saving",cachefile,"\n");
save(list=ls(dataenv),
     file=cachefile,envir=dataenv);
rm(dataenv);
cat("Saving",metafile,"\n");
save(list=c("nData","dataTypes","dataTypesVoom","cfgTable"),
     file=metafile);
cat("Step 1 - read data - done.\n");
