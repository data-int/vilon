##########################
## -- load libraries -- ##

## differential expression
library("limma")
library("edgeR")

## parallel execution
library("foreach")
library("doMC")
ncores <- 8 # choose the number of available threads
registerDoMC(ncores) # change to your number of CPU cores

## set the main working directory
md = 'tmp'
## can <- Sys.getpid();
can <- "vilon.online";
wd <- file.path(md, can);
dir.create(wd,recursive=T);
setwd(wd)

###########################
## -- load functioons -- ##

## DE analysis, RNA-Seq
voomAvsB<-function(counts,adj.method='BY',eB=TRUE,proportion=0.01) {
    
    m <- counts
    m[is.na(m)] <- 0
    
    samples<-sub('[0-9]+','',colnames(m));
    design<-data.frame(As=ifelse(samples=='A',1,0), Bs=ifelse(samples=='B',1,0));
    rownames(design)<-colnames(m);
    design0<-design;
    use<-apply(design0!=0,1,any);
    design<-design0[use,];
    m<-m[,use];
    cm<-makeContrasts(AvsB=As-Bs, # disease-control
                      levels=design);
    m2 <- m

    dge <- DGEList(counts=m2)
    dge <- calcNormFactors(dge)
    v <- voom(dge,design,plot=FALSE)
    f <- lmFit(v,design, method='robust',maxit=99999); # increase until convergence
    if (!any(names(f)=='Amean')) {
        f$Amean<-rowMeans(m);
    }
    cf<-contrasts.fit(f,cm);
    
    if(eB) {    
        fe <- eBayes(cf, proportion); # prior belief: what proportion of genes differential; default is 0.01
        tf.temp <- topTable(fe,number=Inf,adjust.method=adj.method,sort.by='none')
    } else {
        tf.temp <- toptable(cf,number=Inf,adjust.method=adj.method,sort.by='none', genelist=cf$genes)
    }

    len <- dim(m2)[1]
    tf <- data.frame(ID=rownames(m2),
                     logFC=    rep(NA,len),
                     AveExpr=  rep(NA,len),
                     t=        rep(NA,len),
                     P.Value=  rep(NA,len),
                     adj.P.Val=rep(NA,len),
                     B=        rep(NA,len),
                     P=        rep(NA,len));
    tf$logFC=    tf.temp$logFC;
    tf$AveExpr=  cf$Amean;
    tf$t=        tf.temp$t;
    tf$P.Value=  tf.temp$P.Value;
    tf$adj.P.Val=tf.temp$adj.P.Val;
    tf$B=        tf.temp$B;
    tf$P=        exp(tf.temp$B)/(1+exp(tf.temp$B)) # P is the posterior probability
    
    if (!all(tf.temp$ID==rownames(cf))) stop('Internal error.\n');
    return(tf);
}

### DE analysis, other data types
limmaAvsB<-function(log2_counts,adj.method='BY',eB=TRUE,proportion=0.01) {
  
  m <- log2_counts
  m[is.na(m)] <- 0

  samples<-sub('[0-9]+','',colnames(m));
  design<-data.frame(As=ifelse(samples=='A',1,0),
                     Bs=ifelse(samples=='B',1,0));
  rownames(design)<-colnames(m);
  design0<-design;
  use<-apply(design0!=0,1,any);
  design<-design0[use,];

  m<-m[,use];

  cm<-makeContrasts(AvsB=As-Bs,
                    levels=design);
  
  m2 <- m
  
  f<-lmFit(m2,design,
           method='robust',maxit=99999  # increase until convergence (no warnings)
           );
  
  if (!any(names(f)=='Amean')) {
    f$Amean<-rowMeans(m);
  }
  cf<-contrasts.fit(f,cm);
  
  if(eB)
    {    
      
      fe<-eBayes(cf,
                 proportion  # prior belief: what proportion of genes differential; default is 0.01
                 );
      tf.temp<-topTable(fe,number=Inf,adjust.method=adj.method,sort.by='none')
      
    }else{
      tf.temp<-toptable(cf,number=Inf,adjust.method=adj.method,sort.by='none', genelist=cf$genes)
    }
  
  len <- dim(m2)[1]
  tf <- data.frame(ID=rownames(m2),
                   logFC=    rep(NA,len),
                   AveExpr=  rep(NA,len),
                   t=        rep(NA,len),
                   P.Value=  rep(NA,len),
                   adj.P.Val=rep(NA,len),
                   B=        rep(NA,len),
                   P=        rep(NA,len));
  tf$logFC=    tf.temp$logFC;
  tf$AveExpr=  cf$Amean;
  tf$t=        tf.temp$t;
  tf$P.Value=  tf.temp$P.Value;
  tf$adj.P.Val=tf.temp$adj.P.Val;
  tf$B=        tf.temp$B;
  tf$P=        exp(tf.temp$B)/(1+exp(tf.temp$B)) # P is the posterior probability
  
  if (!all(tf.temp$ID==rownames(cf))) stop('Internal error.\n');
  return(tf);
}

#######################
## -- main script -- ##

## collect info
datDir <- 'intermediate_files/'

files <- dir(datDir,ignore.case=T,pattern="*.csv");
fileTypes <- sub(".*[.]","",sub("[.]csv","",files,ignore.case=T));
isClinical <- grepl("clinical",fileTypes,ignore.case=T);
dataTypes <- fileTypes[!isClinical];
nData <- 1:length(dataTypes) # location of each new data type data
dataTypesVoom <- ifelse(grepl("rna-?seq",dataTypes,ignore.case=T),1,0);
names(dataTypesVoom) <- dataTypes;
## which data should be Voom-TMM normalized

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

cfgTable <- as.data.frame(t(read.table("config.txt",sep="\t",row.names=1,
                                       col.names=c("NULL","value"))));

cachefile=paste0(datDir,
                 paste(dataTypes, collapse='.'),
                 ".T.comPatients.",can,".rda");
dataenv <- as.environment(Data);
save(list=ls(dataenv),file=cachefile,envir=dataenv);
rm(dataenv);


## create output directory
outDir  <- file.path(datDir, 'normalized/')
dir.create(outDir, showWarnings=TRUE)

## loop through all data types
for (n in nData) {
    
    i <- which(n==nData)
    type <- dataTypes[i]
    useVoom <- dataTypesVoom[i]

    datA.All <- Data[[n]]   # cancer data
    datB.All <- Data[[n]]   # control data [CANCER]
    
    ## remove duplicated names (unmappable to PWs); can happen for some data types
    datA.All <- datA.All[!duplicated(rownames(datA.All)), ]
    datB.All <- datB.All[!duplicated(rownames(datB.All)), ]

    ## loop through all the cancer donors each
    foreach (nam=colnames(datA.All)) %dopar% {
        
        Donors1 <- nam                                            # 1 cancer sample each, vs
        Donors2 <- colnames(datB.All)[colnames(datB.All) != nam]  # all controls [other CANCERs]

        ## controls data
        datB <- datB.All[, Donors2]
        a <- "B"; b <- 1:length(colnames(datB)); 
        colnames(datB) <- sprintf(paste0(a, "%d"), b)
        
        ## cancer data
        datA <- data.frame(A1=data.frame(datA.All[,Donors1]))
        colnames(datA) <- "A"
        
        ## merge the data
        merged <- cbind(datA, datB)

        ## perform de analysis
        if (useVoom == 1) {
            ## CHECK: input is unloged data!?
            ## set a prior for the calculations
#            testPrior <- voomAvsB(merged)
#            priorTested <- convest(testPrior$P.Val)
            priorTested <-  0.1 # use default, auto sometimes unexpectedly high
            ## calculate statistics
            stats <- voomAvsB(merged, proportion=priorTested)
            method <- "limmaTMMVoom"
        } else {
            ## CHECK: input is loged data!?
            merged <- scale(merged) # normalize
            method <- "scaleLimma"
            ## set a prior for the calculations
#            testPrior <- limmaAvsB(merged)
#            priorTested <- convest(testPrior$P.Val)
            priorTested <-  0.1 # use default, auto sometimes unexpectedly high
            ## calculate statistics
            stats <- limmaAvsB(merged, proportion=priorTested)
        }
        
        ## merge relevant data (if necessary)
        final <- stats

        ## rank the data for the next step
        # 'NA' in P come from high Bs - meaning that P would be very high: exp(B)/(1+exp(B))
        final$P[is.na(final$P)] <- 1
        final.rank <- data.frame(ID=final$ID, P=final$P)
        final.rank <- final.rank[order(final.rank$P, decreasing=TRUE),]
        write.table(final.rank, file = paste0(outDir, can, ".", type, ".", method, ".1t_nt.", nam, ".prior01.P.rnk"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    }
}

