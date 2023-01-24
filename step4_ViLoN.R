
### load essential packages ###
library(igraph)                # for graph preprocessing, vi, etc.
library(blockmodels)           # for SBM methods (instead of Wmixnet; same author)

## loadRData <- function(fileName){
##     ##loads an RData file, and returns it
##     env <- ls()
##     load(fileName)
##     data <- ls()[!(ls() %in% env)]
##}


################################
### --- start --- CONFIG --- ###
################################

### Clinical information
can <- "vilon.online"
datDir <- paste0('tmp/',can,'/intermediate_files/')

files <- dir(datDir,ignore.case=T,pattern="*.csv");
fileTypes <- sub(".*[.]","",sub("[.]csv","",files,ignore.case=T));
isClinical <- grepl("clinical",fileTypes,ignore.case=T);
dataTypes <- fileTypes[!isClinical];
nData <- 1:length(dataTypes) # location of each new data type data

file=paste0(datDir, paste(dataTypes, collapse='.'),".T.comPatients.",can,".rda")
load(file,Data <- new.env())
Data <- as.list(Data)

if ("clinicalData" %in% names(Data)){
    ## patients who were alive at time of last check have "daystolastfollowup", dead patients have "daystodeath", "vitalstatus" tells us which is which so we can combine both times
    clinical <- Data$clinicalData
    clinical[is.na(clinical[,"daystolastfollowup"]),"daystolastfollowup"] = clinical[is.na(clinical[,"daystolastfollowup"]),"daystodeath"]
    clinical[clinical[,"gender"]=="male","gender"] = clinical[clinical[,"gender"]=="male","gender"]=1
    clinical[clinical[,"gender"]=="female","gender"] = clinical[clinical[,"gender"]=="female","gender"]=2
    clinical <- data.frame(sample=rownames(clinical),
                           time=as.numeric(clinical[,"daystolastfollowup"]),
                           censor=as.numeric(clinical[,"vitalstatus"]),
                           status=as.numeric(clinical[,"vitalstatus"])+1, # for library(survival) this means: 1=censored, 2=dead
                           sex=as.numeric(clinical[,"gender"]))
    clinical <- clinical[!is.na(clinical$time),]
    clinical <- clinical[!is.na(clinical$censor),]
    clinical <- clinical[!clinical$time <= 0,]

    rownames(clinical) <- clinical$sample
  ##  samClinical <- rownames(clinical)
}
samClinical <- rownames(t(Data[[paste0(dataTypes[1],"T")]]))


## General params, bipartite graphs ##
netDir <- paste0('tmp/', can, '/intermediate_files/bpgs/')
netSuffix <- ".prior01.bpg"            # bipartite graphs of KEGG-GO-KEGG connected via genes
nodeNames <- read.table("data/c2.cp.kegg.v6.0.symbols.ids", stringsAsFactors=FALSE)[,1]

## Prepare data ##
# list all the input files (per patient, per data type)
ver <- ".1t_nt" # 1 tumor patient vs N tumor patients
allData <- list()
for (d in dataTypes) {

    i <- which(dataTypes == d)
    
    files <- list.files(paste0(netDir, "/"), pattern = paste0(can, ".", d, ".*", ver, ".*", netSuffix)) # for each data type the same patients were used
    files <- sapply(strsplit(basename(files),netSuffix), function(x) paste(x[1:(length(x)-1)]))  
    complete <- c()
    for (nam in samClinical) {
        complete <- c(complete, files[grepl(nam, files)])
    }
    files <- complete
    ##patients <- 5
    ##sel <- sample(length(files), patients)  # in case we want to make a smaller network (e.g., for testing)
    sel <- 1:length(files) 
    allData[[i]] <- files[sel] # data types need to be matched via patients

}
dataTypes <- dataTypes
names(allData) <- dataTypes

## Set thresholds [OPTIONAL] ##
# thr for which edges in bpgs to keep
thrs <- c(rep(list(0), length(dataTypes))) # a thr per each data type

## BM params ##
threads <- 1 # more threads do not seem to make the method run faster

## Parallel execution params ##
library(foreach)
library(doMC)                   # backend for running 'foreach' in parallel
cores <- 8/threads              # change to number of available CPU cores
registerDoMC(cores)             # register the doMC parallel backend

##############################
### --- end --- CONFIG --- ###
##############################

#####################################
### --- start --- MAIN METHOD --- ###
#####################################

### --- Clustering patient networks internally --- ###

## initialize an object with all patients' clusterings
ppClust <- rep(list(c()), length(dataTypes))
names(ppClust) <- c(dataTypes)

## force a specific # of clusterings to make the analysis computable in a reasonable time
# VI can compare groups with different # of clusterings so we can just take the best ICL per patient
cfgTable <- as.data.frame(t(read.table(paste0('tmp/',can,"/config.txt"),sep="\t",row.names=1,
                                       col.names=c("NULL","value"))));

minC <- as.numeric(levels(cfgTable$online.kmin))
maxC <- as.numeric(levels(cfgTable$online.kmax)) #Inf # to allow for automated selection

for (d in 1:length(dataTypes)) {
    
  # select data type
  data <- allData[[d]]
  
  # initiate a df (once - outside of the parallel loop, as it is not saved between parallel instances)
  ppList <- list()
  
  # set the number of patients that need internal clustering
  iterations <- 1:(length(data))
  
  ppClust[[d]] <- foreach(i=iterations) %dopar% { #  for (i in 1:(length(data)-1)) {

    ### first patient ###
    file <- data[i]
    p1 <- tail(unlist(strsplit(file, "[.]")), n=1) # id of a patient
    
    # create an adjacency matrix from patient's bipartite graph
    dat1 <- as.matrix(read.table(paste0(netDir, "/", file, netSuffix), sep="\t", head=F, fill=TRUE, stringsAsFactors=FALSE))
    dat1 <- t(dat1) %*% dat1
      
    # remove an edge filtered on a preselected threshold [optional]    
    thr <- thrs[[d]][i]
    dat1[dat1 <= thr] <- 0
    diag(dat1) <- 0 # get rid of self-loops

    # logit transform the data to make it analysable under Gaussian distribution
    # only Gaussian distribution can potentially be applied to our data (wrt its properties)
    dat1 <- plogis(dat1)
    diag(dat1) <- 0 # get rid of self-loops

    # find communities in the first patient
    res1 <- BM_gaussian(membership_type='SBM', adj=dat1, verbosity=0, autosave='', plotting='', exploration_factor=1.5, explore_min=minC, explore_max=maxC, ncores=threads)
    res1$estimate()

    # select the optimal number of clusters
    # if we set a hard minC and maxC, the optimal will be this number
    selOpt <- which.max(res1$ICL)
   
    # select data with the clusters
    clustering1 <- as.data.frame(res1$memberships[[selOpt]]$Z, row.names=nodeNames)
    clustering1$cluster <- apply(clustering1, 1, function(x) as.numeric(which.max(x))) # label data point's cluster
    comm <- clustering1$cluster
    
    # save results into a patient-patient network
    ppList[[p1]] <- comm
        
    # make the result final - necessary because of the parallel runs
    ppList <- ppList
    
  } 
}
#save.image()


### --- Comparing clustered patients with each other --- ###

## initialize patient-patient network
df <- data.frame(p1=character(), p2=character(), vi=as.numeric())
ppNet <- list()
for (i in 1:length(dataTypes)) {
    ppNet[[i]] <- df
}
names(ppNet) <- dataTypes

## run the comparison
for (d in 1:length(ppClust)) {
    data <- ppClust[[d]]
    
    # initiate a df (once - outside of the parallel loop, as it is not saved between parallel instances)
    ppDf <- df
    
    # set the number of patients (-1) for first loop
    iterations <- 1:(length(data)-1)
    
    ppNet[[names(ppClust[d])]] <- foreach(i=iterations, .combine=rbind) %dopar% {
        
### first patient ###
        p1 <- names(data[i][[1]]) # id of a patient
        
        # get the communities
        comm1 <- unlist(data[i][[1]]) ## compare() requires vectors
        
        js <- c((i+1):length(data))
        for (j in js) {
            
### next patient to be compared with the previous one ###
            p2 <- names(data[j][[1]]) # id of a patient
            
            # create an adjacency matrix from patient's bipartite graph; filter on a preselected threshold [optionally]
            comm2 <- unlist(data[j][[1]])
            
            # get a metric (distance!) of similarity between the two community structures
            vi <- compare(comm1, comm2, method = "vi")
            
            # save results into a patient-patient network
            ppRow <- data.frame(p1=p1, p2=p2, vi=vi)
            ppDf <- rbind(ppDf, ppRow)
            
        }
        
        # make the result final outside of the loop - because of the parallel runs
        ppDf <- rbind(ppDf)
        
    }
}

nPatients <- length(allData[[1]])


### --- Integrating patient-patient multi-valued edge --- ###

## average vi over data types to be integrated
dataInt <- dataTypes #c("rnaseq", "meth")

# initialize integrated data entry with vi from the 1st data type
ppNet$int <- ppNet[[dataInt[1]]]

# loop through the other data types
for (d in dataInt[-1]) { 
    ppNet$int$vi <- ppNet$int$vi + ppNet[[d]]$vi
}

# average the distance based on multiple data types
ppNet$int$vi <- ppNet$int$vi/length(dataInt)

# from VI calculate a weight, create an adjacency matrix
#dataSetups <- c(dataTypes, "int")
dataSetups <- c("int")
adj <- list()
for (d in 1:length(dataSetups)) {
    nam <- dataSetups[d]
    dat <- ppNet[[nam]]
    
    # normalize vi (as suggested in the original publication); distance - the smaller the more similar
    dat$nvi <- 1/log(nPatients)*dat$vi
    
    # weight with vi (for clustering with SBM); weight - the larger the more similar
    dat$wnvi <- 1-dat$nvi
    
    # from edge list create undirected graph
    sel <- c('p1','p2','wnvi')
    edge_list <- dat[,sel]
    G <- graph.data.frame(edge_list, directed=FALSE)
    attribute <- 'wnvi'

    # set edge weights to values from 'vi' attribute
    A <- as_adjacency_matrix(G, type="both", names=TRUE, sparse=FALSE, attr=attribute)
    adj[[d]] <- A    
}
names(adj) <- dataSetups

# only use the selected (common) patients for survival analysis
samData <- colnames(adj$int)
##samAvail <- intersect(samData, samClinical)
##adj <- lapply(adj, function(x) x[samAvail, samAvail])
adj <- lapply(adj, function(x) x[samData, samData])

### --- Clustering patients in the patient-patient network --- ###


clustMet <- c("spectrClust", "SBM")[2]

dataSel <- dataTypes 

## don't trust Firth cox when there are warnings
opt.prev<-options(warn=2);
fits <- list(int=list());
cns<-minC:maxC;
for (cn in cns) {
    patientGroups <- list()
    nodeNames <- rownames(adj$int)
    minC <- cn
    maxC <- cn
    ##for (d in 1:length(dataSetups)) {
    d <- 1; ## there is only 1 data setup - int
    nam <- dataSetups[d]
    dat <- adj[[nam]]        

    if (clustMet == "SBM") {
        ## group patients using SBM
        res <- BM_gaussian(membership_type='SBM', adj=dat, verbosity=0, autosave='', plotting='', exploration_factor=1.5, explore_min=minC, explore_max=maxC, ncores=threads)
        res$estimate()            
        ## select the optimal number of clusters
        selOpt <- cn           
        ## select data with the clusters
        clustering <- as.data.frame(res$memberships[[selOpt]]$Z, row.names=nodeNames)
        clustering$cluster <- apply(clustering, 1, function(x) as.numeric(which.max(x))) # label data point's cluster
        comm <- clustering$cluster
    }
    if (clustMet == "spectrClust") {
        ## group patients using spectral clustering
        C <- cn
        comm = spectralClustering(dat, C); # the final subtypes information
    }

    if ("clinicalData" %in% names(Data)){        
        library("survival")
        library("coxphf")
        ## methods: SNF (some functions)
        library("SNFtool")
        
        ## nodes are now the patients (previously pathways)
        clinTemplate <- clinical[nodeNames,] 
        nevent <- sum(clinTemplate$censor) 

        ## save results into a patient-patient network
        clinTemplate$group <- comm  
        patientGroups[[d]] <- clinTemplate
        ##}
        names(patientGroups) <- dataSetups    
        ## --- Survival analysis --- ##
        ## compare survival times for the different groups obtained from the clustering    
        ## add groupings from all data types
        clinicalAllData <- patientGroups[[1]][,-dim(patientGroups[[1]])[2]]
        clinicalAllData[,dataSetups[1]] <- patientGroups[[1]]$group
                                        #for (d in 2:length(dataSetups)) {
                                        #    nam <- dataSetups[d]
                                        #    clinicalAllData[,nam] <- patientGroups[[nam]]$group
                                        #}
        
        ## get data for survival analysis
        dat <- clinicalAllData[clinicalAllData$time>1,] # should remove earlier
        
        ## refactor to the factor with most patients    
        dat$intF <- as.factor(dat$int);
        maxSam <- which(table(dat$int) == max(table(dat$int)))
        dat$intF <- relevel(dat$intF, ref=maxSam[1])

                                        #    save(dat, file=paste0("patientGroups/VLN.", can, ".", paste(dataInt, collapse='.'), ".clusters", cn, ".Rda"))
        
        ## run the analyses
        fits[['int'   ]][[cn]] <- list(try(coxphf(Surv(time, status)~intF, data=dat)));
        fits[['int'   ]][[cn]][[2]] <- table(dat$intF)

        write.table(clustering,file="./integrated_stratification_vilon.csv",sep=",", dec=".", quote=T)
    } else{
        write.table(clustering,file="./integrated_stratification_vilon.csv",sep=",", dec=".", quote=T)
    }
    
}


if ("clinicalData" %in% names(Data)){
    options(opt.prev);  ## revert to as before
    ##pvalThr <- 1
    ##probThr <- 0.05 # for selecting factors

    pvalThr <- as.numeric(levels(cfgTable$online.pmodel))
    probThr <- as.numeric(levels(cfgTable$online.ppairwise))

    stats<-NULL;
    for (fn in names(fits)) {
        myfits<-fits[[fn]][cns];
        ## calculate the overall pvalue
        ps<-sapply(myfits,function(fit){
            f<-fit[[1]];
            ifelse(inherits(f,"try-error"), NA,
            (1 - pchisq(2 * diff(f$loglik), f$df)));
        });
        ## store the raw p-values
        pr <- ps
        ## adjust p-values for multiple testing
        ps<-p.adjust(ps,method="BH");
        ## select only significant results
        sel<-!is.na(ps) & ps<pvalThr;
        for (fit.n in which(sel)) {
            fit<-myfits[[fit.n]][[1]];
            cn<-cns[fit.n];
            coef.sel<-fit$prob<probThr;
            coefs<-exp(mean(abs(fit$coef[coef.sel])));
            gr.sizes <- myfits[[fit.n]][[2]][-1]; # remove the base group (doesn't count towards HRs)
            gr.base <- myfits[[fit.n]][[2]][1]    # get the size of the base group
            gr.sig <- as.character(paste(c(paste0("[", gr.base, "]"), gr.sizes[coef.sel]), collapse=', ' ));
            gr.ave <- round(mean(gr.sizes[coef.sel]), digit=0) # average size of the significant groups, without the base
            coefs.max <- NA
            if (length(fit$coef[coef.sel]) > 0) { coefs.max <- exp(max(abs(fit$coef[coef.sel]))) };
            coefs.p.all <- as.character(paste(fit$prob[coef.sel], collapse=', ' ))
            coefs.all <- as.character(paste(round(log2(exp((fit$coef[coef.sel]))), digits=2), collapse=', ' ));
            my.stats<-c(model=     fn,
                        AIC=       extractAIC(fit)[2],
                        BIC=       extractAIC(fit,k=log(nevent))[2],
                        cn=        cn,
                        facs=      paste(names(fit$coefficients)[coef.sel],collapse=' '),
                        p=         pr[fit.n],
                        q=         ps[fit.n],
                        avgHR=     coefs,
                        meanLog2HR=log2(coefs),
                        maxHR=     coefs.max,
                        maxLog2HR= log2(coefs.max),
                        log2HRs=   coefs.all,       # need to be log2'ed before as we store them as characters (see above)
                        grSize=    gr.sig,
                        grAve=     gr.ave,
                        nSig=      sum(coef.sel),
                        pctSig=    sum(coef.sel)/(cn-1),
                        pFacs=     coefs.p.all
                        );
            cat(paste(my.stats,collapse="\t"),"\n");
            if (is.null(stats)) stats<-data.frame(t(my.stats),stringsAsFactors=FALSE) else stats<-rbind(stats,my.stats);
        }
    }
    stats$fullSize <- nPatients
    met <- 'ViLoN'
    stats$meth <- met
    write.table(stats,file="./significant_model_results_vilon.csv",sep=",", dec=".", quote=T)

    dTypes <- paste(dataInt, collapse = ".")
    rdaNam <- paste0(can, "_", dTypes, "-int.", met)
    save(stats, file=paste0(rdaNam, ".Rda"))
}
