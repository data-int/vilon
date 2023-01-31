# NOTE: only works from this folder
# NOTE2: only can be run 1 line at a time, or returns some errors
##setwd("/bi/ala/scratch/bibot/ViLoN/7272197557c4a9d1e47eaa473f7cf628")
can <- "vilon.online"
datDir <- paste0('./')
metafile <- paste0("./","meta",".T.comPatients.",can,".rda");
cat("Loading",metafile,"\n");
load(metafile); ## cfgTable, dataTypes, nData, dataTypesVoom


if(cfgTable$online.lra=="lra"){
### --- Load the molecular data

results <- file.path(datDir, "results")
dir.create(results,recursive=T)

file=paste0(datDir, paste(dataTypes, collapse='.'),".T.comPatients.",can,".rda")
load(file,Data <- new.env())
Data <- as.list(Data)

Data1 <- Data[[1]]

### --- Clinical information
if ("clinicalData" %in% names(Data)){
    clinical <- Data$clinicalData
    clinical[is.na(clinical[,"daystolastfollowup"]),"daystolastfollowup"] = clinical[is.na(clinical[,"daystolastfollowup"]),"daystodeath"]
    clinical[clinical[,"gender"]=="male","gender"] = clinical[clinical[,"gender"]=="male","gender"]=1
    clinical[clinical[,"gender"]=="female","gender"] = clinical[clinical[,"gender"]=="female","gender"]=2
    clinical[clinical[,"gender"]=="male","gender"] = clinical[clinical[,"gender"]=="male","gender"]=1
    clinical <- data.frame(sample=rownames(clinical),
                           time=as.numeric(clinical[,"daystolastfollowup"]),
                           censor=as.numeric(clinical[,"vitalstatus"]),
                           status=as.numeric(clinical[,"vitalstatus"])+1, # for library(survival) this means: 1=censored, 2=dead
                           sex=as.numeric(clinical[,"gender"]),
                           age=as.numeric(clinical[,"yearstobirth"]),
                           stage=(clinical[,"pathologicstage"]))
    clinical <- clinical[!is.na(clinical$time),]
    rownames(clinical) <-sapply(clinical$sample,function(n){gsub("[-]",".",n)},USE.NAMES=F)
    samClinical <- rownames(clinical)
    samAll <- colnames(Data[[1]])
    samClinical <- intersect(samClinical, samAll)
    complete <- samClinical
    clinTemplate <- clinical[complete,]
}else{
    samClinical <-colnames(Data[[1]])
    complete <- samClinical
}

## select data to be integrated, and for common patients
dataInt  <- dataTypes
    
### --- Preprocess the data
library("SNFtool")
library("edgeR")

## as for VLN
library("preprocessCore")
data <- sapply(dataTypes,function(d){
    if (d=="rnaseq"){
        Data1 <- DGEList(counts=Data[[paste0(d,"T")]][,complete]);
        Data1 <- calcNormFactors(Data1);
        v <- voom(Data1,design=NULL,plot=FALSE);
        Data1 <- v$E;
        Data1 <- 2^Data1 # positives
    }else{
        Data1 <- scale(Data[[paste0(d,"T")]][,complete])  # log2 values, negatives included
    }
    return(Data1)
},USE.NAMES=T)

### --- Run LRAcluster
source('R/LRAcluster.R')
source('R/poisson.R')
source('R/gaussian.R')
source('R/binary.R')
source('R/na.R')

### --- survival analysis
library("survival")
library("coxphf")

opt.prev<-options(warn=2);
## test for the optimal #C
fits<-list(int=list());
minC <- as.numeric(cfgTable$online.kmin)
maxC <- as.numeric(cfgTable$online.kmax)
cns<-minC:maxC;

## run LRAcluster
types <- c("rnaseq"="poisson", "meth"="gaussian", "cnvsnp"="gaussian")
types <- unname(types[dataTypes])
rlist<-LRAcluster(data=data,types=types,names=names(data),dimension=2)

## test multiple cluster numbers
for (cn in cns) {
    ## cluster patients based on the LRAcluster results
    rclust<-kmeans(t(rlist$coordinate),cn)
    ##write.table(rclust,file=paste0(results,"/integrated_stratification_LRA.",cn,".clusters.csv"),row.names=F,sep=",", dec=".",quote=F)
    if ("clinicalData" %in% names(Data)){
        dat <- clinTemplate[names(rclust$cluster),] # make sure both objects have the same patient order
        dat$int <- as.numeric(rclust$cluster)
        
        ## refactor to the factor with most patients
        dat$intF<-as.factor(dat$int);
        maxSam <- which(table(dat$int) == max(table(dat$int)))
        dat$intF <- relevel(dat$intF, ref=maxSam[1])
        
        ## run the analyses
        ## results for v1.V3.V1
        fits[['int'   ]][[cn]] <- list(try(coxphf(Surv(time, status)~intF, data=dat)));
        fits[['int'   ]][[cn]][[2]] <- table(dat$intF)
    }
}

if ("clinicalData" %in% names(Data)){
    options(opt.prev);  ## revert to as before
    nevent <- sum(clinTemplate$censor)
    pvalThr <- as.numeric(cfgTable$online.pmodel)
    probThr <- as.numeric(cfgTable$online.ppairwise)
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
    stats$fullSize <- dim(Data1)[2]
    met <- 'LRAcluster'
    stats$meth <- met

    dTypes <- paste(dataInt, collapse = ".")
    rdaNam <- paste0(can, "_", dTypes, "-int.", met)
    save(stats, file=paste0("../", rdaNam, ".Rda"))
}

}
