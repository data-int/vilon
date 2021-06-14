# NOTE: only works from this folder
# NOTE2: only can be run 1 line at a time, or returns some errors

### --- Load the molecular data
can <- "THCA"
dataTypes <- c("rnaseq", "meth", "cnvsnp")
datDir <- paste0('../tmp/',can,'/intermediate_files/')
file=paste0(datDir, paste(dataTypes, collapse='.'),".T.comPatients.",can,".rda")
load(file)

Data1 <- rnaseqT

### --- Clinical information
clinical <- clinicalData
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
rownames(clinical) <- clinical$sample
samClinical <- rownames(clinical)
samAll <- colnames(Data1)
complete <- intersect(samClinical, samAll)
clinTemplate <- clinical[complete,]

## select data to be integrated, and for common patients
dataInt  <- c("rnaseq", "meth")
com1 <- rnaseqT[,complete]
com2 <- methT[,complete]
    
### --- Preprocess the data
library("SNFtool")
library("edgeR")

## as for VLN
library("preprocessCore")
Data1 <- DGEList(counts=com1); Data1 <- calcNormFactors(Data1); v <- voom(Data1,design=NULL,plot=FALSE); Data1 <- v$E; Data1 <- 2^Data1 # positives
Data2 <- scale(com2)  # log2 values, negatives included

### --- Run LRAcluster
source('R/LRAcluster.R')
source('R/poisson.R')
source('R/gaussian.R')
source('R/binary.R')
source('R/na.R')

data <- list(Data1, Data2)
names <- dataTypes

### --- survival analysis
library("survival")
library("coxphf")

opt.prev<-options(warn=2);
## test for the optimal #C
fits<-list(int=list());
cns<-2:9;

## run LRAcluster
types <- c("poisson", "gaussian")
rlist<-LRAcluster(data=data,types=types,names=names,dimension=2)

## test multiple cluster numbers
for (cn in cns) {

    ## cluster patients based on the LRAcluster results
    rclust<-kmeans(t(rlist$coordinate),cn)
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
options(opt.prev);  ## revert to as before

nevent <- sum(clinTemplate$censor)

pvalThr <- 1 #0.05
probThr <- 0.05 # for selecting factors
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
