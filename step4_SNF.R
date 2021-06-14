#install.packages("SNFtool")
library("SNFtool")
library("edgeR")

## Perform the below analysis as suggested in the original SNF publication
# we played with some parameters, normalizations, etc. but haven't seen improvements (data not shown)

## First, set all the parameters:
K = 20;		# number of neighbors, usually (10~30)
alpha = 0.5;  	# hyperparameter, usually (0.3~0.8)
T = 10; 	# Number of Iterations, usually (10~20)

## Data is of size n x d_i, where n is the number of patients, d_i is the number of genes / methylated points

## load the relevant data
can <- "THCA"
dataTypes <- c("rnaseq", "meth", "cnvsnp")
datDir <- paste0('tmp/',can,'/intermediate_files/')
file=paste0(datDir, paste(dataTypes, collapse='.'),".T.comPatients.",can,".rda")
load(file)

Data1 <- t(rnaseqT)

## Clinical information
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
samAll <- rownames(Data1)
complete <- intersect(samClinical, samAll)

dataInt <- c("rnaseq", "meth")

Data1 <- t(rnaseqT[,complete])
Data2 <- t(methT[,complete])

## Calculate distance matrices(here we calculate Euclidean Distance, you can use other distance, e.g,correlation)

## If the data are all continuous values, we recommend the users to perform standard normalization before using SNF, though it is optional depending on the data the users want to use.  

Data1 = standardNormalization(Data1);
Data2 = standardNormalization(Data2);

## Calculate the pair-wise distance; If the data is continuous, we recommend to use the function "dist2" as follows; if the data is discrete, we recommend the users to use ""
Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));

## next, construct similarity graphs
W1 = affinityMatrix(Dist1, K, alpha)
W2 = affinityMatrix(Dist2, K, alpha)

## next, we fuse all the graphs
## then the overall matrix can be computed by similarity network fusion(SNF):
W = SNF(list(W1, W2), K, T)

## With this unified graph W of size n x n, you can do clustering


### --- survival analysis --- ###
library("survival")
library("coxphf")
    
clinTemplate <- clinical[complete,]

## don't trust Firth cox when there are warnings
opt.prev<-options(warn=2);
## test for the optimal #C
fits<-list(int=list());
cns<-2:9;
for (cn in cns) {
    comm.SNF.int = spectralClustering(W, cn);
    comm.SNF <- data.frame(sample=colnames(W1), int=comm.SNF.int)
    rownames(comm.SNF) <- comm.SNF$sample
    
    ## make sure both objects have the same aptient order
    dat <- clinTemplate[rownames(comm.SNF),]
    dat <- cbind(dat, comm.SNF$int)
    colnames(dat)[length(colnames(dat))]  <- "int"
    nevent <- sum(clinTemplate$censor)

    ## refactor to the factor with most patients   
    dat$intF <- as.factor(dat$int);
    maxSam <- which(table(dat$int) == max(table(dat$int)))
    dat$intF <- relevel(dat$intF, ref=maxSam[1])

#    save(dat, file=paste0("patientGroups/SNF.", can, ".", paste(dataInt, collapse='.'), ".clusters", cn, ".Rda"))
    
    ## run the analyses
    fits[['int'   ]][[cn]] <- list(try(coxphf(Surv(time, status)~intF, data=dat)));
    fits[['int'   ]][[cn]][[2]] <- table(dat$intF)
    
}
options(opt.prev);  ## revert to as before

pvalThr <- 1
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
stats$fullSize <- dim(Data1)[1]
met <- 'SNF'
stats$meth <- met

dTypes <- paste(dataInt, collapse = ".")
rdaNam <- paste0(can, "_", dTypes, "-int.", met)
save(stats, file=paste0(rdaNam, ".Rda"))
