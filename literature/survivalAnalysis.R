
### Risk labels developed in the main TCGA-THCA Cell 2014 paper
risk <- read.csv('mmc3_THCA-TP.csv', sep='\t', stringsAsFactors=FALSE, na.strings=c("", "NA"))
rownames(risk) <- gsub("-01.*", "", risk$sample) # to match the samples to our dataset
sel <- c("Risk", "GDAC_MACIS_category")
risk <- risk[,sel]
#risk <- risk[complete.cases(risk),] # rm NA labels

### Clinical information, available on TCGA
can <- "THCA"
dataTypes <- c("rnaseq", "meth", "cnvsnp")
file=paste0("/bi/aim/scratch/mkandula/TCGA/TCGA2STAT/", can,"/",paste(dataTypes, collapse = "."),".T.comPatients.", can, ".rda")
load(file)
## patients who were alive at time of last check have "daystolastfollowup", dead patients have "daystodeath", "vitalstatus" tells us which is which so we can combine both times
clinical <- clinicalData
clinical[is.na(clinical[,"daystolastfollowup"]),"daystolastfollowup"] = clinical[is.na(clinical[,"daystolastfollowup"]),"daystodeath"]
clinical <- data.frame(sample=rownames(clinical),
                       time=as.numeric(clinical[,"daystolastfollowup"]),
                       censor=as.numeric(clinical[,"vitalstatus"]),
                       status=as.numeric(clinical[,"vitalstatus"])+1, # for library(survival) this means: 1=censored, 2=dead
                       sex=as.numeric(clinical[,"gender"]))
clinical <- clinical[!is.na(clinical$time),]
clinical <- clinical[!is.na(clinical$censor),]
clinical <- clinical[!clinical$time <= 0,]

rownames(clinical) <- clinical$sample
samClinical <- rownames(clinical)

samRisk <- rownames(risk)

clin.lit <- cbind(clinical[samRisk,], risk[samRisk,])

### survival analysis
library("survival")
library("coxphf")

## prep dataset
dat <- clin.lit[clin.lit$time>1,] # all good
dat <- dat[!is.na(clin.lit$time),]

res <- c()
sel <- c('Risk', "GDAC_MACIS_category")[1]

## refactor to the factor with most patients
dat <- dat[!is.na(dat[,sel]),]
dat[,sel] <- as.factor(dat[,sel]);
maxSam <- which(table(dat[,sel]) == max(table(dat[,sel])))
dat$riskF <- relevel(dat[,sel], ref=maxSam[1])

## run Cox Regression
fit <- coxphf(Surv(time, status)~riskF, data=dat)

res <- data.frame(
    model=c(rep('tcga14_risk',2), rep('tcga14_macis',3)),
    met=rep('clinical',5),
    cn=c(rep(3,2),rep(4,3)),
    HR=c(7.875000,1.085256,5.77739,24.40269,50.88727),
    p=c(0.003447496,0.905320434,1.514900e-01,8.427458e-05,2.073818e-06),
    p.model=c(rep(0.01146057,2),rep(1.277974e-06,3)),
    grSize=c(24,171,63,38,27),
    refSize=c(rep(259,2),rep(316,3)),
    fullSize=c(rep(454,2),rep(444,3))
)

res$q <- p.adjust(res$p, method="BH")
res$q.model <- c(rep(p.adjust(unique(res$p.model), method="BH")[1],2), rep(p.adjust(unique(res$p.model), method="BH")[2],3))
res$log2HR=log2(res$HR)
res$absHR=round(2^abs(res$log2HR), 2)
res$AIC=extractAIC(fit)[2]
nevent <- sum(dat$censor)
res$BIC=extractAIC(fit,k=log(nevent))[2]
res$prod<-round((res$absHR-1)*res$grSize, 2) / res$fullSize * 1000; # bring to 100 patients

write.table(res, 'tcga_cell14_clinRiskRes.txt', quote=F, row.names=F, sep='\t')
