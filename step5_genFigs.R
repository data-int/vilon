### --- TCGA KIRC  --- ###

can <- "THCA"
dataTypes <- c("rnaseq", "meth")
datDir <- ''

## ViLoN
met <- 'ViLoN'
nam <- paste0(datDir, can, '_', paste(dataTypes, collapse='.'),'-int.',met,".Rda")
load(nam); d1 <- stats

## SNF
met <- 'SNF'
nam <- paste0(datDir, can, '_', paste(dataTypes, collapse='.'),'-int.',met,".Rda")
load(nam); d2 <- stats

## LRAcluster
met <- 'LRAcluster'
nam <- paste0(datDir, can, '_', paste(dataTypes, collapse='.'),'-int.',met,".Rda")
load(nam); d3 <- stats

d <- rbind(d1, d2, d3)

## select only integration results
d <- d[d$model == 'int',]

d$id<-seq(NROW(d));
dd<-data.frame(log2HR=as.numeric(unlist(strsplit(as.character(d$log2HRs),", "))));
dd$grSize<-as.numeric(unlist(sapply(strsplit(d$grSize,", "), function(r){r[-1]})));
dd$fullSize<-rep(d$fullSize,sapply(strsplit(as.character(d$log2HRs),", "),length));
dd$met<-rep(d$met,sapply(strsplit(as.character(d$log2HRs),", "),length));
dd$model<-rep(d$model,sapply(strsplit(as.character(d$log2HRs),", "),length));
dd$cn<-rep(d$cn,sapply(strsplit(as.character(d$log2HRs),", "),length));
dd$AIC<-rep(d$AIC,sapply(strsplit(as.character(d$log2HRs),", "),length));
dd$id<-rep(d$id,sapply(strsplit(as.character(d$log2HRs),", "),length));
dd$absHR<-round(2^abs(dd$log2HR), 2);
dd$prod<-round((dd$absHR-1)*dd$grSize, 2) / dd$fullSize * 1000; # bring to 100 patients
dd$refSize<-as.numeric(unlist(sapply(strsplit(d$grSize,", "), function(r){len <- length(r); rep(gsub("]", "", gsub("\\[", "", r[1])),len-1)})));
dd$pctSig<-rep(d$pctSig,sapply(strsplit(as.character(d$log2HRs),", "),length));
dd$q.model <- rep(d$q,sapply(strsplit(as.character(d$log2HRs),", "),length));

# we look at cn=2..9, integrated, adjusted p-value relevant
dd$p<-as.numeric(unlist(strsplit(d$pFacs,", ")))

# adjust p-values per tool
metNam <- 'ViLoN'
dat <- dd[dd$met == metNam, ]
dat$q <- p.adjust(as.numeric(dat$p), method='BH')
cc <- dat

metNam <- 'SNF'
dat <- dd[dd$met == metNam, ]
dat$q <- p.adjust(as.numeric(dat$p), method='BH')
cc <- rbind(cc, dat)

metNam <- 'LRAcluster'
dat <- dd[dd$met == metNam, ]
dat$q <- p.adjust(as.numeric(dat$p), method='BH')
dd <- rbind(cc, dat)

# select only significant factors
dd <- dd[as.numeric(dd$q.model) < 0.05,]
dd <- dd[order(-dd$prod),]

## Published results from TCGA, Nature 2013 ccRCC paper [https://www.nature.com/articles/nature12222]

litProd <- read.table("literature/tcga_cell14_clinRiskRes.txt", sep="\t", header=T, stringsAsFactors=F)
litProd$id <- c((max(dd$id)+1):((max(dd$id))+length(litProd$model)))
litProd$AIC <- NA
litProd$q.model <- litProd$q
litProd <- litProd[litProd$q<.05,]

## save the table
tabSel <- c("met", "model", "cn", "q", "q.model", "AIC", "grSize", "refSize", "fullSize", "absHR", "log2HR", "prod")
tabSav <- rbind(dd[,tabSel], litProd[,tabSel])
write.table(tabSav[order(-tabSav$prod),], 'THCA.all.tab', sep = "\t", row.names=F, quote=F)

sel <- c("met", "model", "cn", "q", "id", "absHR", "log2HR", "grSize", "prod", "refSize", "fullSize")
ddd <- rbind(dd[,sel], litProd[,sel])
ddd <- ddd[order(-ddd$prod),]

mets<-sort(unique(ddd$met));
mods<-sort(unique(ddd$model));

combs <- data.frame(met=factor(ddd$met,levels=c('ViLoN','SNF', 'LRAcluster', 'clinical')),
                    mod=factor(ddd$model,levels=c('int', 'tcga14_macis', 'tcga14_risk')),cn=ddd$cn)
my.order<-order(combs$met,combs$mod,-as.numeric(combs$cn));
combs<-combs[my.order,];
ddd<-ddd[my.order,];

col.idx<-apply(combs[,c('mod','cn')],1,paste,collapse='.');

allcols<-c(int.9='black',
           int.8='darkred',
           int.7='brown4',
           int.6='brown2',
           int.5='darkorange4',
           int.4='darkorange1',
           int.3='red',
           int.2='magenta',
           tcga14_macis.4='cyan',
           tcga14_risk.3='plum2'
           );

msel<-T;  # ddd$met %in% c('ViLoN','Kocak','Theissen');
#msel <- ddd$met %in% c('clinical')

pchMap<-c(
    ViLoN=6,
    LRAcluster=1,
    SNF=2,
    clinical=0
)

pchs<-pchMap[ddd$met];
ymax<-100*max(ddd$absHR)+1;


## main plot
X11(width=3*2.1,height=3*2.1)
plot((ddd$grSize/ddd$fullSize*100)[msel],
     (100*(ddd$absHR-1))[msel],
     ylim=c(0,ymax),
     xlim=c(0,50),
     pch=pchs,
     cex=sqrt(ddd$prod/100)[msel],
     col='black',
     main='Patient-stratification trade offs\nPatient numbers vs Risk prediction',
     sub='symbol size = product-score',
     xlab='affected patient numbers [%]',
     ylab='relative risk change [%]',
     cex.axis=1.3,
     cex.lab=1.3,
     cex.sub=1.3);

## mark the highest products per each method
ddd.best <- c()
for (met in unique(ddd$met)) {
    d <- ddd[ddd$met==met,]
    ddd.best <- rbind(ddd.best, d[d$prod == max(d$prod),])
}
points(x=(ddd.best$grSize/ddd.best$fullSize*100), y=(100*(ddd.best$absHR-1)),
       cex=sqrt(ddd.best$prod/100), pch= c(25, 24, 15, 21), 
       col=c('magenta','darkorange1','cyan', 'black'), #allcols[paste(ddd.best$model,ddd.best$cn,sep='.')],
       bg=c('magenta','darkorange1','cyan', 'black') #allcols[paste(ddd.best$model,ddd.best$cn,sep='.')]
           )

points(0,0,col='black',pch=20)

labelsAll <- unique(unlist(apply(combs[msel,],1,paste,collapse='.')));
colSel <- unlist(lapply(labelsAll, function(x) paste(unlist(strsplit(x, '\\.'))[2:3], collapse='.')))
pchAll <- unique(data.frame(met=ddd$met[msel], pch=pchs))
rownames(pchAll) <- pchAll$met
pchSel <- unlist(lapply(labelsAll, function(x) paste(unlist(strsplit(x, '\\.'))[1], collapse='.')))

labels.uq<-sub('[.][0-9]','',labelsAll);
sel.uq<-!duplicated(labels.uq);

legend("topright", c('ViLoN', 'SNF', 'LRAcluster', 'reference', 'best model'),
       pch=c(6,2,1,0,20), # as.numeric(factor(combs[,1],levels=mets)),
       col=c(rep('black',4),'white'), cex=1.1)

abline(h=0,lty='dotted');
dev.copy2pdf(file=paste0('THCA.int+lit.black.v3.pdf'))

