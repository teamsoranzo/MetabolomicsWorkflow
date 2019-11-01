library("lattice")
library("gdata")
umethods<-c("singlevar")
mycols<- c("pale green","dark sea green","sea green","dark cyan","steel blue","dark turquoise","sky blue","royal blue","dark blue","dark slate blue")
  p <- NULL
  obs.max <- exp.max <- NULL
getwd()
fn<-"/lustre/scratch115/realdata/mdt0/projects/int_wes_metabol/raremetal/AllMetabolites/raremetalworker_inverse_norm/ALLChr_Trial_M53106.singlevar.score.txt"
 f <- read.csv(fn, sep="", stringsAsFactors=FALSE, header=T)
f<-f[f$ALL_AF>=0.001,]
Pval<-"PVALUE"
   p <- f[,Pval]
    p <- p[!is.na(p)]
    p <- as.numeric(gsub("<","",p))
    ord <- order(p, na.last=NA)
    p <- p[ord]
    ## allocate min to p==0
    p[p==0] <- min(p[p>0])
    obs.p <- p
    log.obs.p <- -(log10(obs.p))
    exp.p <- c(1:length(obs.p))
    log.exp.p <- -(log10( (exp.p-0.5)/length(exp.p)))
    obs.max <- round(max(log.obs.p),2)
    exp.max <- round(max(log.exp.p),2)
xlim <- c(0,max(exp.max))
ylim <- c(0,max(obs.max)+0.3)
## for confidence interval
n <- ceiling(10^(max(exp.max)))
e <- -log10(1:n/n)
c95 <- sapply(1:n, function(i) qbeta(0.95,i,n-i+1))
c05 <- sapply(1:n, function(i) qbeta(0.05,i,n-i+1))
dat <- cbind(c05, c95)
getwd()
png("/lustre/scratch115/realdata/mdt0/projects/int_wes_metabol/raremetal/AllMetabolites/raremetalworker_inverse_norm/qq-plot_TrialFromraremetalworker_MAF0.001_M53106.png")
plot(xlim, ylim, type="n", xlab=expression(paste(Expected -log[10], italic(p),sep="")), ylab=expression(paste(Observed -log[10], italic(p),sep="")), main="qq-plot_TrialFromraremetalworker_M53106",  las=1, axes=F, cex.lab=1.2)
polygon(c(e,exp.max), c(-log10(c95),exp.max), col="lightgrey", border="lightgrey")
polygon(c(e,exp.max), c(-log10(c05),exp.max), col="lightgrey", border="lightgrey")
lines(c(0,exp.max), c(0,exp.max), col="darkgray", lwd=1)
axis(1)
axis(2, las=1)
lambda <- NULL
  points(log.exp.p, log.obs.p, pch=20, cex=1.2, col=mycols[1])
  ind <- which(p<1e-50)
  if (length(ind)!=0)  {
    for (j in 1:length(ind))  {
      text(log.exp.p[ind[j]], log.obs.p[ind[j]],"M53106", pos=2, col=mycols[1])
    }
  }
X <- qchisq(1 - obs.p, df=1)
lambda <- round(median(X) / qchisq(0.5, df=1),2)
legend.text <- paste("M53106"," ", lambda, sep="")
legend(0, ylim[2], legend=legend.text, fill=mycols, box.col="white", border=mycols, cex=0.6)
dev.off()



f <- read.csv(fn, sep="", stringsAsFactors=FALSE, header=T)
f$LogPVALUE<-(-log10(f$PVALUE))

pdf("AlleleFreq.pdf")
par(mfrow=c(2,1))
hist(f$ALL_AF[f$LogPVALUE>=2.5 & f$LogPVALUE<=4.2], 500)
hist(f$ALL_AF, 500)
dev.off()


pdf("Pvals.pdf")
par(mfrow=c(2,1))
hist(f$LogPVALUE[f$LogPVALUE>=2.5 & f$LogPVALUE<=4.2], 500)
abline(h=100)
hist(f$LogPVALUE, 500)
dev.off()

tbl<-table(f$LogPVALUE[f$LogPVALUE>=2.5 & f$LogPVALUE<=5])

PValblock<-names(tbl)[as.vector(tbl) >=100]
AF_PValblock <- NULL
for (i in PValblock){
AF_PValblock[[i]]<-f$ALL_AF[which(f$LogPVALUE==i)]
}


pdf("BoxplotPValblock.pdf")
par(mar=c(12,10,1,1)+0.1,mgp=c(8,1,0))
boxplot(AF_PValblock, las=2, xlab="Logpvalue", ylab="AF")
dev.off()

## Look if the same sample has a singletons in that pvalue and if there another sample has a singletons for the other pvalue.
SNPselection<-NULL
for (y in names(tbl)[as.vector(tbl) >=100][1:2]){

Chr<-f[which(f$LogPVALUE==y),1]
Pos<-f[which(f$LogPVALUE==y),2]
SNPselection[[y]]<-paste0(Chr,":",Pos,"-",Pos)
write.table(paste0(Chr,":",Pos,"-",Pos), file=paste0("SNPselection",y), quote=F, col.names=F, row.names=F)
}





library(VariantAnnotation)

#######2.5447697794195

#~~~To be Run in Bash
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# generating the header
tabix -fh /lustre/scratch115/projects/int_wes_metabol/WES_QCeddata_37/INTERVAL_QCv1.0.AC-AN.vcf.gz chr1:1 > genotypes_SNPselection2.5447697794195.vcf
# then you can call all your positions
xargs -a SNPselection2.5447697794195  -I {} tabix -f /lustre/scratch115/projects/int_wes_metabol/WES_QCeddata_37/INTERVAL_QCv1.0.AC-AN.vcf.gz  {} >> genotypes_SNPselection2.5447697794195.vcf
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fl <- "/lustre/scratch115/realdata/mdt0/projects/int_wes_metabol/raremetal/AllMetabolites/raremetalworker_inverse_norm/genotypes_SNPselection2.5447697794195.vcf"
vcf <- readVcf(fl, "hg19")


listIndsKeepNoREl<-read.table("/lustre/scratch115/realdata/mdt0/projects/int_wes_metabol/Metabolon/CleanedData_160414/Pheno/listIndsKeepNoREl.txt", header=F, stringsAsFactors=F)
listIndsKeepNoREl<-as.vector(as.matrix(listIndsKeepNoREl))
GenotypesIndsKeepNoREl<-subset(GenotypesALL, select = colnames(GenotypesALL)[colnames(GenotypesALL) %in% listIndsKeepNoREl])


Saple<-NULL
for (z in 1:dim(GenotypesIndsKeepNoREl)[1]){
Saple[z]<-names(which(GenotypesIndsKeepNoREl[z,]=="0/1"))
}

GenotypesIndsKeepNoREl2.5447697794195<-GenotypesIndsKeepNoREl
Saple2.5447697794195<-Saple

########2.55615677599041

#~~~To be Run in Bash
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# generating the header
tabix -fh /lustre/scratch115/projects/int_wes_metabol/WES_QCeddata_37/INTERVAL_QCv1.0.AC-AN.vcf.gz chr1:1 > genotypes_SNPselection2.55615677599041.vcf
# then you can call all your positions
xargs -a SNPselection2.55615677599041  -I {} tabix -f /lustre/scratch115/projects/int_wes_metabol/WES_QCeddata_37/INTERVAL_QCv1.0.AC-AN.vcf.gz  {} >> genotypes_SNPselection2.55615677599041.vcf
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


fl <- "/lustre/scratch115/realdata/mdt0/projects/int_wes_metabol/raremetal/AllMetabolites/raremetalworker_inverse_norm/genotypes_SNPselection2.55615677599041.vcf"
vcf <- readVcf(fl, "hg19")


listIndsKeepNoREl<-read.table("/lustre/scratch115/realdata/mdt0/projects/int_wes_metabol/Metabolon/CleanedData_160414/Pheno/listIndsKeepNoREl.txt", header=F, stringsAsFactors=F)
listIndsKeepNoREl<-as.vector(as.matrix(listIndsKeepNoREl))
GenotypesIndsKeepNoREl<-subset(GenotypesALL, select = colnames(GenotypesALL)[colnames(GenotypesALL) %in% listIndsKeepNoREl])


Saple<-NULL
for (z in 1:dim(GenotypesIndsKeepNoREl)[1]){
if(GenotypesIndsKeepNoREl[z,]=="0/1"){
Saple[z]<-names(which(GenotypesIndsKeepNoREl[z,]=="0/1"))
print(names(which(GenotypesIndsKeepNoREl[z,]=="0/1")))
#}else if(){
#There are tripletons#### Check this out
#}else{
#}
}




#idx <- "/lustre/scratch115/projects/int_wes_metabol/WES_QCeddata_37/Chr20/Chr20_INTERVAL_QCv1.0.AC-AN_SelectInds.vcf.gz.tbi"
#rngs <- GRanges("20", IRanges(c(251967,251967), c(3650380, 3650380)))
#       names(rngs) <- c("20:251967", "20:3650380")
#       param <- ScanVcfParam(which=rngs)
#       tab <- TabixFile(fl, idx)
#       vcf <- readVcf(tab, "hg19", param)

