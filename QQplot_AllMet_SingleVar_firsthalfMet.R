#!/usr/bin/Rscript
library("lattice")
library("gdata")

##### Number of variants per approach
# DistrVarALL <- read.table("~/Work/INTERVAL/SKAT/RES_AllchrsAllMethods_SKAT_ResStand_NewFAM/DistrVarALL.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)
# colnames(DistrVarALL)<-c("Chr","naive","functional","LoF")
# ymax=max(as.vector(DistrVarALL[,-1]))*1.1
# png(paste("/Users/lb17/Work/INTERVAL/SKAT/RES_AllchrsAllMethods_SKAT_ResStand_NewFAM/DistributionVarAll.png",sep = ""), pointsize=18, width=700, height=600)
# par(mai=c(1.3,1.3,0.3,0.3))
# barplot(t(DistrVarALL[,-1]), names.arg = t(DistrVarALL[,1]) ,las=2, beside = T, ylim =c(0,ymax) , col = c("dark turquoise","sky blue","royal blue"))
# legend("topright",legend = c("naive","functional","LoF"),fill=c("dark turquoise","sky blue","royal blue"))
# dev.off()
# colSums(DistrVarALL[,-1])
##### QQ plots for exome-wide SKAT with multiple traits ##### NAIVE

print(commandArgs(trailingOnly=TRUE))
## Reads ${LSB_JOBINDEX}
get.Args <- function() {
   ## Extract argument names and corresponding values
   my_args_list   <- strsplit(grep("=", gsub("--", "", commandArgs(trailingOnly=TRUE)), value = T), "=")
   my_args        <- lapply(my_args_list, function(x) x[2])
   names(my_args) <- lapply(my_args_list, function(x) x[1])
   ## Assign argument values to global variables in R
   for (arg_name in names(my_args)) {
     eval(parse(text=paste(arg_name, " <<- my_args[[arg_name]]", sep="")))
   }
}
## ID is the global parameter from ${LSB_JOBINDEX}
get.Args()
a <- as.numeric(ID)
m <- as.numeric(Meth)

Approaches <- c("naive") # kept like this just for convinience, no meaning for naive in singlevar
print(Approaches[a])

if(Approaches[a]!="LoF"){
	MinNumVarWin=5
}else{
	MinNumVarWin=2
}

umethods<-c("singlevar")
print(umethods[m])

#dir<-list.dirs(path = paste0("/lustre/scratch115/projects/int_wes_metabol/raremetal/AllMetabolites/raremetal_inverse_norm_AllMet/",Approaches[a]),full.names=F,recursive=F)
#TT<-matrix(unlist(strsplit(dir, "/")),ncol=10, byrow=TRUE)[,10]
#traits<-grep("^M[0-9]+", TT, value=T, perl=T)

traits<-readLines("/lustre/scratch115/projects/int_wes_metabol/Metabolon/CleanedData_160414/Pheno/AllListMet_No11Fail_No2SegFault.txt")

PathF<-paste0("/lustre/scratch115/realdata/mdt0/projects/int_wes_metabol/raremetal/AllMetabolites/raremetal_inverse_norm_AllMet/","Singlevar")

Filein<-list.files(path = paste0(PathF,"/"), full.names = F,  recursive = F, pattern =paste0("Combined_ALL_INTERVAL_QCv1.0_inverse_norm_",Approaches[a],".*.",traits[1],".meta.*.results"))
ptmA <- proc.time()
fspl<-strsplit(Filein,"\\.")
#me<-NULL
ma<-NULL
for(y in 1:length(fspl)){
#  	me[y]<-fspl[[y]][7]
  	ma[y]<-paste0(fspl[[y]][3],".",fspl[[y]][4])
}
#umethods<-unique(me)
maf<-unique(ma)
lab.traits<-traits

#mycol <- c("chocolate","dark orange","orange","gold","yellow green","green yellow","olive drab","dark green","lime green","green","spring green","pale green","dark sea green","sea green","dark cyan","steel blue","dark turquoise","sky blue","royal blue","blue","dark blue","dark slate blue","slate blue","medium purple","blue violet","dark orchid","dark magenta","medium orchid","orchid","plum","medium violet red")
#mycols<-mycol[1:length(traits)]
mycols<- c("pale green","dark sea green","sea green","dark cyan","steel blue","dark turquoise","sky blue","royal blue","dark blue","dark slate blue")
for (mf in 1:length(maf)){
  print(maf[mf])
  #m=1
  p <- NULL
  obs.max <- exp.max <- NULL
  for (i in 1:length(traits)){
    #i=1
    print(traits[i])
    fn <- paste(PathF,"/Combined_ALL_INTERVAL_QCv1.0_inverse_norm_",Approaches[a],".",maf[mf],".",traits[i],".meta.",umethods[m],".results", sep="")
    f <- read.csv(fn, sep="", stringsAsFactors=FALSE, header=T)
    f<-f[f$POOLED_ALT_AF>=0.001,]
    #head(f)
    if(umethods[m] == "SKAT_"){
      Pval<-"PVALUE_LIU"
      f<-f[which(f$NUM_VAR>=MinNumVarWin),]
    } else if (umethods[m] == "singlevar"){
      Pval<-"PVALUE"
    } else {
      Pval<-"PVALUE"
      f<-f[which(f$NUM_VAR>=MinNumVarWin),]
    }
    if (dim(f)[1]==0) {
    next
    }
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
    obs.max[i] <- round(max(log.obs.p),2)    
    exp.max[i] <- round(max(log.exp.p),2)
  }
  Ords<-order(obs.max)
  traits<-traits[Ords]
  lab.traits<-traits
  nby=10
  sequi<-c(seq(1,length(traits),by = nby),length(traits)) 
  mycols<-rep(mycols,length(sequi))
  for (interval in 2:length(sequi)){
    #interval=2
    if (interval != length(sequi)){
      sta=sequi[interval-1]
      endi=(sequi[interval]-1)
      print(paste0(sta,"-----",endi))
      #traits<-traits[sta:endi]
      #lab.traits<-traits
    } else {
      sta=sequi[interval-1]
      endi=(sequi[interval])
      #traits<-traits[sta:endi]
      #lab.traits<-traits
      print(paste0(sta,"-----",endi))
    }
## get complete distribution for maximum observed and expected
#Plotting limits if plot several traits in teh same plot
#m=1
p <- NULL
obs.max <- exp.max <- NULL
ptm <- proc.time()
for (i in sta:endi){
  #i=sta
  print(i)
  fn <- paste(PathF,"/Combined_ALL_INTERVAL_QCv1.0_inverse_norm_",Approaches[a],".",maf[mf],".",traits[i],".meta.",umethods[m],".results", sep="")  
  f <- read.csv(fn, sep="", stringsAsFactors=FALSE, header=T)
  f<-f[f$POOLED_ALT_AF>=0.001,]  
  #head(f)
    if(umethods[m] == "SKAT_"){
      Pval<-"PVALUE_LIU"
      f<-f[which(f$NUM_VAR>=MinNumVarWin),]
    } else if (umethods[m] == "singlevar"){
      Pval<-"PVALUE"
    } else {
      Pval<-"PVALUE"
      f<-f[which(f$NUM_VAR>=MinNumVarWin),]
    }
    if (dim(f)[1]==0) {
    next
    }
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
  obs.max[i] <- round(max(log.obs.p),2)    
  exp.max[i] <- round(max(log.exp.p),2)
}
obs.max<-obs.max[which(is.na(obs.max)==FALSE)]
exp.max<-exp.max[which(is.na(exp.max)==FALSE)]
xlim <- c(0,max(exp.max))
ylim <- c(0,max(obs.max)+0.3)
## for confidence interval
n <- ceiling(10^(max(exp.max)))
e <- -log10(1:n/n)
## Calculate 95% confidence intervals. The jth order statistic from a
## uniform(0,1) sample has a beta(j,n-j+1) distribution (Casella &
## Berger, 2002, 2nd edition, pg 230, Duxbury)
c95 <- sapply(1:n, function(i) qbeta(0.95,i,n-i+1))
c05 <- sapply(1:n, function(i) qbeta(0.05,i,n-i+1))
dat <- cbind(c05, c95)
#save(dat, file="/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/rare_variants_analysis/core_traits_exome-wide/analysis/p-values_confidence_interval_unique_gene_names.RData")
#load("/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/rare_variants_analysis/core_traits_exome-wide/analysis/p-values_confidence_interval_unique_gene_names.RData")  ## dat
#c05 <- dat[,1]
#c95 <- dat[,2]
#png("~/uk10k/manuscripts/main/for_resubmission/plots/qq-plot_exome-wide_MetaSKAT_unique_gene_names.png", pointsize=16, width=600, height=600)
#png("qq-plot_INTERVAL_Met_exome-wide_Trail_Multitraits_SKAT-O.png", pointsize=18, width=700, height=600)
#png(paste("/Users/lb17/Work/INTERVAL/SKAT/Skat_on_PCA/RES_All/qq-plot_INTERVAL_Met_exome-wide_Trail_Multitraits_SKAT-O_PCA_",umethods[m],".png",sep = ""), pointsize=18, width=700, height=600)
png(paste(PathF,"/QQplot/qq-plot_INTERVAL_Met_Exome_Multitraits_inverse_norm_",Approaches[a],"_Meth",umethods[m],"_MAF_",maf[mf],"_MetOrd_",sta,"_",endi,".png",sep = ""), pointsize=18, width=700, height=600)
par(mai=c(1.3,1.3,0.3,0.3))
#xlim <- c(0,5.01)   ## calculated from above
#ylim <- c(0,7.51)   ## calculated from above
plot(xlim, ylim, type="n", xlab=expression(paste(Expected -log[10], italic(p),sep="")), ylab=expression(paste(Observed -log[10], italic(p),sep="")), main=paste0(Approaches[a],"_",umethods[m],"_",maf[mf],"_Met",sta,"_",endi), las=1, axes=F, cex.lab=1.2)
polygon(c(e,exp.max), c(-log10(c95),exp.max), col="lightgrey", border="lightgrey")
polygon(c(e,exp.max), c(-log10(c05),exp.max), col="lightgrey", border="lightgrey")
lines(c(0,exp.max), c(0,exp.max), col="darkgray", lwd=1)
axis(1)
axis(2, las=1)      
## loop through traits
lambda <- NULL
for (i in sta:endi)  {
  #i=2
  print(paste0("plot",i))
  #fn <- paste("/lustre/scratch115/projects/int_wes_metabol/RunSKAT/results/Approaches[a]/skat/",traits[i],"/results_SKAT_chr13_",traits[i],"_midpoint.txt", sep="")
  fn <- paste(PathF,"/Combined_ALL_INTERVAL_QCv1.0_inverse_norm_",Approaches[a],".",maf[mf],".",traits[i],".meta.",umethods[m],".results", sep="")
  f <- read.csv(fn, sep="", stringsAsFactors=FALSE, header=T)
  f<-f[f$POOLED_ALT_AF>=0.001,]  
  #head(f)
    if(umethods[m] == "SKAT_"){
      Pval<-"PVALUE_LIU"
      f<-f[which(f$NUM_VAR>=MinNumVarWin),]
    } else if (umethods[m] == "singlevar"){
      Pval<-"PVALUE"
    } else {
      Pval<-"PVALUE"
      f<-f[which(f$NUM_VAR>=MinNumVarWin),]
    }
    if (dim(f)[1]==0) {
    next
    }
    p <- f[,Pval]
    p <- p[!is.na(p)]
    p <- as.numeric(gsub("<","",p))
  assign(paste("f",umethods[m],traits[i],sep="."),f)  
  ord <- order(p, na.last=NA)
  p <- p[ord]
  ## allocate min to p==0
  p[p==0] <- min(p)
  obs.p <- p
  log.obs.p <- -(log10(obs.p))
  exp.p <- c(1:length(obs.p))
  log.exp.p <- -(log10( (exp.p-0.5)/length(exp.p)))
  ## QQ-plot for each trait
  points(log.exp.p, log.obs.p, pch=20, cex=1.2, col=mycols[i])
  ind <- which(p<1e-500)
  if (length(ind)!=0)  {
    for (j in 1:length(ind))  {
      text(log.exp.p[ind[j]], log.obs.p[ind[j]], lab.traits[i], pos=2, col=mycols[i])
    }
  }
  ## calculate lambda
  X <- qchisq(1 - obs.p, df=1)
  lambda[i] <- round(median(X) / qchisq(0.5, df=1),2)
}
 lambda<-lambda[which(is.na(lambda)==FALSE)]
 traitst<-traits[sta:endi]
 lab.traitst<-traitst
 legend.text <- paste(lab.traitst," ", lambda, sep="")
 cat("Trait lambda",file=paste0(PathF,"/QQplot/lambda_MetOrd_",Approaches[a],"_",umethods[m],"_",sta,"_",endi,".txt"),sep="\n")
 #legend.text<-paste0("M323223"," 1.1")
 cat(legend.text,file=paste0(PathF,"/QQplot/lambda_MetOrd_",Approaches[a],"_",umethods[m],"_",sta,"_",endi,".txt"),sep="\n", append=TRUE)
 print(legend.text)
 legend(0, ylim[2], legend=legend.text, fill=mycols, box.col="white", border=mycols, cex=0.6)
dev.off()
Tt<-proc.time() - ptm
print(Tt)
}
}

Ft<-proc.time() - ptmA
print(Ft)

