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
#a=2
Approaches <- c("naive","functional", "LoF")
print(Approaches[a])

path1="/lustre/scratch115/projects/int_wes_metabol/raremetal/AllMetabolites/raremetal_inverse_norm_AllMet/Singlevar"
path2=paste0("/lustre/scratch115/projects/int_wes_metabol/raremetal/AllMetabolites/raremetal_inverse_norm_AllMet/",Approaches[a],"/AllMerged")
path3="/lustre/scratch115/realdata/mdt0/projects/int_wes_metabol/Metabolon/CleanedData_160414"


## Combined results only significants MAF > 0.001 SINGLEVAR
S<-read.table(paste0(path1,"/ALL_Metabolites_INTERVAL_QCv1.0_inverse_norm_SingleVar_Significant1.6e7_MAF0.001above_GenesAnnot.results"), header=TRUE, stringsAsFactors=F)
## Combined results only significants all variants SINGLEVAR
Sa<-read.table(paste0(path1,"/ALL_Metabolites_INTERVAL_QCv1.0_inverse_norm_SingleVar_Significant1.6e7_CommonANDrare_GenesAnnot.results"), header=TRUE, stringsAsFactors=F)
## Combined results only significants MAF > 0.001  GENE-BASED
G<-read.table(paste0(path2,"/ALL_Metabolites_INTERVAL_QCv1.0_inverse_norm_Genebased_",Approaches[a],".MAF0.001_Significant.results"), header=TRUE, stringsAsFactors=F)
## Metabolite summary table
Mi<-read.delim(paste0(path3,"/Metabolite_summary_data_all.txt"), header=TRUE, stringsAsFactors=F)
for  (i in 1:dim(Mi)[1]){
	if(nchar(Mi$COMP_ID[i])==2){
	Mi$COMP_ID[i]<-paste0("M000",Mi$COMP_ID[i])
	}
	if(nchar(Mi$COMP_ID[i])==3){
        Mi$COMP_ID[i]<-paste0("M00",Mi$COMP_ID[i])
        }
	if(nchar(Mi$COMP_ID[i])==4){
        Mi$COMP_ID[i]<-paste0("M0",Mi$COMP_ID[i])
        }
	if(nchar(Mi$COMP_ID[i])==5){
        Mi$COMP_ID[i]<-paste0("M",Mi$COMP_ID[i])
        }
}

##MERGING

#### Based on Varinats
# Splits colmuns per variants in gene based and then merging with single var.
library(splitstackshape)
GVar<-data.frame(cSplit(G, splitCols=c('VARs','MAFs','SINGLEVAR_EFFECTs','SINGLEVAR_PVALUEs'),sep=c(';',',',',',','), 'long', drop = FALSE), stringsAsFactors=F)
GVars<-data.frame(na.omit(GVar),stringsAsFactors=F)
GVars$VARs<-as.character(GVars$VARs)
str(GVars)
##Create a unique identifier for variants in the Single variant dataset
S$VARsS<-paste(S$CHROM.1,S$POS,S$REF,S$ALT, sep=":")
Sa$VARsS<-paste(Sa$CHROM.1,Sa$POS,Sa$REF,Sa$ALT, sep=":")

########################### CHANGED FOR UNIQUE VARS

### Include Biochemical
G$BIOCHEMICAL<-rep(NA,dim(G)[1])
G$SUPER_PATHWAY<-rep(NA,dim(G)[1])
G$SUB_PATHWAY<-rep(NA,dim(G)[1])

for (i in 1:dim(G)[1]){
        G$BIOCHEMICAL[i]<-Mi[which(Mi$COMP_ID== G$METABOLITE[i]),2]
        G$SUPER_PATHWAY[i]<-Mi[which(Mi$COMP_ID== G$METABOLITE[i]),3]
        G$SUB_PATHWAY[i]<-Mi[which(Mi$COMP_ID== G$METABOLITE[i]),4]
}


for (i in 1:dim(Sa)[1]){
        Sa$BIOCHEMICAL[i]<-Mi[which(Mi$COMP_ID== Sa$METABOLITE[i]),2]
        Sa$SUPER_PATHWAY[i]<-Mi[which(Mi$COMP_ID== Sa$METABOLITE[i]),3]
        Sa$SUB_PATHWAY[i]<-Mi[which(Mi$COMP_ID== Sa$METABOLITE[i]),4]
}

GVars$BIOCHEMICAL<-rep(NA,dim(GVars)[1])
GVars$SUPER_PATHWAY<-rep(NA,dim(GVars)[1])
GVars$SUB_PATHWAY<-rep(NA,dim(GVars)[1])

for (i in 1:dim(GVars)[1]){
        GVars$BIOCHEMICAL[i]<-Mi[which(Mi$COMP_ID== GVars$METABOLITE[i]),2]
        GVars$SUPER_PATHWAY[i]<-Mi[which(Mi$COMP_ID== GVars$METABOLITE[i]),3]
        GVars$SUB_PATHWAY[i]<-Mi[which(Mi$COMP_ID== GVars$METABOLITE[i]),4]
}




## Identify variables that makes the Variants variable unique
length(unique(paste0(GVars$VARs, GVars$METABOLITE, GVars$TEST,GVars$WINDOW_NAME)))/length(paste0(GVars$VARs, GVars$METABOLITE, GVars$TEST,GVars$WINDOW_NAME))
#1
length(unique(paste0(Sa$VARs,Sa$METABOLITE)))/length(paste0(Sa$VARs,Sa$METABOLITE))
#1
length(unique(paste0(S$VARs,S$METABOLITE)))/length(paste0(S$VARs,S$METABOLITE))
#1
## aggregate the variables to have a unique variant table - SINGLEVAR
SaUnique<-as.data.frame(aggregate(cbind(Sa$CHROM, Sa$START, Sa$END, Sa$REGION, Sa$WINDOW, Sa$GENE, Sa$NVARinWIN, Sa$POS, Sa$REF, Sa$ALT, Sa$N, Sa$POOLED_ALT_AF, Sa$DIRECTION_BY_STUDY, Sa$EFFECT_SIZE, Sa$EFFECT_SIZE_SD, Sa$H2, Sa$PVALUE, Sa$METABOLITE, Sa$BIOCHEMICAL, Sa$SUPER_PATHWAY, Sa$SUB_PATHWAY),list(Sa$VARsS), paste,  collapse="||"), stringsAsFactors=F)
apply(SaUnique, 2, function(y) length(which(sapply(strsplit(y, split="\\|\\|"), function(x) length(unique(x)) == 1) == FALSE)))

SaUniqueU<-apply(SaUnique[,-c(12:22)],2,function(y) sapply(strsplit(y, split="\\|\\|"), function(x) x[1]))
SaUniqueU<-as.data.frame(SaUniqueU, stringsAsFactors=F)
SaUniqueU<-cbind(SaUniqueU,SaUnique[,c(12:22)])
colnames(SaUniqueU)<-c(colnames(Sa)[20],colnames(Sa)[c(1:7,9:19,21:23)])


## aggregate the variable to have a unique variant table - GENEBASED
GVarsUnique<-as.data.frame(aggregate(cbind(GVars$WINDOW_NAME, GVars$CHR, GVars$START, GVars$END, GVars$REGION, GVars$GENE_NAME, GVars$NVARWin, GVars$NUM_VAR, GVars$MAFs, GVars$SINGLEVAR_EFFECTs, GVars$SINGLEVAR_PVALUEs, GVars$AVG_AF, GVars$MIN_AF, GVars$MAX_AF, GVars$EFFECT_SIZE, GVars$PVALUE, GVars$METABOLITE, GVars$TEST,  GVars$BIOCHEMICAL, GVars$SUPER_PATHWAY, GVars$SUB_PATHWAY),list(GVars$VARs), paste,  collapse="||"), stringsAsFactors=F)
# check  if Windows are unique and there are no variants ovelappin gin different windows
apply(GVarsUnique, 2, function(y) length(which(sapply(strsplit(y, split="\\|\\|"), function(x) length(unique(x)) == 1) == FALSE)))

#Group.1      V1      V2      V3      V4      V5      V6      V7      V8      V9
#      0      19       0      13       6      19      19      19     259     207
#    V10     V11     V12     V13     V14     V15     V16     V17     V18
#    286     286     376     207     376     841     900     286     834

GVarsUniqueU<-apply(GVarsUnique[,c(1,3)],2,function(y) sapply(strsplit(y, split="\\|\\|"), function(x) x[1]))
GVarsUniqueU<-as.data.frame(GVarsUniqueU, stringsAsFactors=F)
GVarsUniqueU<-cbind(GVarsUniqueU,GVarsUnique[,-c(1,3)])
colnames(GVarsUniqueU)<-c(colnames(GVars)[9],colnames(GVars)[2], colnames(GVars)[c(1,3:8,10:22)])


## Merge Gene based with singlevar with MAF > 0.001 -- NO OVERLAP (appending)
M_SaG<-merge(SaUniqueU,GVarsUniqueU,by.x=1 , by.y=1 , all=TRUE, suffixes=c(".SingleVar",".GeneBased"))


write.table(M_SaG, file=paste0("/lustre/scratch115/projects/int_wes_metabol/raremetal/AllMetabolites/raremetal_inverse_norm_AllMet/Table_UniqueVar_ALL_Metabolites_INTERVAL_QCv1.0_inverse_norm_Genebased_",Approaches[a],".MAF0.001_Significant_MergedSingleVar_CommonANDrare_GenesAnnot.txt"), quote=F, row.names=FALSE, col.names=TRUE, sep="\t")


### OLD
## Merge Gene based with singlevar with MAF > 0.001 -- NO OVERLAP (appending)
# M_SG<-merge(S,GVars,by.x=20 , by.y=9 , all=TRUE, suffixes=c(".SingleVar",".GeneBased"))
## Merge Gene based with singlevar with Common and rare -- Possible overlap
# M_SaG<-merge(Sa,GVars,by.x=20 , by.y=9 , all=TRUE, suffixes=c(".SingleVar",".GeneBased"))


#M_SG$BIOCHEMICAL.SingleVar<-rep(NA,dim(M_SG)[1])
#M_SG$SUPER_PATHWAY.SingleVar<-rep(NA,dim(M_SG)[1])
#M_SG$SUB_PATHWAY.SingleVar<-rep(NA,dim(M_SG)[1])
#M_SG$BIOCHEMICAL.GeneBased<-rep(NA,dim(M_SG)[1])
#M_SG$SUPER_PATHWAY.GeneBased<-rep(NA,dim(M_SG)[1])
#M_SG$SUB_PATHWAY.GeneBased<-rep(NA,dim(M_SG)[1])

#for (i in 1:dim(M_SG)[1]){
#	if(!is.na(M_SG$METABOLITE.SingleVar[i])){
#	M_SG$BIOCHEMICAL.SingleVar[i]<-Mi[which(Mi$COMP_ID== M_SG$METABOLITE.SingleVar[i]),2]
#	M_SG$SUPER_PATHWAY.SingleVar[i]<-Mi[which(Mi$COMP_ID== M_SG$METABOLITE.SingleVar[i]),3]
#	M_SG$SUB_PATHWAY.SingleVar[i]<-Mi[which(Mi$COMP_ID== M_SG$METABOLITE.SingleVar[i]),4]
#	}
#	if(!is.na(M_SG$METABOLITE.GeneBased[i])){
#	M_SG$BIOCHEMICAL.GeneBased[i]<-Mi[which(Mi$COMP_ID== M_SG$METABOLITE.GeneBased[i]),2]
#        M_SG$SUPER_PATHWAY.GeneBased[i]<-Mi[which(Mi$COMP_ID== M_SG$METABOLITE.GeneBased[i]),3]
#        M_SG$SUB_PATHWAY.GeneBased[i]<-Mi[which(Mi$COMP_ID== M_SG$METABOLITE.GeneBased[i]),4]
#	}
#}


#M_SaG$BIOCHEMICAL.SingleVar<-rep(NA,dim(M_SaG)[1])
#M_SaG$SUPER_PATHWAY.SingleVar<-rep(NA,dim(M_SaG)[1])
#M_SaG$SUB_PATHWAY.SingleVar<-rep(NA,dim(M_SaG)[1])
#M_SaG$BIOCHEMICAL.GeneBased<-rep(NA,dim(M_SaG)[1])
#M_SaG$SUPER_PATHWAY.GeneBased<-rep(NA,dim(M_SaG)[1])
#M_SaG$SUB_PATHWAY.GeneBased<-rep(NA,dim(M_SaG)[1])

#for (i in 1:dim(M_SaG)[1]){
#        if(!is.na(M_SaG$METABOLITE.SingleVar[i])){
#        M_SaG$BIOCHEMICAL.SingleVar[i]<-Mi[which(Mi$COMP_ID== M_SaG$METABOLITE.SingleVar[i]),2]
#        M_SaG$SUPER_PATHWAY.SingleVar[i]<-Mi[which(Mi$COMP_ID== M_SaG$METABOLITE.SingleVar[i]),3]
#        M_SaG$SUB_PATHWAY.SingleVar[i]<-Mi[which(Mi$COMP_ID== M_SaG$METABOLITE.SingleVar[i]),4]
#        }
#        if(!is.na(M_SaG$METABOLITE.GeneBased[i])){
#        M_SaG$BIOCHEMICAL.GeneBased[i]<-Mi[which(Mi$COMP_ID== M_SaG$METABOLITE.GeneBased[i]),2]
#        M_SaG$SUPER_PATHWAY.GeneBased[i]<-Mi[which(Mi$COMP_ID== M_SaG$METABOLITE.GeneBased[i]),3]
#        M_SaG$SUB_PATHWAY.GeneBased[i]<-Mi[which(Mi$COMP_ID== M_SaG$METABOLITE.GeneBased[i]),4]
#        }
#}




#VarsUnique<-as.data.frame(aggregate(cbind(M_SG$CHROM ,M_SG$START.SingleVar ,M_SG$END.SingleVar ,M_SG$REGION.SingleVar ,M_SG$WINDOW ,M_SG$GENE ,M_SG$NVARinWIN ,M_SG$CHROM.1 ,M_SG$POS ,M_SG$REF ,M_SG$ALT ,M_SG$N ,M_SG$POOLED_ALT_AF ,M_SG$DIRECTION_BY_STUDY ,M_SG$EFFECT_SIZE.SingleVar ,M_SG$EFFECT_SIZE_SD ,M_SG$H2 ,M_SG$PVALUE.SingleVar ,M_SG$METABOLITE.SingleVar ,M_SG$WINDOW_NAME ,M_SG$REGION.GeneBased ,M_SG$GENE_NAME ,M_SG$NVARWin ,M_SG$NUM_VAR ,M_SG$MAFs ,M_SG$SINGLEVAR_EFFECTs ,M_SG$SINGLEVAR_PVALUEs ,M_SG$AVG_AF ,M_SG$MIN_AF ,M_SG$MAX_AF ,M_SG$EFFECT_SIZE.GeneBased ,M_SG$PVALUE.GeneBased ,M_SG$METABOLITE.GeneBased ,M_SG$TEST),list(M_SG$VARsS) , paste,  collapse="||"), stringsAsFactors=F)


#write.table(M_SG, file=paste0("/lustre/scratch115/projects/int_wes_metabol/raremetal/AllMetabolites/raremetal_inverse_norm_AllMet/Table_ALL_Metabolites_INTERVAL_QCv1.0_inverse_norm_Genebased_",Approaches[a],".MAF0.001_Significant_MergedSingleVar_MAF0.001above_GenesAnnot.txt"), quote=F, row.names=FALSE, col.names=TRUE)

#write.table(M_SaG, file=paste0("/lustre/scratch115/projects/int_wes_metabol/raremetal/AllMetabolites/raremetal_inverse_norm_AllMet/Table_ALL_Metabolites_INTERVAL_QCv1.0_inverse_norm_Genebased_",Approaches[a],".MAF0.001_Significant_MergedSingleVar_CommonANDrare_GenesAnnot.txt"), quote=F, row.names=FALSE, col.names=TRUE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### Based on Genes
#Use only Gene based and collapse columns with windows into genes


AGGALL<-as.data.frame(aggregate(cbind(G$WINDOW_NAME,G$CHR,G$START,G$END,G$REGION,G$NVARWin,G$NUM_VAR,G$VARs,G$MAFs,G$SINGLEVAR_EFFECTs,G$SINGLEVAR_PVALUEs,G$AVG_AF,G$MIN_AF,G$MAX_AF,G$EFFECT_SIZE,G$PVALUE,G$METABOLITE, G$TEST, G$BIOCHEMICAL,G$SUPER_PATHWAY,G$SUB_PATHWAY),list(G$GENE_NAME), paste,  collapse="||"), stringsAsFactors=F)

apply(AGGALL, 2, function(y) length(which(sapply(strsplit(y, split="\\|\\|"), function(x) length(unique(x)) == 1) == FALSE)))

AGGALLU<-apply(AGGALL[,c(1,3)],2,function(y) sapply(strsplit(y, split="\\|\\|"), function(x) x[1]))
AGGALLU<-as.data.frame(AGGALLU, stringsAsFactors=F)
AGGALLU<-cbind(AGGALLU,AGGALL[,-c(1,3)])
colnames(AGGALLU)<-c(colnames(G)[6],colnames(GVars)[c(1:5,7:22)])


write.table(AGGALLU, file=paste0("/lustre/scratch115/projects/int_wes_metabol/raremetal/AllMetabolites/raremetal_inverse_norm_AllMet/Table_UniqueGene_ALL_Metabolites_INTERVAL_QCv1.0_inverse_norm_Genebased_",Approaches[a],".MAF0.001_Significant.txt"), quote=F, row.names=FALSE, col.names=TRUE, sep="\t")

#OLD

#AGGALL<-as.data.frame(aggregate(cbind(G$WINDOW_NAME,G$CHR,G$START,G$END,G$REGION,G$NVARWin,G$NUM_VAR,G$VARs,G$MAFs,G$SINGLEVAR_EFFECTs,G$SINGLEVAR_PVALUEs,G$AVG_AF,G$MIN_AF,G$MAX_AF,G$EFFECT_SIZE,G$PVALUE),list(G$GENE_NAME,G$METABOLITE, G$TEST), paste,  collapse="||"), stringsAsFactors=F)

#AGGALL[,5]<-apply(AGGALL,1,function(x) {unlist(strsplit(x[5], "\\|\\|"))[1]})
#AGGALL[,7]<-apply(AGGALL,1,function(x) {max(unlist(strsplit(x[7], "\\|\\|")))})
#AGGALL[,6]<-apply(AGGALL,1,function(x) {min(unlist(strsplit(x[6], "\\|\\|")))})

#colnames(AGGALL)<-colnames(G)[c(6,18:19,1:5,7:17)]

#AGGALL$BIOCHEMICAL<-rep(NA,dim(AGGALL)[1])
#AGGALL$SUPER_PATHWAY<-rep(NA,dim(AGGALL)[1])
#AGGALL$SUB_PATHWAY<-rep(NA,dim(AGGALL)[1])

#for (i in 1:dim(AGGALL)[1]){
#        if(!is.na(AGGALL$METABOLITE[i])){
#        AGGALL$BIOCHEMICAL[i]<-Mi[which(Mi$COMP_ID== AGGALL$METABOLITE[i]),2]
#        AGGALL$SUPER_PATHWAY[i]<-Mi[which(Mi$COMP_ID== AGGALL$METABOLITE[i]),3]
#        AGGALL$SUB_PATHWAY[i]<-Mi[which(Mi$COMP_ID== AGGALL$METABOLITE[i]),4]
#        }
#}


#write.table(AGGALL, file=paste0("/lustre/scratch115/projects/int_wes_metabol/raremetal/AllMetabolites/raremetal_inverse_norm_AllMet/Table_ALL_Metabolites_INTERVAL_QCv1.0_inverse_norm_Genebased_",Approaches[a],".MAF0.001_Significant.txt"), quote=F, row.names=FALSE, col.names=TRUE)


#AGGALL_UNIQG<-as.data.frame(aggregate(cbind(G$WINDOW_NAME,G$CHR,G$START,G$END,G$REGION,G$NVARWin,G$NUM_VAR,G$VARs,G$MAFs,G$SINGLEVAR_EFFECTs,G$SINGLEVAR_PVALUEs,G$AVG_AF,G$MIN_AF,G$MAX_AF,G$EFFECT_SIZE,G$PVALUE,G$METABOLITE,G$TEST),list(G$GENE_NAME), paste,  collapse="||"), striingsAsFactors=F)

#AGGALL_UNIQG[,3]<-apply(AGGALL_UNIQG,1,function(x) {unlist(strsplit(x[3], "\\|\\|"))[1]})
#AGGALL_UNIQG[,4]<-apply(AGGALL_UNIQG,1,function(x) {max(unlist(strsplit(x[4], "\\|\\|")))})
#AGGALL_UNIQG[,5]<-apply(AGGALL_UNIQG,1,function(x) {min(unlist(strsplit(x[4], "\\|\\|")))})

#colnames(AGGALL_UNIQG)<-colnames(G)[c(6,1:5,7:19)]





