#!/usr/bin/Rscript

library(plyr)
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
c <- as.numeric(ID)

#c=22
Chr<-c(1:22,"X")
print(Chr[c])

#R<-read.table("/lustre/scratch115/realdata/mdt0/projects/int_wes_metabol/raremetal/AllMetabolites/raremetal_inverse_norm_AllMet/Singlevar/ALL_Metabolites_INTERVAL_QCv1.0_inverse_norm_SingleVar_Significant1.6e7_MAF0.001above.results", stringsAsFactors=F , header=T)
Ra<-read.table("/lustre/scratch115/realdata/mdt0/projects/int_wes_metabol/raremetal/AllMetabolites/raremetal_inverse_norm_AllMet/Singlevar/ALL_Metabolites_INTERVAL_QCv1.0_inverse_norm_SingleVar_1e-4_CommonAndRare.results", stringsAsFactors=F , header=T)
W<-read.table(paste0("/lustre/scratch115/projects/int_wes_metabol/VEP_85_GENCODE24_GFF3/Chr",Chr[c],"/naive/Chr",Chr[c],"_windows_naive_1.txt"), stringsAsFactors=F, header=F)

R<-Ra

HR<-colnames(R)
HW<-c("CHROM", "START", "END", "REGION", "WINDOW", "GENE", "NVARinWIN")
H<-c(HW,HR)
W<-W[order(W[,1], W[,2]),]
R<-subset(R, R[,1]==Chr[c])
R<-R[order(R[,1], R[,2]),]

#count<-0
#for(i in 1:dim(W)[1]){
#	if(W[i,3]>W[i+1,2]){
#	count<-count+1
#	print(paste0("N:",count," Overlap",W[i,2],"-",W[i,3]," | ",W[i+1,2],"-",W[i+1,3]))
#	}
#}



WR<-NULL
for (y in 1:dim(R)[1]){
	id<-0
	for (i in 1:dim(W)[1]){
		int<-seq(W[i,2],W[i,3])
		if(R[y,2] %in% int){
		WR[[y]]<-data.frame(W[i,],R[y,],stringsAsFactors=F)
		print(paste0("GENE: ",y))
		#print(W[i,],R[y,])
		id<-1
		}
	
	}
		if (id==0){
			WR[[y]]<-data.frame(matrix(rep(NA,7),1,7),R[y,],stringsAsFactors=F)
			
	}
}


WR<-as.data.frame(matrix(unlist(WR), byrow=T,nrow=length(WR)))

#WR<-rbind.fill(lapply(WR, function(f) {as.data.frame(Filter(Negate(is.null), f))}))
colnames(WR)<-H
#write.table(WR, paste0("/lustre/scratch115/realdata/mdt0/projects/int_wes_metabol/raremetal/AllMetabolites/raremetal_inverse_norm_AllMet/Singlevar/Chr",Chr[c],"_GenesAnnot_ALL_Metabolites_INTERVAL_QCv1.0_inverse_norm_SingleVar_Significant1.6e7_MAF0.001above.results"), quote=F, col.names=T, row.names=F)

write.table(WR, paste0("/lustre/scratch115/realdata/mdt0/projects/int_wes_metabol/raremetal/AllMetabolites/raremetal_inverse_norm_AllMet/Singlevar/Chr",Chr[c],"_GenesAnnot_ALL_Metabolites_INTERVAL_QCv1.0_inverse_norm_SingleVar_1e-4_CommonANDrare.results"), quote=F, col.names=T, row.names=F)

#tmp <- with(WR, by(WR$PVALUE, WR$GENE, function(x) min(x)))
#sapply(tmp, min)
#WRMinP<-t(sapply(with(WR, by(WR, WR$GENE, function(x) x[which.min(x$PVALUE),])),function(x) x[which.min(x$PVALUE),]))

#WRMinP<-ddply(WR, "GENE", function(x) x[which.min(x$PVALUE),])

#write.table(WRMinP, paste0("/lustre/scratch115/realdata/mdt0/projects/int_wes_metabol/raremetal/AllMetabolites/raremetal_inverse_norm_AllMet/Singlevar/Chr",Chr[c],"_GenesAnnot_MinPVALUE_ALL_Metabolites_INTERVAL_QCv1.0_inverse_norm_SingleVar_Significant1.6e7_MAF0.001above.results"), quote=F, col.names=T, row.names=F)

