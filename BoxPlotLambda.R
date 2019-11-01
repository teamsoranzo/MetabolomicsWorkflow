

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

Approaches <- c("naive","functional","LoF")
print(Approaches[a])
#a=1

Path="/lustre/scratch115/realdata/mdt0/projects/int_wes_metabol/raremetal/AllMetabolites/raremetal_inverse_norm_AllMet/"
filess <-list.files(path=paste0(Path,Approaches[a],"/QQplot"), pattern="*ALL.txt")
files <- list.files(path=paste0(Path,Approaches[a],"/QQplot"), pattern="*ALL.txt", full.names=T)
Met <- read.table(files[1], header=FALSE, stringsAsFactors=F)[,1]
df <- do.call(cbind,lapply(files,function(fn){read.table(fn,header=FALSE, stringsAsFactors=F)[,2]}))
rownames(df) <- Met
colnames(df)<-unlist(lapply(filess, function(x) strsplit(x, "_")[[1]][4]))

pdf(paste0(Path,Approaches[a],"/QQplot/Boxplotlambda_",Approaches[a],".pdf"))
boxplot(df, col=topo.colors(4), main=Approaches[a])
abline(h=1, col="blue")
dev.off()



