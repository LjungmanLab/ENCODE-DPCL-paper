suppressMessages(library(tidyverse))
suppressMessages(library(scales))
suppressMessages(library(tools))

args <- commandArgs(trailingOnly=TRUE)

all_peak_file <- args[1]
dels_bipeak_file <- args[2]
dels_unipeak_file <- args[3]

all_peaks <- read.delim(all_peak_file,sep="\t",header=TRUE,check.names=FALSE)
dels_bipeaks <- read.delim(dels_bipeak_file,sep="\t",header=TRUE,check.names=FALSE)
dels_unipeaks <- read.delim(dels_unipeak_file,sep="\t",header=TRUE,check.names=FALSE)


cmd1 <- paste0("dels_bipeaks_unique <- unique(dels_bipeaks$bipeak_id)")
write(cmd1,stderr())
eval(parse(text=cmd1))

cmd2 <- paste0("dels_unipeaks_unique <- unique(dels_unipeaks$unipeak_id)")
write(cmd2,stderr())
eval(parse(text=cmd2))

cmd3 <- paste0("erna_peaks <- all_peaks[which(all_peaks$class_id %in% dels_bipeaks_unique | all_peaks$class_id %in% dels_unipeaks_unique),]")
write(cmd3,stderr())
eval(parse(text=cmd3))


write.table(x=erna_peaks,file=paste(file_path_sans_ext(all_peak_file),".dELS_intersect.bed",sep=""),sep="\t",quote=FALSE,na="NA",row.names=FALSE,col.names=TRUE)
