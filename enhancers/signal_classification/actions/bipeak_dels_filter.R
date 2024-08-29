suppressMessages(library(tidyverse))
suppressMessages(library(scales))
suppressMessages(library(tools))

args <- commandArgs(trailingOnly=TRUE)

bipeak_file <- args[1]
dels_peak_file <- args[2]

bipeaks <- read.delim(bipeak_file,sep="\t",header=TRUE,check.names=FALSE)
dels_peaks <- read.delim(dels_peak_file,sep="\t",header=FALSE,check.names=FALSE)


cmd1 <- paste0("dels_peaks_unique <- unique(dels_peaks$V4)")
write(cmd1,stderr())
eval(parse(text=cmd1))


cmd2 <- paste0("erna_bipeaks <- bipeaks[which(bipeaks$peak1_id %in% dels_peaks_unique & bipeaks$peak2_id %in% dels_peaks_unique),]")
write(cmd2,stderr())
eval(parse(text=cmd2))


write.table(x=erna_bipeaks,file=paste(file_path_sans_ext(bipeak_file),".dELS_intersect.bed",sep=""),sep="\t",quote=FALSE,na="NA",row.names=FALSE,col.names=TRUE)
