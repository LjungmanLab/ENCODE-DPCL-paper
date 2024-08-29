suppressMessages(library(tidyverse))

# # read in file as trailing argument to use script outside of count_bed (or can cat file and pipe to Rscript)
# args <- commandArgs(trailingOnly=TRUE)
# bed_file <- args[1]
# bed_df <- read.delim(bed_file,sep="\t",header=FALSE)

bed_df <- suppressWarnings(read.delim("/dev/stdin",sep="\t",header=FALSE))

# chroms <- c("chr1")

chroms <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13",
            "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")

cmd1 <- paste0("filtered_df = bed_df[bed_df$V1 %in% chroms,]")
# write(cmd1,stderr())
eval(parse(text=cmd1))

tryCatch(write.table(x=filtered_df,file="",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE),error=function(e) if(!grepl("ignoring SIGPIPE signal",e$message))stop(e))
