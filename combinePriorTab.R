args=(commandArgs(TRUE))
library("plyr"
#        ,lib.loc=loc3_3
)
library("readbulk"
#         ,lib.loc=loc3_3
)

inp_dir = args[[1]]
out_file = args[[2]]


raw_data <- read_bulk(directory=inp_dir,subdirectories=TRUE, extension=".prior.tab",stringsAsFactor=FALSE, fun=read.delim)
raw_data$CHROM <- factor(raw_data$CHROM, levels=c(seq(1,22),"X","Y","MT"))
raw_data$POS <- as.numeric(raw_data$POS)
data <- raw_data[order(raw_data$CHROM, raw_data$POS),]
write.table(data,file=paste(inp_dir,"/",out_file,sep=""),sep="\t",row.names=F,quote=F)

