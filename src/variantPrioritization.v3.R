#!/usr/bin/env Rscript
#
#
# Author: Alba Sanchis Juan (as2635@cam.ac.uk)
# Project: NGC
# Description: This script performs the variant prioritization
# of the filtered by frequency tab file
# v2: new prioritization method and compatible with Cellbase

#####################
##Loading libraries##
#####################

#loc3_3="/home/aak64/R/x86_64-pc-linux-gnu-library/3.3"

library("optparse"
#	, lib.loc=loc3_3
	)
library("plyr"
#	, lib.loc=loc3_3
	)
library("data.table"
#	, lib.loc=loc3_3
	)
library("parallel")
library("tools")

option_list = list(
  make_option(c("-v", "--vars"), type="character", default=NULL,
              help="variants file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.filt.txt",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-l", "--lengths"), type="character", default=NULL,
              help="Ensembl transcript lengths", metavar="character"),
  make_option(c("-c", "--cores"), type="character", default=1,
              help="number of cores used", metavar="integer")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$vars) || is.null(opt$lengths)){
  print_help(opt_parser)
  stop("At least two arguments must be supplied (variants file and transcript lengths).n", call.=FALSE)
}

#####################
##Parsing arguments##
#####################
variants_path <- opt$vars
out_path <- opt$out
length_path <- opt$lengths
n_cores <- opt$cores

#################
##Reading files##
#################
raw_variants <- fread(variants_path, header = T, stringsAsFactors = F, quote="")
raw_variants <- data.frame(raw_variants)

ENST <- fread(length_path, skip = 1, col.names = c("ENSG", "ANNO_Feature", "CDS_length", "ANNO_Symbol", "transcript_type", "gene_type"), stringsAsFactors = F)
ENST <- data.frame(ENST)

#############
##Functions##
#############
transcriptPrior <- function(df){
#  PAR1start <- 60001
#  PAR1stop <- 2699520
#  PAR2start <- 154931044
#  PAR2stop <- 155270560

#  df["PAR"] <- ifelse((df$CHROM=="X" & df$POS>PAR1start & df$POS<PAR1stop)
#                      | (df$CHROM=="X" & df$POS>PAR2start & df$POS<PAR2stop), T, F)
#  df["AUT_OR_PAR"] <- ifelse(df$CHROM=="MT"
#                             | (df$CHROM=="X" & df$PAR==F)
#                             | df$CHROM=="Y", F, T)
  df["VAR"] <- paste(df$CHROM, df$POS, df$REF, df$ALT, df$ANNO_Gene)
  df["Effect_impact_rankings"] <- ifelse(df$ANNO_IMPACT=="HIGH", 1,
                                         ifelse(df$ANNO_IMPACT=="MODERATE", 2,
                                                ifelse(grepl("splice_donor_5th_base_variant",df$ANNO_SpliceRegion), 3,
                                                       ifelse(grepl("splice_polypyrimidine_tract_variant",df$ANNO_SpliceRegion), 4,
                                                              ifelse(df$ANNO_IMPACT=="LOW", 5, 6)))))
  df <- join(df, ENST, by = "ANNO_Feature")
  df["biotype_rank"] <- ifelse(df$ANNO_BIOTYPE=="protein_coding", 1, 2)
  df["Canonical_rank"] <- ifelse(df$ANNO_CANONICAL=="YES", 1, 2)

  df_prior <- do.call("rbind", mclapply(unique(df$VAR), function(variant){
    x <- df[df$VAR == variant,]
    x <- x[with(x, order(biotype_rank, Effect_impact_rankings, Canonical_rank, -CDS_length)), ]
    x <- x[1,]
    x
  }, mc.cores=n_cores))

  return(df_prior[,!names(df_prior) %in% c("Effect_impact_rankings"
                                           , "biotype_rank"
                                           , "Canonical_rank"
                                           , "CDS_length")])
}

##########################
##Variant prioritization##
##########################
variants <- transcriptPrior(raw_variants)

###############
##Saving file##
###############
print(dim(variants))
root_file = basename(out_path)
path_dir = dirname(out_path)

#save(variants,file=paste(path_dir,"/",root_file,".dat",sep=""))
write.table(variants,file=out_path,sep ="\t",quote=F,row.names=F)
