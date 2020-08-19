#!/usr/bin/env Rscript
#
#author: Alba Sanchis Juan (as2635@cam.ac.uk)
#
# last edited: 11 Sep 2019 - CEF

#############
##Functions##
#############

##freq filtering functions##

dominant_filter <- function(df){
  if(nrow(df) >=1){
      x <- subset(df, ((is.na(GNOMADg_AC) | GNOMADg_AC <= 5)
		    & (is.na(GNOMADgV3_AC) | GNOMADgV3_AC <= 5)
		    & (is.na(GNOMADe_AC) | GNOMADe_AC <= 5)
                    & (is.na(EXAC_AC) | EXAC_AC <= 5) 
		    & (is.na(TOPMED_AC) | TOPMED_AC <= 10)
		    & (is.na(grch37_NGC_AC) | grch37_NGC_AC <= 5)
		    & (is.na(NGC_AC) | NGC_AC <= 5)
                    ))
  }else{
      x <- df
  }
  return(x)
}

recessive_filter <- function(df){
  if(nrow(df) >=1){
      x <- subset(df, ((is.na(GNOMADg_AF) | GNOMADg_AF <= 0.01)
		    & (is.na(GNOMADgV3_AF) | GNOMADgV3_AF <= 0.01)
		    & (is.na(GNOMADe_AF) | GNOMADe_AF <= 0.01)
                    & (is.na(EXAC_AF) | EXAC_AF <= 0.01)
                    & (is.na(TOPMED_AF) | TOPMED_AF <= 0.01)
                    & (is.na(grch37_NGC_AF) | grch37_NGC_AF <= 0.10)
                    & (is.na(NGC_AF) | NGC_AF <= 0.10)
		    & (is.na(EXAC_AC_Hom) | EXAC_AC_Hom <= 5)
		    & (is.na(GNOMADe_nhomalt) | GNOMADe_nhomalt <= 5)
		    & (is.na(GNOMADg_nhomalt) | GNOMADg_nhomalt <= 5)
		    & (is.na(GNOMADgV3_nhomalt) | GNOMADgV3_nhomalt <= 5)
		    & (is.na(TOPMED_Hom) | TOPMED_Hom <= 10)
                    ))
  }else{
      x <- df
  }
  return(x)
}

hemi_filter <- function(df){
  if(nrow(df) >=1){
      x <- subset(df, ((is.na(EXAC_AC_Hemi) | EXAC_AC_Hemi <= 5)
		    & (is.na(GNOMADe_AC_male) | GNOMADe_AC_male <= 5)
		    & (is.na(GNOMADg_AC_male) | GNOMADg_AC_male <= 5)
		    & (is.na(GNOMADgV3_AC_male) | GNOMADgV3_AC_male <= 5)
		    & (is.na(TOPMED_Hom) | TOPMED_Hom <= 3)
                    ))
  }else{
      x <- df
  }
  return(x)
}

red_pen_filter <- function(df){
  if(nrow(df) >=1){
      x <- subset(df, ((is.na(GNOMADg_AF) | GNOMADg_AF <= 0.001)
		    & (is.na(GNOMADgV3_AF) | GNOMADgV3_AF <= 0.001)
		    & (is.na(GNOMADe_AF) | GNOMADe_AF <= 0.001)
                    & (is.na(EXAC_AF) | EXAC_AF <= 0.001)
                    & (is.na(TOPMED_AF) | TOPMED_AF <= 0.001)
                    & (is.na(grch37_NGC_AF) | grch37_NGC_AF <= 0.05)
                    & (is.na(NGC_AF) | NGC_AF <= 0.05)
                    ))
  }else{
      x <- df
  }
  return(x)
}



