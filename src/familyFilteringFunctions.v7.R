#!/usr/bin/env Rscript
#
#author: Alba Sanchis Juan (as2635@cam.ac.uk)
#
# last edited: 11 Sep 2019 - CEF
# Last edited: 30 Jul 2020 - AAK
#####################
##Loading libraries##
#####################

library("optparse")
library("plyr")
library("dplyr")
library("data.table")

`%ni%` <- Negate(`%in%`)

#############
##Functions##
#############

getGT <- function(df, samples, ped){
  
    for (sample in ped$ngc_id){
        df[, paste0(sample, "_AD")] <- sapply(df[, sample], function(x)strsplit(x, ":")[[1]][2])
        df[, paste0(sample, "_PL")] <- sapply(df[, sample], function(x)strsplit(x, ":")[[1]][5])
        df[, sample] <- sapply(df[, sample], function(x)strsplit(x, ":")[[1]][1])
    }
    return(df)
}

addgeneList <- function(df, gl){
    df$in_gene_list <- df$ANNO_SYMBOL %in% gl
    return(df)
}
addgeneList2 <- function(df, gl){
    df$in_haem_list <- df$ANNO_SYMBOL %in% gl
    return(df)
}

addPolyCols <- function(df, gl){
  if(gl[1, "GENE"] == "all_genes"){
    df$poly_hom <- TRUE
    df$poly_comphet <- TRUE
    return(df)
  }
  df$poly_hom <- df$ANNO_SYMBOL %in% subset(gl, "hom" %in% strsplit(INHERITANCE, ",")[[1]])$GENE
  df$poly_comphet <- df$ANNO_SYMBOL %in% subset(gl, "comphet" %in% strsplit(INHERITANCE, ",")[[1]])$GENE
  return(df)
}

addOMIM <- function(df, om){
    df$OMIM_disorder <- omim$disorder[match(df$ANNO_SYMBOL, omim$gene_symbol)]
    df$OMIM_gene_moi <- omim$moi[match(df$ANNO_SYMBOL, omim$gene_symbol)]
    return(df)
}

readDepthFilter <- function(df, affectedH) {
    if(nrow(df) >=1){
        first_flag = TRUE
        for (affected1 in affectedH){
          affected_AD <- sapply(affected1, function(x)paste0(x, "_AD"))
          tmp1 <- df[apply(data.frame(df[,affected_AD])
                       , MARGIN = 1
                       , function(x) all(as.numeric(unlist(strsplit(x, "\\,"))[2]) > 11))
                  ,]
          tmp2 <- df[apply(data.frame(df[,affected_AD])
                       , MARGIN = 1
                       , function(x) all(as.numeric(unlist(strsplit(x, "\\,"))[2]) > 2))
                  ,]
          tmp2B <-  subset(tmp2,  FILTER == 'PASS')
          if(nrow(tmp2) >=1){
            tmp2C <- tmp2[apply(data.frame(tmp2[,affected_AD])
                            , MARGIN = 1
                            , function(x) is.character(x) & ratioAl(x) > 0.2 & ratioAl(x) < 0.8  )
                  ,]
          }else{
            tmp2C <- tmp2
          }
          if(first_flag){
            x <- unique(rbind(tmp1, tmp2B, tmp2C))
            first_flag <- FALSE
          } else{
            x <- inner_join(x, unique(rbind(tmp1, tmp2B, tmp2C)))
          } 
        }
    }else{
        x <- df
    }
    return(x)
}

ratioAl <- function(x){
  
    a <- as.numeric(unlist(strsplit(x, "\\,"))[1])
    b <- as.numeric(unlist(strsplit(x, "\\,"))[2])
    r <- b / (a+b)
    return(r)
}

dominant_inheritance_filter <- function(df){
    if(nrow(df) >=1){
        tmp1 <- subset(df, is.na(df$OMIM_gene_moi)
                       | df$OMIM_gene_moi == "NA"
                       | df$OMIM_gene_moi == ""
                       )
        tmp2 <- df[apply(data.frame(df[,"OMIM_gene_moi"])
                       , MARGIN = 1
                       , function(x) FALSE %in% paste(unlist(strsplit(x, " | ", fixed = TRUE)) %in% c("Autosomal recessive", "X-linked recessive")))
                  ,]
        x <- unique(rbind(tmp1, tmp2))
    }else{
        x <- df
    }
    return(x)
}

impact_filter <- function(df){

  x <- subset(df, ANNO_IMPACT %in% c("HIGH","MODERATE") 
              | HGMD_ID != "" 
              | grepl("splice_donor_5th_base_variant",ANNO_SpliceRegion)
              | (grepl("splice_donor_0th_base_variant",ANNO_SpliceRegion) & grepl("splice_region_variant",ANNO_Consequence))
              | (grepl("splice_donor_region_variant",ANNO_SpliceRegion) & grepl("splice_region_variant",ANNO_Consequence))
              | (grepl("splice_polypyrimidine_tract_variant",ANNO_SpliceRegion) & (nchar(REF) != 1 | nchar(ALT) != 1 | (REF %in% c("T","C") & ALT %in% c("A","G")) | (REF %in% c("A","G") & ALT %in% c("T","C"))))
              | grepl("Pathogenic|Likely_pathogenic|Uncertain_significance", CLNSIG)
              | grepl("Pathogenic|Likely_pathogenic", CLNSIGCONF))

  return(x)
}

##MOI functions##

denovo <- function(df, affected, non_affected, mother, father, imprintedgenes, name){

    ##Variants where the non-affected relatives are all ref
    ##Het if no non_affected; het+hom+hemi for others
    ##Reports for autosomal, X-female, and X-male
    
    tmp <- readDepthFilter(df, affected)
    
    if(length(non_affected)==0){ #Singletons
        tmp2 <- tmp[apply(data.frame(tmp[,affected])
                        , MARGIN = 1
                        , function(x) all(x %in% c("0/1", "1/0")))
                   ,]

        tmp3 <- tmp2

    }else{ # Duos/ Trios
        tmp2 <- tmp[apply(data.frame(tmp[,affected])
                        , MARGIN = 1
                        , function(x) all(x %in% c("0/1", "1/0", "1", "1/1")))
                   ,]
        
        tmp3 <- tmp2[apply(data.frame(tmp2[,non_affected])
                       , MARGIN = 1
                       , function(x) all(x %in% c("0/0", "0", "./.", ".")))
                  ,]
    }
    
    tmp4 <- dominant_filter(tmp3)

    tmp5 <- subset(subset(tmp3, in_haem_list), ANNO_IMPACT %in% c("HIGH") 
                   | grepl("DM", HGMD_CLASS)
                   | grepl("Pathogenic|Likely_pathogenic", CLNSIG)
                   | grepl("Pathogenic|Likely_pathogenic", CLNSIGCONF))

    tmp6 <- imprinted(subset(tmp, ANNO_SYMBOL %in% imprintedgenes$GENE), non_affected, affected, mother, father, imprintedgenes)

    tmp7 <- parental_mosaic(subset(tmp,  FILTER == 'PASS'), non_affected, affected, mother, father)
    
    x <- unique(rbind(tmp4, tmp5, tmp6, tmp7))
    if(nrow(x) >=1){x$MOI <- name}
    return(x)

}


AR_hom <- function(df, affected, non_affected, parents, name){

    ##Variant is homozygous in affected
    ##Variant is not homozygous in non-affected
    ##At least one allele is inherited if have both parents
    ##Reports for autosomal, and X-female where inherited from mother
  
    tmp <- readDepthFilter(df, affected)
    
    tmp2 <- tmp[apply(data.frame(tmp[,affected])
                    , MARGIN = 1
                    , function(x) all(x %in% c("1/1")))
               ,]
    
    if(length(non_affected)==0){

        tmp3 <- tmp2

    }else{
        
        tmp2A <- tmp2[apply(data.frame(tmp2[,non_affected])
                         , MARGIN = 1
                         , function(x) all(x %ni% c("1/1", "1")))
                    ,]

        if(length(parents)==2){
            tmp3 <- tmp2A[apply(data.frame(tmp2A[,parents])
                            , MARGIN = 1
                            , function(x) any(x %in% c("0/1", "1/0")))
                       ,]
        }else{
            tmp3 <- tmp2A
        }
    }
    
    tmpA1 <- recessive_filter(tmp3)
    tmpA2 <- subset(subset(tmp3, poly_hom), (grepl("Pathogenic|Likely_pathogenic", CLNSIG) | grepl("Pathogenic|Likely_pathogenic", CLNSIGCONF) | grepl("DM", HGMD_CLASS)) 
                    & !grepl("Benign|Likely_benign", CLNSIG) 
                    & !grepl("Benign|Likely_benign", CLNSIGCONF))
    x <- unique(rbind(tmpA1, tmpA2))
    if(nrow(x) >=1){x$MOI <- name}
    return(x)
  
}


XL <- function(df, affected, non_affected, name){

    ##Hemizygous variants inherited from mother
    ##Reports for X-male, also reports de novo

    tp <- subset(df, CHROM == "X")
    tmp <- readDepthFilter(tp, affected)

    tmp2 <- tmp[apply(data.frame(tmp[,affected])
                  , MARGIN = 1
                  , function(x) all(x %in% c("1")))
            ,]
  
    if(length(non_affected)==0){

        tmp3 <- tmp2

    }else{
        tmp3 <- tmp2[apply(data.frame(tmp2[,non_affected])
                         , MARGIN = 1
                         , function(x) all(x %ni% c("1/1", "1")))
                    ,]
    }

    tmp4 <- subset(tmp, ANNO_SYMBOL %in% c("PCDH19"))
    tmp5 <- tmp4[apply(data.frame(tmp4[,affected])
                  , MARGIN = 1
                  , function(x) all(x %in% c("0/1", "1/0")))
                ,]
    
    if(length(non_affected)==0){

        tmp6 <- tmp5
        
    }else{
        tmp6 <- tmp5[apply(data.frame(tmp5[,non_affected])
                         , MARGIN = 1
                         , function(x) all(x %ni% c("0/1", "1/0")))
                    ,]
    }
                  

    x <- unique(rbind(hemi_filter(tmp3), red_pen_filter(tmp6)))
        
    if(nrow(x) >=1){x$MOI <- name}
    return(x)
}


AR_comphet <- function(df, affected, non_affected, parents, mother, father, name){

    ##AR_comphet for single/duo/trio, any number affected:
    ##Genes with 2 het variants: at least one denovo (if duo/trio) and one inherited (if trio) OR
    ##at least one from each parent (if trio)   
    ##Reports for autosomal and X-female where inherited is from mother (denovo_inh only)
    ##Doesnt check if CH in unaffected 
    ##Note: dont know if they are in trans unless trio and one from each parent
    ##Note: doesnt work when affected have different parents (NGC00117) - overcalls


    tmp <- readDepthFilter(df, affected)
    tmpA1 <- recessive_filter(tmp)
    tmpA2 <- subset(subset(tmp, poly_hom | poly_comphet), (grepl("Pathogenic|Likely_pathogenic", CLNSIG) | grepl("Pathogenic|Likely_pathogenic", CLNSIGCONF) | grepl("DM", HGMD_CLASS)) 
                    & !grepl("Benign|Likely_benign", CLNSIG) 
                    & !grepl("Benign|Likely_benign", CLNSIGCONF))
    tmpA <- unique(rbind(tmpA1, tmpA2))
    
    affected_hets <- tmpA[apply(data.frame(tmpA[,affected])
                            , MARGIN = 1
                            , function(x) all(x %in% c("0/1", "1/0")))
                       ,]
    
    if(length(parents) == 0){ #singleton
        x1 <- affected_hets[apply(data.frame(affected_hets[,"ANNO_SYMBOL"])
                  , MARGIN = 1
                  , function(g) length(subset(affected_hets, ANNO_SYMBOL == g)$ANNO_SYMBOL) >1)
                   ,]
        x1a <- subset(x1, !poly_comphet | poly_hom)
        x1b <- subset(subset(x1, poly_comphet & !poly_hom), ANNO_SYMBOL %in% subset(x1, POS %in% tmpA1$POS)$ANNO_SYMBOL)
        x <- unique(rbind(x1a, x1b))
    }else{
        nonaffected_nonhom <- affected_hets[apply(data.frame(affected_hets[,non_affected])
                                                , MARGIN = 1
                                                , function(x) all(x %ni% c("1/1", "1")))
                                           ,]
        denovo <- nonaffected_nonhom[apply(data.frame(nonaffected_nonhom[,non_affected])
                                         , MARGIN = 1
                                         , function(x) all(x %in% c("0/0", "0", "./.", ".")))
                                    ,]
        denovo_genes <- unique(denovo$ANNO_SYMBOL)

        if(length(parents) == 1){ #duo
            tmp2 <- nonaffected_nonhom[apply(data.frame(nonaffected_nonhom[,"ANNO_SYMBOL"])
                                           , MARGIN = 1
                                           , function(g) length(subset(nonaffected_nonhom, ANNO_SYMBOL == g)$ANNO_SYMBOL) >1)
                                      ,]
            x1 <- subset(tmp2, ANNO_SYMBOL %in% denovo_genes)
            x1a <- subset(x1, !poly_comphet | poly_hom)
            x1b <- subset(subset(x1, poly_comphet & !poly_hom), ANNO_SYMBOL %in% subset(x1, POS %in% tmpA1$POS)$ANNO_SYMBOL)
            x <- unique(rbind(x1a, x1b))
            
        }else{ #trio
            genes_denovo_comphet <- sapply(denovo_genes,
                                    function(g) any(apply(data.frame(nonaffected_nonhom[nonaffected_nonhom$ANNO_SYMBOL == g, parents]),
                                                          MARGIN = 1, function(x) any(x %in% c("0/1", "1/0")) & any(x %in% c("0/0", "0", "./.", ".")))))
  
            tmp3 <- subset(nonaffected_nonhom, ANNO_SYMBOL %in% attributes(genes_denovo_comphet[genes_denovo_comphet == TRUE])[[1]])

            tmp4  <- AR_comphet_true(nonaffected_nonhom, affected, non_affected, parents, mother, father)
            
            x1 <- unique(rbind(tmp3, tmp4))
            x1a <- subset(x1, !poly_comphet | poly_hom)
            x1b <- subset(subset(x1, poly_comphet & !poly_hom), ANNO_SYMBOL %in% subset(x1, POS %in% tmpA1$POS)$ANNO_SYMBOL)
            x <- unique(rbind(x1a, x1b))
            
        }
    }
  
  if(nrow(x) >=1){x$MOI <- name}
  return(x)
  
}


AR_comphet_true <- function(nonaffected_nonhom, affected, non_affected, parents, mother, father){

    ##Genes with 2 variants, one from each parent
    ##Reports for autosomal only
    ##Needs both parents
    ##Doesnt check if CH in unaffected 
      
  nonaffected_nonhom$motherMask <- nonaffected_nonhom[,mother] %in% c("0/1", "1/0")
  nonaffected_nonhom$fatherMask <- nonaffected_nonhom[,father] %in% c("0/1", "1/0")
  
  nonaffected_nonhom$compHetFatherMask <- !nonaffected_nonhom$motherMask & nonaffected_nonhom$fatherMask
  nonaffected_nonhom$compHetMotherMask <- nonaffected_nonhom$motherMask & !nonaffected_nonhom$fatherMask
  
  
  ##?tapply
  symbolList <- sapply(unique(nonaffected_nonhom$ANNO_SYMBOL), function(ANNO_SYMBOL){
    any(nonaffected_nonhom[nonaffected_nonhom$ANNO_SYMBOL == ANNO_SYMBOL,]$compHetMotherMask) & 
      any(nonaffected_nonhom[nonaffected_nonhom$ANNO_SYMBOL == ANNO_SYMBOL,]$compHetFatherMask)
  })
  
  if(length(symbolList[symbolList == TRUE]) > 0){
    y <- nonaffected_nonhom[nonaffected_nonhom$ANNO_SYMBOL %in% attributes(which(symbolList) == TRUE)$names, 
            !names(nonaffected_nonhom) %in% c("fatherMask", "motherMask", "compHetFatherMask", "compHetMotherMask")]
  }else{
    y <- NULL
  }
  
  return(y)
  
}


MT <- function(df, affected){
  x <- subset(df, CHROM == "nothing")
  if(nrow(x) >=1){x$MOI <- "MT"}
  return(x)
}


inherited <- function(df, affected, parents, name){

    ##Variants where the affected relatives are all heterozygous
    ##Reports for autosomal, X-female (and incidentally for male XLR)
    ##Works for single, duo, trio, quad (both affected have to have variant)
    
  tmp2A <- df[apply(data.frame(df[,affected])
                  , MARGIN = 1
                  , function(x) all(x %in% c("0/1", "1/0", "1")))
            ,]

  if(length(parents) == 2){ #trio
      tmp2 <- tmp2A[apply(data.frame(tmp2A[,parents])
                       , MARGIN = 1
                       , function(x) any(x %in% c("0/1", "1/0", "1")))
                  ,]
  }else{ #single/duo
      tmp2 <- tmp2A
  }
  
  tmp3 <- red_pen_filter(tmp2)
  tmp4 <- readDepthFilter(tmp3, affected)
  x <- dominant_inheritance_filter(tmp4)
  
  if(nrow(x) >=1){x$MOI <- name}
  return(x)
}


imprinted <- function(df, non_affected, affected, mother, father, imprintedgenes){

    tmp <- df[apply(data.frame(df[,affected])
                  , MARGIN = 1
                  , function(x) all(x %in% c("0/1", "1/0", "1/1")))
            ,]
    tmp2 <- tmp[apply(data.frame(tmp[,non_affected])
                    , MARGIN = 1
                    , function(x) all(x %ni% c("1/1")))
              ,]
    
    if(nrow(tmp2) ==0){return(tmp2)}
    if(mother == "None" & father == "None"){
          tmp3 <- tmp2
      }else if(mother == "None" | father == "None"){
        if(mother != "None"){
            tmp3 <- tmp2[apply(data.frame(tmp2[,c("ANNO_SYMBOL", mother)])
                          , MARGIN = 1
                          , function(x) checkImprint(x, "mother", imprintedgenes))
                    ,]
        }else if(father != "None"){
            tmp3 <- tmp2[apply(data.frame(tmp2[,c("ANNO_SYMBOL", father)])
                          , MARGIN = 1
                          , function(x) checkImprint(x, "father", imprintedgenes))
                    ,]
        }
      }else{
          tmp3 <- tmp2[apply(data.frame(tmp2[,c("ANNO_SYMBOL", mother, father)])
                        , MARGIN = 1
                        , function(x) checkImprint(x, "trio", imprintedgenes))
                  ,]
      }

    y <- red_pen_filter(tmp3)
    
    return(y)

}

checkImprint <- function(t, tp, imprintedgenes){
  g <- imprintedgenes[imprintedgenes$GENE == t[1],]$EXPRESSION
  if(tp == "trio"){
    if(g == "Both"){
      return(any(t[c(2,3)] %in% c("0/1", "1/0")))
    }else if(g == "Maternal"){
      return(t[2] %in% c("0/1", "1/0"))
    }else if(g == "Paternal"){
      return(t[3] %in% c("0/1", "1/0"))
    }
  }else if(tp == "mother"){
    if(g == "Maternal"){
      return(t[2] %in% c("0/1", "1/0"))
    }else{
      return(TRUE)
    }
  }else if(tp == "father"){
    if(g == "Paternal"){
      return(t[2] %in% c("0/1", "1/0"))
    }else{
      return(TRUE)
    }
  }
  return(FALSE)
}


parental_mosaic <- function(df, non_affected, affected, mother, father){
  
  tmpA <- df[apply(data.frame(df[,affected])
                  , MARGIN = 1
                  , function(x) all(x %in% c("0/1", "1/0")))
            ,]
  tmp <- tmpA[apply(data.frame(tmpA[,sapply(affected, function(x)paste0(x, "_AD"))])
                   , MARGIN = 1
                   , function(x) all(sapply(x, function(y) is.character(y) & ratioAl(y) > 0.4)))
             ,]

  tmp2A <- tmp[apply(data.frame(tmp[,non_affected])
                    , MARGIN = 1
                    , function(x) all(x %ni% c("1/1")))
              ,]
  tmp2 <- tmp2A[apply(data.frame(tmp2A[,non_affected[non_affected %ni% c(mother, father)]])
                     , MARGIN = 1
                     , function(x) all(x %in% c("0/0", "0", "./.", ".")))
               ,]
  tmp3m <- NULL
  tmp3f <- NULL
  if(mother != "None" | father != "None"){
    if(mother != "None"){
      tmp3A <- tmp2[apply(data.frame(tmp2[,mother])
                         , MARGIN = 1
                         , function(x) x %in% c("0/1", "1/0"))
                   ,]
      tmp3B <- tmp3A[apply(data.frame(tmp3A[,paste0(mother, "_AD")])
                           , MARGIN = 1
                           , function(x) is.character(x) & ratioAl(x) < 0.4)
                     ,]
      if(nrow(tmp3B) >=1){
      
        tmp3C <- tmp3B[apply(data.frame(tmp3B[,paste0(mother, "_PL")])
                           , MARGIN = 1
                           , function(x) is.character(x) & as.numeric(unlist(strsplit(x, "\\,"))[1]) < 20)
                     ,]
      }else{
        tmp3C <- tmp3B
      }
      if(father != "None"){
        tmp3m <- tmp3C[apply(data.frame(tmp3C[,father])
                           , MARGIN = 1
                           , function(x) x %in% c("0/0", "0", "./.", "."))
                     ,]
      }else{
        tmp3m <- tmp3C
      }
    }
    if(father != "None"){
      tmp3A <- tmp2[apply(data.frame(tmp2[,father])
                          , MARGIN = 1
                          , function(x) x %in% c("0/1", "1/0"))
                    ,]
      tmp3B <- tmp3A[apply(data.frame(tmp3A[,paste0(father, "_AD")])
                           , MARGIN = 1
                           , function(x) is.character(x) & ratioAl(x) < 0.4)
                     ,]
      if(nrow(tmp3B) >=1){
        tmp3C <- tmp3B[apply(data.frame(tmp3B[,paste0(father, "_PL")])
                           , MARGIN = 1
                           , function(x) is.character(x) & as.numeric(unlist(strsplit(x, "\\,"))[1]) < 20)
                     ,]
      }else{
        tmp3C <- tmp3B
      }
      if(mother != "None"){
        tmp3f <- tmp3C[apply(data.frame(tmp3C[,mother])
                            , MARGIN = 1
                            , function(x) x %in% c("0/0", "0", "./.", "."))
                      ,]
      } else{
        tmp3f <- tmp3C
      }
    }
  } else{
    return(NULL)
  }

  y <- dominant_filter(unique(rbind(tmp3m, tmp3f)))
  return(y)
  
}

