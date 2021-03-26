#!/usr/bin/env Rscript
#
#author: Alba Sanchis Juan (as2635@cam.ac.uk)
#
# last edited: 19 May 2020 - CEF

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

source(file.path(script.basename,"/familyFilteringFunctions.v7.R"))

option_list = list(
  make_option(c("-p", "--ped"), type="character", default=NULL,
              help="pedigree file", metavar="character"),
  make_option(c("-f", "--family"), type="character", default=NULL, 
              help="family id", metavar="character"),
  make_option(c("-v", "--vars"), type="character", default=NULL, 
              help="variants file name", metavar="character"),
  make_option(c("-g", "--genelist"), type="character", default=NULL, 
              help="variants file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.filt.txt", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-c", "--columns"), type="character", default=NULL, 
              help="columns.txt", metavar="character"),
  make_option(c("-n", "--pheno"), type="character", default=NULL, 
              help="sample_phenotypes.txt", metavar="character"),
  make_option(c("-d", "--hpodb"), type="character", default=NULL, 
              help="phenotype_to_genes.txt", metavar="character"),
  make_option(c("-m", "--omim"), type="character", default=NULL, 
              help="_geneInfoBase.txt", metavar="character"),
  make_option(c("-s", "--haem_genes"), type="character", default=NULL, 
              help="haem_somatic_mosaicism_genes.txt", metavar="character"),
  make_option(c("-r", "--imprinted_genes"), type="character", default=NULL, 
              help="imprinted_genes.txt", metavar="character"),
  make_option(c("-l", "--poly_genes"), type="character", default=NULL, 
              help="polymorphic_genes.txt", metavar="character"),
  make_option(c("-i", "--inherited_gene_list"), type="character", default=NULL, 
              help="gene list for inherited var filtering [default: use HPO terms only]", metavar="character"),
  make_option(c("-b", "--freq_script"), type="character", default="", 
              help="freq functions script", metavar="character")

); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$vars)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (variants file).n", call.=FALSE)
}


#####################
##Parsing arguments##
#####################

fam <- opt$family
variants_path <- opt$vars
out_path <- opt$out
genelist_path <- opt$genelist
pheno_path <- opt$pheno
hpodb_path <- opt$hpodb
omim_path <- opt$omim
ped_path <- opt$ped
cols_path <- opt$columns
haemgenes_path <- opt$haem_genes
imprintedgenes_path <- opt$imprinted_genes
polygenes_path <- opt$poly_genes
inhgenes_path <- opt$inherited_gene_list
freq_script <- opt$freq_script

source(paste0(script.basename,"/familyFilteringFunctions_freqs",freq_script,".v7.R"))

#################
##Reading files##
#################

cols <- fread(cols_path, header = F, stringsAsFactors = F)
cols <- data.frame(cols)

raw_variants <- fread(variants_path, header = T, stringsAsFactors = F)
raw_variants <- data.frame(raw_variants)

if (pheno_path != '.') {
    hpo <- fread(pheno_path, stringsAsFactors = F, header = T)
    hpo <- data.frame(hpo)

}else {
    hpo <- data.frame()
}


hpodb <- fread(hpodb_path, skip = 1, header = F, stringsAsFactors = F, col.names = c("HPO_ID", "HPO_DESCRIPTION", "ID", "GENE_SYMBOL", "INFO", "SOURCE", "DISEASE"))
hpodb <- data.frame(hpodb)

genelist <- fread(genelist_path, header = F, stringsAsFactors = F)
genelist <- data.frame(genelist)

haemgenes <- fread(haemgenes_path, header = F, stringsAsFactors = F)
haemgenes <- data.frame(haemgenes)

imprintedgenes <- fread(imprintedgenes_path, header = F, stringsAsFactors = F, col.names = c("GENE", "EXPRESSION"))
imprintedgenes <- data.frame(imprintedgenes)

polygenes <- fread(polygenes_path, header = F, stringsAsFactors = F, col.names = c("GENE", "INHERITANCE"))
polygenes <- data.frame(polygenes)

if(!is.null(inhgenes_path)){
  inhgenes <- fread(inhgenes_path, header = F, stringsAsFactors = F)
  inhgenes <- inhgenes$V1
}else{
  inhgenes <- c()
}

omim <- fread(omim_path, header = T, stringsAsFactors = F)
omim <- data.frame(omim)

ped <- fread(ped_path, header = T, stringsAsFactors = F)
ped <- data.frame(ped)

#####################
##Variant filtering##
#####################

ped_fam <- subset(ped, family_id == fam)
samples <- ped_fam$ilmn_id
# family_structure <- length(samples)

affected <- subset(ped_fam, affected == 2)$ngc_id
non_affected <- subset(ped_fam, affected == 1)$ngc_id
#unknown_pheno <- subset(ped_fam,family_id == fam & affected == 0)$ngc_id

variants <- getGT(raw_variants, samples, ped_fam) # 25944 * 1750; 2nd Exeter run:748359*1759 


##Getting only variants where affected is/are not REF
variants_case <- variants[apply(data.frame(variants[,affected])
                                , MARGIN = 1
                                , function(x) any(x %in% c("0/1", "1/1", "1/0", "1")))
                          ,] # 5 * 1750

#Adding column for MOI
variants_case$MOI <- NA

#Getting hpo genes for specific phenotype
if (pheno_path !='.'){
    sample_hpo <- subset(hpo, ngc_id %in% subset(ped_fam, affected == 2)$ngc_id)

    hpo_genes <- unique(subset(
                         hpodb, HPO_ID %in% unique(
                                                    unlist(
                                                           sapply(
                                                                  sample_hpo$observed_HPO_terms, 
                                                                  function(x) strsplit(x, "::")[[1]]
                                                                  )
                                                           )
                                                    )
                            )$GENE
                       )
}else{
    hpo_genes = as.character()

}

genelist_complete <- unique(c(genelist$V1, hpo_genes))

#Adding gene list info to table
variants_case_gene <- addgeneList(variants_case, genelist_complete)
variants_case_gene <- addgeneList2(variants_case_gene, haemgenes$V1)
variants_case_gene <- addPolyCols(variants_case_gene, polygenes)

#Adding omim info to table
variants_case_gene_omim <- addOMIM(variants_case_gene, omim)


#Subset for only variant with high/moderate impact OR in hgmd OR in clinvar
variants_case_gene_imp <- impact_filter(variants_case_gene_omim)
#Filtering by MOI
final <- c()

first_flag = TRUE
for (affected1 in affected){
    
      variants_case_gene_imp1 <- variants_case_gene_imp[apply(data.frame(variants_case_gene_imp[,affected1])
                                                                      , MARGIN = 1
                                                                      , function(x) any(x %in% c("0/1", "1/1", "1/0", "1"))
                                                              ),]
      
      denovo1 <- NULL
      AR_hom1 <- NULL
      AR_comphet1 <- NULL
      MT1 <- NULL
      XL1 <- NULL
      IH1 <- NULL
      
      mother <- subset(ped_fam, ngc_id == affected1)$mother_id
      father <- subset(ped_fam, ngc_id == affected1)$father_id
   
      if(mother == "None" & father == "None"){ #Singeleton
          parents <- c()
          variants_case_gene_imp_filt <- subset(variants_case_gene_imp1, in_gene_list)
      }else if(mother == "None" | father == "None"){
        if(mother != "None"){
          parents <- c(mother)
        }else if(father != "None"){
          parents <- c(father)
        }
        variants_case_gene_imp_filt <- subset(variants_case_gene_imp1, in_gene_list)
      }else{
        parents <- c(mother, father)
        variants_case_gene_imp_filt <- variants_case_gene_imp1
      }
      
      denovo1 <- denovo(variants_case_gene_imp_filt, affected1, non_affected, mother, father, imprintedgenes, "denovo")
      AR_hom1 <- AR_hom(variants_case_gene_imp_filt, affected1, non_affected, parents, "AR_hom")
      XL1 <- XL(variants_case_gene_imp_filt, affected1, non_affected, "XLR")
      AR_comphet1 <- AR_comphet(variants_case_gene_imp_filt, affected1, non_affected, parents, mother, father, "AR_comphet")
      IH1 <- inherited(subset(variants_case_gene_imp_filt, ANNO_SYMBOL %in% hpo_genes | ANNO_SYMBOL %in% inhgenes), affected1, parents, "inherited")
      MT1 <- MT(variants_case_gene_imp1, affected1) # currently doesn't return anything
      
      if(first_flag){
        final <- rbind(denovo1, AR_hom1, AR_comphet1, XL1, IH1, MT1)
        first_flag <- FALSE
      } else{
        final <- rbind(final, denovo1, AR_hom1, AR_comphet1, XL1, IH1, MT1)
      } 
}

if(length(affected)>1){
      variants_case_gene_imp1 <- variants_case_gene_imp[apply(data.frame(variants_case_gene_imp[,affected])
                                                                      , MARGIN = 1
                                                                      , function(x) all(x %in% c("0/1", "1/1", "1/0", "1"))
                                                              ),]    
      denovo1 <- NULL
      AR_hom1 <- NULL
      AR_comphet1 <- NULL
      XL1 <- NULL
      IH1 <- NULL
      
      mother <- subset(ped_fam, ngc_id == affected1)$mother_id
      father <- subset(ped_fam, ngc_id == affected1)$father_id
    
      if(mother != "None" & father != "None"){
          parents <- c(mother, father)
      }else if(mother != "None"){
          parents <- c(mother)
      }else if(father != "None"){
          parents <- c(father)
      }else{
          parents <- c()
      }

      denovo1 <- denovo(variants_case_gene_imp1, affected, non_affected, mother, father, imprintedgenes, "denovo_shared")
      AR_hom1 <- AR_hom(variants_case_gene_imp1, affected, non_affected, parents, "AR_hom_shared")
      XL1 <- XL(variants_case_gene_imp1, affected, non_affected, "XLR_shared")
      AR_comphet1 <- AR_comphet(variants_case_gene_imp1, affected, non_affected, parents, mother, father, "AR_comphet_shared")
      IH1 <- inherited(subset(variants_case_gene_imp1, ANNO_SYMBOL %in% hpo_genes | ANNO_SYMBOL %in% inhgenes), affected, parents, "inherited_shared")
         
      final <- rbind(final, denovo1, AR_hom1, AR_comphet1, XL1, IH1)
} 

####################
##Defining columns##
####################

pre_cols <- cols$V1[1:7]
ngc_ids <- unique(ped_fam$ngc_id)
samples_cols <- c(ngc_ids, paste0(ngc_ids, "_AD"))
post_cols <- cols$V1[8:length(cols$V1)]
columns <- c(pre_cols, samples_cols, post_cols)

##Writing output
columns_rev = columns[columns %in% colnames(final)]
#out_r_data = save.image(file='/home/aak64/rds/rds-flr24-wgs10k-ngc/WGS10K/data/NGC/us/misc/TEST_SNV/XTR/example/20200802/tmp_binaries/tmp.RData')
message('-- Writing the filtered list of variants')
#write.table(subset(final, select = columns), out_path, sep = "\t", quote = F, row.names = F)
write.table(subset(final, select = columns_rev), out_path, sep = "\t", quote = F, row.names = F)
