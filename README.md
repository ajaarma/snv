# SNV Pipeline

# Installation 

git clone https://github.com/ajaarma/snv.git


#							#
# Required Installation packages			#
#							#

#Install anaconda v2.0
Follow this link for installation: https://docs.anaconda.com/anaconda/install/linux/


##### Conda environment commands ##########
	$ conda create --name snv
	$ source activate snv
		$ cond env list
	$ conda install python=2.7.16
		$ pip install xmltodict
		$ pip install dicttoxml
	
	$ conda install samtools=1.3
	$ conda install vcftools=0.1.14
	$ conda install bcftools=1.9
	$ conda install gcc #(OSX)
		$ conda install gcc_linux-64 #(Linux)
	$ conda install parallel
	$ conda install -c mvdbeek ucsc_tools
	** conda-develop -n <snv-environment-name> <path-to-vep-plugin/vep/Plugins/ **
	  -- Example: $ conda-develop -n snv <path-prefix>/XTR/softwares/vep/Plugins/

	$ conda install -c r r-optparse
	$ conda install -c r r-dplyr
	$ conda install -c r r-plyr
	$ conda install -c r r-data.table
	$ conda install -c aakumar r-readbulk
	$ conda install -c bioconda ensembl-vep=100.4
		$ vep_install -a cf -s homo_sapiens -y GRCh37 -c <path-prefix>/XTR/softwares/vep/grch37 --CONVERT

#						#	
# Data directory and datasets			#
#						#

##### Default datasets provided #####
	1. exac_pli: XTR/resources/gnomad/grch37/gnomad.v2.1.1.lof_metrics.by_transcript_forVEP.txt
	2. ensembl: XTR/resources/ensembl/grch37/ensBioMart_grch37_v98_ENST_lengths_191208.txt
	3. region-exons: XTR/resources/regions/grch37/hg19_refseq_ensembl_exons_50bp_allMT_hgmd_clinvar_20200519.bed
	4. regio-pseudo-autosomal-region: XTR/resources/regions/grch37/hg19_non_pseudoautosomal_regions_X.txt
	5. HPO: XTR/resources/hpo/phenotype_to_genes.tar.gz
	
##### Other datasets that require no entry to user-configuration file #####
	6. Curated: 
		6.1. Genelist: XTR/resources/curated/NGC_genelist_allNamesOnly-20200519.txt
		6.2. Somatic mosaicism genes: XTR/resources/curated/haem_somatic_mosaicism_genes_20191015.txt
		6.3. Imprinted gene list: XTR/resources/curated/imprinted_genes_20200424.txt
		6.4. Polymorphic gene list: XTR/resources/curated/polymorphic_genes_20200509.txt
	7. OMIM: XTR/resources/omim/omim_20200421_geneInfoBase.txt

# Download link for following dataset and place them in corresponding directories as shown
   
   1. HPO: Extract HPO phenotypes mapping:
   	$ cd path-prefix/XTR/resources/hpo/
   	$ tar -zxvf phenotypes_to_genes.tar.gz 

   2. REFERENCE SEQUENCE GENOME (FASTA file alongwith Index)
   	Download link: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/
	Put this in folder: XTR/resources/genomes/grch37/Homo_sapiens.GRCh37.74.dna.fasta
		Also get the corresponding index file for the above '.fasta' file

   3. GNOMAD
   	Download link: https://gnomad.broadinstitute.org/downloads
	Put it in this folder: XTR/resources/gnomad/grch37/gnomad.genomes.r2.1.1.sites.vcf.bgz
	Edit User config flat file CONFIG/UserConfig.txt : gnomad_g=gnomad/grch37/gnomad.genomes.r2.1.1.sites.vcf.bgz
  
   4. ExAC:
  	Download Link: https://gnomad.broadinstitute.org/downloads (ExAC tab)
	Put it in this folder: XTR/resources/exac/grch37/ExAC.r0.3.1.sites.vep.decompose.norm.prefixed_PASS-only.vcf.gz
	Edit User config flat file CONFIG/UserConfig.txt : 
		exac=exac/grch37/ExAC.r0.3.1.sites.vep.decompose.norm.prefixed_PASS-only.vcf.gz
		exac_t=exac/grch37/ExAC.r0.3.1.sites.vep.decompose.norm.prefixed_PASS-only.vcf.gz
  
   5. CADD:
  	Download link: https://cadd.gs.washington.edu/download
	Put it in this directory: XTR/resources/cadd/grch37/whole_genome_SNVs.tsv.gz
				  XTR/resource/cadd/grch37/InDels.tsv.gz
	Edit the user config flat file CONFIG/UserConfig.txt :
				cadd_snv=cadd/grch37/whole_genome_SNVs.tsv.gz
				cadd_indel=cadd/grch37/InDels.tsv.gz

   6. REVEL:
  	Download link: https://sites.google.com/site/revelgenomics/downloads
	Put it in this directory: XTR/resources/revel/grch37/new_tabbed_revel.tsv.gz
	Edit the user config flat file CONFIG/UserConfig.txt : 
				revel=revel/grch37/new_tabbed_revel.tsv.gz

   7. HGMD:
  	Download link: http://www.hgmd.cf.ac.uk/ac/index.php (Require personal access login)
	Put it in this directory: XTR/resources/hgmd/grch37/hgmd_pro_2019.4_hg19_wID.vcf.gz
	Edit the user config flat file CONFIG/UserConfig.txt :
				hgmd=hgmd/grch37/hgmd_pro_2019.4_hg19_wID.vcf.gz

   8. CLINVAR:
  	Download link: https://www.ncbi.nlm.nih.gov/variation/docs/ClinVar_vcf_files/
	Put it in this directory: XTR/resources/clinvar/grch37/clinvar_20200506.vcf.gz
	Edit the user config flat file CONFIG/UserConfig.txt :
				clinvar=clinvar/grch37/clinvar_20200506.vcf.gz

 # Customized Curated Annotation sets	  #

 Deafult present with this distribution. Can be found in XML file with these tags:
    (1) GeneList: <genelist>
    (2) Somatic mosaicism genes: <haemGenesFile>
    (3) Imprinted genes: <imprintedGenesFile>
    (4) Polymorphic genes: <polymorphicGenesFile>
    (5) HPO terms: <hpo>
    (6) OMIM: <omim>
 
#									#
# Example run for 3 families						#
	(1) TwEx_EX1910699-TwEx_EX1910700-TwEx_EX1910701 
	(2) WE0282-WE0283-WE0323			 
	(3) WE0361-WE0362-WE0363			 

### Activate the conda environment ###
	$ source activate snv

### Step - 1: ###

	1. Create user defined XML file from input User Configuration flat file
      and Base-XML file
	2.  UserConfig.txt: 
   		(a) Add the absolute path prefix for the resources directory with tag: resourceDir. 
		    An example can be seen in CONFIG/Example-UserConfig.txt file.
		(b) Manually check if datasets corresponding to other field tags are correctly downloaded and 
		    put in respective folders.

##### Command: #####

	$ python createAnalysisXML.py -u <user-config-flat-file> 
			    	      -b <existing-base-xml-file> 
			              -o <out-xml-file-name>

##### Example: #####

	$ python createAnalysisXML.py -u CONFIG/UserConfig.txt 
			       	      -b CONFIG/Analysis_base_grch37.xml 
			              -o CONFIG/Analysis_user_grch37.xml

Outputs: => CONFIG/Analysis_user_grch37.xml

### Step-2: ###
   1. Put the respective vcf files in the directory: XTR/example/vcf/ 
   2. Create manifest file from input directory or list of manifest file
   3. Assign gender to each family members
   4. Normalize original vcf files
   5. List of all the family ids that needed to be analyzed. See example sample_pedigree.txt file.

##### Command: #####
	$ python processXTR.py -a <user-config-analysis-xml-file>
		     	-i <input-family-trio-mapping> 
		     	-d <path-to-trio-vcf-files>
		     	-o <path-to-output-directory>
		     	-e <analysis-type: norm>

##### Example: #####
	$ python processXTR.py -a CONFIG/Analysis_user_grch37.xml \ 
			-i <path-prefix>/XTR/example/sample_pedigree.txt \
			-d <path-prefix>/XTR/example/vcf \
			-o <path-prefix>/XTR/example/out \
			-e norm

#### List of all the family ids that needed to be analyzed #####
	$ cut -f 1 <path-prefix>/XTR/example/out/manifest.txt | grep -v '^family' > <path-prefix>/XTR/example/out/extrFamAll.txt
 

### Step -3: ###

   Generate all the shell scripts that can be incorporated into any HPC cluster network.

##### Command: #####
	$ python processSNV.py -a <USER-XML-FILE>
		    	-p <Project-date>
		      	-m <manifest-file>
		     	-e <analysis-type: exter_merge>
			-w <work-dir>
		     	-g <gender-file>
		     	-d <samples-fof>
		     	-f <List-of-family-ids>
		     
##### Example: #####
	$ python processSNV.py -a CONFIG/Analysis_user_grch37.xml \ 
			-p 20200726 \
			-m <path-prefix>/XTR/example/out/manifest/manifest.txt \
			-e exeter \
			-w <path-prefix>/XTR/example/ \
			-g <path-prefix>/XTR/example/out/manifest/gender_genomic.txt \
			-d <path-prefix>/XTR/example/exeter_samples_norm.fof \
			-f <path-prefix>/XTR/example/out/manifest/extrFamALL.txt 

##### Outputs: #####
	Three scripts in the directory: <path-prefix>/XTR/example/20200726/tmp_binaries/
	Launch the scripts in these 3 stages sequentially after each of them gets finished.

	   (1) mergeVCF.sh : For merging all normalized VCF files. 
	   (2) splitAndAnnotate_chr*.sh where *=chr1..chr22,chrX and chrY
		scatter the annotation and frequency filtering per chromosome for all families.
	   (3) mergeAndFilter.sh:
		Merge all the chromosome and apply inheritance filtering.


### Step-4 ###

	Final output of list of filtered variant is present in:
	<path-prefix>/XTR/example/20200726/fam_filter/<family-id>/<fam-id>.filt_<project-date>.txt


#									  #
# Replicating for full cohort of 514 families (n=1543) 			  #
#									  #

    (1) Replace the example 3 families file: <path-prefix>/XTR/example/sample_pedigree.txt 
	with 514 families file : <path-prefix>/XTR/example/exeter_pedigree.txt  
	in Step-2.

    (2) Re-run the analysis Step3-4

    N.B: Average run time for full cohort of 514 families in our cluster: 5-6 hours. The split out is:
	(1) MergeVCF step: approx 40 mins
	(2) splitAndAnnotate step: 7-8 hours in HPC cluster
	(3) mergeAndFilter step: approx 4 hours.


************************* The End ***************************
