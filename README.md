# SNV Pipeline

SNV calling pipeline developed explicitly to process individual or trio vcf files obtained from Illumina based pipeline (grch37/grch38).
The pipeline requires user defined datasets & annotation sources, available tools and input set of vcf files. It generates analysis scripts that can be incorporated into high performance cluster (HPC) computing to process the samples. This results in list of filtered variants per family that can be used for interpreation, reporting and further downstream analysis.

# Installation 

	git clone https://github.com/ajaarma/snv.git


#							#
# Required Installation packages			#
#							#

##### Install anaconda v2.0 #####
	Follow this link for installation: https://docs.anaconda.com/anaconda/install/linux/


##### Conda environment commands ##########
	$ conda create --name snv
	$ source activate snv
	$ conda install python=2.7.16
	$ pip install xmltodict
	$ pip install dicttoxml
	
	$ conda install -c bioconda gvcfgenotyper
	$ conda install -c anaconda gawk	
	$ conda install samtools=1.3
	$ conda install vcftools=0.1.14
	$ conda install bcftools=1.9
	$ conda install gcc #(OSX)
	$ conda install gcc_linux-64 #(Linux)
	$ conda install parallel
	$ conda install -c mvdbeek ucsc_tools
	** conda-develop -n <snv-environment-name> <path-to-vep-plugin/vep/Plugins/ **
	  -- Example: $ conda-develop -n snv <path-prefix>/demo/softwares/vep/Plugins/

	$ conda install -c r r-optparse
	$ conda install -c r r-dplyr
	$ conda install -c r r-plyr
	$ conda install -c r r-data.table
	$ conda install -c aakumar r-readbulk
	$ conda install -c bioconda ensembl-vep=100.4
	$ vep_install -a cf -s homo_sapiens -y GRCh37 -c <path-prefix>/demo/softwares/vep/grch37 --CONVERT
	$ vep_install -a cf -s homo_sapiens -y GRCh38 -c <path-prefix>/demo/softwares/vep/grch38 --CONVERT
#						#	
# Data directory and datasets			#
#						#

##### Default datasets provided #####
	1. exac_pli: demo/resources/gnomad/grch37/gnomad.v2.1.1.lof_metrics.by_transcript_forVEP.txt
	2. ensembl: demo/resources/ensembl/grch37/ensBioMart_grch37_v98_ENST_lengths_191208.txt
	3. region-exons: demo/resources/regions/grch37/hg19_refseq_ensembl_exons_50bp_allMT_hgmd_clinvar_20200519.txt
	4. region-pseudo-autosomal: demo/resources/regions/grch37/hg19_non_pseudoautosomal_regions_X.txt
	5. HPO: demo/resources/hpo/phenotype_to_genes.tar.gz
	
##### Other datasets that require no entry to user-configuration file #####
	6. Curated: 
		6.1. Genelist: demo/resources/curated/NGC_genelist_allNamesOnly-20200519.txt
		6.2. Somatic mosaicism genes: demo/resources/curated/haem_somatic_mosaicism_genes_20191015.txt
		6.3. Imprinted gene list: demo/resources/curated/imprinted_genes_20200424.txt
		6.4. Polymorphic gene list: demo/resources/curated/polymorphic_genes_20200509.txt
	7. OMIM: demo/resources/omim/omim_20200421_geneInfoBase.txt

# Download link for following dataset and place them in corresponding directories as shown
   
	1. HPO: Extract HPO phenotypes mapping:
   		$ cd <path-prefix>/demo/resources/hpo/
   		$ tar -zxvf phenotypes_to_genes.tar.gz 

	2. REFERENCE SEQUENCE GENOME (FASTA file alongwith Index)
   		Download link: https://drive.google.com/drive/folders/1Ro3pEYhVdYkMmteSr8YRPFeTvb_K0lVf?usp=sharing
		Download file: Homo_sapiens.GRCh37.74.dna.fasta
			Get the corresponding index and dict files: *.fai and *.dict
		Put this in folder: demo/resources/genomes/grch37/Homo_sapiens.GRCh37.74.dna.fasta

	3. GNOMAD
   		Download link (use wget): 
		Genomes: https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz
		Exomes: https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz
		Put it in this folder: 
			demo/resources/gnomad/grch37/gnomad.genomes.r2.1.1.sites.vcf.bgz
			demo/resources/gnomad/grch37/gnomad.exomes.r2.1.1.sites.vcf.bgz
		Edit User config flat file CONFIG/UserConfig.txt : 
			gnomad_g=gnomad/grch37/gnomad.genomes.r2.1.1.sites.vcf.bgz
			gnomad_e=gnomad/grch37/gnomad.exomes.r2.1.1.sites.vcf.bgz
  
	4. ExAC:
  		Download Link: https://drive.google.com/drive/folders/11Ya8XfAxOYmlKZ9mN8A16IDTLHdHba_0?usp=sharing
		Download file: ExAC.r0.3.1.sites.vep.decompose.norm.prefixed_PASS-only.vcf.gz
			also the index files (*.csi and *.tbi)
		Put it in this folder as: 
			demo/resources/exac/grch37/ExAC.r0.3.1.sites.vep.decompose.norm.prefixed_PASS-only.vcf.gz
		Edit User config flat file CONFIG/UserConfig.txt : 
			exac=exac/grch37/ExAC.r0.3.1.sites.vep.decompose.norm.prefixed_PASS-only.vcf.gz
			exac_t=exac/grch37/ExAC.r0.3.1.sites.vep.decompose.norm.prefixed_PASS-only.vcf.gz
  
	5. CADD:
  		Download link (use wget):
			https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz
			https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/InDels.tsv.gz
			(Also download the corresponding tabix index files as well)
		Put it in this directory: 
			demo/resources/cadd/grch37/whole_genome_SNVs.tsv.gz
			demo/resource/cadd/grch37/InDels.tsv.gz
		Edit the user config flat file CONFIG/UserConfig.txt :
			cadd_snv=cadd/grch37/whole_genome_SNVs.tsv.gz
			cadd_indel=cadd/grch37/InDels.tsv.gz

	6. REVEL:
  		Download link: https://drive.google.com/drive/folders/12Tl1YU5bI-By_VawTPVWHef7AXzn4LuP?usp=sharing
		Download file: new_tabbed_revel.tsv.gz
		         Also the index file: *.tbi
		Put it in this directory: demo/resources/revel/grch37/new_tabbed_revel.tsv.gz
		Edit the user config flat file CONFIG/UserConfig.txt : 
			revel=revel/grch37/new_tabbed_revel.tsv.gz

	7. HGMD:
  		Download link: http://www.hgmd.cf.ac.uk/ac/index.php (Require personal access login)
		Put it in this directory: demo/resources/hgmd/grch37/hgmd_pro_2019.4_hg19.vcf

		Use this command to process HGMD file inside this directory:
			$ cat hgmd_pro_2019.4_hg19.vcf | sed '/##comment=.*\"/a  ##INFO=<ID=ID,Number=1,Type=String,Description=\"HGMD ID\">' | awk -v OFS="\t" '{ if(/^#/){ print }else{ print $1,$2,$3,$4,$5,$6,$7,"ID="$3";"$8 } }' | bgzip -c  > hgmd_pro_2019.4_hg19_wID.vcf.gz
			$ bcftools index -t hgmd_pro_2019.4_hg19_wID.vcf.gz
			$ bcftools index hgmd_pro_2019.4_hg19_wID.vcf.gz	

		Put it in this directory: demo/resources/hgmd/grch37/hgmd_pro_2019.4_hg19_wID.vcf.gz
		Edit the user config flat file CONFIG/UserConfig.txt :
			hgmd=hgmd/grch37/hgmd_pro_2019.4_hg19_wID.vcf.gz

	8. CLINVAR:
  		Download link: 
			https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/weekly/clinvar_20200506.vcf.gz
			https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/weekly/clinvar_20200506.vcf.gz.tbi
		Put it in this directory: demo/resources/clinvar/grch37/clinvar_20200506.vcf.gz
		Edit the user config flat file CONFIG/UserConfig.txt :
			clinvar=clinvar/grch37/clinvar_20200506.vcf.gz

##### Customized Curated Annotation sets	#####

	Default present with this distribution. Can be found in XML file with these tags:
		(1) GeneList: <genelist>
		(2) Somatic mosaicism genes: <haemGenesFile>
		(3) Imprinted genes: <imprintedGenesFile>
		(4) Polymorphic genes: <polymorphicGenesFile>
		(5) HPO terms: <hpo>
		(6) OMIM: <omim>
 
#									#
### Activate the conda environment ###
	$ source activate snv

### Step - 1: ###

	1. Edit CONFIG/UserConfig.txt: 
   		(a) Add the absolute path prefix for the resources directory with tag: resourceDir. 
		    An example can be seen in CONFIG/Example-UserConfig.txt file.
		(b) Manually check if datasets corresponding to other field tags are correctly downloaded and 
		    put in respective folders.
	2. Create user defined XML file from input User Configuration flat file and Base-XML file

##### Command: #####

	$ python createAnalysisXML.py -u <user-config-flat-file> 
			    	      -b <existing-base-xml-file> 
			              -o <out-xml-file-name>

##### Example: #####

	$ python createAnalysisXML.py -u CONFIG/UserConfig.txt 
			       	      -b CONFIG/Analysis_base_grch37.xml 
			              -o CONFIG/Analysis_user_grch37.xml

##### Outputs: ##### 
	CONFIG/Analysis_user_grch37.xml

### Step-2: ###
	1. Put the respective vcf files in the directory. For example: demo/example/vcf/ 
	2. Create manifest file in same format as shown in demo/example/example_manifest.txt
	3. Assign gender to each family members (illumina or sample id). For example: demo/example/example_genders.txt
	4. List of all the family ids that needed to be analyzed. See example sample_pedigree.txt file.
	         For e.g: demo/example/example_family_analysis.txt

### Step -3: ###

   Generate all the shell scripts that can be incorporated into user specific HPC cluster network.

##### Command: #####
	$ python processSNV.py 	-a <USER-XML-FILE>
		    		-p <Project-date>
		      		-m <manifest-file>
		     		-e <analysis-type: gvcfGT>
				-w <work-dir>
		     		-g <gender-file>
		     		-d <samples-fof>
		     		-f <List-of-family-ids-tobe-analyzed>
				-s <List-of-all-family-members-in-manifest>
				-r <family-header-mapping-file (optional)> 
		     
##### Example: #####
	$ python processSNV.py 	-a CONFIG/Analysis_user_grch37.xml \ 
				-p 20210326 \
				-m <path-prefix>/demo/example/example_manifest.txt \
				-e gvcfGT \
				-w <path-prefix>/demo/example/ \
				-g <path-prefix>/demo/example/example_genders.txt \
				-d <path-prefix>/demo/example/exeter_samples_norm.fof \
				-f <path-prefix>/demo/example/manifest/example_family_analysis.txt \
				-s <path-prefix>/demo/example/manifest/example_family.fof \
				-r <path-prefix>/demo/example/manifest/example_family_header.txt (optional)


##### Outputs: #####
	Two scripts in the directory: <path-prefix>/demo/example/20210326/tmp_binaries/
	Launch the scripts in these 2 stages sequentially after each of them gets finished.

	   (1) genotypeAndAnnotate_chr%.sh where %=1..22,X,Y and MT
		scatter the annotation and frequency filtering per chromosome for all families.
	   (2) mergeAndFilter.sh:
		Merge all the chromosome and apply inheritance filtering.


### Step-4 ###

	Final output of list of filtered variant is present in:
	<path-prefix>/demo/example/20210326/fam_filter/<family-id>/<fam-id>.filt_<project-date>.txt


##### For any questions/issues/bugs please mail us at: #####
	Ajay: aak@ebi.ac.uk
	Courtney: cf458@cam.ac.uk
	Alba: as2635@cam.ac.uk
