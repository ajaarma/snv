#!/usr/bin/python

#########################################################################################
# Program to launch scripts for running WGS-SNV pipeline on Exeter cohort (n=514)       #
#                                                                                       #
# =>Input: User designed Configuration-XML file and other parameters (see documentation)#
# => Output: Shell scripts for                                                          #
#        (1) Merging up of normlaized VCF(.gz) files.                                   #
#        (2) Annotation of sample specific vcf files split by chromosome                #
#        (3) Merge by chomosome and family inheritance filtering                        #
#                                                                                       #
#########################################################################################

vers="1.0"

import re,sys,getopt,os,datetime,subprocess
from collections import OrderedDict

script_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(script_path+'/config/')
sys.path.append(script_path+'/cluster/')

from config import *
from cluster import *

##################
####   MAIN   ####
##################


if __name__=="__main__":

    chr_list = [str(x) for x in range(1,23)]+['X','Y','MT'] 
    
    #################################################################
    #  Get command line arguments                                   #
    #################################################################

    objC = CONFIG()
    cmd_dict = objC.processSNVArguments()
    
    #sys.exit()
    xml_file = cmd_dict['xml']; proj_date = cmd_dict['proj'];
    exp_type = cmd_dict['analType']; data_file = cmd_dict['allSample']
    work_dir = cmd_dict['workDir'];manifest = cmd_dict['manifest']
    gender_file = cmd_dict['gender'];family_list = cmd_dict['fam']
    email_id = cmd_dict['email'];flag_samples = cmd_dict['flagSample']
    reheader = cmd_dict['header'];launch_flag = cmd_dict['launch']


    ################################################################
    # Initializing all temporary output and analysis directories   #
    ################################################################
    
    tmp_dict = objC.processInit(work_dir,sys.argv,proj_date)
    tmp_dir = tmp_dict['tmpDir'];tmp_bin = tmp_dict['tmpBin']; 
    tmp_data = tmp_dict['tmpData']; tmp_status = tmp_dict['tmpStat']
    tmp_log = tmp_dict['tmpLog']; sb_log = tmp_dict['sbLog']

    
    objC = CONFIG()
    config_dict = objC.getConfigDict(xml_file)
    objS = CLUSTER()

    '''
    ################################################################
    # STEP-1: Merge VCF samples                                    #
    ################################################################
     
    cluster_file, swh = objS.getClusterWriteHandle(tmp_bin,'mergeVCF')

    # Remove the quotes/comments for adding SLURM-Cluster specific details
      
    swh = objS.writeClusterTop(config_dict,proj_date,vers,swh)
    script_out, script_err, swh = objS.writeClusterInit(config_dict,sb_log,
                                                    'mergeVCF',email_id,
                                                            exp_type,swh)

    if email_id and launch_flag:
        swh = objS.writeClusterSpecific(swh)
        swh = objS.writeClusterModule(config_dict,swh)

    

    tmp_stat_file = tmp_status+"/Job_status_merge.txt"

    # Function for merging VCF files. 
    tmp_merge_file, swh = objS.writeClusterMergeVCF(config_dict,data_file,
                                                        tmp_dir,tmp_log,
                                                      tmp_stat_file,swh)

    swh.close()

    '''

    ################################################################################
    # STEP-1: Steps and functions for launching SLURM/SHELL scripts per chromosome #
    ################################################################################
    
    for chr_num in chr_list:
        chr_dir = tmp_data
        tmp_data = tmp_data+"/chr"+str(chr_num)
        os.system('mkdir -p '+tmp_data)
        tmp_stat_file = tmp_status+"/Job_status_"+str(chr_num)+".txt"

        # Write the header and job specific name
        if exp_type == 'gvcfGT':
            cluster_file, swh = objS.getClusterWriteHandle(tmp_bin,
                                                       "genotypeAndAnnotate",
                                                       chr_num
                                                      )

        elif re.search('exeter',exp_type,re.IGNORECASE):
            cluster_file, swh = objS.getClusterWriteHandle(tmp_bin,
                                                       'splitAndAnnotate',
                                                       chr_num
                                                      )
        
        # Remove the quotes/comments for adding SLURM-Cluster specific details
        ''' 
        swh = objS.writeClusterTop(config_dict,proj_date,vers,swh)
        script_out, script_err, swh = objS.writeClusterInit(config_dict,
                                                          sb_log,chr_num,
                                                          email_id,exp_type,swh
                                                         )
        '''
        if email_id and launch_flag: 
            swh = objS.writeClusterTop(config_dict,proj_date,vers,swh)
            script_out, script_err, swh = objS.writeClusterInit(config_dict,
                                                          sb_log,chr_num,
                                                          email_id,exp_type,swh
                                                         )
            swh = objS.writeClusterSpecific(swh)
            swh = objS.writeClusterModule(config_dict,swh)
        
        # If the analysis task is joint genotyping using gvcfgenotyper (Illumina) 
        if exp_type =="gvcfGT":
            tmp_out_file, swh = objS.writeClusterGvcfGT(config_dict,data_file,
                                                      tmp_data,tmp_log,
                                                      tmp_stat_file,exp_type,
                                                      swh,chr_num,reheader)

        # Explicitly for joint merging of Exter samples (n=514). As a replacement 
        # for above joint genptyping using gvcfGenotype (Illumina).
        elif re.search('exeter',exp_type,re.IGNORECASE):
            tmp_merge_file = chr_dir+'/exeter.merged.bcf'
            tmp_out_file,swh = objS.writeClusterSplitByChrom(tmp_merge_file,
                                                           tmp_data,tmp_log,
                                                           tmp_stat_file,swh,
                                                           chr_num
                                                          )
            # Excludes gfvcfGenotyper step. Only uses left normalization steps            
            tmp_out_file, swh = objS.writeClusterGvcfGT(config_dict,data_file,
                                                            tmp_data,tmp_log,
                                                      tmp_stat_file,exp_type,
                                                      swh,chr_num,reheader
                                                     )

        # Extract Exonic region
        tmp_out_file, swh = objS.writeClusterRegion(config_dict,tmp_data,tmp_out_file,
                                                            tmp_stat_file,chr_num,swh
                                                   )

        # Annotate the extracted region using VEP
        tmp_out_file, swh = objS.writeClusterVEP(config_dict,tmp_data,tmp_out_file,
                                                           tmp_stat_file,chr_num,
                                                                 script_path,swh
                                              )

        # Add custom annotations
        tmp_out_file, swh = objS.writeClusterCustomAnnot(config_dict,tmp_data,
                                                       tmp_out_file,chr_num,swh
                                                      )

        # Fix ploidy
        #tmp_out_file = tmp_data+'/tmpFile_'+str(chr_num)+'_AC0.exon.vep.anno.bcf'
        tmp_out_file, swh = objS.writeClusterPloidy(config_dict,tmp_out_file,
                                                  chr_num,gender_file,tmp_stat_file,
                                                                exp_type,script_path,
                                                                    swh,flag_samples
                                                 )
        # Frequency filter
        tmp_out_file, swh = objS.writeClusterFreqFilter(config_dict,tmp_out_file,
                                                         chr_num,tmp_stat_file,
                                                               script_path,swh
                                                     )
        # Impact filter
        tmp_out_file, tmp_out_file_tab,\
                swh = objS.writeClusterImpactFilter(config_dict,tmp_out_file,
                                                     chr_num,tmp_stat_file,
                                                           script_path,swh
                                                 )
        # Variant Prioritization
        tmp_out_file_tab, swh = objS.writeClusterVariantPrioritization(config_dict,
                                                                     tmp_out_file_tab,
                                                                     tmp_stat_file,
                                                                     script_path,swh
                                                                    )
        # End status
        swh = objS.writeClusterEndStatus(tmp_stat_file,swh)
        swh.close()
      
        ''' 
        if to_launch:
            out = subprocess.check_output("sbatch "+cluster_file,shell=True)
            job_strs = re.split("\W|\s",out)
            job_id = job_strs[3]
            jobDict[job_id] = chr_num
        '''
        tmp_data = chr_dir
        #sys.exit()

    ########################################################################
    # STEP-2: Combining all chromosome output and family filtering         # 
    ########################################################################
    
    # combining all the chromosome output
    out_merge_file,swh = objS.writeClusterCombineFiles(config_dict,proj_date,
                                                     vers,tmp_dir,tmp_status,
                                                            "mergeAndFilter",
                                                            exp_type,email_id,
                                                                script_path)
    
    # Family filtering
    swh = objS.writeClusterFamilyFilter(config_dict,work_dir,tmp_dir,manifest,
                                          family_list,proj_date,script_path,
                                                        out_merge_file,swh)
    swh.close()




