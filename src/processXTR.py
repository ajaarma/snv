#!/usr/bin/python
#####################################################################################
#                                                                                   #
# Description: Script to process Exeter pedigree mapping file, VCF files. Outputs   #
#              VCF files that are normalized                                        #
#              Header of AD field in VCF is changed to numberic for ease of merging #
#              Creates a manifest file                                              #
#                                                                                   #
#####################################################################################

import re,sys,os

script_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(script_path+'/config/')

from config import *


if __name__=='__main__':

    # Process Command line arguments
    objC = CONFIG()
    cmd_dict = objC.processXTRArguments() 

    xml_file = cmd_dict['xml']
    xtr_map_file = cmd_dict['xtrMap']
    vcf_dir = cmd_dict['vcfDir']
    out_dir = cmd_dict['outDir']
    exp_type = cmd_dict['expType']

    # Read in the Configuration XML file
    config_dict = objC.getConfigDict(xml_file)

    # Reference genome file
    grch37_ref = '/'.join([ config_dict['general']['resourceDir'],
                            config_dict['general']['refFasta']
                          ]
                         )

    '''
    grch37_ref = '/rds/project/flr24/rds-flr24-wgs10k-ngc/WGS10K/data/NGC/us/'+\
                 'resources//genomes/ensembl_b37_74/Homo_sapiens.'+\
                 'GRCh37.74.dna.fasta'
    '''
    
    # Create internal subdirectories to store normalized,compressed 
    # and manifest file
    bgzip_dir = os.path.abspath(out_dir+'/bgzip_comp')
    bgzip_norm_dir = os.path.abspath(out_dir+'/bgzip_norm')
    manifest_dir = os.path.abspath(out_dir+'/manifest')
    log_dir = os.path.abspath(out_dir+'/logs')

    if exp_type == 'norm':
        norm_dir = os.path.abspath(out_dir+'/bgzip_norm')
        os.system('mkdir -p '+norm_dir)

    os.system('mkdir -p '+manifest_dir) 
    os.system('mkdir -p '+bgzip_dir)
    os.system('mkdir -p '+log_dir)
    #os.system('module load samtools/1.3')

    lh = open(log_dir+'/Commands_log.txt','w')
    mh = open(manifest_dir+'/manifest.txt','w')
    manifest_head = ['family_id','ngc_id','father_id','mother_id','gender_declared',
                    'affected','ilmn_id','family_size','load_hpc_date','bam_path',
                    'vcf_path']
    
    print >>mh,'\t'.join(manifest_head)

    genderDict,vcfDict,mh = objC.createManifestFile(xtr_map_file,vcf_dir,out_dir,mh)
    mh.close()

    # Gender information
    gh = open(manifest_dir+'/gender_genomic.txt','w')

    for keys,values in genderDict.items():
        print >>gh,keys+' '+values

    gh.close()

    # Creating all sample vcf fof
    fof = open(out_dir+'/exeter_samples.fof','w')
    
    count = 1
    for keys in vcfDict.keys():
        
        bgzip_file = bgzip_dir+'/'+os.path.basename(keys)+'.gz'
        bgzip_cmd = 'bgzip -c '+keys+' > '+bgzip_file
        tabix_cmd = 'tabix -p vcf '+bgzip_file

        if exp_type =='norm':
            
            print >>lh,'######### bgzip compress the vcf file: '+keys+' #########'
            print >>lh, bgzip_cmd+' ; '+tabix_cmd
            os.system(bgzip_cmd+' ; '+tabix_cmd)
            tmp_file = bgzip_norm_dir+'/tmp.'+os.path.basename(keys)
            
            # Assigne numberic variable type to Allelic depth field (AD)
            # for exeter samples
            head_tmp_file_cmd = 'bcftools view -h '+bgzip_file+' | '+\
                    'sed \'s/##FORMAT=<ID=AD,Number=./'+\
                             '##FORMAT=<ID=AD,Number=R/\''+' > '+tmp_file
            
            rest_concat_cmd = 'bcftools view -H '+bgzip_file +' >> '+tmp_file
           
            print >>lh,'######### Edit the header information #########'
            print >>lh, head_tmp_file_cmd+' ; '+rest_concat_cmd
            os.system(head_tmp_file_cmd+' ; '+rest_concat_cmd)

            bgzip_norm_file = bgzip_norm_dir+'/'+os.path.basename(keys)+'.gz'
            norm_cmd = 'bcftools norm -m - -f '+grch37_ref+' '+tmp_file+\
                        ' | bgzip -c > '+bgzip_norm_file
            norm_cmd_index = 'bcftools index '+bgzip_norm_file
          
            print >>lh,'########## Nornalize the VCF for bi/multi-allelic sites #########'
            print >>lh, norm_cmd+' ; '+norm_cmd_index+' ; rm '+tmp_file
            print >>lh,'\n'
            os.system(norm_cmd+' ; '+norm_cmd_index+' ; rm '+tmp_file)
            print >>fof,bgzip_norm_file
        
        count = count+1

    fof.close()
    lh.close()
    gh.close()
   
    
    print '\n\n'

