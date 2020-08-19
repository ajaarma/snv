#!/usr/bin/python

#########################################################################################
#                                                                                       # 
# Description: Generic class for creating the Slurm or Shell script                     #
#                                                                                       #
#########################################################################################

import re,sys,os,subprocess,datetime,time

class SLURM:

    def __init__(self, elements=[]):
        self.__elements={}
        for e in elements:
            self.__elements[e]=1

    def display(self):
        print 'Inside SLURM class. Creating SLURM/SHELL script for launching jobs'+\
               'in the cluster'

    def getSlurmWriteHandle(self,bin_dir,slurm_name,chr_num=[]):
        ''' Subroutine to generate SLURM/SHELL script file name '''

        if len(chr_num)!=0:
            slurm_file = os.path.abspath(bin_dir+"/"+slurm_name+"_chr"+chr_num+".sh")
            wh = open(slurm_file,"w")
        else:
            slurm_file = os.path.abspath(bin_dir+"/"+slurm_name+".sh")
            wh = open(slurm_file,"w")

        return slurm_file, wh

    def writeSlurmTop(self,config_dict,proj_name,vers,wh):
        ''' Subroutine to write the top header part of the SLURM/SHELL batch 
            script'''

        now = datetime.datetime.now()
        print >>wh,"#!/bin/bash"
        print >>wh,"#!"
        print >>wh,"## date: "+"-".join([str(now.year),str(now.month),str(now.day)])
        print >>wh,"## project name: "+proj_name
        print >>wh,"## pipeline version: "+vers
        print >>wh,"## genome version: "+config_dict["general"]["genomeBuild"]
        print >>wh,"##\n\n"

        return wh

    def writeSlurmInit(self,config_dict,sb_log,chr_num,email_id,gt_type,wh):
        ''' Subroutine to provide initial details for SLURM batch script. 
            Automatically extracted from user provided XML file '''
        
        slurm_dict = config_dict["slurm"]
        script_out = sb_log+"/"+gt_type+"_"+str(chr_num)+".o.txt"
        script_err = sb_log+"/"+gt_type+"_"+str(chr_num)+".e.txt"
        
        for slr in slurm_dict:
            dash = slurm_dict["dash"]
            doubDash = slurm_dict["doubDash"]
            if slr=="params1":

                par1_dict = slurm_dict[slr]
                for par1 in par1_dict:
                    if par1=="J":
                        print >>wh,"#SBATCH "+dash+par1+" "+str(chr_num)+"_"+par1_dict[par1]
                    else:
                        print >>wh,"#SBATCH "+dash+par1+" "+par1_dict[par1]

            if slr=="params2":
                par2_dict = slurm_dict[slr]
                for par2 in par2_dict:
                    try:
                        if par2=="mail-user":
                            print >>wh,"#SBATCH "+doubDash+par2+"="+email_id
                        elif par2=="output":
                            print >>wh,"#SBATCH "+doubDash+par2+"="+script_out
                        elif par2=="error":
                            print >>wh,"#SBATCH "+doubDash+par2+"="+script_err
                        else:
                            print >>wh,"#SBATCH "+doubDash+par2+par2_dict[par2]
                    except:
                        print >>wh,"#SBATCH "+doubDash+par2

                print >>wh,"\n"
        return script_out, script_err, wh

    def writeSlurmSpecific(self,wh):
        
        ''' Subroutine for writing SLURM specific details in the batch script '''

        print >>wh,"#! Number of nodes and tasks per node allocated by SLURM (do not change):"
        print >>wh,"numnodes=$SLURM_JOB_NUM_NODES"
        print >>wh,"numtasks=$SLURM_NTASKS"
        print >>wh,"mpi_tasks_per_node=$(echo \"$SLURM_TASKS_PER_NODE\" | sed -e  's/^\([0-9][0-9]*\).*$/\\1/')"
        print >>wh,"\n"
        
        return wh

    def writeSlurmModule(self,config_dict,wh):
        
        ''' Subroutine for writing SLURM specific module details. Automatically 
            extracted from user provided XML file. '''

        mod_list = config_dict["module"]["value"]
        print >>wh,". /etc/profile.d/modules.sh"
        for modKey in mod_list:
            print >>wh,"module "+modKey
        print >>wh,"\n"

        print >>wh, 'source activate snv'
        
        '''
        if config_dict["export"]:
            exp = config_dict["export"]["value"]
            if isinstance(exp, basestring):
                exp_list = [exp]
            else:
                exp_list = exp
            for expKey in exp_list:
                print >>wh, expKey
        ''' 

        print >>wh,"\n\n"
        print >>wh,"JOBID=$SLURM_JOB_ID"
        print >>wh,"echo -e \"JobID: $JOBID\n======\""
        print >>wh,"echo \"Time: `date`\""
        print >>wh,"echo \"Running on master node: `hostname`\""
        print >>wh,"echo \"Current directory: `pwd`\n\""
        #print >>wh,"echo -e \"numtasks=$numtasks, numnodes=$numnodes, mpi_tasks_per_node=$mpi_tasks_per_node (OMP_NUM_THREADS=$OMP_NUM_THREADS)\""
        print >>wh,"\n\n"

        return wh

    def writeSlurmNGCvsXTRCompare(self,config_dict,inp_data,tmp_data,tmp_log,
                                        tmp_stat_file,exp_type,swh,chr_num):

        ''' Subroutine to count variants at each stage of the snv step '''

        os.chdir(inp_data)
        file_list = sorted(os.listdir(inp_data),key=os.path.getmtime)

        for ele in file_list:
            print ele
            if ele.endswith('.bcf'):
                path_file = os.path.join(inp_data,ele)
                cmd = 'a=`bcftools view -H '+path_file+' | wc -l` '
                cmd_out = 'echo \"'+ele+' : \" $a >> '+tmp_data
                print >>swh,cmd
                print >>swh,cmd_out
            elif ele.endswith('.tab'):
                path_file = os.path.join(inp_data,ele)
                cmd = 'a=`wc -l '+path_file+' `'
                cmd_out = 'echo \"'+ele+' : \" $a >> '+tmp_data
                print >>swh,cmd
                print >>swh,cmd_out
            elif ele.endswith('.gz'):
                path_file = os.path.join(inp_data,ele)
                cmd = 'a=`bcftools view -H '+path_file+' | wc -l `'
                cmd_out = 'echo \"'+ele+' : \" $a >> '+tmp_data
                print >>swh,cmd
                print >>swh,cmd_out 
       
        return swh


    def writeSlurmMergeVCF(self,config_dict,data_file,tmp_dir,tmp_log,
                                                    tmp_stat_file,wh):

        ''' Subroutine to merge the bgzip compressed VCF files '''

        merge_out_file = tmp_dir+'/tmp_data/exeter.merged.bcf'
        
        print >>wh,'############# START: Merging Exeter VCF files ###########'

        fh = open(data_file)
        sid_list = []

        for lines in fh:
            lines = lines.strip()
            sid_list.append(lines)

        fh.close()

        cmd_merge = 'bcftools merge -m all '+' '.join(sid_list)+' -Ob -o '+\
                                                                    merge_out_file
        cmd_index = 'bcftools index '+merge_out_file
        print >>wh,cmd_merge
        print >>wh,'\n',cmd_index,'\n'

        wh = self.checkErrorFile(merge_out_file,wh,'error',tmp_stat_file)

        print >>wh,'############# END: Merging Exeter VCF files #############'

        return merge_out_file,wh

    def writeSlurmSplitByChrom(self,tmp_merge_file,tmp_data,tmp_log,
                                                    tmp_stat_file,wh,chr_num):
        
        ''' Subroutine to split the Merged bcf file into individual 
            chromosomes '''

        inputD = tmp_merge_file
        outD = tmp_data
        tmp_out_file = outD+'/exeter_tmp_chr'+str(chr_num)+'.bcf'

        if re.search('chr',chr_num,re.IGNORECASE):
            chr_num = chr_num
        else:
            chr_num = chr_num

        print >>wh,'############# START: Split by Chromosome ################ '
        cmd_split_chr = 'bcftools view -r '+str(chr_num)+' '+tmp_merge_file+' -Ob -o '+tmp_out_file
        print >>wh,cmd_split_chr

        cmd_split_chr_index = 'bcftools index '+tmp_out_file
        print >>wh,cmd_split_chr_index

        wh = self.checkErrorFile(tmp_out_file,wh,'error',tmp_stat_file)

        return tmp_out_file,wh


    def writeSlurmGvcfGT(self,config_dict,data_file,tmp_data,tmp_log,tmp_stat_file,
                                                     exp_type,wh,chr_num,reheader):

        ''' Subroutine to implement running up of joint genotyping using gvcfGenotyper.
            Also incorporates normlaization using bcftools 
        '''
    
        tmp_out_fileA = tmp_data+'/'+exp_type+"_tmp_chr"+str(chr_num)+".bcf"
        tmp_log_file = tmp_log+"/"+exp_type+"_tmp_chr"+str(chr_num)+".log"
        tmp_out_file = tmp_data+"/tmpFile_"+str(chr_num)+"_AC0.bcf"

        #ref_seq = config_dict["general"]["refDir"]+"/"+config_dict["general"]["refGen"]
        ref_seq = '/'.join([
                            config_dict["general"]["resourceDir"],
                            config_dict["general"]["refFasta"]
                           ])
                          
        
        gvcfDict = config_dict["gvcfGenotyper"]
        par_dict = gvcfDict["params"]
        
        if exp_type == 'gvcfGT':
            gvcfDict = config_dict["gvcfGenotyper"]
            command = []

            print >>wh,"############# START: Illumina gvcfGenotyper Step ###############\n"
            print >>wh,"echo \"[`date`] Running Illumina gvcfGenoTyper\"\n"

            binary = gvcfDict["binary"]
            dash = gvcfDict["dash"]
            par_dict = gvcfDict["params"]
            command.append(binary)
                    
            for pd in par_dict:
                if pd=="l":
                    cmd = dash+pd+" "+data_file
                    command.append(cmd)
                if pd=="O":
                    cmd = dash+pd+par_dict[pd]
                    command.append(cmd)
                if pd=="o":
                    cmd = dash+pd+" "+tmp_out_fileA
                    command.append(cmd)
                if pd=="f":
                    cmd = dash+pd+" "+ref_seq
                    command.append(cmd)
                if pd=="r":
                    chrom = str(chr_num)
                    if "C" in par_dict and par_dict["C"] != None and par_dict["C"] != "":
                        chrom = par_dict["C"]+str(chr_num)
                    if chr_num == "MT" and "M" in par_dict and par_dict["M"] != None and par_dict["M"] != "":
                        chrom = par_dict["M"]
                    cmd = dash+pd+" "+chrom
                    command.append(cmd)
                if pd=="L":
                    cmd = dash+pd+" "+tmp_log_file
                    command.append(cmd)

            print >>wh," ".join(command)

            cmd_bcf_1 = "\nbcftools index "+tmp_out_fileA
            print >>wh,cmd_bcf_1

        print >>wh,"\necho \"[`date`] Normalize indels and filter for AC>0\""
        cmd1 = "\nbcftools norm -m - -f "+ref_seq+" "+tmp_out_fileA
        cmd2, cmd3 = "", ""
        if "M" in par_dict and par_dict["M"] != None and par_dict["M"] != "":
            cmd2 = "|sed 's/"+par_dict["M"]+"/MT/'"
        if "C" in par_dict and par_dict["C"] != None and par_dict["C"] != "":
            cmd3 = "|sed 's/="+par_dict["C"]+"/=/'|sed 's/^"+par_dict["C"]+"//'"
        cmd4 = "|bcftools view -e 'AC==0' -Ob -o "+tmp_out_file

        print >>wh,cmd1+cmd2+cmd3+cmd4
        print >>wh,"bcftools index -f "+tmp_out_file

        if reheader:
            print >>wh,"\necho \"[`date`] Reheader step\""
            print >>wh,"\nbcftools reheader -s "+reheader+" "+tmp_out_file+" > "+tmp_out_file+".tmp"
            print >>wh,"mv "+tmp_out_file+".tmp "+tmp_out_file            
            print >>wh,"bcftools index -f "+tmp_out_file          

        wh = self.checkErrorFile(tmp_out_file,wh,"error",tmp_stat_file)
        print >>wh,"\n########## END: Joint genotyping (gfvcfHenotyper)/Merging Step ##########"
        return tmp_out_file, wh

    
    def writeSlurmRegion(self,config_dict,tmp_data,tmp_out_file,tmp_stat_file,
                                                          chr_num,exp_type,wh):

        ''' Subroutine to extract the exonic regions '''

        prefix_path = re.split("\.bcf",tmp_out_file)[0]
        exon_file = prefix_path+".exon.bcf"
        sorted_file = prefix_path+".exon.sorted.bcf"

        regFile = '/'.join([ 
                            config_dict["general"]["resourceDir"],
                            config_dict['regions']['exons']
                           ])

        print >>wh,"\n\n############## START: Extract Exonic Region ################\n"
        
        ''' if exp_type=='exeter':
            norm_file = prefix_path+'.norm.vcf'
            print >>wh,'echo \"[`date`] normalizing vcf for multi-allelic site\" \n'
            cmd_norm = 'bcftools norm -m - '+tmp_out_file+' -Ob -o '+norm_file
            cmd_norm_index = 'bcftools index '+norm_file
            print >>wh, cmd_norm
            print >>wh, cmd_norm_index
            wh = self.checkErrorFile(norm_file,wh,"error",tmp_stat_file)
            tmp_out_file = norm_file
            exon_file = prefix_path+'.norm.exon.bcf'
            sorted_file = prefix_path+".norm.exon.sorted.bcf"
        ''' 

        print >>wh,"echo \"[`date`] Getting exonic regions\" \n"
        print >>wh,"bcftools view -R "+regFile+" "+tmp_out_file+" -Ob -o "+exon_file
        print >>wh,"bcftools view "+exon_file+" | bcftools sort | bcftools view -Ob -o "+sorted_file
        print >>wh,"mv "+sorted_file+" "+exon_file
        print >>wh,"bcftools index "+exon_file
      
        wh = self.checkErrorFile(exon_file,wh,"error",tmp_stat_file)
        print >>wh,"\n############ END: Extract Exonic Region  ####################"

        return exon_file, wh

    
    def writeSlurmPloidy(self,config_dict,tmp_out_file,chr_num,gender_file,
                                tmp_stat_file,gt_type,script_path,wh,flag_samples):

        ''' Subroutine for implementing fix ploidy steps '''

        prefix_path = re.split("\.bcf",tmp_out_file)[0]

        pf_in = tmp_out_file
        pf_out = prefix_path+".sites.bcf"
        pl_tmp = prefix_path+".sites.ploidy.tmp"
        pl_out = prefix_path+".sites.ploidy.bcf"
        fr_tmp = prefix_path+".sites.ploidy.ngcfreq.tmp"
        fr_out = prefix_path+".sites.ploidy.ngcfreq.bcf"

        print >>wh,"\n\n########## START: Add Flags, Fix Ploidy, and add internal cohort frequency steps ##########\n "
        print >>wh,"\necho \"[`date`] Adding pass flags\"\n"
        if gt_type=="gvcfGT":
            cmd1 = "bcftools annotate -x FILTER -Ou "+pf_in+" | bcftools +fill-AN-AC | bcftools filter -e \'AN<(1*N_SAMPLES)\' -s \'LOWCALL\' -m + | bcftools filter -i \'FT==\"PASS\"\' -s \'LOWPF\' -m + -Ob > "+pf_out
        
        elif gt_type == 'exeter':
            cmd1 = "bcftools annotate -x FILTER -Ou "+pf_in+" | bcftools +fill-AN-AC | bcftools filter -e \'AN<(1*N_SAMPLES)\' -s \'LOWCALL\' -m 'LOWPF\' -m + -Ob > "+pf_out
        cmd1_index = "\nbcftools index "+pf_out
        print >>wh,cmd1
        print >>wh,cmd1_index
        wh = self.checkErrorFile(pf_out,wh,"error",tmp_stat_file)

        print >>wh,"\necho \"[`date`] Fixploidy step\""

        if chr_num == "X":
            cmd1, cmd2, cmd3, cmd4, cmd5, cmd6, cmd7 = [], [], [], [], [], [], []

            #cmd1.append("s=`grep M$ "+gender_file+" | cut -f1 -d' ' | tr '\\n' ',' | sed 's/,$//'` && bcftools view -Ob -s $s "+pf_out+" >"+pl_tmp+".males.bcf") # general
            cmd1.append("s=`grep M$ "+gender_file+" | cut -f1 | tr '\\n' ',' | sed 's/,$//'` && bcftools view -Ob -s $s "+pf_out+" >"+pl_tmp+".males.bcf") #Exeter
            #cmd1.append("s=`grep M$ "+gender_file+" | cut -f1 -d' ' | tr '\\n' ',' | sed 's/,$//'` && bcftools view -Ob -s ^${s} "+pf_out+" >"+pl_tmp+".females.bcf") #general
            cmd1.append("s=`grep M$ "+gender_file+" | cut -f1 | tr '\\n' ',' | sed 's/,$//'` && bcftools view -Ob -s ^${s} "+pf_out+" >"+pl_tmp+".females.bcf") #Exeter
            cmd1.append("bcftools index "+pl_tmp+".males.bcf")
            cmd1.append("bcftools index "+pl_tmp+".females.bcf")

            cmd2.append("bcftools view -Ob -g het -T "+config_dict["general"]["resourceDir"]+"/"+config_dict['regions']['PAR']+" "+pl_tmp+".males.bcf >"+pl_tmp+".males.nonPARhet.bcf")
            cmd2.append("bcftools index "+pl_tmp+".males.nonPARhet.bcf")

            cmd3.append("bcftools isec -Ob -n~10 "+pl_tmp+".males.bcf "+pl_tmp+".males.nonPARhet.bcf -w 1 > "+pl_tmp+".males.NOTnonPARhet.bcf")
            cmd3.append("bcftools index "+pl_tmp+".males.NOTnonPARhet.bcf")

            cmd4.append("bcftools view "+pl_tmp+".males.nonPARhet.bcf | python "+script_path+"/makeHaploid.py --in_file - | bcftools view -Ob >"+pl_tmp+".males.nonPARhet.fixed.bcf")
            cmd4.append("bcftools index "+pl_tmp+".males.nonPARhet.fixed.bcf")

            cmd5.append("bcftools concat -Ob "+pl_tmp+".males.nonPARhet.fixed.bcf "+pl_tmp+".males.NOTnonPARhet.bcf > "+pl_tmp+".males.fixed.bcf")
            cmd5.append("bcftools sort -Ob "+pl_tmp+".males.fixed.bcf > "+pl_tmp+".males.fixed.sorted.bcf")
            cmd5.append("bcftools index "+pl_tmp+".males.fixed.sorted.bcf")

            cmd6.append("bcftools merge -Ob "+pl_tmp+".males.fixed.sorted.bcf "+pl_tmp+".females.bcf >"+pl_tmp+".fixed.bcf")
            cmd6.append("bcftools index "+pl_tmp+".fixed.bcf")
        
            cmd7.append("bcftools view "+pl_tmp+".fixed.bcf -s `bcftools query -l "+pf_out+" | xargs | tr \' \' \',\'` -Ob -o "+pl_out)
            cmd7.append("bcftools index "+pl_out)

            print >>wh,"\n\n".join(["\n".join(cmd1), "\n".join(cmd2), "\n".join(cmd3), "\n".join(cmd4), "\n".join(cmd5), "\n".join(cmd6), "\n".join(cmd7)])
            print >>wh,"rm "+pl_tmp+"*"
        else:
            print >>wh,"cp "+pf_out+" "+pl_out
            print >>wh,"bcftools index "+pl_out
            
      
        wh = self.checkErrorFile(pl_out,wh,"error",tmp_stat_file)
        print >>wh,"\necho \"[`date`] Adding Internal cohort (NGC) AFs, excluding flagged samples\"\n"

        if flag_samples:
            cmd4 = "bcftools view -S ^"+flag_samples+" "+pl_out
        else:
            cmd4 = "bcftools view "+pl_out
        cmd4_sed1 = "sed '/^#CHROM/i ##INFO=<ID=NGC_AF,Number=A,Type=Float,Description=\"NGC alternate allele frequencies\">'"
        cmd4_sed2 = "sed '/^#CHROM/i ##INFO=<ID=NGC_AC,Number=A,Type=Integer,Description=\"NGC allele count in genotypes\">'"
        cmd4_sed3 = "sed '/^#CHROM/i ##INFO=<ID=NGC_AN,Number=1,Type=Integer,Description=\"NGC total number of alleles in called genotypes\">'"
        cmd4_sed4 = "perl -F\\\\t -ane 'print && next if /^#/; @in = split /;/, $F[7]; for $i (@in) {$an = $1 if $i =~ /^AN=([\d]+)/; $ac = $1 if $i =~ /^AC=([\d]+)/;} $af = 0; $af =$ac/$an if $an > 0;$F[7] .= \";NGC_AN=$an;NGC_AC=$ac;NGC_AF=$af\"; print join(\"\\t\", @F)'"
        cmd4_sed5 = "bcftools view - -Ob -o "+fr_tmp
        cmd4 = " | ".join([cmd4,cmd4_sed1,cmd4_sed2,cmd4_sed3,cmd4_sed4,cmd4_sed5])

        cmd4_index = "\nbcftools index "+fr_tmp

        print >>wh,cmd4
        print >>wh,cmd4_index

        print >>wh, "\necho \"[`date`] Annotating original file with Internal Cohort AFs\""

        cmd5 = "\nbcftools annotate -Ob -o "+fr_out+" -c INFO/NGC_AC,INFO/NGC_AN,INFO/NGC_AF -a "+fr_tmp+" "+pl_out
        cmd5_index = "\nbcftools index "+fr_out
        print >>wh,cmd5
        print >>wh,cmd5_index
        print >>wh,"rm "+fr_tmp+"*"

        wh = self.checkErrorFile(fr_out,wh,"error",tmp_stat_file)
        
        #print >>wh,"a=`diff <(bcftools view -H "+fr_out +" | cut -f1,2,4,5 | sort | uniq) <(bcftools view -H -r "+chr_num+" "+config_dict['general']['prevVars']+" | cut -f1,2,4,5 | sort | uniq) | grep \"^>\" | wc -l`"
        #print >>wh,"if [  $a != 0 ]; then"
        #print >>wh,"  echo -e \"Warning: $a missing variants in "+fr_out+"\""
        #print >>wh,"fi"
        print >>wh,"\n\n########## END: Add Flags, Fix Ploidy, and add internal cohort frequency steps ##########\n "

        return fr_out, wh        
        

    def writeSlurmCustomAnnot(self,config_dict,tmp_data,tmp_out_file,chr_num,wh):
        
        ''' Subroutine to add custom annotations : gnomAD,ExAC,HGMD, Clinvar '''

        refDir = config_dict["general"]["resourceDir"]
        prefix_path = re.split("\.bcf",tmp_out_file)[0]
        tmp_file = prefix_path+".tmp.bcf"
        exon_tmp_anno_file = prefix_path+".anno.bcf"
        
        cust_annot = config_dict["annotation"]["cust-annot"]
        cmd_all = []

        first = 1
        for ann in cust_annot:
            if first:
                inp = tmp_out_file
                first = 0
            else:
                inp = exon_tmp_anno_file
                
            if "tag" in cust_annot[ann]:
                c = ",".join(["INFO/" + cust_annot[ann]["tag"] + "_" + i + ":=INFO/" + i for i in cust_annot[ann]["fields"].split(',')])
            elif cust_annot[ann]["fields"] == "INFO":
                c = cust_annot[ann]["fields"]
            else:
                c = ",".join(["INFO/" + i for i in cust_annot[ann]["fields"].split(',')])

            if "vcfDir" in cust_annot[ann]:
                vcfDir = cust_annot[ann]["vcfDir"]
            else:
                vcfDir = refDir

            vcfFile = cust_annot[ann]["vcf"]
            if "chr" in cust_annot[ann]:
                vcfFile = cust_annot[ann]["vcf"] % chr_num
                if "excl_chr" in cust_annot[ann] and chr_num in cust_annot[ann]["excl_chr"].split(","):
                    if "repl_chr" in cust_annot[ann]:
                        vcfFile = cust_annot[ann]["vcf"] % cust_annot[ann]["repl_chr"]
                    else:
                        continue                    

            cmd = "bcftools annotate -c " + c + " -a " + vcfDir + "/" + vcfFile + " -r " + chr_num + " " + inp + " -Ob -o " + tmp_file
            cmd_all.append(cmd)

        print >>wh,"\n\n########### START: Add Custom Annotations ############### \n"
        print >>wh,"\necho \"[`date`] Adding custom annotations\"\n"

        for cmd in cmd_all:
            print >>wh, cmd
            print >>wh, "mv -f %s %s" % (tmp_file, exon_tmp_anno_file)
            print >>wh, "bcftools index -f %s\n" % exon_tmp_anno_file

        print >>wh,"\n\n########### END: Add Custom Annotations ############### \n"
        
        return exon_tmp_anno_file, wh
 

    def writeSlurmFreqFilter(self,config_dict,tmp_out_file,chr_num,tmp_stat_file,
                                                                  script_path,wh):
        
        ''' Subroutine for implementing frequency filtering steps '''

        prefix_path = re.split("\.bcf",tmp_out_file)[0]
        exon_tmp_anno_file_001_tmp = prefix_path+".0_01.tmp.bcf"
        exon_tmp_anno_file_001_tmp2 = prefix_path+".0_01.tmp2.bcf"
        exon_tmp_anno_file_001_tmp3 = prefix_path+".0_01.tmp3.bcf"
        exon_tmp_anno_file_001 = prefix_path+".0_01.bcf"

        refDir = config_dict['general']['resourceDir']
        bcf_bin = config_dict["freqImpact"]["bcfBinary"]
        dash = config_dict["freqImpact"]["dash"]

        exclDict = config_dict["freqImpact"]["exclude"]
        cmd_all = []

        first = 1
        for keys in exclDict:
            exclStr = "\'"+keys+exclDict[keys]+"\'"
            if first:
                cmd = bcf_bin+" view -Ob -e "+exclStr+" "+tmp_out_file
                cmd_all.append(cmd)
                first = 0
            else:
                cmd = bcf_bin+" view -Ob -e "+exclStr
                cmd_all.append(cmd)

        cmd = " | ".join(cmd_all)+ " >"+exon_tmp_anno_file_001_tmp
           
        if "keep" in config_dict["freqImpact"] and config_dict["freqImpact"]["keep"] != "":

            cmd2 = bcf_bin+" view -Ob -R "+'/'.join([refDir,config_dict["freqImpact"]["keep"]])+\
                                    " "+tmp_out_file+" > "+exon_tmp_anno_file_001_tmp2
            cmd2 = cmd2+"\n"+bcf_bin+" concat -Ob "+exon_tmp_anno_file_001_tmp+" "+\
                              exon_tmp_anno_file_001_tmp2+" > "+exon_tmp_anno_file_001_tmp3
            cmd2 = cmd2+"\n"+bcf_bin+" view "+exon_tmp_anno_file_001_tmp3+\
                                " | bcftools sort | bcftools view -Ob -o "+exon_tmp_anno_file_001
        else:
            cmd2 = "cp %s %s" % (exon_tmp_anno_file_001_tmp, exon_tmp_anno_file_001)

    
        print >>wh,"\n\n########## START: Frequency Filtering Step #############"
        print >>wh,"\necho \"[`date`] Frequency filtering \"\n"
        print >>wh, cmd
        print >>wh, "bcftools index -f %s\n" % exon_tmp_anno_file_001_tmp            
        print >>wh, cmd2
        print >>wh, "bcftools index -f %s\n" % exon_tmp_anno_file_001            

        wh = self.checkErrorFile(exon_tmp_anno_file_001,wh,"error",tmp_stat_file)
        
        print >>wh,"\n########### END: Frequency Filtering Step #############"

        return exon_tmp_anno_file_001, wh

                
    def writeSlurmVEP(self,config_dict,tmp_data,tmp_out_file,tmp_stat_file,
                                                        chr_num,script_path,wh):

        ''' Subroutine to implementing VEP annotation '''
        
        prefix_path = re.split("\.bcf",tmp_out_file)[0]
        tmp1 = prefix_path+".noGT.vcf.gz"
        tmp2 = prefix_path+".noGT.new.vcf.gz"
        vep_file_new = prefix_path+".noGT.new.vep.vcf.gz"
        stats_file = prefix_path+".noGT.vep_stats.txt"
        tmp3 = prefix_path+".vep.tmp.bcf"
        exon_anno_file_001 = prefix_path+".vep.bcf"

        print >>wh,"\n\n########### START: VEP Annotation #################### \n"

        print >>wh,"echo \"[`date`] Getting vcf without GT and pulling new variants\n\""

        print >>wh,"bcftools view -G "+tmp_out_file+" -Oz -o "+tmp1
        print >>wh,"bcftools index "+tmp1

        if "prevVEP" in config_dict["annotation"] and config_dict["annotation"]["prevVEP"] != None and config_dict["annotation"]["prevVEP"] != "":
            vep_file = config_dict["annotation"]["prevVEP"]
            print >>wh, "bcftools isec -Oz -n~10 "+tmp1+" "+config_dict["annotation"]["prevVEP"]+" -w 1 > "+tmp2
            print >>wh, "bcftools index "+tmp2
        else:
            vep_file = tmp1
            tmp2 = tmp1

        print >>wh,"\necho \"[`date`] Running VEP\"\n"

        ref_seq = '/'.join([
                            config_dict["general"]["resourceDir"],
                            config_dict["general"]["refFasta"]
                           ])
        
        refGenome = config_dict['general']['genomeBuild']
        defPar = []
        allPar = []
        plugIn = []

        refDir = config_dict["general"]["resourceDir"]
        outformat = config_dict["annotation"]["outFormat"] 
        defDict = config_dict["annotation"]["defPar"]
        allDict = config_dict["annotation"]["all"]
        plugInDict = config_dict["annotation"]["plugIn"]
        doubDash = config_dict["annotation"]["doubDash"]
        vep_bin = config_dict["annotation"]["vep-bin"]
        
        for keys,values in defDict.items():
            if re.search('^dir\_cache',keys):
                defPar.append(doubDash+keys)
                defPar.append('/'.join([os.path.dirname(refDir),values,refGenome]))
            elif keys=="fasta":
                defPar.append(doubDash+keys)
                defPar.append(ref_seq)
            elif defDict[keys]==None:
                defPar.append(doubDash+keys)
            else:
                defPar.append(doubDash+keys)
                defPar.append(defDict[keys])

        for keys in allDict:
            if allDict[keys]==None:
                allPar.append(doubDash+keys)
            else:
                allPar.append(doubDash+keys)
                allPar.append(allDict[keys])

        for keys,values in plugInDict.items():
            if re.search('CADD',keys):
                plugin_str = ' '.join(['--plugin',','.join([keys,'/'.join([
                                                    refDir,values['cadd_snv']
                                                    ]),
                                                                '/'.join([
                                                    refDir,values['cadd_indel']
                                                    ])
                                            ])
                                      ])
            elif values !='.':
                plugin_str = ' '.join(['--plugin',','.join([keys,'/'.join([
                                                        refDir,values
                                                    ])
                                                ]) 
                                      ]) 
            else:
                plugin_str = ' '.join(['--plugin',keys])

            plugIn.append(plugin_str)

        ''' 
        for keys in plugInDict:
            c = [keys]
            if plugInDict[keys]!=None:
                pre = refDir
                if "dir" in plugInDict[keys]:
                    pre = plugInDict[keys]["dir"]
                for k in plugInDict[keys]:
                    if k == "dir":
                        continue
                    if k.startswith("file"):
                        c.append(pre+plugInDict[keys][k])
                    else:
                        c.append(plugInDict[keys][k])                        
            plugIn.append(doubDash+"plugin "+",".join(c))
            
        '''

        cmd = "bcftools view -r "+str(chr_num)+" "+tmp2+" | "+vep_bin+" --"+\
                outformat+" -o STDOUT "+" ".join(defPar)+" "+" ".join(allPar)+" "+\
                " ".join(plugIn)
        cmd = cmd+" --stats_text --stats_file "+stats_file
        cmd = cmd+" | bgzip -c > "+vep_file_new
        cmd_index = "bcftools index -f "+vep_file_new
        print >>wh,cmd
        print >>wh,cmd_index
        
        print >>wh,"\necho \"[`date`] Annotating bcf with VEP results\"\n"
        if "prevVEP" in config_dict["annotation"] and config_dict["annotation"]["prevVEP"] != None and config_dict["annotation"]["prevVEP"] != "":
            print >>wh,"bcftools annotate -c INFO/ANN -a "+vep_file+" "+tmp_out_file+" -Ob -o "+tmp3
            print >>wh,"bcftools index "+tmp3
            print >>wh,"\nbcftools annotate -c INFO/ANN -a "+vep_file_new+" "+tmp3+" -Ob -o "+exon_anno_file_001
            print >>wh,"bcftools index "+exon_anno_file_001
        else:
            print >>wh,"\nbcftools annotate -c INFO/ANN -a "+vep_file_new+" "+tmp_out_file+" -Ob -o "+exon_anno_file_001
            print >>wh,"bcftools index "+exon_anno_file_001            
        
        wh = self.checkErrorFile(exon_anno_file_001,wh,"error",tmp_stat_file)
        print >>wh,"\n############### END: VEP Annotation  #################### \n"

        return exon_anno_file_001, wh
    

    def writeSlurmImpactFilter(self,config_dict,tmp_out_file,chr_num,tmp_stat_file,
                                                                    script_path,wh):

        ''' Subroutine for implementing impact filtering steps '''

        prefix_path = re.split("\.bcf",tmp_out_file)[0]
        out1 = prefix_path+".imp.bcf"
        out2 = prefix_path+".imp.tab"

        v2t_bin = script_path+"/"+config_dict["freqImpact"]["vcf2tab"]
        bcf_bin = config_dict["freqImpact"]["bcfBinary"]
        dash = config_dict["freqImpact"]["dash"]
    
        extRegList = config_dict["freqImpact"]["extReg"]["value"]['item']
        extRegStr = "\'"+"|".join(extRegList)+"\'"
        
        cmd_bcf = bcf_bin+" view "+tmp_out_file+" | grep -E "+extRegStr+" > "+out1
        cmd_v2t = "python "+v2t_bin +" -t "+out1+" -o "+out2

        print >>wh,"\n\n########## START: Impact Filtering and vcf2tab Step ##############"
        print >>wh,"\necho \"[`date`] Impact filtering \"\n"
        print >>wh,cmd_bcf
        wh = self.checkErrorFile(out1,wh,"error",tmp_stat_file)
        print >>wh,"\necho \"[`date`] Converting VCF to TAB \"\n"
        print >>wh,cmd_v2t
        wh = self.checkErrorFile(out2,wh,"error",tmp_stat_file)
        print >>wh, "\na=`grep -v \"^#\" "+ out1 +" | cut -f1,2,4,5 | sort | uniq | wc -l`"
        print >>wh, "b=`grep -v \"^CHROM\" "+ out2 +" | cut -f1,2,4,5 | sort | uniq | wc -l`"
        
        print >>wh,"if [  $a != $b ]; then"
        print >>wh,"  echo -e \"\\nWarning: missing variants in "+out2+" : only ${b} of ${a}\""
        #print >>wh,"  exit 1"
        print >>wh,"fi"
        print >>wh,"\n\n########## END: Impact Filtering and vcf2tab Step ##############"

        return out1, out2, wh

    
    def writeSlurmVariantPrioritization(self,config_dict,tmp_out_file_tab,
                                                 tmp_stat_file,script_path,wh):

        ''' Subroutine for implementing variant prioritization '''

        prefix_path = re.split("\.tab",tmp_out_file_tab)[0]
        out_file = prefix_path+".prior.tab"

        print >>wh,"\n\n######## START: Variant Prioritization Step ############"
        print >>wh,"echo \"[`date`] Performing variant prioritization\"\n"

        R_script = script_path+"/"+config_dict["varPrior"]["vpBin"]
        varPriorDict = config_dict["varPrior"]["params"]
        dash = config_dict["varPrior"]["dash"]

        params_list = []

        for keys in varPriorDict:
            if keys=="v":
                params_list.append(dash+keys+" "+tmp_out_file_tab)
            elif keys=="o":
                params_list.append(dash+keys+" "+out_file)
            elif keys=="l":
                params_list.append(dash+keys+" "+config_dict["general"]["resourceDir"]+"/"+varPriorDict[keys])
            else:
                params_list.append(dash+keys+" "+varPriorDict[keys])

        print >>wh,"\nRscript "+R_script+" "+" ".join(params_list)

        wh = self.checkErrorFile(out_file,wh,"error",tmp_stat_file)

        print >>wh,"c=`awk -F'\\t' ' { for (i = 1; i <= NF; ++i) print i, $i; exit } ' "+tmp_out_file_tab+" | grep Gene | cut -f1 -d' '`"
        print >>wh,"a=`cut -f1,2,4,5,$c "+tmp_out_file_tab +" | sort | uniq | wc -l`"
        print >>wh,"b=`wc -l "+out_file+" | cut -f1 -d' '`"
        print >>wh,"if [  $a != $b ]; then"
        print >>wh,"  echo -e \"\\nWarning: missing variants in "+out_file+" : only ${b} of ${a}\""
        print >>wh,"  echo \"0\" | cat > "+tmp_stat_file
        print >>wh,"  exit 1"
        print >>wh,"  else "
        print >>wh,"        echo \"1\" | cat > "+tmp_stat_file
        print >>wh,"fi"
        print >>wh,"\necho \"[`date`] Variant Prioritization Completed\""

        print >>wh,"\n\n######## END: Variant Prioritization Step ############"

        return out_file, wh

    def writeSlurmCombineFiles(self,config_dict,proj_name,vers,tmp_dir,
                               tmp_status,slurm_name,gt_type,email_id,script_path):

        ''' Subroutine to combine all the output files split across chromosomes'''

        objS = SLURM()
        tmp_bin = tmp_dir+"/tmp_binaries"
        tmp_data = tmp_dir+"/tmp_data"
        sb_log =  tmp_dir+"/sb_log"
        tmp_stat_file = tmp_status+"/Job_status_"+slurm_name+".txt"

        slurm_file, wh = objS.getSlurmWriteHandle(tmp_bin,slurm_name)
        wh = objS.writeSlurmTop(config_dict,proj_name,vers,wh)
        script_out, script_err, wh = objS.writeSlurmInit(config_dict,sb_log,
                                                    "CombChr",email_id,gt_type,wh)
        wh = objS.writeSlurmSpecific(wh)
        wh = objS.writeSlurmModule(config_dict,wh)
        
        print >>wh,"\n\n######## START: Combining Variants for All Chromsomes ###########"
        
        print >>wh,"bcftools concat `ls -v "+tmp_data+"/chr*/tmpFile_*_"+config_dict["combChromTab"]["bcfSuffix"]+".bcf` -Ob -o "+tmp_data+"/merged.all."+config_dict["combChromTab"]["bcfSuffix"]+".bcf"
        print >>wh,"bcftools index "+tmp_data+"/merged.all."+config_dict["combChromTab"]["bcfSuffix"]+".bcf"

        R_script = script_path+"/"+config_dict["combChromTab"]["combBin"]
        out_file = "merged.all."+config_dict["combChromTab"]["tabSuffix"]+".tab"

        print >>wh,"\n\n####### END: Combining Prioritized Variants for All Chromsomes #######"
        print >>wh,"Rscript "+R_script+" "+tmp_data+" "+out_file

        out_file_path = os.path.abspath(tmp_data+"/"+out_file)

        wh = self.checkErrorFile(out_file_path, wh,"error",tmp_stat_file)

        return out_file_path,wh


    def writeSlurmFamilyFilter(self,config_dict,work_dir,tmp_dir,manifest,
                                        family_list,proj_name,script_path,
                                                        out_merge_file,wh):

        ''' Subroutine to implement family inheritance filtering steps '''

        print >>wh,"\n########## START: Family Filtering #############\n"
        print >>wh,"echo \"Family Filtering started\"\n"
        famDict = config_dict["famFilter"]
        famScript = script_path+"/"+famDict["famBin"]
        liftoverScript = script_path+"/"+famDict["liftoverBin"]
        crossmap_path = '$CONDA_PREFIX/bin/'+famDict["crossmap"]
        chain_path = config_dict["general"]["resourceDir"]+"/"+famDict["chain"]
        new_col =famDict["newCol"]
        var_file = out_merge_file
        ped_file = manifest
        col_file = script_path+"/"+famDict["columns"]
        fam_file = family_list
        gene_list = config_dict['general']['resourceDir']+'/'+famDict["genelist"]
        haem_file = config_dict['general']['resourceDir']+'/'+famDict["haemGenesFile"]
        imprint_file = config_dict['general']['resourceDir']+'/'+famDict["imprintedGenesFile"]
        poly_file = config_dict['general']['resourceDir']+'/'+famDict["polymorphicGenesFile"]
        pheno_file = famDict["phenotypes"]
        hpo = config_dict["general"]["resourceDir"]+"/"+famDict["hpo"]
        omim = config_dict["general"]["resourceDir"]+"/"+famDict["omim"]
        fam_dir = tmp_dir+"/fam_filter"
         
        fh = open(fam_file)
        for lines in fh:
            lines = lines.strip()
            strs = re.split("\t|\s",lines)
            fam_id = strs[0].strip()
            fam_id_path = os.path.abspath(fam_dir+"/"+fam_id+"/filtering")
            os.system("mkdir -p "+fam_id_path)

        freqScript = ""
        if "freqScript" in famDict and famDict["freqScript"] != None and famDict["freqScript"] != "":
            freqScript = " --freq_script "+famDict["freqScript"]
        inhGenes = ""
        if "inheritedGenelist" in famDict and famDict["inheritedGenelist"] != None and famDict["inheritedGenelist"] != "":
            inhGenes = " --imprinted_genes "+famDict["inheritedGenelist"]

        cmd = "for i in `cut -f1 -d' ' "+fam_file+"`;do echo \"Rscript "+famScript+" --vars "+var_file+" --ped "+ped_file+" --family $i --genelist "+gene_list+" --pheno "+pheno_file+" --hpodb "+hpo+" --omim "+omim+" --columns "+col_file+" --haem_genes "+haem_file+" --imprinted_genes "+imprint_file+" --poly_genes "+poly_file+inhGenes+freqScript+" --out "+fam_dir+"/$i/filtering/${i}.filt_"+proj_name+".txt 2>&1 | tee "+tmp_dir+"/tmp_log/Fam_filter.$i.log\";done | parallel -j15 --no-notice "

        print >>wh,cmd
        
        cmd1 = "for i in `cut -f1 -d' ' "+fam_file+"`;do echo \"sh "+liftoverScript+" "+fam_dir+"/$i/filtering/${i}.filt_"+proj_name+".txt "+script_path+" "+crossmap_path+" "+chain_path+" "+new_col+"\";done | parallel -j15 --no-notice "
        
        print >>wh, "\n"
        print >>wh, cmd1 + "\n"
        

        print >>wh,"\n\n############# END: Family Filtering ###############\n\n"
        return wh

        fh.close()


    def checkErrorFile(self,out_file,wh,flag,status_file=[]):

        ''' Subroutine to implement Error/Warning status '''
        
        if_cmd = "\nif [ ! -e "+out_file+" ]; then \n"
        if flag=="error":
            echo_1 = "   echo -e \"ERROR: "+out_file+" doesnot exist\" \n"
            echo_2 = "   echo \"0\" | cat > \""+status_file+"\"\n"
            echo_3 = "   exit 1\n"
            echo_4 = "   else\n"
            echo_5 = "          echo \"1\" | cat > \""+status_file+"\""
        end_cmd = "\nfi"
        print >>wh,if_cmd+echo_1+echo_2+echo_3+echo_4+echo_5+end_cmd

        return wh

    def writeSlurmEndStatus(self,tmp_stat_file,wh):

        ''' Subroutine to implement End status of the script'''

        print >>wh, "echo \"2\" | cat > "+tmp_stat_file

        return wh

    
