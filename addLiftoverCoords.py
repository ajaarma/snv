#!usr/bin/env python

#
# this script adds annotation information from a file 
#

import sys
import optparse 
import os
import pdb
import re
import shutil

#############
# CONSTANTS #
#############


#################
# END CONSTANTS #
#################


class OptionParser(optparse.OptionParser):
    """
    Adding a method for required arguments.
    Taken from:
    http://www.python.org/doc/2.3/lib/optparse-extending-examples.html
    """
    def check_required(self, opt):
        option = self.get_option(opt)

        # Assumes the option's 'default' is set to None!
        if getattr(self.values, option.dest) is None:
            print "%s option not supplied" % option
            self.print_help()
            sys.exit(1)

def main():
	
    opt_parser = OptionParser()
   
    opt_parser.add_option("--annot_file",
                          dest="annot_file",
                          type="string",
                          help="annotation file",
                          default=None)
    opt_parser.add_option("--in_file",
                          dest="in_file",
                          type="string",
                          help="file to add annot to",
                          default=None)
    opt_parser.add_option("--new_col",
                          dest="new_col",
                          type="string",
                          help="new column header",
                          default=None)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("--annot_file")
    opt_parser.check_required("--in_file")
    opt_parser.check_required("--new_col")

    annot_file = open(options.annot_file)
    in_file = open(options.in_file)
    new_col = options.new_col
  
    h_fl = True   
    annot = {}                             # annot[key] = [annots]

    for line in annot_file:
        line=formatLine(line)
        line = line.split("\t")
        n = line[0] + ":" + line[1]
        if line[3] not in annot:
            annot[line[3]] = n
        else:
            a = annot[line[3]].split(",")
            f = 0
            for aa in a:
                if n == aa:
                    f = 1
            if not f:
                annot[line[3]] = annot[line[3]] + "," + n
    annot_file.close()

    for ln in in_file:
        ln=formatLine(ln)
        ln = ln.split("\t")
        if h_fl:
            ln.append(new_col)
            print "\t".join(ln)
            h_fl = False
            continue

        k = ln[0] + ":" + ln[1]
        
        if k in annot:
            ln.append(annot[k])
        else:
            ln.append("NA")
                
        print "\t".join(ln)

    in_file.close()


    
    sys.exit(0)



############
# END_MAIN #
############

def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line
                             
#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
