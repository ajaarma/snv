#!usr/bin/env python

#
# this script edits all genotypes to hemizygous and
# adds MOS to GT filter field if they were het
# [0/1, 1/0, 1/1] > 1, [0/0] > 0, [./., 1/., ./1] > .
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
   
    opt_parser.add_option("--in_file",
                          dest="in_file",
                          type="string",
                          help="in file (- for stdin)",
                          default=None)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("--in_file")

    if options.in_file == "-":
        in_file = sys.stdin
    else:
        in_file = open(options.in_file)
  
    for line in in_file:
        if line.startswith('#'):
            print formatLine(line)
            continue
        line=formatLine(line)
        line = line.split("\t")
        nline = line[0:10]
        for i in range(10, len(line)):
            nline.append(makeHaploid(line[i]))
        print "\t".join(nline)

    if options.in_file != "-":
        in_file.close() 
   
def makeHaploid(gt):
    g = gt.split(':')
    if len(g[0].split('/')) == 1:
        return gt
    g[0], g[1] = newGTFT(g[0], g[1])
    return ":".join(g)

def newGTFT(gg, ff):
    ngts = {'0/1':'1',
                '1/0':'1',
                '1/1':'1',
                '1/.':'1',
                './1':'1',
                '0/0':'0',
                './.':'.'
                }
    mos = ['0/1', '1/0', '1/.', './1']
    if gg not in ngts:
        print >>sys.stderr, "WARNING: %s not changed" % gg
        return gg, ff
    
    ## Commented out for Exeter cohort; Should be modified 
    ## for generalized version
    
    #if gg in mos:
    #    if ff == '.':
    #        nf = 'MOS'
    #    else:
    #        nf = ff + ";MOS"
    #else:
    nf = ff
    return ngts[gg], nf
    
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
