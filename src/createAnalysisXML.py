#!/usr/bin/python

#################################################################################
#                                                                               #
# Description: Script to create user specific XML file                          #
#                                                                               #
# Input: (a) User specified config flat file.( See README for the description   #
#        (b) Base XML file provided default with this package                   #
#                                                                               #
# Output: User defined XML file                                                 #
#                                                                               #
# Usage: Please see the README file for usage.                                  #
#                                                                               #
#################################################################################

import re,sys,os

script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(script_path+'/CONFIG/')

from CONFIG import *
import dicttoxml
from xml.dom.minidom import parseString

if __name__=='__main__':

    # Process command line arguments
    objC = CONFIG()
    cmd_dict = objC.processXMLArguments()

    user_config = cmd_dict['userConfig']
    base_xml = cmd_dict['baseXml']
    out_xml = cmd_dict['outXml']

    # Read and update the Base Template XML file with user config file
    base_config_dict = objC.getConfigDict(base_xml)
    user_config_dict = objC.getUserConfigDict(user_config)
    # Update 
    base_config_dict = objC.updateBaseXML(base_config_dict,user_config_dict)

    # Convert upadted XML to user XML and print it out
    xml = dicttoxml.dicttoxml(base_config_dict,custom_root='config',attr_type=False)
    dom = parseString(xml)
    wh = open(out_xml,'w')
    print >>wh, dom.toprettyxml()
    wh.close()
    
    print '\n\n'


    
