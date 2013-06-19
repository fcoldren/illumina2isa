#!/usr/bin/env python

# makes a layout from an ISA-Tab configuration file in txt

import xml.etree.ElementTree as ET
import sys

# "/Users/fcoldren/src/isatab-templates/assay_configuration_files/transcription_seq.xml"
# use file above to generate RNA-seq template

# first command line input is the absolute path to ISA-Tab configuration file
config_file = sys.argv[1]

def create_ISATAB_txt(ISATABconfig):
    """Creates a list of column headers for the ISA-Tab txt file"""
    a_config_tree = ET.parse(ISATABconfig)
    root = a_config_tree.getroot()
    prefix = "{http://www.ebi.ac.uk/bii/isatab_configuration#}"
    field_ = prefix + "field"
    unit_ = prefix + "unit-field"
    column_headers = []

    for element in root.iter():
        if "header" in element.attrib:
            column_headers.append(element.get('header'))
            if element.get('is-forced-ontology') == 'true':
                column_headers.append("Term Source REF")
                column_headers.append("Term Accession Number")
        if element.tag == unit_:
            column_headers.append("Unit")
            column_headers.append("Term Source REF")
            column_headers.append("Term Accession Number")
        elif "protocol-type" in element.attrib:
            column_headers.append("Protocol REF")
    return column_headers

file_column_headers = create_ISATAB_txt(config_file)

# still need to clean up the output name and location option
outfile = open("/Users/fcoldren/Desktop/a_studyName_transcription_seq.txt","w")
outfile.write("\t".join(file_column_headers))
outfile.close()
