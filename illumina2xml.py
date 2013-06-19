#!/usr/bin/env python

import xml.etree.ElementTree as ET
import sys

config_file = "/Users/fcoldren/src/isatab-templates/assay_configuration_files/transcription_seq.xml"
a_config_tree = ET.parse(config_file)
root = a_config_tree.getroot()
prefix = "{http://www.ebi.ac.uk/bii/isatab_configuration#}"
field_ = prefix + "field"
unit_ = prefix + "unit-field"
experiment_root = ET.Element("Experiment")

def make_new_sample_entry(root_config,root_experiment,sample_name):
    """Builds simplified xml based on an ISA-Tab xml config file"""
    sample_root = ET.SubElement(root_experiment,"Sample",name=sample_name)

    # element attribute value holds the column header text
    # element.get("value") returns the column header

    #makes the xml structure for a sample
    for element in root.iter():
        x = str(element.get("header"))
        if "header" in element.attrib and (x.endswith("Name") or x.endswith("File")):
            last_node = ET.SubElement(sample_root,"column-header",value=element.get("header"))
        elif "protocol-type" in element.attrib:
            last_node = ET.SubElement(sample_root,"column-header",value="Protocol REF",protocol=element.get("protocol-type"))
        else:
            if "header" in element.attrib:
                ET.SubElement(last_node,"column-header",value=element.get("header"))
            if element.tag == unit_:
                ET.SubElement(last_node,"column-header",value='Unit', unit_attr1 = 'Term Source REF', unit_attr2 = 'Term Accession Number')
    return root_experiment

# "/Users/fcoldren/scripts_tools/my_scripts/Illumina_files/DemultiplexConfig.xml"
# this is the DemultiplexConfig.xml
f1 = "/Users/fcoldren/scripts_tools/my_scripts/Illumina_files/DemultiplexConfig.xml"
demultiplex_config_tree = ET.parse(f1)
demultiplex_root = demultiplex_config_tree.getroot()
list_of_samples = []
for element in demultiplex_root.iter('Lane'):
    for child in element:
        if child.get('Index') != "Undetermined":
            child.set('lane',element.get('Number'))
            sample_name = child.get('SampleId')
            assay = make_new_sample_entry(root,experiment_root,sample_name)


for element in demultiplex_root.iter('Sample'):
    for item in assay.iter('Sample'):
        if item.get('name') == str(element.get('SampleId')):
            for child in item.iter():
                if child.get('value') == "Comment[PathBio submission]":
                    child.text = element.get('ProjectId')
                if child.get('value') == "Parameter Value[mid]":
                    child.text = element.get('Index')
                if child.get('value') == "Parameter Value[lane]":
                    child.text = element.get('lane')


ET.dump(assay)

        



