#!/usr/bin/env python

import xml.etree.ElementTree as ET
import sys

# files required for this script to run
# 1: ISA-Tab configuration file
# 2: DemultiplexConfig.xml
# 3: runParameters.xml

config_file = "/Users/fcoldren/src/isatab-templates/assay_configuration_files/transcription_seq.xml"
f1 = "/Users/fcoldren/scripts_tools/my_scripts/Illumina_files/DemultiplexConfig.xml"
f2 = "/Users/fcoldren/scripts_tools/my_scripts/Illumina_files/runParameters.xml"

a_config_tree = ET.parse(config_file)
root = a_config_tree.getroot()
prefix = "{http://www.ebi.ac.uk/bii/isatab_configuration#}"
field_ = prefix + "field"
unit_ = prefix + "unit-field"
experiment_root = ET.Element("Experiment")

# "/Users/fcoldren/scripts_tools/my_scripts/Illumina_files/DemultiplexConfig.xml"
# this is the DemultiplexConfig.xml
parameters = {}
demultiplex_config_tree = ET.parse(f1)
demultiplex_root = demultiplex_config_tree.getroot()


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
                unit_node = ET.SubElement(last_node,"column-header",value='Unit', unit_attr1 = 'Term Source REF', unit_attr2 = 'Term Accession Number')
                ET.SubElement(unit_node,"column-header",value='Term Source REF'
    return root_experiment

def parse_runParameters4isa(file,run_parameters_to_write):
    """ Parses some info from the Illumina file runParameters.xml """
    runparameters_config_tree = ET.parse(file)
    runparameters_root = runparameters_config_tree.getroot()
    p = []
    for element in runparameters_root.iter():
        if element.tag == "Pe":
            p.append(((element.text).strip()).split('\n')) # to clean up weird whitespace in this text
            run_parameters_to_write["Parameter Value[cBot clustering kit]"] = p[0][0]
        if element.tag == "RunStartDate":
            d = element.text
            year = d[:2]
            month = d[2:4]
            day = d[4:]
            date = "20" + year + "-" + month + "-" + day
            run_parameters_to_write["Date"] = date # date for sequencing
        if element.tag == "ScannerID":
            run_parameters_to_write["Parameter Value[sequencing instrument ID]"] = element.text
        if element.tag == "ApplicationName":
            software = element.text
        if element.tag == "ApplicationVersion":
            software_version = str(element.text)
        if element.tag == "RTAVersion":
            rtaversion = element.text
        if element.tag == "FCPosition":
            run_parameters_to_write["Parameter Value[slot]"] = element.text
    software_full = software + " " + software_version
    rta = "RTA" + " " + rtaversion
    run_parameters_to_write["Parameter Value[sequencing instrument control software]"] = software_full
    run_parameters_to_write["Parameter Value[base caller]"] = rta
    run_parameters_to_write["Parameter Value[quality score]"] = rta
    return run_parameters_to_write

# puts lane information into samples' tag
for element in demultiplex_root.iter('Lane'):
    for child in element:
        if child.get('Index') != "Undetermined":
            child.set('lane',element.get('Number'))
            sample_name = child.get('SampleId')
            assay = make_new_sample_entry(root,experiment_root,sample_name)

# get non-sample specific info from DemultiplexConfig.xml
for element in demultiplex_root:
    if element.tag == "FlowcellInfo":
        parameters["Parameter Value[flow cell]"] = element.get('ID')
    if element.tag == "Software":
        parameters["Parameter Value[demultiplex]"] = element.get('Version')

# this will not fill in the field Parameter Value[sample1]
# parses details from DemultiplexConfig.xml into xml structure
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

run = parse_runParameters4isa(f2,parameters)

# fill in experiment wide parameters parsed from runParameters.xml
# and DemultiplexConfig.xml to the xml structure for the samples

for element in assay.iter():
    for key in run:
        if key != "Date" and element.get('value') == key:
            element.text = run[key]

for element in assay.iter('column-header'):
    if element.get('protocol') == "nucleic acid sequencing":
        for child in element:
            if child.get('value') == "Date":
                child.text = run["Date"]

count = 0
# write headers to file
for element in assay.iter("Sample"):
    j = element.get('name')
    if count == 0:
        headers = []
        for child in element.iter():
            headers.append(child.get('value'))
        



