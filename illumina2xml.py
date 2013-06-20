#!/usr/bin/env python

import xml.etree.ElementTree as ET
import sys
import os

# on command line first input is directory where Illumina files live
# files required for this script to run
# 1: ISA-Tab configuration files
# 2: DemultiplexConfig.xml
# 3: runParameters.xml

# directory where DemultiplexConfig.xml and runParameters.xml live
# test locally using "/Users/fcoldren/scripts_tools/my_scripts/Illumina_files"
d = sys.argv[1]


assay_config_file = "/Users/fcoldren/src/isatab-templates/ISATab_configuration_files_2012-06-20/transcription_seq.xml"
#study_config_file = "/Users/fcoldren/src/isatab-templates//Users/fcoldren/src/isatab-templates/ISATab_configuration_files_2012-06-20/studySample.xml"
#investigation_config_file ="/Users/fcoldren/src/isatab-templates//Users/fcoldren/src/isatab-templates/ISATab_configuration_files_2012-06-20/investigation.xml"

def check(dir):
    """ Checks that specified directory exists and is specified correctly """
    if os.path.exists(dir) is False:
        print dir," directory does not exist"
        quit()
    if dir.endswith("/") is False:
        dir = dir + "/"
    demultiplex_file_opt1 = dir + "Unaligned/DemultiplexConfig.xml"
    demultiplex_file_opt2 = dir + "DemultiplexConfig.xml"
    run_parameters_file = dir + "runParameters.xml"
    if (os.path.exists(demultiplex_file_opt1) is False) or (os.path.exists(run_parameters_file) is False):
        if (os.path.exists(demultiplex_file_opt2) is False) or (os.path.exists(run_parameters_file) is False):
            print "A required Illumina file is missing"
            quit()
    if (os.path.exists(demultiplex_file_opt1) is True):
        xml_files_to_parse = (demultiplex_file_opt1,run_parameters_file)
    elif (os.path.exists(demultiplex_file_opt2) is True):
        xml_files_to_parse = (demultiplex_file_opt2,run_parameters_file)
    print xml_files_to_parse
    return xml_files_to_parse


def make_ISATab_dir(dir):
    ISATab_dir = d + "/ISATAB/"
    if os.path.exists(ISATab_dir) is False:
        os.makedirs(ISATab_dir)
    return ISATab_dir

def make_new_sample_entry(root_config,root_experiment,sample_name):
    """Builds simplified xml structure based on an ISA-Tab xml config file"""
    sample_root = ET.SubElement(root_experiment,"Sample",name=sample_name)

    # element attribute value holds the column header text
    # element.get("value") returns the column header

    prefix = "{http://www.ebi.ac.uk/bii/isatab_configuration#}" #prefix for ISA tags
    unit_ = prefix + "unit-field"
    
    #makes the xml structure for a sample
    for element in root_config.iter():
        x = str(element.get("header"))
        if "header" in element.attrib and (x.endswith("Name") or x.endswith("File")):
            last_node = ET.SubElement(sample_root,"column-header",value=element.get("header"))
        elif "protocol-type" in element.attrib:
            last_node = ET.SubElement(sample_root,"column-header",value="Protocol REF",protocol=element.get("protocol-type"))
        else:
            if "header" in element.attrib:
                non_node_column = ET.SubElement(last_node,"column-header",value=element.get("header"))
                if element.get('is-forced-ontology') == 'true':
                    ET.SubElement(non_node_column,"column-header",value='Term Source REF')
                    ET.SubElement(non_node_column,"column-header",value='Term Accession Number')
            if element.tag == unit_:
                unit_node = ET.SubElement(last_node,"column-header",value='Unit')
                ET.SubElement(unit_node,"column-header",value='Term Source REF')
                ET.SubElement(unit_node,"column-header",value='Term Accession Number')
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

def make_sample_row_as_list(sample_element):
    sample_info = []
    for child in sample_element.iter("column-header"):
            if child.get('value') == "Sample Name":
                sample_info.append(j)
            elif child.get('value') == "Protocol REF":
                sample_info.append(child.get('protocol'))
            elif child.text is None:
                sample_info.append("")
            else:
                sample_info.append(child.text)
    return sample_info

illumina_xml_files = check(d) #check to see if the files exist and return file locations
output_directory = make_ISATab_dir(d)

# parse the assay config file
a_config_tree = ET.parse(assay_config_file)
ISA_assay_config_root = a_config_tree.getroot()

demultiplex_config_tree = ET.parse(illumina_xml_files[0])
demultiplex_root = demultiplex_config_tree.getroot()


# dictionary will store experiment wide parameters
parameters = {}

# starting the xml container for holding parsed data
experiment_root = ET.Element("Experiment")

# puts lane information into samples' tag
# and constructs the xml container for each sample
for element in demultiplex_root.iter('Lane'):
    for child in element:
        if child.get('Index') != "Undetermined":
            child.set('lane',element.get('Number'))
            sample_name = child.get('SampleId')
            assay = make_new_sample_entry(ISA_assay_config_root,experiment_root,sample_name)

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

# get experiment wide parameters from runParameters.xml
run = parse_runParameters4isa(illumina_xml_files[1],parameters)

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
all = []

output_assay_txt_file = output_directory + "/a_studyID_transcription profiling_nucleotide sequencing.txt"
outfile = open(output_assay_txt_file, "w")
for element in assay.iter("Sample"):
    j = element.get('name')
    sample_row = []
    if count == 0:
        headers = []
        for child in element.iter("column-header"):
            headers.append(child.get('value'))
        outfile.write("\t".join(headers) + "\n")
        count = 1
        sample_row = make_sample_row_as_list(element)
    else:
        sample_row = make_sample_row_as_list(element)
    all.append(sample_row)

for i in all:
    outfile.write("\t".join(i) + "\n")
outfile.close()
