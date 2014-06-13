#! /usr/bin/env python

import sys
import os
import glob
import re
import yaml
from collections import namedtuple

def expandOsPath(path):
    """
    To expand the path with shell variables.
    Arguments:
    - `path`: path string
    """
    return os.path.expanduser(os.path.expandvars(path))

def genFilesWithPattern(pathList, Pattern):
    """
    To generate files list on the fly.
    Arguments:
    - `pathList`: the path of the files
    - `Pattern`: pattern like config["input_files"]
    """
    pathList.append(Pattern)
    Files = glob.glob(expandOsPath(os.path.join(*pathList)))
    return Files

def parse_bowtie1_log(s):
    total_pattern = re.compile(r"""\#\s+reads\s+processed:\s(?P<total_reads>.+)\s*""",  # total_reads
                        re.VERBOSE)
    unique_mapped_pattern = re.compile("""\#\s+reads\s+with\s+at\s+least\s+one\s+reported\s+alignment:\s+(?P<unique_mapped_reads>\S+)\s+\(\S+\)""", # unique_mapped_reads
                        re.VERBOSE)
    for line in s:
        match = total_pattern.match(line)
        if match:
            total_reads = match.group("total_reads")
        match = unique_mapped_pattern.match(line)
        if match:
            unique_mapped_reads = match.group("unique_mapped_reads")
    res = namedtuple('res', ['total_reads', 'unique_mapped_reads'])
    r = res(total_reads=total_reads, unique_mapped_reads=unique_mapped_reads)
    return r

def parse_bowtie2_log(s):
    total_pattern = re.compile(r"""(?P<total_reads>\d+)\s+reads;\s+of\s+these:""",  # total_reads
                        re.VERBOSE)
    unique_mapped_pattern = re.compile("""\s*(?P<unique_mapped_reads>\d+)\s+\(\S+\).+exactly\s+1\s+time""", # unique_mapped_reads
                        re.VERBOSE)
    multiple_mapped_pattern = re.compile("""\s+(?P<multiple_mapped_reads>\d+)\s+\(\S+\).+aligned\s+>1\s+times""", # unique_mapped_reads
                        re.VERBOSE)
    for line in s:
        match = total_pattern.match(line)
        if match:
            total_reads = match.group("total_reads")
        match = unique_mapped_pattern.match(line)
        if match:
            unique_mapped_reads = match.group("unique_mapped_reads")
        match = multiple_mapped_pattern.match(line)
        if match:
            multiple_mapped_reads = match.group("multiple_mapped_reads")
    res = namedtuple('res', ['total_reads', 'unique_mapped_reads', 'multiple_mapped_reads'])
    r = res(total_reads=total_reads, 
        unique_mapped_reads=unique_mapped_reads, 
        multiple_mapped_reads=multiple_mapped_reads)
    return r

def parse_rmdup_log(s):
    pattern = re.compile(r'\[bam_rmdupse_core\]\s+(?P<dup_reads>\d+)\s/\s\d+', re.VERBOSE)
    for line in s:
        match = pattern.match(line)
        if match:
            dup_reads = match.group("dup_reads")
    res = namedtuple('res', ['dup_reads'])
    r = res(dup_reads=dup_reads)
    return r

def parse_phantomPeak_log(s):
    NSC_pattern = re.compile(r'.*\(NSC\)\s*(?P<NSC>\d*\.\d*).+', re.VERBOSE)
    RSC_pattern = re.compile(r'.*\(RSC\)\s*(?P<RSC>\d*\.\d*).+', re.VERBOSE)
    for line in s:
        match = NSC_pattern.match(line)
        if match:
            NSC = match.group("NSC")
        match = RSC_pattern.match(line)
        if match:
            RSC = match.group("RSC")
    res = namedtuple('res', ['NSC', 'RSC'])
    r = res(NSC=NSC, RSC=RSC)
    return r

def getSummaryFiles(input_type, config, search_paths):
    """
    Get all summary files under the folders.
    input_type: file types.
    config: config loaded from yaml.
    """
    input_type = "*" + input_type
    files = genFilesWithPattern([config["project_dir"], config["data_dir"]], input_type)
    for search_path in search_paths:
        files.extend(genFilesWithPattern([config["project_dir"], config["data_dir"], search_path],
            input_type))
    return files

def getFileId(file_basename):
    """
    Remove suffix of the summary file to get file id.
    """
    suffixes = ['.fastq.alignment.log', '.fq.alignment.log', '.gz.alignment.log', '.bam.rmdup.log', '.bam.phantomPeak.log']
    for suffix in suffixes:
        file_basename = file_basename.replace(suffix, '')
    return file_basename

## Search subdirectories under data folder.
search_paths = ["fastq", "rmdup"]

## Used for final results.
summary_dict = {}

## Load the same config yaml file of the pipeline.
config_name = sys.argv[1]
config_f = open(config_name, "r")
config = yaml.load(config_f)
config_f.close()

if config["aligner"] == "bowtie":
    ## To be used in debug
    # input_files = {".alignment.log":("total_reads", "unique_mapped_reads")}

    ## Summary files used for summarizing.
    input_files = {
    ".alignment.log":("total_reads", "unique_mapped_reads"), 
    ".rmdup.log":("dup_reads"), 
    ".phantomPeak.log":("NSC", "RSC")
    }

    ## Decide the parser here by a dict.
    parser_dict = {
    ".alignment.log": parse_bowtie1_log,
    ".rmdup.log": parse_rmdup_log,
    ".phantomPeak.log": parse_phantomPeak_log
    }

    ## Used to assign the output field in output file.
    output_header = ["sample", "total_reads", "unique_mapped_reads", "dup_reads", "NSC", "RSC"]
elif config["aligner"] == "bowtie2":
    ## to be used in debug
    # input_files = {".alignment.log":("total_reads", "unique_mapped_reads", "multiple_mapped_reads")}

    ## Summary files used for summarizing.
    input_files = {
    ".alignment.log":("total_reads", "unique_mapped_reads", "multiple_mapped_reads"), 
    ".rmdup.log":("dup_reads"), 
    ".phantomPeak.log":("NSC", "RSC")
    }

    ## Decide the parser here by a dict.
    parser_dict = {
    ".alignment.log": parse_bowtie2_log,
    ".rmdup.log": parse_rmdup_log,
    ".phantomPeak.log": parse_phantomPeak_log
    }

    ## Used to assign the output field in output file.
    output_header = ["sample", "total_reads", "unique_mapped_reads", "multiple_mapped_reads", "dup_reads", "NSC", "RSC"]

## Scan the files to summarize the pipeline.
for input_type, summary_types in input_files.items():
    summary_files = getSummaryFiles(input_type, config, search_paths)
    if len(summary_files) != 0:
        for summary_file in summary_files:
            file_id = getFileId(os.path.basename(summary_file))
            if file_id not in summary_dict:
                summary_dict[file_id] = {'sample':file_id}
            input_file = file(summary_file)
            lines = input_file.readlines()
            input_file.close()
            ## Here the value of the dict is the parser function!
            res = parser_dict[input_type](lines)
            ## Unpack the results into dict.
            for i in range(len(res._fields)):
                if res._fields[i] not in output_header:
                    output_header.append(res._fields[i])
                summary_dict[file_id][res._fields[i]] = res[i]

## Output to file, and the columns order is decided by output_header.
output_file = file("summary_stats.txt", "w")
header_line = "\t".join(output_header) + "\n"
output_file.write(header_line)
for sample in summary_dict.keys():
    output_list = []
    for stat in output_header:
        if stat in summary_dict[sample]:
            output_list.append(summary_dict[sample][stat])
        else:
            output_list.append("NA")
    line = "\t".join(output_list) + "\n"
    output_file.write(line)
output_file.close()
