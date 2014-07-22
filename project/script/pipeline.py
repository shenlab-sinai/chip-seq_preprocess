#! /usr/bin/env python

import os
import sys
import yaml
from ruffus import *
import glob
import subprocess
import string

def expandOsPath(path):
    """
    To expand the path with shell variables.
    Arguments:
    - `path`: path string
    """
    return os.path.expanduser(os.path.expandvars(path))

def genFilesWithPattern(pathList, Pattern):
    """
    To generate files list on the fly based on wildcards pattern.
    Arguments:
    - `pathList`: the path of the files
    - `Pattern`: pattern like config["input_files"]
    """
    pathList.append(Pattern)
    Files = expandOsPath(os.path.join(*pathList))
    return Files

config_name = sys.argv[1]
config_f = open(config_name, "r")
config = yaml.load(config_f)
config_f.close()
inputfiles = expandOsPath(os.path.join(config["project_dir"], config["data_dir"], "fastq", config["input_files"]))
FqFiles = [x for x in glob.glob(inputfiles)]
fq_name, fq_ext = os.path.splitext(config["input_files"])
fq_ext_suffix = ".alignment.log"
Bam_path = expandOsPath(os.path.join(config["project_dir"], config["data_dir"])) + "/"
FastQC_path = expandOsPath(os.path.join(config["project_dir"], config["data_dir"], "FastQC"))
rmdup_path = expandOsPath(os.path.join(config["project_dir"], config["data_dir"], "rmdup"))

@transform(FqFiles, formatter(fq_ext), os.path.join(Bam_path, "{basename[0]}.bam"), config)
def alignFastqByBowtie(FqFileName, OutputBamFileName, config):
    """
    To align '.fastq' to genome.
    Arguments:
    - `FqFileName`: file to be processed
    """
    if "aligner" in config:
        if config["aligner"] == "bowtie":
            cmds = ['fastq2bam_by_bowtie.sh']
            cmds.append(FqFileName)
            cmds.append(expandOsPath(config['bowtie_index']))
        elif config["aligner"] == "bowtie2":
            cmds = ['fastq2bam_by_bowtie2.sh']
            cmds.append(FqFileName)
            cmds.append(config['bowtie_index'])
        else:
            raise KeyError
    else:
        cmds = ['fastq2bam_by_bowtie.sh']
        cmds.append(FqFileName)
        cmds.append(expandOsPath(config['bowtie_index']))

    target = expandOsPath(os.path.join(config["project_dir"], config["data_dir"]))
    cmds.append(target)
    cmds.append(config["pair_end"])
    cores = int(int(config['cores'])/len(FqFiles))
    if cores == 0:
        cores = 1
    cmds.append(str(cores))
    logfile = FqFileName + ".alignment.log"
    p = subprocess.Popen(
        cmds, stdout=open(logfile, "w"), stderr=open(logfile, "w"),
        bufsize=1)
    stdout, stderr = p.communicate()
    return stdout

@follows(alignFastqByBowtie, mkdir(FastQC_path))
@transform(alignFastqByBowtie, suffix(".bam"), ".bam.fastqc.log", config)
def runFastqc(BamFileName, fastqcLog, config):
    """
    To run FastQC
    Arguments:
    - `BamFileName`: bam file
    - `config`: config
    """
    cmds = ['fastqc']
    cmds.append("-o")
    cmds.append(expandOsPath(os.path.join(config["project_dir"], config["data_dir"], "FastQC")))
    if "fastqc_threads" in config:
        cmds.append("-t")
        cmds.append(str(config["fastqc_threads"]))
    else:
        cmds.append("-t")
        cmds.append("2")
    cmds.append(BamFileName)
    logfile = BamFileName + ".fastqc.log"
    p = subprocess.Popen(
        cmds, stdout=open(logfile, "w"), stderr=open(logfile, "w"),
        bufsize=1)
    stdout, stderr = p.communicate()
    return stdout

@follows(runFastqc, mkdir(rmdup_path))
@transform(alignFastqByBowtie, formatter(".bam"), os.path.join(rmdup_path, "{basename[0]}.bam"), config)
def rmdupBam(BamFileName, rmdupFile, config):
    """
    To run FastQC
    Arguments:
    - `BamFileName`: bam file
    - `config`: config
    """
    if config["pair_end"]=="no":
        cmds = ['rmdup.bam.sh']
    else:
        cmds = ['rmdup_PE.bam.sh']
    cmds.append(BamFileName)
    cmds.append(rmdup_path)
    if "bam_sort_buff" in config:
        cmds.append(config["bam_sort_buff"])
    logfile = BamFileName + ".rmdup.log"
    p = subprocess.Popen(
        cmds, stdout=open(logfile, "w"), stderr=open(logfile, "w"),
        bufsize=1)
    stdout, stderr = p.communicate()
    return stdout

@follows(rmdupBam, mkdir(expandOsPath(os.path.join(rmdup_path, "tdf"))))
@transform(rmdupBam, suffix(".bam"), ".bam.tdf.log", config)
def genTDF(BamFileName, tdfLog, config):
    """
    To generate TDF files for IGV
    Arguments:
    - `BamFileName`: bam file
    - `config`: config
    """
    cmds = ['igvtools']
    cmds.append("count")
    cmds.append(BamFileName)
    TDFPath = expandOsPath(os.path.join(rmdup_path, "tdf"))
    baseName = os.path.basename(BamFileName)
    cmds.append(os.path.join(TDFPath, baseName.replace(".bam", ".tdf")))
    cmds.append(config["IGV_genome"])
    logfile = BamFileName + ".tdf.log"
    p = subprocess.Popen(
        cmds, stdout=open(logfile, "w"), stderr=open(logfile, "w"),
        bufsize=1)
    stdout, stderr = p.communicate()
    return stdout

@follows(rmdupBam)
@transform(rmdupBam, suffix(".bam"), ".bam.phantomPeak.log", config)
def runPhantomPeak(BamFileName, Log, config):
    """
    To check data with phantomPeak
    Arguments:
    - `BamFileName`: bam file
    - `config`: config
    """
    cmds = ['runPhantomPeak.sh']
    cmds.append(BamFileName)
    logfile = BamFileName + ".phantomPeak.log"
    p = subprocess.Popen(
        cmds, stdout=open(logfile, "w"), stderr=open(logfile, "w"),
        bufsize=1)
    stdout, stderr = p.communicate()
    return stdout

@follows(runPhantomPeak, genTDF)
@merge(rmdupBam, expandOsPath(os.path.join(rmdup_path, config["project_name"]+".ngsplot.all.log")), config)
def runNgsplotAll(BamFileNames, Log, config):
    """
    To check data with phantomPeak
    Arguments:
    - `BamFileName`: bam file
    - `config`: config
    """
    cmds = ['ngsplot_all.sh']
    cmds.append(rmdup_path)
    cmds.append(config["ngsplot_genome"])
    cmds.append(config["project_name"])
    cmds.append(str(config["cores"]))
    cmds.append(str(config["ngsplot_fraglen"]))
    logfile = expandOsPath(os.path.join(rmdup_path, config["project_name"]+".ngsplot.all.log"))
    p = subprocess.Popen(
        cmds, stdout=open(logfile, "w"), stderr=open(logfile, "w"),
        bufsize=1)
    stdout, stderr = p.communicate()
    return stdout

## Run to FastQC step
# pipeline_run([runFastqc], multiprocess=config["cores"])

## Run the whole pipeline
pipeline_run([runNgsplotAll], multiprocess=config["cores"])

## Plot the pipeline flowchart
# pipeline_printout_graph("all_flowchart.png", "png", [runNgsplotAll], pipeline_name="Preprocessing of ChIP-seq")