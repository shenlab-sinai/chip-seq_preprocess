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
alignment_path = expandOsPath(os.path.join(config["project_dir"], config["data_dir"],"Alignment"))
fastqc_path = expandOsPath(os.path.join(config["project_dir"], config["data_dir"], "FastQC"))
rmdup_path = expandOsPath(os.path.join(config["project_dir"], config["data_dir"], "Alignment", "rmdup"))
tdf_path = expandOsPath(os.path.join(config["project_dir"], config["data_dir"], "Alignment", "rmdup", "tdf"))
phantompeak_path = expandOsPath(os.path.join(config["project_dir"], config["data_dir"], "PhantomPeak"))
diffrepeat_path = expandOsPath(os.path.join(config["project_dir"], config["data_dir"], "DiffRepeat"))
log_path = expandOsPath(os.path.join(config["project_dir"], config["data_dir"], "Log"))
ngsplot_path = expandOsPath(os.path.join(config["project_dir"], config["data_dir"], "NgsPlot"))

@merge(FqFiles, expandOsPath(os.path.join(log_path,"mkdir.log")), config)
def makeDir(input_files,output_file,config):
    """
    make directories
    """
    cmds = ['mkdir']
    cmds.append(log_path)
    cmds.append(fastqc_path)
    cmds.append(alignment_path)
    cmds.append(rmdup_path)
    cmds.append(tdf_path)
    cmds.append(phantompeak_path)
    cmds.append(diffrepeat_path)
    cmds.append(ngsplot_path)
    cores = 1
    logfile = expandOsPath(os.path.join(config["project_dir"], config["data_dir"])) + "/" + "mkdir.log"

    p = subprocess.Popen(
        cmds, stdout=open(logfile, "w"), stderr=open(logfile, "w"),
        bufsize=1)
    stdout, stderr = p.communicate()
    return stdout

@follows(makeDir)
@transform(FqFiles, formatter(fq_ext), os.path.join(alignment_path, "{basename[0]}.uniqmapped.bam"), config)
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
            cmds.append(config['parseAln_mapq'])
            cmds.append(config['parseAln_diff'])
        else:
            raise KeyError
    else:
        cmds = ['fastq2bam_by_bowtie.sh']
        cmds.append(FqFileName)
        cmds.append(expandOsPath(config['bowtie_index']))

    cmds.append(alignment_path)
    cmds.append(config["pair_end"])
    cores = int(config['cores'])
    if cores == 0:
        cores = 1
    cmds.append(str(cores))
    logfile = log_path + "/" +  os.path.basename(FqFileName) + ".alignment.log"

    p = subprocess.Popen(
        cmds, stdout=open(logfile, "w"), stderr=open(logfile, "w"),
        bufsize=1)
    stdout, stderr = p.communicate()
    return stdout

@follows(makeDir)
@transform(FqFiles, formatter(fq_ext), os.path.join(log_path, "{basename[0]}.fastq.fastqc.log"), config)
def runFastqc(FqFileName, fastqcLog, config):
    """
    To run FastQC
    Arguments:
    - `FqFileName`: fastq file
    - `config`: config
    """
    cmds = ['runFastQC.sh']
    #cmds.append("-o")
    cmds.append(fastqc_path)
    cores = int(config['cores'])
    if cores == 0:
        cores = 1
    #cmds.append("-t")
    cmds.append(str(cores))
    cmds.append(FqFileName)

    cmds.append(config["pair_end"])

    logfile = fastqcLog

    p = subprocess.Popen(
        cmds, stdout=open(logfile, "w"), stderr=open(logfile, "w"),
        bufsize=1)
    stdout, stderr = p.communicate()
    return stdout

@transform(alignFastqByBowtie, formatter(".bam"), os.path.join(rmdup_path, "{basename[0]}_rmdup.bam"), config)
def rmdupBam(BamFileName, rmdupFile, config):
    """
    To remove duplicates
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
    #if "bam_sort_buff" in config:
    #    cmds.append(config["bam_sort_buff"])
    logfile = log_path + "/" +  os.path.basename(BamFileName) + ".rmdup.log"

    cores = 1

    p = subprocess.Popen(
        cmds, stdout=open(logfile, "w"), stderr=open(logfile, "w"),
        bufsize=1)
    stdout, stderr = p.communicate()
    return stdout

@transform(rmdupBam, formatter(".bam"), os.path.join(log_path, "{basename[0]}.bam.tdf.log"), config)
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
    baseName = os.path.basename(BamFileName)
    cmds.append(os.path.join(tdf_path, baseName.replace(".bam", ".tdf")))
    cmds.append(config["IGV_genome"])
    logfile = tdfLog

    cores = 1

    p = subprocess.Popen(
        cmds, stdout=open(logfile, "w"), stderr=open(logfile, "w"),
        bufsize=1)
    stdout, stderr = p.communicate()
    return stdout

@transform(rmdupBam, formatter(".bam"), os.path.join(log_path, "{basename[0]}.bam.phantomPeak.log"), config)
def runPhantomPeak(BamFileName, PPLog, config):
    """
    To check data with phantomPeak
    Arguments:
    - `BamFileName`: bam file
    - `config`: config
    """
    cmds = ['runPhantomPeak.sh']
    cmds.append(BamFileName)
    cmds.append(str(config["cores"]))
    cmds.append(phantompeak_path)
    logfile = PPLog

    cores = int(config['cores'])
    if cores == 0:
        cores = 1

    p = subprocess.Popen(
        cmds, stdout=open(logfile, "w"), stderr=open(logfile, "w"),
        bufsize=1)
    stdout, stderr = p.communicate()
    return stdout

@merge(alignFastqByBowtie, expandOsPath(os.path.join(diffrepeat_path, config["project_name"]+".diffrepeat.resulttable")), config)
def runDiffrepeat(BamFileNames, ResultFile, config):
    """
    To run diffrepeats
    Arguments:
    - `BamFileNames`: bam files
    - `config`: config
    """
    cmds = ['runDiffrepeat.sh']
    cmds.append(diffrepeat_path)
    cmds.append(alignment_path)
    cmds.append(config["repbase_db"])
    cmds.append(ResultFile)
    cmds.append(config["diffrepeat_editdist"])
    cmds.append(config["diffrepeat_mapq"])
    logfile = expandOsPath(os.path.join(log_path, config["project_name"]+".diffrepeat.log"))

    cores = int(config['cores'])
    if cores == 0:
        cores = 1

    p = subprocess.Popen(
        cmds, stdout=open(logfile, "w"), stderr=open(logfile, "w"),
        bufsize=1)
    stdout, stderr = p.communicate()
    return stdout

@merge(rmdupBam, expandOsPath(os.path.join(ngsplot_path, config["project_name"]+".ngsplot.all.log")), config)
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
    logfile1 = expandOsPath(os.path.join(log_path, config["project_name"]+".ngsplot.all.log1"))
    logfile2 = expandOsPath(os.path.join(log_path, config["project_name"]+".ngsplot.all.log2"))
    p = subprocess.Popen(
        cmds, stdout=open(logfile1, "w"), stderr=open(logfile2, "w"),
        bufsize=1)
    stdout, stderr = p.communicate()
    return stdout

if config["run_phantompeak"]=="yes":
    @follows(runFastqc, genTDF, runPhantomPeak, runNgsplotAll, runDiffrepeat)
    def complete():
        """
        dummy step to consolidate pipeline
        """
        return 0
else:
    @follows(runFastqc, genTDF, runNgsplotAll, runDiffrepeat)
    def complete():
        """
        dummy step to consolidate pipeline
        """
        return 0

## Run pipeline
pipeline_printout_graph("all_flowchart.png", "png", [complete], pipeline_name="Preprocessing of ChIP-seq")
pipeline_run([complete], multiprocess=config["cores"])

