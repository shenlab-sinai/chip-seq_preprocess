#! /usr/bin/env python

import os
import sys
import yaml
from ruffus import *
import glob
import subprocess
import string
import drmaa
from ruffus.drmaa_wrapper import run_job, error_drmaa_job

my_drmaa_session = drmaa.Session()
my_drmaa_session.initialize()

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

def cluster_options(config, task_name, cores, logfile):
    """
    Generate a string of cluster options to feed an LSF job.
    Arguments:
    - `config`: configuration as associative array from the YAML file.
    - `task_name`: the specific task name, such as runPhantomPeak.
    - `cores`: number of cores to use for this task.
    - `logfile`: log file name.
    """

    str_options = "-W %s -n %d -o %s -e %s -q %s -R span[hosts=1] -L /bin/bash" % \
        (config["wall_time"][task_name], cores, logfile, logfile, config["queue"]) 
    # Name of the partition where you want to run your job. By default, the cluster will assign one for you.
    if "m" in config:
        str_options = str_options + " -m %s" % (config["m"])
    # Allocation account. By default, will use the free account (i.e. scavenger on Minerva).
    if "alloc" in config:
        str_options = str_options + " -P %s" % (config["alloc"])
    return str_options


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

script_path = os.path.dirname(os.path.realpath(__file__))

@follows(mkdir(alignment_path), mkdir(log_path))
@transform(FqFiles, formatter(fq_ext), os.path.join(alignment_path, "{basename[0]}.bam"), config)
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

    cmds.append(alignment_path)
    cmds.append(config["pair_end"])
    cores = int(config['cores'])
    if cores == 0:
        cores = 1
    cmds.append(str(cores))
    logfile = log_path + "/" +  os.path.basename(FqFileName) + ".alignment.log"

    run_job(" ".join(cmds),
        job_name = os.path.basename(FqFileName) + "_alignFastqByBowtie",
        job_other_options = cluster_options(config, "alignFastqByBowtie", cores, logfile),
        job_script_directory = os.path.dirname(os.path.realpath(__file__)),
        job_environment={ 'BASH_ENV' : '~/.bash_profile' },
        retain_job_scripts = True, drmaa_session=my_drmaa_session)

    return 0

@follows(mkdir(fastqc_path))
@transform(FqFiles, formatter(fq_ext), os.path.join(log_path, "{basename[0]}.bam.fastqc.log"), config)
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
    
    run_job(" ".join(cmds),
        job_name = os.path.basename(FqFileName) + "_fastqc",
        job_other_options = cluster_options(config, "runFastqc", cores, logfile),
        job_script_directory = os.path.dirname(os.path.realpath(__file__)),
        job_environment={ 'BASH_ENV' : '~/.bash_profile' },
        retain_job_scripts = True, drmaa_session=my_drmaa_session)

    return 0

@follows(alignFastqByBowtie, mkdir(rmdup_path))
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

    run_job(" ".join(cmds),
        job_name = os.path.basename(BamFileName) + "_rmdup",
        job_other_options = cluster_options(config, "rmdupBam", cores, logfile),
        job_script_directory = os.path.dirname(os.path.realpath(__file__)),
        job_environment={ 'BASH_ENV' : '~/.bash_profile' },
        retain_job_scripts = True, drmaa_session=my_drmaa_session)

    return 0

@follows(rmdupBam, mkdir(tdf_path))
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

    run_job(" ".join(cmds),
        job_name = os.path.basename(BamFileName) + "_genTDF",
        job_other_options = cluster_options(config, "genTDF", cores, logfile),
        job_script_directory = os.path.dirname(os.path.realpath(__file__)),
        job_environment={ 'BASH_ENV' : '~/.bash_profile' },
        retain_job_scripts = True, drmaa_session=my_drmaa_session)

    return 0

@follows(rmdupBam, mkdir(phantompeak_path))
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

    run_job(" ".join(cmds),
        job_name = os.path.basename(BamFileName) + "_runPhantomPeak",
        job_other_options = cluster_options(config, "runPhantomPeak", cores, logfile),
        job_script_directory = os.path.dirname(os.path.realpath(__file__)),
        job_environment={ 'BASH_ENV' : '~/.bash_profile' },
        retain_job_scripts = True, drmaa_session=my_drmaa_session)

    return 0

@follows(alignFastqByBowtie, mkdir(diffrepeat_path))
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
    logfile = expandOsPath(os.path.join(log_path, config["project_name"]+".diffrepeat.log"))

    cores = int(config['cores'])
    if cores == 0:
        cores = 1

    run_job(" ".join(cmds),
        job_name = "runDiffRepeat",
        job_other_options = cluster_options(config, "runDiffrepeat", cores, logfile),
        job_script_directory = os.path.dirname(os.path.realpath(__file__)),
        job_environment={ 'BASH_ENV' : '~/.bash_profile' },
        retain_job_scripts = True, drmaa_session=my_drmaa_session)

    return 0

if config["run_phantompeak"]=="yes":
    @follows(runFastqc, genTDF, runPhantomPeak, runDiffrepeat)
    def complete():
        """
        dummy step to consolidate pipeline
        """
        return 0
else:
    @follows(runFastqc, genTDF, runDiffrepeat)
    def complete():
        """
        dummy step to consolidate pipeline
        """
        return 0

if __name__ == '__main__':
    pipeline_printout_graph("all_flowchart.png", "png", [complete], pipeline_name="Preprocessing of ChIP-seq on LSF")

    ## run to step of PhantomPeak
    ## multithread number need to be changed!
    pipeline_run([complete], multithread=6)
    my_drmaa_session.exit()
