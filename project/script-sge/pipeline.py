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

def writeLog(logfile, stdout_res, stderr_res):
    """
    To write the log information from drmaa to log file.
    Arguments:
    - `logfile`: file name of target logfile.
    - `stdout_res`: list, return from drmaa, stdout of the job.
    - `stderr_res`: list, return from drmaa, stderr of the job.
    """
    log_f = file(logfile, "w")
    log_f.write("Stdout:\n")
    log_f.write("".join(stdout_res))
    log_f.write("\nStderr:\n")
    log_f.write("".join(stderr_res))
    log_f.close()
    return 0

def cluster_options(config, task_name, cores, logfile):
    """
    Generate a string of cluster options to feed an LSF job.
    Arguments:
    - `config`: configuration as associative array from the YAML file.
    - `task_name`: the specific task name, such as runPhantomPeak.
    - `cores`: number of cores to use for this task.
    - `logfile`: log file name.
    """
    ## Here are the paramters for SGE.
    str_options = "-cwd -V -pe shm %d -q %s -j y -o %s" % \
        (cores, config["queue"], logfile)
    if "h_vmem" in config:
        str_options = str_options + " -l h_vmem=%s" % (config["h_vmem"])
    return str_options


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

scipt_path = os.path.dirname(os.path.realpath(__file__))

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
    cores = int(config['cores'])
    if cores == 0:
        cores = 1
    cmds.append(str(cores))
    logfile = FqFileName + ".alignment.log"

    stdout_res, stderr_res = run_job(" ".join(cmds),
        job_name = "alignFastqByBowtie_" + os.path.basename(FqFileName),
        job_other_options = cluster_options(config, "alignFastqByBowtie", cores, logfile),
        job_script_directory = os.path.dirname(os.path.realpath(__file__)),
        job_environment={ 'BASH_ENV' : '~/.bashrc' },
        retain_job_scripts = True, drmaa_session=my_drmaa_session)

    writeLog(logfile, stdout_res, stderr_res)

    return 0

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
    cores = int(config['cores'])
    if cores == 0:
        cores = 1
    cmds.append("-t")
    cmds.append(str(cores))
    cmds.append(BamFileName)
    logfile = BamFileName + ".fastqc.log"

    stdout_res, stderr_res = run_job(" ".join(cmds),
        job_name = "fastqc_" + os.path.basename(BamFileName),
        job_other_options = cluster_options(config, "runFastqc", cores, logfile),
        job_script_directory = os.path.dirname(os.path.realpath(__file__)),
        job_environment={ 'BASH_ENV' : '~/.bashrc' },
        retain_job_scripts = True, drmaa_session=my_drmaa_session)

    writeLog(logfile, stdout_res, stderr_res)

    return 0

@follows(runFastqc, mkdir(rmdup_path))
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
    logfile = BamFileName + ".rmdup.log"

    cores = 1

    stdout_res, stderr_res = run_job(" ".join(cmds),
        job_name = "rmdup_" + os.path.basename(BamFileName),
        job_other_options = cluster_options(config, "rmdupBam", cores, logfile),
        job_script_directory = os.path.dirname(os.path.realpath(__file__)),
        job_environment={ 'BASH_ENV' : '~/.bashrc' },
        retain_job_scripts = True, drmaa_session=my_drmaa_session)

    writeLog(logfile, stdout_res, stderr_res)

    return 0

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

    cores = 1

    stdout_res, stderr_res = run_job(" ".join(cmds),
        job_name = "genTDF_" + os.path.basename(BamFileName),
        job_other_options = cluster_options(config, "genTDF", cores, logfile),
        job_script_directory = os.path.dirname(os.path.realpath(__file__)),
        job_environment={ 'BASH_ENV' : '~/.bashrc' },
        retain_job_scripts = True, drmaa_session=my_drmaa_session)

    writeLog(logfile, stdout_res, stderr_res)

    return 0

@follows(genTDF)
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
    cmds.append(str(config["cores"]))
    logfile = BamFileName + ".phantomPeak.log"

    cores = int(config['cores'])
    if cores == 0:
        cores = 1

    stdout_res, stderr_res = run_job(" ".join(cmds),
        job_name = "runPhantomPeak_" + os.path.basename(BamFileName),
        job_other_options = cluster_options(config, "runPhantomPeak", cores, logfile),
        job_script_directory = os.path.dirname(os.path.realpath(__file__)),
        job_environment={ 'BASH_ENV' : '~/.bashrc' },
        retain_job_scripts = True, drmaa_session=my_drmaa_session)

    writeLog(logfile, stdout_res, stderr_res)

    return 0

if __name__ == '__main__':
    ## run to step of PhantomPeak
    ## multithread number need to be changed!
    pipeline_run([runPhantomPeak], multithread=200)

    my_drmaa_session.exit()
