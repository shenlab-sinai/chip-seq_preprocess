#! /usr/bin/env python
import difflib
import os
import string
import glob
import sys

def expandOsPath(path):
    """
    To expand the path with shell variables.
    Arguments:
    - `path`: path string list
    """
    return os.path.expanduser(os.path.expandvars(os.path.join(*path)))

def main():
    """
    To generate config.norm.txt for ngs.plot to draw plot with input as background.
    sys.argv[1]: folder holds bam files
    """
    bams = [ os.path.basename(x) for x in glob.glob(expandOsPath([sys.argv[1], "*.bam"]))]
    bams.sort()
    inputs = [ x for x in bams if ("input" in x) or ("Input" in x)]
    chips = [ x for x in bams if ("input" not in x) and ("Input" not in x)]
    for chip in chips:
        if len(inputs) > 1:
            bk = difflib.get_close_matches(chip, inputs, 1, 0.01)[0]
        else:
            bk = inputs[0]
        print("%s:%s\t-1\t%s:%s"%(chip, bk, string.replace(chip, ".bam", ""), string.replace(bk, ".bam", "")))


if __name__ == '__main__':
    main()