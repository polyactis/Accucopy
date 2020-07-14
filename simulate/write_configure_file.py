#!/usr/bin/python
from subprocess import *
import sys
import os


def run(cmd, wkdir=None):
    sys.stderr.write("Running %s ...\n" % cmd)
    p = Popen(cmd, shell=True, cwd=wkdir)
    p.wait()
    return p.returncode

def touchFile(outputDir = "", AccurityPath = None):
    conf = outputDir + "configure"
    of = open(conf, "w")
    stri = """reference_genome_name	hs37d5
read_length	101
window_size	500
reference_index_folder_path	/y/hlab/AccurityTestData/refData/
reference_genome_fasta_path	/y/hlab/AccurityTestData/refData/hs37d5.fa
samtools_path	/y/program/bin/samtools
freebayes_path	/y/program/bin/freebayes
accurity_path	%s""" %(AccurityPath)
    of.write(stri)
    of.close()

touchFile(AccurityPath=sys.argv[1])