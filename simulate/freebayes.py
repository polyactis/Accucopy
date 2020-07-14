#!/usr/bin/python
from subprocess import *
import sys
import os


def run(cmd, wkdir=None):
    sys.stderr.write("Running %s ...\n" % cmd)
    p = Popen(cmd, shell=True, cwd=wkdir)
    p.wait()
    return p.returncode

def runFreeBayes(freebayes=None, ref_genome_fasta_file=None, oneThousandSNPFilepath=None, tumor_bam=None, vcf_file_path=None):
    cmd = "%s -f %s --no-indels --no-mnps --no-complex -m 20 -q 20 -t %s %s|gzip > %s" % (
        freebayes, ref_genome_fasta_file, oneThousandSNPFilepath, tumor_bam, vcf_file_path)
    run(cmd)

runFreeBayes(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
