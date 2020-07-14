#!/usr/bin/python
from subprocess import *
import sys
import os

dire = sys.argv[1]
normalBam = os.path.dirname(dire) + "/normalSample/normal.aligned.sort.bam"

def run(cmd, wkdir=None):
    sys.stderr.write("Running %s ...\n" % cmd)
    p = Popen(cmd, shell=True, cwd=wkdir)
    p.wait()
    return p.returncode



cmd = "/bin/cp /y/home/luozhihui/src/normalize_cpp/workflow/sbin/configure ./"
run(cmd, dire)
cmd = "/bin/ln -s %s %s"%(os.path.abspath(normalBam), os.path.abspath(dire) + "/normal.bam")
run(cmd)
cmd = "/y/home/luozhihui/src/normalize_cpp/Accurity/main.py -c configure -t tumor.aligned.sort.bam  -n normal.bam -o infer_result --debug 1"

a = run(cmd, dire)
assert not a, "error in purity estimate"