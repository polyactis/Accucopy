#!/usr/bin/python
from subprocess import *
import sys
import os

def run(cmd, wkdir=None):
    sys.stderr.write("Running %s ...\n" % cmd)
    p = Popen(cmd, shell=True, cwd=wkdir)
    p.wait()
    return p.returncode

dire = sys.argv[-1]
atr = " ".join(sys.argv[1:-1])
cmd = "/y/Sunset/tcga/gdc-tools/Puricise/simulation_pipeline/sbin/calculate.pl " + atr
run(cmd, dire)
