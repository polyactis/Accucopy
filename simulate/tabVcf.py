#!/usr/bin/python
from subprocess import *
import sys
import os

def run(cmd, wkdir=None):
    sys.stderr.write("Running %s ...\n" % cmd)
    p = Popen(cmd, shell=True, cwd=wkdir)
    p.wait()
    return p.returncode

upDir = os.path.dirname(sys.argv[2])
if not os.path.exists(upDir):
    cmd = "mkdirhier " + upDir
    a = run(cmd)
    assert not a, "error in make dir " + upDir
cmd = "/bin/cp " + sys.argv[1] + " " + upDir  + "/"
run(cmd)
cmd = "/y/home/luozhihui/program/libexec/EAGLE/applyCopyNumber.pl -i " + sys.argv[1] + " -o " + sys.argv[2]
run(cmd)
