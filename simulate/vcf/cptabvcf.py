#!/usr/bin/python
from subprocess import *
import sys
import os

def run(cmd, wkdir=None):
    sys.stderr.write("Running %s ...\n" % cmd)
    p = Popen(cmd, shell=True, cwd=wkdir)
    p.wait()
    return p.returncode

num = len(sys.argv)
print num
assert not num < 3, "the argument number is wrong!!! "
Dir = sys.argv[1]
for i in sys.argv[2:]:
    cmd = "/bin/cp " + i + " " + Dir + "/"
    b = run(cmd)
    assert not b, "error in copy the file: " + i
