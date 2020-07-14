#!/usr/bin/python
from subprocess import *
import sys
import os


def run(cmd, wkdir=None):
    sys.stderr.write("Running %s ...\n" % cmd)
    p = Popen(cmd, shell=True, cwd=wkdir)
    p.wait()
    return p.returncode

def eagle(coverage, vcf_1, outdir , mk_j = 8, mkbam_j=8):
    coverage_n = str(coverage)
    updir = os.path.dirname(outdir)
    if not os.path.exists(updir):
        cmd = "mkdirhier " + updir
	a = run(cmd)
	assert not a, "error in make dir " + updir
    if os.path.exists(outdir):
        cmd = "rm -r " + outdir
        a = run(cmd)
        assert not a, "error in delete dir " + outdir
    cmd = "/y/home/luozhihui/my_source_code/EAGLE/packages/bin/configureEAGLE.pl" + "\
    --run-info=/y/home/luozhihui/EAGLE/simulation/RunInfo_PairedReadsBarcode8x32Tiles.xml \
    --reference-genome=/y/Sunset/tcga/gdc-tools/Puricise/bwa_aln/hs37d5_namechr.fa  \
    --variant-list=/y/Sunset/tcga/gdc-tools/Puricise/1000_genomes/biallel_snp_vcf/total_vcf.vcf \
    --variant-list=" + vcf_1  + "\
    --coverage-depth=" + coverage_n + "\
    --genome-mutator-options='--organism-ploidy=2 --ploidy-chromosome=chrY --ploidy-level=0' \
    " + outdir
    dirname = os.path.dirname(outdir)
    a = run(cmd)
#    assert not a ,"error in configureEAGLE.pl"
    cmd = "make -j " + str(mk_j)
    a = run(cmd, outdir)
#    assert not a, "error in make -j" + str(mk_j)
    cmd = "make bam -j " + str(mkbam_j)
    a = run(cmd, outdir)
#    assert not a, "error in make bam -j" + str(mkbam_j)


eagle(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
