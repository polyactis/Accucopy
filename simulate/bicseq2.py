import sys,os
from subprocess import *

def run(cmd, wkdir=None):
    sys.stderr.write("Running %s ...\n" % cmd)
    if wkdir is None:
        p = Popen(cmd, shell=True)
    else:
        p = Popen(cmd, shell=True, cwd=wkdir)
    p.wait()
    return p.returncode

class BicSeq:
    def __init__(self, tumorBam=None, normalBam=None):
        self.normpath = "/y/home/luozhihui/my_source_code/NBICseq-norm_v0.2.4/NBICseq-norm.pl"
        self.segPath = '/y/home/luozhihui/my_source_code/NBICseq-seg_v0.7.2/NBICseq-seg.pl'
        """
        chromName  	faFile 	MapFile 	readPosFile 	binFileNorm
        chr5	/y/home/luozhihui/my_source_code/NBICseq-norm_v0.2.4/test/ref/chr5.fa	map/hg19CRG.100bp/hg19.CRC.100mer.chr5.txt	myfolder/chr5.seq	chr5.norm.bin
        """
        self.normconfig = "/y/home/luozhihui/src/cancer_purity/workflow/sbin/configFile.norm"
        """
        chromNamei      binFileNorm
        chr5    /y/home/luozhihui/my_source_code/NBICseq-norm_v0.2.4/test/chr5.norm.bin
        """
        self.segconfig = "/y/home/luozhihui/src/cancer_purity/workflow/sbin/configFile.seg"
        self.newSamtools = "/y/home/luozhihui/my_source_code/NBICseq-norm_v0.2.4/test/samtools-0.1.7a_getUnique-0.1.3/samtools"
        self.tumorBam = tumorBam
        self.normalBam = normalBam
        self.mapPath = "/y/home/luozhihui/my_source_code/NBICseq-norm_v0.2.4/test/map/hg19CRG.100bp"


#step two: setup configure file for norm
    def seqNorm(self, workdir=None, bamfile=None):
        if not os.path.exists(workdir):
            os.makedirs(workdir)
        if not os.path.exists(workdir + "/SeqDir"):
            os.makedirs(workdir + "/SeqDir")
        cmd = self.newSamtools + " view " + " -U BWA," + workdir +"/SeqDir/,N,N  -q 20 " + bamfile
        run(cmd)
        confFile = open(workdir + "/configFile", "w")
        confFile.write("chromName\tfaFile\tMapFile\treadPosFile\tbinFileNorm\n")
        for i in range(1, 23):
            chr = "chr" + str(i)
            ref = "/y/home/luozhihui/my_source_code/NBICseq-norm_v0.2.4/test/ref/" + chr + ".fa"
            map = self.mapPath + "/hg19.CRC.100mer." + chr + ".txt"
            posi = workdir + "/SeqDir/" + chr + ".seq"
            binFile = workdir + "/" + chr +".norm.bin"
            confFile.write(chr + "\t" + ref + "\t" + map + "\t" + posi + "\t" + binFile + "\n")
        confFile.close()
        cmd = self.normpath + " -l=100 -s=300 -b=500 --fig=" + workdir + "/normalize.pdf " + workdir + "/configFile " + workdir +"/output.file"
        run(cmd)

    def seqSeg(self, workdir=None, tumorDir=None, normDir=None):
        if not os.path.exists(workdir):
            os.makedirs(workdir)
        confFile = open(workdir + "/configFile", "w")
        confFile.write("chromName\tbinFileNorm.Case\tbinFileNorm.Control\n")
        for i in range(1,23):
            chr = "chr" + str(i)
            caseBinFile = tumorDir + "/" + chr + ".norm.bin"
            controlBinFile = normDir + "/" + chr + ".norm.bin"
            confFile.write(chr + "\t" + caseBinFile + "\t" + controlBinFile + "\n")
        confFile.close()
        cmd = self.segPath + " --control  --nrm --strict --lambda=100  --fig=" + workdir + "/cnv.pdf " + workdir + "/configFile " + workdir + "/output.file"
        run(cmd)

def runbic():
    bic = BicSeq(tumorBam="/y/Sunset/db/individual_alignment/398_572_TCGA-HI-7169_Broad_GSC_vs_1_by_method5_realigned0_reduced0.bam",\
                 normalBam= "/y/Sunset/db/individual_alignment/383_557_TCGA-HI-7169_Broad_GSC_vs_1_by_method5_realigned0_reduced0.bam")
    tumorDir = "./tumorBicSeq"
    bic.seqNorm(workdir=tumorDir, bamfile=bic.tumorBam)
    normalDir = "./normalBicSeq"
    bic.seqNorm(workdir=normalDir, bamfile=bic.normalBam)
    bic.seqSeg(workdir="seg", tumorDir=tumorDir, normDir=normalDir)

runbic()
