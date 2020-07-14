#!/usr/bin/python
from pegaflow.DAX3 import *
import os
import sys


workflow_AC = ADAG("SilicoValidation")

def registerFile(workflow, filename):
    file = File(os.path.basename(filename))
    file.addPFN(PFN("file://" + os.path.abspath(filename), "local"))
    workflow.addFile(file)
    return file

def registerExecutefile(workflow, executeFile):
    architecture = "x86"
    operatingSystem = "linux"
    namespace = "pegasus"
    version = "1.0"
    executeName = os.path.basename(executeFile)
    execute = Executable(namespace=namespace, name=executeName, os=operatingSystem, arch=architecture, installed=True)
    execute.addPFN(PFN("file://" + os.path.abspath(executeFile), "local"))
    workflow.addExecutable(execute)
    return executeName

def addEagleJobToWorkflow(workflow, excute, argv):
    step_1 = Job(namespace="pegasus", name=excute)
    step_1.addArguments(argv[0], argv[1], argv[2])
    workflow.addJob(step_1)

def addPicardtoWorkflow(workflow , excute, picard,  eagledir):
#    eagle_bam = registerFile(workflow_AC, eagledir + "/eagle.bam")
    eagle_bam = eagledir + "/eagle.bam"
    fastq_1 = File(eagledir + "/eagle_1.fastq")
    fastq_2 = File(eagledir + "/eagle_2.fastq")
    step_1 = Job(namespace="pegasus", name=excute)
    step_1.addArguments("-jar ", picard, "SamToFastq", "I=", eagle_bam, "FASTQ=", fastq_1, "SECOND_END_FASTQ=", fastq_2, "VALIDATION_STRINGENCY=LENIENT")
    step_1.uses(picard , link=Link.INPUT)
#    step_1.uses(eagle_bam, link=Link.INPUT)
    step_1.uses(fastq_1, link=Link.OUTPUT)
    step_1.uses(fastq_2, link=Link.OUTPUT)
    workflow.addJob(step_1)
    return fastq_1, fastq_2

def addPrefxToWorkflow(workflow, excute, eagledir, fastq, part , num):
    step_1 = Job(namespace="pegasus", name=excute)
    fastq_add = File(eagledir + "/eagle_"+ part +"_add.fastq")
    step_1.addArguments(fastq, fastq_add, num)
    step_1.uses(fastq, link=Link.INPUT)
    step_1.uses(fastq_add, link=Link.OUTPUT)
    workflow.addJob(step_1)
    return fastq_add

def addMixToWorkflow(workflow, excute, mixdir, fastq, part):
    step_1 = Job(namespace="pegasus", name=excute)
    fastq_mix = File(mixdir + "/mix_" + part + ".fastq")
    for i in fastq:
        step_1.uses(i, link=Link.INPUT)
    fastq.append(fastq_mix)
    step_1.addArguments(*fastq)
    step_1.uses(fastq_mix, link=Link.OUTPUT)
    workflow.addJob(step_1)
    return fastq_mix

def addBwaToWorkflow(workflow, excute, argv):
    step_1 = Job(namespace="pegasus", name=excute)
    sam = File(argv[4])
    step_1.addArguments(*argv)
    step_1.uses(sam, link=Link.OUTPUT)
    step_1.uses(argv[2], link=Link.INPUT)
    step_1.uses(argv[3], link=Link.INPUT)
    workflow.addJob(step_1)
    return sam


def addSamtoolsViewToWorkflow(workflow, excute, argv):
    step_1 = Job(namespace="pegasus", name=excute)
    bam = File(argv[3])
    step_1.addArguments(*argv)
    step_1.uses(bam, link=Link.OUTPUT)
    step_1.uses(argv[1], link=Link.INPUT)
    workflow.addJob(step_1)
    return bam

def addSamtoolsSortToWorkflow(workflow, excute, argv):
    step_1 = Job(namespace="pegasus", name=excute)
    sortBam = File(argv[1])
    step_1.addArguments(*argv)
    step_1.uses(sortBam, link=Link.OUTPUT)
    step_1.uses(argv[2], link=Link.INPUT)
    workflow.addJob(step_1)
    return sortBam

def addPurityToWorkflow(workflow, excute, dir):
    step_1 = Job(namespace="pegasus", name=excute)
    step_1.addArguments(dir)
    workflow.addJob(step_1)


##########################execute register################################################################
eagle_exe = registerExecutefile(workflow_AC, "/y/hlab/tcga/gdc-tools/Puricise/simulation_pipeline/sbin/eagle.py")
java_exe = registerExecutefile(workflow_AC, "/y/home/luozhihui/program/java/jdk1.8.0_101/bin/java")
picard = registerFile(workflow_AC, "/y/home/luozhihui/my_source_code/picard.jar")
add_prefix = registerExecutefile(workflow_AC, "/y/hlab/tcga/gdc-tools/Puricise/simulation_pipeline/sbin/addPrefix.pl")
mix_exe = registerExecutefile(workflow_AC, "/y/hlab/tcga/gdc-tools/Puricise/simulation_pipeline/sbin/mix.pl")
bwa_exe = registerExecutefile(workflow_AC, "/y/hlab/tcga/gdc-tools/Puricise/simulation_pipeline/sbin/wrap.sh")
samtools_exe = registerExecutefile(workflow_AC, "/y/home/luozhihui/program/bin/samtools")
purity_exe = registerExecutefile(workflow_AC, "/y/hlab/tcga/gdc-tools/Puricise/simulation_pipeline/sbin/purity.py")


###########################the workflow#####################################################################
curdir = os.getcwd()
vcf = "/y/hlab/tcga/gdc-tools/Puricise/simulation_pipeline/vcf/None.vcf"
coverage, step = 5, 0.2
project = "covergae-"+ str(coverage) +"-purity-step-"+str(step)
simulationDir = curdir + "/" + project
if not os.path.exists(simulationDir):
    os.mkdir(simulationDir)


normalSample = simulationDir + "/normalSample"
addEagleJobToWorkflow(workflow_AC, eagle_exe, [str(coverage), vcf, normalSample])
fq = addPicardtoWorkflow(workflow_AC, java_exe, picard, normalSample)
normalSamOutput = addBwaToWorkflow(workflow_AC, bwa_exe, ["10", "/y/hlab/tcga/gdc-tools/Puricise/1000_genomes/hs37d5.fa", fq[0], fq[1], normalSample + "/normal.aligned.sam"])
normalBamOutput = addSamtoolsViewToWorkflow(workflow_AC, samtools_exe, ["view -b -S", normalSamOutput, "-o", normalSample + "/normal.aligned.bam"])
normalSortBamOutput = addSamtoolsSortToWorkflow(workflow_AC, samtools_exe, ["sort -@ 10 -o", normalSample + "/normal.aligned.sort.bam", normalBamOutput])

#for i in range(0, int(1/step)):
i =0
purity = (i+1) * step
dirname = simulationDir + "/" + "purity." + str(purity)
if not os.path.exists(dirname):
    os.mkdir(dirname)

normalCoverage = coverage * (1-purity)
tumorCoverage = coverage * purity
eagleList = (dirname + "/eagle_tumor_sample_normal_part_" + str(normalCoverage), dirname + "eagle_tumor_sample_tumor_part_" + str(tumorCoverage))
addPrefixOrder = 0
eagleDict, mixDict = {}, {}
fq_1, fq_2 = [], []
#    for subdir in eagleList:
subdir = eagleList[0]
addEagleJobToWorkflow(workflow_AC, eagle_exe, [str(normalCoverage), vcf, subdir])
eagleDict[subdir] = addPicardtoWorkflow(workflow_AC, java_exe, picard, subdir)
addPrefixOrder += 1
eagleDict[subdir][2] = addPrefxToWorkflow(workflow_AC, add_prefix, subdir, eagleDict[subdir][0], "1", str(addPrefixOrder))
eagleDict[subdir][3] = addPrefxToWorkflow(workflow_AC, add_prefix, subdir, eagleDict[subdir][1], "2", str(addPrefixOrder))
fq_1.append(eagleDict[subdir][2])
fq_2.append(eagleDict[subdir][3])

subdir = eagleList[1]
addEagleJobToWorkflow(workflow_AC, eagle_exe, [str(normalCoverage), vcf, subdir])
eagleDict[subdir] = addPicardtoWorkflow(workflow_AC, java_exe, picard, subdir)
addPrefixOrder += 1
eagleDict[subdir][2] = addPrefxToWorkflow(workflow_AC, add_prefix, subdir, eagleDict[subdir][0], "1", str(addPrefixOrder))
eagleDict[subdir][3] = addPrefxToWorkflow(workflow_AC, add_prefix, subdir, eagleDict[subdir][1], "2", str(addPrefixOrder))
fq_1.append(eagleDict[subdir][2])
fq_2.append(eagleDict[subdir][3])





    mixDict[eagleList[0]] = addMixToWorkflow(workflow_AC, mix_exe, dirname, fq_1, "1")
    mixDict[eagleList[1]] = addMixToWorkflow(workflow_AC, mix_exe, dirname, fq_2, "2")
    tumorSamOutput = addBwaToWorkflow(workflow_AC, bwa_exe, ["10", "/y/hlab/tcga/gdc-tools/Puricise/1000_genomes/hs37d5.fa", mixDict[eagleList[0]], mixDict[eagleList[1]], dirname + "/tumor.aligned.sam"])
    tumorBamOutput = addSamtoolsViewToWorkflow(workflow_AC, samtools_exe, ["view -b -S", tumorSamOutput, "-o", dirname + "/tumor.aligned.bam"])
    tumorSortBamOutput = addSamtoolsSortToWorkflow(workflow_AC, samtools_exe, ["sort -@ 10 -o", dirname + "/tumor.aligned.sort.bam", tumorBamOutput])
    addPurityToWorkflow(workflow_AC, purity_exe, dirname)

workflow_AC.writeXML(sys.stdout)
