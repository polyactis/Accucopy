import os,sys

def splitFile(inputF=None):
    data = open(inputF, "r")
    header = data.readline()
    for line in data:
        array = line.split("\t")
        chrom = array[0]
        start = array[1]
        end = array[2]
        segfile = chrom + ".newseg"
        if not os.path.exists(segfile):
            if 'OP' in locals().keys():
                OP.close()
            OP = open(segfile, "w")
        OP.write("1\t1\t1\t" + start + "\t" + end + "\n")
    OP.close()

splitFile("/y/home/luozhihui/try/accurity/bicseq2/seg/output.file")