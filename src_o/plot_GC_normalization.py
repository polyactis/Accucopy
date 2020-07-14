#!/usr/bin/env python
import matplotlib
matplotlib.use("Agg")
import os, sys
import matplotlib.pyplot as plt
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-a', '--adjFactor_filepath', help='the adjust factor file ')
parser.add_argument('-r', '--regin_filepath', help='the regression input file')
parser.add_argument('-o', '--output_filepath', default='None', help='the output png full path, i.e. output/cnv.png')
args = parser.parse_args()


def draw_normalization_plot(inFile_adj_factor=None, reginFile=None, out_png=None):
    plt.figure(1, figsize=(10, 15))
    plt.subplots_adjust(right=0.95, left=0.1, top=0.95, bottom=0.05)
    if os.path.isfile(inFile_adj_factor):
        data = pd.read_table(inFile_adj_factor, header=None)
        if data.size>0:
            # plot the distribution
            ax = plt.subplot(211)
            ax.set_title("adjust factor")
            plt.xlabel("GC-ratio*1000")
            plt.ylabel("adj-factor")
            plt.plot(range(len(data)), data[0], ".", alpha=0.3)
    else:
        sys.stderr.write("WARNING: %s does not exist.\n"%inFile_adj_factor)

    if os.path.isfile(reginFile):
        data_reg = pd.read_table(reginFile, sep="\t", header=0)
        if data_reg.size>0:
            ax = plt.subplot(212)
            ax.set_title("plot the regression input file")
            ax.set_xlabel("GC Ratio")
            ax.set_ylabel("read count)")
            ax.set_yscale('log')
            ax.plot(data_reg.gcRatio, data_reg.readCount, '.', alpha=0.1)
    else:
        sys.stderr.write("WARNING: %s does not exist.\n"%(reginFile))

    plt.savefig(out_png, dpi=100)

draw_normalization_plot(args.adjFactor_filepath, args.regin_filepath, args.output_filepath)