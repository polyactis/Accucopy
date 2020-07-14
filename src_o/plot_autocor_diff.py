#!/usr/bin/env python
import matplotlib
matplotlib.use("Agg")
import re, os, sys
import matplotlib.pyplot as plt
import pandas as pd
import argparse
from matplotlib.ticker import MultipleLocator
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_filepath', help='the file holding the period detection GADA input of infer')
parser.add_argument('-a', '--auto_filepath', help='the file holding the auto-correlation result from infer')
parser.add_argument('-s', '--segment_filepath', help='the file holding the period detection GADA segmentation output')
parser.add_argument('-o', '--output_filepath', default='None', help='the output png full path, i.e. output/cnv.png')
args = parser.parse_args()


def get_diff_data(inFile=None):
    infile = open(inFile, "r")
    header = infile.readline().strip("\n\r").split("\t")
    data = pd.DataFrame(columns=header)
    comment_re = re.compile("^#")
    i = 0
    cutoff = -1
    for line in infile:
        if comment_re.search(line):
            cutoff = line.strip("\n\r")
        else:
            oneline = line.strip("\n\r").split("\t")
            data.loc[i] = oneline
            i = i + 1
    return data, cutoff

def get_auto_data(inFile=None):
    #names=["read_count_ratio",	"correlation"]
    data = pd.read_table(inFile, header=0)
    #data.window_count_smoothed = np.round(data.window_count_smoothed, decimals=0)
    data["correlation_log10"] = np.log10(data.correlation)
    return data

def draw_diff_plot(data=None, cutoff=None, segFile=None, out_png=None, auto_data=None):
    shif_diff_thresold_pair = cutoff.strip("\n\r").split(" : ")[1].split()
    #0 is the firm separator
    shif_diff_thresold_pair[0] = min(0, float(shif_diff_thresold_pair[0]))
    shif_diff_thresold_pair[1] = max(0, float(shif_diff_thresold_pair[1]))

    data.cor_shift_diff = data.cor_shift_diff.astype('float64')
    data.round_int = data.round_int.astype('int')
    seg = pd.read_table(segFile, header=0)

    plt.figure(1, figsize=(12, 17))
    #figure, ax_array = plt.subplots(2, sharex=True)
    plt.subplots_adjust(right=0.95, left=0.07, top=0.95, bottom=0.05)
    # plot the distribution
    ax = plt.subplot(321)
    ax.set_xlabel("shift-1 diff of log10(auto-cor)")
    ax.set_ylabel("count")
    ax.xaxis.grid(True, which='major')
    ax.yaxis.grid(True, which='major')
    ax.hist(data.cor_shift_diff, bins=200)

    #plot the tre
    ax = plt.subplot(322)
    ax.cla()
    ax.set_xlabel("period gap")
    ax.set_ylabel("auto-cor")
    ax.set_xlim([0, 1])
    #locs, labels = ax.xticks()
    #ax.set_xticks(np.arange(0, 3, 0.25).tolist(), np.arange(0, 3, 0.25).tolist())
    #ax.tick_params(axis="x", which="both", width=0.05)
    xmajorLocator = MultipleLocator(0.2)
    xminorLocator = MultipleLocator(0.05)
    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_minor_locator(xminorLocator)
    ax.xaxis.grid(True, which='major')
    #ax.yaxis.grid(True, which='major')
    ax.plot(auto_data.read_count_ratio, auto_data.correlation, ".", markerfacecolor="None")

    #plot autocorrelation log 10
    ax = plt.subplot(323)
    ax.cla()
    ax.set_xlabel("period gap")
    ax.set_ylabel("log10(auto-cor)")
    ax.set_yscale('log')
    ax.set_xlim([0, 1])
    #locs, labels = ax.xticks()
    #ax.set_xticks(np.arange(0, 3, 0.25).tolist(), np.arange(0, 3, 0.25).tolist())
    #ax.tick_params(axis="x", which="both", width=0.05)
    xmajorLocator = MultipleLocator(0.2)
    xminorLocator = MultipleLocator(0.05)
    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_minor_locator(xminorLocator)
    ax.xaxis.grid(True, which='major')
    #ax.yaxis.grid(True, which='major')
    ax.plot(auto_data.read_count_ratio, auto_data.correlation, ".", markerfacecolor="None")

    # plot the dot dot graph of diff
    ax = plt.subplot(324)
    ax.set_xlabel("period gap * 1000")
    ax.set_ylabel("shift-1 diff of log10(auto-cor)")
    ax.set_xlim([0, 1000])
    ax.plot(data.cor_shift_diff, '.', markerfacecolor="None")
    ax.set_xlabel(cutoff)
    xmajorLocator = MultipleLocator(200)
    xminorLocator = MultipleLocator(50)
    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_minor_locator(xminorLocator)
    ax.xaxis.grid(True, which='major')
    ax.hlines(shif_diff_thresold_pair[0], 0, 1000)
    ax.hlines(shif_diff_thresold_pair[1], 0, 1000)

    # plot the cutoff plot
    ax = plt.subplot(326)
    ax.set_title("black line is GADA output, green dot is cutoff value")
    ax.set_ylim([-1.1, 1.1])
    ax.set_xlim([0, 1000])
    ax.set_xlabel("period gap * 1000")
    ax.set_ylabel("shift-1 diff of log10(auto-cor)")
    xmajorLocator = MultipleLocator(200)
    xminorLocator = MultipleLocator(50)
    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_minor_locator(xminorLocator)
    ax.xaxis.grid(True, which='major')
    ax.plot(data.round_int[(data.cor_shift_diff < shif_diff_thresold_pair[0]) | (data.cor_shift_diff > shif_diff_thresold_pair[1])],
            '.', markerfacecolor="None", alpha=0.1)
    for i in range(len(seg)):
        #print seg.iloc[i, 3], seg.iloc[i, 0], seg.iloc[i, 1]
        ax.hlines(seg.iloc[i, 3], seg.iloc[i, 0], seg.iloc[i, 1])

    # save png
    plt.savefig(out_png, dpi=250)


if os.path.isfile(args.input_filepath) and os.path.isfile(args.auto_filepath) and os.path.isfile(args.segment_filepath):
    data_tuple = get_diff_data(args.input_filepath)
    auto_daraframe = get_auto_data(args.auto_filepath)
    draw_diff_plot(data=data_tuple[0], cutoff=data_tuple[1], segFile=args.segment_filepath,
                   out_png=args.output_filepath, auto_data=auto_daraframe)
else:
    sys.stderr.write("ERROR: At least one of the three input files does not exist! Exit normally.\n")
