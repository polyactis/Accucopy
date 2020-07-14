#!/usr/bin/env python
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_file_path', help='the file holding coverage information ,chr22.ratio.w500.csv.gz')
parser.add_argument('-o', '--output_file_path', default='None', help='the output png full path, i.e. coverage_compare_by_adjust.png')
args = parser.parse_args()

def read_count_change(inFile=None, output_file_path=None):
    data = pd.read_table(inFile, comment="#", compression="infer", header=0, sep=",")

    figure, ax_array = plt.subplots(1, 2, figsize=(15, 8))
    plt.subplots_adjust(right=0.95, left=0.07, top=0.95, bottom=0.07)
    ax1 = ax_array[0]
    ax1.set_title("Raw")
    ax1.set_xlabel("coverage normal")
    ax1.set_ylabel("coverage tumor")
    ax1.plot(data.coverage_normal, data.coverage_tumor, ".", alpha=0.01)

    #after adjust
    ax2 = ax_array[1]
    ax2.set_title("Post adjusment")
    ax2.set_xlabel("coverage normal adj")
    ax2.set_ylabel("coverage tumor adj")
    ax2.plot(data.coverage_normal_adj, data.coverage_tumor_adj, ".", alpha=0.01)

    #save ong
    # Fine-tune figure; make subplots close to each other and hide x ticks for
    # all but bottom plot.
    #figure.subplots_adjust(hspace=0)
    #plt.setp([a.get_xticklabels() for a in figure.axes[:-1]], visible=False)
    plt.savefig(output_file_path, dpi=200)

read_count_change(args.input_file_path, args.output_file_path)