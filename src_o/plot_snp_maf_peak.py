#!/usr/bin/env python
import matplotlib
matplotlib.use("Agg")
import re, os, sys
import matplotlib.pyplot as plt
import pandas as pd
import argparse
from matplotlib.ticker import MultipleLocator
import numpy as np
import matplotlib.ticker as ticker

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_folder_path', help='the folder holding snp_maf_pdf_of_peak_xx.tsv files')
parser.add_argument('-o', '--output_file_path', default='None', help='the output figure full path, i.e. output.png or output.jpg')
args = parser.parse_args()


if os.path.isdir(args.input_folder_path):
	figure, ax_array = plt.subplots(2, sharex=True)
	#, figsize=(15, 10)
	plt.subplots_adjust(right=0.97, left=0.1, top=0.95, bottom=0.07)
	ax_upper = ax_array[0]
	ax_lower = ax_array[1]
	xmajorLocator = MultipleLocator(0.1)
	xminorLocator = MultipleLocator(0.02)

	ax_upper.set_title("SNP MAF histogram")
	ax_upper.set_ylabel("Smoothed SNP count")
	ax_upper.xaxis.set_major_locator(xmajorLocator)
	ax_upper.xaxis.set_minor_locator(xminorLocator)
	ax_upper.xaxis.grid(True, which='major')

	ax_lower.set_ylabel("log10(SNP count)")
	ax_lower.set_yscale("log")
	ax_lower.set_xlabel("SNP MAF")
	ax_lower.set_xscale("log")
	ax_lower.xaxis.set_major_locator(xmajorLocator)
	ax_lower.xaxis.set_minor_locator(xminorLocator)
	ax_lower.xaxis.grid(True, which='major')
	#ax_lower.xaxis.set_ticks(np.arange(0.45, 1.0, 0.1))    #set_ticks() will mess up xtick labels
	ax_lower.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
	ax_lower.set_xlim([0.45, 1.0])
	re_peak_index = re.compile("snp_maf_pdf_of_peak_(?P<peak_index>[0-9]*).tsv")
	plot_list = []
	label_list = []
	for peak_index in xrange(9):
		#max 9 peaks. stop plotting the rest (insignificant)
		abs_path = os.path.join(args.input_folder_path, "snp_maf_pdf_of_peak_%s.tsv"%peak_index)
		sys.stderr.write("Making plot for %s ... "%abs_path)
		if os.path.isfile(abs_path):
			data_label = "peak %s"%peak_index
			df = pd.read_table(abs_path, comment="#", header=0)
			if (df["count"] > 0).any(axis=0):
				plot_upper, = ax_upper.plot(df["maf"], df["count"], ".", markerfacecolor="None", alpha=0.3,
											label=data_label, markeredgewidth=1.5)
				plot_lower, = ax_lower.plot(df["maf"], df["count"], ".", markerfacecolor="None", alpha=0.3,
											label=data_label, markeredgewidth=1.5)
				label_list.append(data_label)
				plot_list.append(plot_lower)
				sys.stderr.write("\n")
			else:
				sys.stderr.write("Empty data. Skipped.\n")
		else:
			sys.stderr.write("Not present.\n")

	ax_upper.legend(plot_list, label_list)
	# Fine-tune figure; make subplots close to each other and hide x ticks for
	# all but bottom plot.
	figure.subplots_adjust(hspace=0)
	#plt.setp([a.get_xticklabels() for a in figure.axes[:-1]], visible=False)
	plt.savefig(args.output_file_path, dpi=200)
else:
	sys.stderr.write("ERROR: %s is not a folder.\n"%args.input_folder_path)
	sys.exit(2)