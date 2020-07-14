#!/usr/bin/env python
import matplotlib
matplotlib.use("Agg")
import os, sys
import matplotlib.pyplot as plt
import pandas as pd
import argparse
from matplotlib.ticker import MultipleLocator
import math

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--coverage_ratio_file_path', help='the file containing the smoothed coverage ratio histogram data ')
parser.add_argument('-p', '--peak_bounds_file_path', help='the file containing coverage ratio peak bounds  ')
parser.add_argument('-o', '--output_file_path', default='None', help='the output figure full path, i.e. output.png, output.jpg')
args = parser.parse_args()


def draw_normalization_plot(input_file_path=None, peak_bounds_file_path=None, output_file_path=None):
	#plt.figure(1, figsize=(20, 15))
	if os.path.isfile(input_file_path):
		data = pd.read_table(input_file_path, sep="\t", comment="#")
		data = data[data.window_count_smoothed != 0]
		if data.size>0:
			# plot the distribution
			figure, ax_array = plt.subplots(2, sharex=True)
			plt.subplots_adjust(right=0.95, left=0.12, top=0.95, bottom=0.07)
			ax1 = ax_array[0]
			ax1.set_title("TRE ratio histogram")
			ax1.set_xlabel("TRE ratio")
			ax1.set_ylabel("Window count")
			ax1.plot(data.read_count_ratio, data.window_count_smoothed, ".", markerfacecolor="None", alpha=0.5)

			data.window_count_smoothed = [math.log10(x) for x in data.window_count_smoothed]
			ax2 = ax_array[1]
			ax2.set_xlabel("TRE ratio")
			ax2.set_ylabel("log10(Window count)")
			# ax2.set_yscale("log")
			xmajorLocator = MultipleLocator(0.5)
			xminorLocator = MultipleLocator(0.1)
			ax2.xaxis.set_major_locator(xmajorLocator)
			ax2.xaxis.set_minor_locator(xminorLocator)
			ax2.plot(data.read_count_ratio, data.window_count_smoothed, ".", markerfacecolor="None", alpha=0.5)
			if os.path.isfile(peak_bounds_file_path):
				data_peak_bounds = pd.read_table(peak_bounds_file_path, sep="\t", comment="#")
				if data_peak_bounds.size>0:
					for i in xrange(data_peak_bounds.shape[0]):
						ax1.axvspan(data_peak_bounds.loc[i].lowerBound, data_peak_bounds.loc[i].upperBound,
						           alpha=0.3, color='red')
						ax2.axvspan(data_peak_bounds.loc[i].lowerBound, data_peak_bounds.loc[i].upperBound,
						            alpha=0.3, color='red')
			else:
				sys.stderr.write("WARNING: %s does not exist.\n"%peak_bounds_file_path)
			# Fine-tune figure; make subplots close to each other and hide x ticks for
			# all but bottom plot.
			figure.subplots_adjust(hspace=0)
			#plt.setp([a.get_xticklabels() for a in figure.axes[:-1]], visible=False)
			plt.savefig(output_file_path, dpi=200)
	else:
		sys.stderr.write("WARNING: %s does not exist.\n"%input_file_path)




draw_normalization_plot(args.coverage_ratio_file_path, args.peak_bounds_file_path, args.output_file_path)