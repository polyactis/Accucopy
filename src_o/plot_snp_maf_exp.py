#!/usr/bin/env python
import matplotlib
matplotlib.use("Agg")
import re, os, sys
import matplotlib.pyplot as plt
import pandas as pd
import argparse
from matplotlib.ticker import MultipleLocator
import numpy as np

def axvlines(ax=None, xs=None, **plot_kwargs):
	"""
	Draw vertical lines on plot
	:param xs: A scalar, list, or 1D array of horizontal offsets
	:param plot_kwargs: Keyword arguments to be passed to plot
	:return: The plot object corresponding to the lines.
	"""
	xs = np.array((xs, ) if np.isscalar(xs) else xs, copy=False)
	lims = ax.get_ylim()
	x_points = np.repeat(xs[:, None], repeats=3, axis=1).flatten()
	y_points = np.repeat(np.array(lims+(np.nan, ))[None, :], repeats=len(xs), axis=0).flatten()
	plot, = ax.plot(x_points, y_points, scaley = False, **plot_kwargs)
	return plot

def axhlines(ys, **plot_kwargs):
	"""
	Draw horizontal lines across plot
	:param ys: A scalar, list, or 1D array of vertical offsets
	:param plot_kwargs: Keyword arguments to be passed to plot
	:return: The plot object corresponding to the lines.
	"""
	ys = np.array((ys, ) if np.isscalar(ys) else ys, copy=False)
	lims = plt.gca().get_xlim()
	y_points = np.repeat(ys[:, None], repeats=3, axis=1).flatten()
	x_points = np.repeat(np.array(lims + (np.nan, ))[None, :], repeats=len(ys), axis=0).flatten()
	plot, = plt.plot(x_points, y_points, scalex = False, **plot_kwargs)
	return plot

def draw_snp_maf_exp_one_period(df, period_int, output_file_path):
	df_of_period = df.loc[df['period_int'] == period_int]
	cps_bf_1st_peak_ls = df_of_period["no_of_copy_nos_bf_1st_peak"].unique()
	cps_bf_1st_peak_ls.sort()
	sys.stderr.write("Drawing for period %s to %s, %s cps_bf_1st_peak ... \n"%
	                 (period_int, output_file_path, len(cps_bf_1st_peak_ls)))

	#maximum number of different no_of_copy_nos_bf_1st_peak hypotheses to be tested
	no_of_hypos = 5
	cps_bf_1st_peak_ls = cps_bf_1st_peak_ls[:no_of_hypos]
	no_of_hypos = len(cps_bf_1st_peak_ls)

	figure, ax_array = plt.subplots(no_of_hypos*2, sharex=True, figsize=(15, 10))
	plt.subplots_adjust(right=0.97, left=0.1, top=0.95, bottom=0.07)
	ax_top = ax_array[0]
	ax_bottom = ax_array[-1]
	xmajorLocator = MultipleLocator(0.1)
	xminorLocator = MultipleLocator(0.02)
	ax_bottom.set_xlim([0.45, 1])

	ax_top.set_title("SNP MAF expected")
	ax_bottom.set_xlabel("SNP MAF")
	for ax in ax_array:
		ax.xaxis.set_major_locator(xmajorLocator)
		ax.xaxis.set_minor_locator(xminorLocator)
		ax.xaxis.grid(True, which='major')

	#add y-axis label
	for i in xrange(len(cps_bf_1st_peak_ls)):
		no_of_copy_nos_bf_1st_peak = cps_bf_1st_peak_ls[i]
		for ax_index in [2*i, 2*i+1]:
			ax = ax_array[ax_index]
			if ax_index%2==0:
				ax.set_ylabel("bf=%s"% no_of_copy_nos_bf_1st_peak)
			else:
				ax.set_ylabel("bf=%s adj"%no_of_copy_nos_bf_1st_peak)


	max_no_of_peaks_to_draw = 5
	for i in xrange(len(cps_bf_1st_peak_ls)):
		no_of_copy_nos_bf_1st_peak = cps_bf_1st_peak_ls[i]
		sys.stderr.write("\tbf=%s\n"%no_of_copy_nos_bf_1st_peak)
		df_of_period_of_cps_bf = df_of_period.loc[df_of_period['no_of_copy_nos_bf_1st_peak'] == no_of_copy_nos_bf_1st_peak]
		peak_index_list = df_of_period_of_cps_bf["peak_index"].unique()
		peak_index_list.sort()
		ax_upper = ax_array[2*i]
		ax_lower = ax_array[2*i+1]
		plot_list = []
		label_list = []
		for peak_index in peak_index_list[:max_no_of_peaks_to_draw]:
			df_of_period_of_cps_bf_of_peak = df_of_period_of_cps_bf.loc[df_of_period_of_cps_bf['peak_index'] == peak_index]
			data_label = "peak %s"%peak_index
			sys.stderr.write("\t%s\n"%data_label)
			plot_upper = axvlines(ax_upper, xs=df_of_period_of_cps_bf_of_peak["major_allele_fraction_exp"],
			                      alpha=0.5, linewidth=2.0)
			plot_lower = axvlines(ax_lower, xs=df_of_period_of_cps_bf_of_peak["maf_exp_adjusted"],
			                      alpha=0.5, linewidth=2.0)
			label_list.append(data_label)
			plot_list.append(plot_upper)
		leg = ax_upper.legend(plot_list, label_list, fancybox=True, prop={'size': 10})
		leg.get_frame().set_alpha(0.5)
	# Fine-tune figure; make subplots close to each other and hide x ticks for
	# all but bottom plot.
	figure.subplots_adjust(hspace=0)
	plt.setp([a.get_yticklabels() for a in figure.axes], visible=False)
	plt.savefig(output_file_path, dpi=200)
	sys.stderr.write("Done.\n")


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input_file_path', help='the path to the snp_maf_exp_vs_adj.tsv file')
	parser.add_argument('-o', '--output_folder', default='None', help='the folder to hold output figures')
	args = parser.parse_args()
	if os.path.isfile(args.input_file_path):
		df = pd.read_table(args.input_file_path, comment="#", header=0)
		period_int_list = df["period_int"].unique()
		period_int_list.sort()
		sys.stderr.write("%s unique periods.\n"%len(period_int_list))
		for period_int in period_int_list:
			output_file_path = os.path.join(args.output_folder, "snp_maf_exp_period_%s.png"%period_int)
			draw_snp_maf_exp_one_period(df, period_int, output_file_path)

	else:
		sys.stderr.write("ERROR: %s does not exist.\n"%args.input_file_path)
		sys.exit(2)