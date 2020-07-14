#!/usr/bin/env python
import matplotlib

matplotlib.use('Agg')
import pandas as pd
from matplotlib import collections as mc
from matplotlib import pyplot as plt
import numpy as np
import argparse
import os, sys
import re

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--truth', default=None, help='the file holding the true CNV profile (simulation)')
parser.add_argument('-i', '--input_filepath', help='the file holding the software CNV result')
parser.add_argument('-o', '--output_filepath', default='None', help='the output png full path, i.e. output/cnv.png')
parser.add_argument('-r', '--ref_dict', default='None', help='the reference genome dictionary store the chrosomes length')
args = parser.parse_args()

def readChrosomeLengthFromDictFile(filename):
# dict file's chromosome id should be well sorted
	chr_len = []
	pattern = re.compile(r'^chr\d+$')
	with open(filename, 'r') as f:
		for line in f:
			if line.startswith('@SQ'):
				records = line.split('\t')
				chr_name = records[1].split(':')[1]
				m = pattern.findall(chr_name)
				if m:
					chr_len.append(int(records[2].split(':')[1]))
					
					# if chr_name.startswith('chr'):
					# 	chr_idx = chr_name[3:]
					# 	if chr_idx in [str(x) for x in range(1,23)]:
					# 		chr_len[int(chr_name[3:]) - 1] = int(records[2].split(':')[1])
	return chr_len



def mylog(x):
	if x > 2:
		x = np.log2(x - 1) + 2
	return x


def readDataFromTruthResult(filename):
	df = pd.read_table(filename, sep='\t', comment='#')
	df = df[~(df['chr'] == 'chrX')]
	df = df[~(df['chr'] == 'chrY')]
	df = df[~(df['chr'] == 'chrMT')]
	df['chr'] = df['chr'].map(lambda x: int(x[3:]))
	return df


def readDataFromAccurityResult(filename):
	df = pd.read_table(filename, sep='\t', comment='#')
	#df = df[[0, 3, 4, 10, 11]]
	#df.columns = ["chr", "cp", "major_allele_cp", "start", "end"]
	# df = df[df['cp']!=2]
	df.loc[np.isnan(df.major_allele_cp), "IsClonal"] = 'F'
	df.loc[~np.isnan(df.major_allele_cp), "IsClonal"] = 'T'
	return df


def genomicSegments(dataframe, chr_len):
	t = dataframe.to_dict('list')
	cnv_list = {}
	for chr_id, start, end in zip(t['chr'], t['start'], t['end']):
		if chr_id not in cnv_list:
			cnv_list[chr_id] = [(start, end)]
		else:
			cnv_list[chr_id].append((start, end))
	chr_list = {}
	for i in range(1, 23):
		chr_list[i] = [1, chr_len[i - 1]]
	for chr in cnv_list:
		for segment in cnv_list[chr]:
			if segment[0] == 1:
				chr_list[chr].remove(1)
				chr_list[chr].append(segment[1] + 1)
			elif segment[0] == chr_len[chr - 1]:
				chr_list[chr].remove(chr_len[chr - 1])
				chr_list[chr].append(segment[0] - 1)
			else:
				chr_list[chr].append(segment[0] - 1)
				chr_list[chr].append(segment[1] + 1)
	for chr in chr_list:
		chr_list[chr].sort()
	data = {'chr': [], 'start': [], 'end': []}
	for chr_id in chr_list:
		position = iter(chr_list[chr_id])
		while True:
			try:
				data['start'].append(position.next())
				data['end'].append(position.next())
			except StopIteration:
				break
			data['chr'].append(chr_id)
	df = pd.DataFrame(data)
	df["cp"] = 2
	df["major_allele_cp"] = 1
	df["IsClonal"] = 'T'
	df = pd.concat([dataframe, df], axis=0).sort_values(by=['chr', 'start'])
	df.index = range(len(df))
	return df


def plotSegment(ax, df, ylabel, ymax, shift):
	lines = []
	c = []
	for ix, col in df.iterrows():
		start = col['start']
		end = col['end']
		cp = col['cp']
		chr_id = int(col['chr'])
		if col['IsClonal'] == 'T':
			c.append((1, 0, 0, 1))
		else:
			c.append((0, 1, 0, 1))
		lines.append([(start + shift[chr_id - 1], mylog(cp)), (end + shift[chr_id - 1], mylog(cp))])
	lc = mc.LineCollection(lines, colors=c, linewidths=2)
	ax.add_collection(lc)
	ax.autoscale()
	for i in range(len(shift)):
		ax.axvline(x=shift[i], color='#999999')
	ax.set_ylim(bottom=-0.1, top=mylog(ymax))
	ax.set_xlim(left=1, right=shift[-1])
	ax.get_xaxis().set_ticks([])
	ax.get_xaxis().set_ticklabels([])
	ax.yaxis.grid(linestyle='-', color='#d9d9d9')
	ax.yaxis.set_label_position('right')
	ax.set_ylabel(ylabel)
	yticks = range(6, int(ymax) + 4, 2)
	ytickslabel = [0, 1, 2, 3, 4, 5] + yticks
	yticks = [mylog(x) for x in ytickslabel]
	ax.set_yticks(yticks)
	ax.set_yticklabels(ytickslabel, fontsize=8)


def plotGenomicSegment(data, output_filepath, chr_len, shift):
	# plot copy number
	ymax = 0
	for key in data:
		for ix, col in data[key].iterrows():
			if col['end'] - col['start'] > 3000000:
				if ymax < col['cp']:
					ymax = col['cp']
	noOfPlotRows = 0
	if data.has_key("truth"):
		noOfPlotRows += 2
	if data.has_key("software"):
		noOfPlotRows += 2
	if noOfPlotRows<=0:
		return

	fig, axes = plt.subplots(nrows=noOfPlotRows, ncols=1)

	plotRowIndex = 0
	if data.has_key("truth"):
		plotSegment(axes[plotRowIndex], data['truth'], "Ground Truth", ymax, shift)
		plotRowIndex += 1
	if data.has_key("software"):
		plotSegment(axes[plotRowIndex], data['software'], 'Accucopy', ymax, shift)
		plotRowIndex += 1

	# plot major allele copy number
	ymax = 0
	for key in data:
		drop_index = []
		for ix, col in data[key].iterrows():
			if col['end'] - col['start'] > 3000000:
				if ymax < col['major_allele_cp']:
					ymax = col['major_allele_cp']
			if col['IsClonal'] == 'F':
				drop_index.append(ix)
		data[key].drop(drop_index, inplace=True)
		data[key].drop('cp', axis=1, inplace=True)
		data[key].rename(columns={'major_allele_cp': 'cp'}, inplace=True)
	if data.has_key("truth"):
		plotSegment(axes[plotRowIndex], data['truth'], "Ground Truth", ymax, shift)
		plotRowIndex += 1
	if data.has_key("software"):
		plotSegment(axes[plotRowIndex], data['software'], "Accucopy", ymax, shift)
		plotRowIndex += 1

	title = os.path.splitext(os.path.split(output_filepath)[1])[0]
	plt.suptitle(title)
	plt.subplots_adjust(hspace=0.1)
	xlabels = [x for x in range(1, 23)]
	x = []
	for i in range(22):
		x.append(shift[i] + chr_len[i] / 2)
	plt.xticks(x, xlabels)
	plt.tick_params(axis='x', which='both', length=0)
	plt.subplots_adjust(right=0.98, left=0.03, top=0.95, bottom=0.08)
	fig.text(0.01, 0.75, 'Absolute copy number', va='center', rotation='vertical')
	fig.text(0.01, 0.3, 'Major allele copy number', va='center', rotation='vertical')
	fig.text(0.5, 0.02, 'Chromosome number', ha='center')
	fig.set_size_inches(17.387, 7.638)
	fig.savefig(output_filepath, format='png', dpi=300)


if __name__ == '__main__':
	chr_len = readChrosomeLengthFromDictFile(args.ref_dict)
	shift = [0]
	for i in range(len(chr_len)):
		shift.append(shift[i] + chr_len[i])
	data = {}
	if args.truth and os.path.isfile(args.truth):
		truth = readDataFromTruthResult(args.truth)
		truth = genomicSegments(truth, chr_len)
		data['truth'] = truth

	if args.input_filepath and os.path.isfile(args.input_filepath):
		accurity = readDataFromAccurityResult(args.input_filepath)
		accurity = genomicSegments(accurity, chr_len)
		data['software'] = accurity

	outputDir = os.path.split(args.output_filepath)[0]

	if not os.path.exists(outputDir):
		os.makedirs(outputDir)
	if len(data)>0:
		plotGenomicSegment(data, args.output_filepath, chr_len, shift)
	else:
		sys.stderr.write("ERROR: Two input files, %s, %s, do not exist! Exit normally.\n"%(args.truth, args.input_filepath))
