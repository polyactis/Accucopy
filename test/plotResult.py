#!/usr/bin/env python
import matplotlib

matplotlib.use('Agg')
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file', help='the file holding result')
parser.add_argument('-o', '--outdir', default='.', help='the output directory')
args = parser.parse_args()

right = 0.99
left = 0.05
top = 0.93
bottom = 0.04
purity_max_diff = 0.05
ploidy_max_diff = 0.2
recall_min = 0.9
precision_min = 0.9


def GroupData(df, class_label,commit_id_ref):
	df_gr = df.groupby(['sample_class', 'truth_purity'])
	data = {}
	commit_id_dict = {}
	for each_record in df_gr[[class_label,'commit_id']]:
		sample_class = each_record[0][0]
		sample_purity = str(each_record[0][1])
		datapoint = []
		commit_id_list = []
		for index,row in each_record[1][['commit_id',class_label]].iterrows():
			datapoint.append(row[class_label])
			commit_id_list.append(commit_id_ref[row['commit_id']])
		if sample_class not in data:
			data[sample_class] = {sample_purity:datapoint}
		else:
			data[sample_class][sample_purity]=datapoint
		if sample_class not in commit_id_dict:
			commit_id_dict[sample_class] = {sample_purity:commit_id_list}
		else:
			commit_id_dict[sample_class][sample_purity]=commit_id_list
	return data, commit_id_dict


def plothline(ax, class_label, data):
	if class_label == 'delta_purity':
		if max(data) > purity_max_diff:
			ax.axhline(y=purity_max_diff, color='#EE1111')
			ax.set_ylim(min(min(data), purity_max_diff) - 0.05, max(max(data), purity_max_diff) + 0.05)
	elif class_label == 'delta_ploidy':
		if max(data) > ploidy_max_diff:
			ax.axhline(y=ploidy_max_diff, color='#EE1111')
			ax.set_ylim(min(min(data), ploidy_max_diff) - 0.05, max(max(data), ploidy_max_diff) + 0.05)
	elif class_label == 'cnv_recall':
		if min(data) < recall_min:
			ax.axhline(y=recall_min, color='#EE1111')
			ax.set_ylim(min(min(data), recall_min) - 0.05, max(max(data), recall_min) + 0.05)
	elif class_label == 'cnv_precision':
		if min(data) < precision_min:
			ax.axhline(y=precision_min, color='#EE1111')
			ax.set_ylim(min(min(data), precision_min) - 0.05, max(max(data), precision_min) + 0.05)


def plotfigure(df, class_label, commit_id_ref, outdir):
	data, commit_id_dict = GroupData(df, class_label, commit_id_ref)
	data_class = df.sample_class.unique()
	originpath = os.path.abspath('.')
	os.chdir(outdir)
	isExistsSubdir = os.path.exists(class_label)
	if not isExistsSubdir:
		os.makedirs(class_label)
	os.chdir(class_label)
	plotsubplot(data, data_class, class_label, commit_id_dict)
	os.chdir(originpath)


def plotsubplot(df, df_class, class_label, commit_id_dict):
	for each_class in df_class:
		data = df[each_class]
		data_size = len(data)
		data_label = data.keys()
		data_label.sort()
		subfigsize = int(np.ceil(np.sqrt(data_size)))
		filename = each_class + '.png'
		fig, axes = plt.subplots(nrows=subfigsize, ncols=subfigsize)
		for i,label in zip(range(data_size),data_label):
			rownum = int(i / subfigsize)
			colnum = int(i - rownum * subfigsize)
			if data_size == 1:
				ax = axes
			else:
				ax = axes[rownum][colnum]
			plot_data = data[label]
			commit_id_list = commit_id_dict[each_class][label]
			title = 'purity' + label
			ax.plot(range(len(commit_id_list)), plot_data, 'bo',markersize=3)
			plothline(ax, class_label, plot_data)
			ax.set_xlim([-0.05, len(commit_id_list) - 1 + 0.05])
			ax.set_title(title)
			ax.set_ylabel(class_label)
		plt.suptitle(each_class)
		fig.set_size_inches(19.96, 9.29)
		plt.subplots_adjust(right=right, left=left, top=top, bottom=bottom)
		fig.savefig(filename, format='png', dpi=300)
		plt.close()


if __name__ == '__main__':
	df = pd.read_table(args.file)
	commit_id = []
	commit_id_ref = dict()
	for commit in df.commit_id:
		if commit not in commit_id:
			commit_id.append(commit)
	path = args.outdir.strip()
	isExists = os.path.exists(path)
	if not isExists:
		os.makedirs(path)
	class_label_list = ['Accurity_purity', 'Accurity_ploidy', 'delta_purity', 'delta_ploidy', 'cnv_recall',
						'cnv_precision', 'period', 'first_peak_int', 'likelihood']

	i = 0
	with open(os.path.join(path, 'commit_id_to_num.tsv'), 'w') as f:
		f.write('commit_id\tnum_in_plot\n')
		for commit in commit_id:
			f.write(commit + '\t' + str(i) + '\n')
			i += 1
			commit_id_ref[commit]=i

	for class_label in class_label_list:
		plotfigure(df, class_label, commit_id_ref, path)
