#!/usr/bin/env python

import os, sys
import argparse


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--main_out_file', required=True, type=str, help='the main output infer.out file')
	parser.add_argument('-u', '--true_purity', required=True, help='true purity')
	parser.add_argument('-d', '--true_ploidy', required=True, help='true ploidy')
	parser.add_argument('-v', '--max_purity_delta', default=0.05, help="deviation between true and estimated purity.")
	parser.add_argument('-e', '--max_ploidy_delta', default=0.3, help="deviation between true and estimated ploidy.")
	args = parser.parse_args()
	commit_id = os.environ.get('CI_COMMIT_SHA')[:8]
	true_purity = float(args.true_purity)
	true_ploidy = float(args.true_ploidy)
	max_purity_delta = float(args.max_purity_delta)
	max_ploidy_delta = float(args.max_ploidy_delta)

	sys.stderr.write("Judge result, %s, of commit ID %s with true purity=%s, ploidy=%s ...\n"%(args.main_out_file, commit_id, true_purity, true_ploidy))

	input_file = open(args.main_out_file, 'r')
	line_1 = input_file.next()
	line_2 = input_file.next()
	result_ls = line_2.strip().split('\t')
	main_purity=float(result_ls[0])
	main_ploidy=float(result_ls[1])

	delta_purity=abs(true_purity - main_purity)
	delta_ploidy=abs(true_ploidy - main_ploidy)
	sys.stderr.write("Estimated purity=%s, ploidy=%s.\n"%(main_purity, main_ploidy))
	exit_code = 0
	if delta_purity>max_purity_delta:
		sys.stderr.write("Purity estimate %s deviated from the truth %s by %s, more than %s.\n"%(main_purity, true_purity, delta_purity, max_purity_delta))
		exit_code=1

	if delta_ploidy>max_ploidy_delta:
		sys.stderr.write("Ploidy estimate %s deviated from the truth %s by %s, more than %s.\n"%(main_ploidy, true_ploidy, delta_ploidy, max_ploidy_delta))
		exit_code=1

	sys.stderr.write("exit code is %s.\n"%exit_code)
	sys.exit(exit_code)
