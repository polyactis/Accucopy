#!/usr/bin/env python
import argparse
import os
import shutil
import sys
from pyflow import WorkflowRunner
sys.path.insert(0, os.getcwd())
from software.main import MainFlow
binary_folder = "software"

class TCGATestFlow(WorkflowRunner):
	def __init__(self, alignment_folder=None, tumor_bam_file_prefix=None, normal_bam_file_prefix=None,
	             call_snps=False,
				 output_dir=None, segment_stddev_divider=10, commit_id=None, nCores=4, ref_version=None):
		self.alignment_folder = alignment_folder
		self.tumor_bam_file_prefix = tumor_bam_file_prefix
		self.normal_bam_file_prefix = normal_bam_file_prefix
		self.call_snps = call_snps
		self.output_dir = output_dir
		self.segment_stddev_divider = segment_stddev_divider
		self.commit_id = commit_id
		self.nCores = nCores
		self.ref_version = ref_version
		self.configure_file_path = "configure"

	def workflow(self):
		############################################
		# STEP 0: Generate the configure file #
		############################################
		sys.stderr.write("Generating the configure file\n")
		with open(self.configure_file_path,'w') as f:
			f.write("read_length\t101\n")
			f.write("window_size\t500\n")
			f.write("reference_folder_path\t/y/hlab/AccurityTestData/%s\n" % self.ref_version)
			f.write("samtools_path\t/y/program/bin/samtools\n")
			f.write("strelka_path\t/y/program/strelka-2.9.10.centos6_x86_64\n")
			f.write("binary_folder\t%s\n" % binary_folder)


		##############################
		# For each sub sample, do....#
		##############################
		tumor_bam_file_source_path = os.path.join(self.alignment_folder, "%s.bam"%self.tumor_bam_file_prefix)
		tumor_bai_file_source_path = os.path.join(self.alignment_folder, "%s.bam.bai"%self.tumor_bam_file_prefix)

		normal_bam_file_source_path = os.path.join(self.alignment_folder, "%s.bam"%self.normal_bam_file_prefix)
		normal_bai_file_source_path = os.path.join(self.alignment_folder, "%s.bam.bai"%self.normal_bam_file_prefix)

		if not os.path.isfile(tumor_bam_file_source_path):
			sys.stderr.write("ERROR: Tumor sample %s does not exist.\n" % tumor_bam_file_source_path)
			sys.exit(4)
		if not os.path.isfile(normal_bam_file_source_path):
			sys.stderr.write("ERROR: Normal sample %s does not exist.\n" % normal_bam_file_source_path)
			sys.exit(4)


		tumor_bam = "tumor.bam"
		normal_bam = "normal.bam"

		os.symlink(tumor_bam_file_source_path, tumor_bam)
		os.symlink(tumor_bai_file_source_path, "tumor.bam.bai")
		os.symlink(normal_bam_file_source_path, normal_bam)
		os.symlink(normal_bai_file_source_path, "normal.bam.bai")

		#default starting from step 1 (call SNPs)
		main_starting_step = 1
		snp_vcf_strelka_file_source_path = os.path.join(self.alignment_folder, "%s_strelka_snp.vcf.gz"%self.tumor_bam_file_prefix)
		if os.path.isfile(snp_vcf_strelka_file_source_path) and not self.call_snps:
			#create dir to put the symbol link to variants.vcf.gz
			os.makedirs(os.path.join("strelka_snp/results/variants"))
			snp_vcf_strelka_file_dst_path = os.path.join("strelka_snp/results/variants/variants.vcf.gz")
			os.symlink(snp_vcf_strelka_file_source_path, snp_vcf_strelka_file_dst_path)
			main_starting_step = 2


		sys.stderr.write("Working on %s vs %s, output to %s ...\n"%
		                 (self.tumor_bam_file_prefix,
		                  self.normal_bam_file_prefix, self.output_dir))

		########################
		# STEP 3: Run the main job #
		########################
		sys.stderr.write("##### Starting step %s.\n"%main_starting_step)
		main_flow = MainFlow(self.configure_file_path, tumor_bam, normal_bam,
		                             output_dir = self.output_dir,
		                             segment_stddev_divider=self.segment_stddev_divider,
		                             snp_coverage_min=2,
		                             snp_coverage_var_vs_mean_ratio=10.0,
		                             clean=0, step=main_starting_step, debug=1, 
									 nCores=self.nCores)
		main_flow.readConfigureFile(self.configure_file_path)
		main_flow.readDictFile()
		main_job = self.addWorkflowTask(("Main.%s"%self.tumor_bam_file_prefix).replace(".", "_"),
		                                   main_flow)

		####################
		# STEP 7: Clean....#
		####################
		sys.stderr.write("Clean ...\n")
		cmd = "mv %s/*ratio.w*.csv.gz %s/../"%(self.output_dir, self.output_dir)
		clean_job = self.addTask(("Clean.%s"%self.tumor_bam_file_prefix).replace(".", "_"), cmd, dependencies=main_job)

	def print_result_to_console(self):
		"""
		print
		    .../infer.out.tsv
            .../infer.out.details.tsv

		:return:
		"""
		result_file = open(os.path.join(self.output_dir, "infer.out.tsv"))
		sys.stderr.write("Content of infer.out.tsv:\n")
		for line in result_file:
			sys.stderr.write(line)
		result_file.close()

		result_file = open(os.path.join(self.output_dir, "infer.out.details.tsv"))
		sys.stderr.write("First 15 lines of infer.out.details.tsv:\n")
		counter = 0
		for line in result_file:
			sys.stderr.write(line)
			counter += 1
			if counter>=15:
				break


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-f', '--alignment_folder', default="/y/Sunset/db/individual_alignment/", type=str, help='individual alignment folder')
	parser.add_argument('-t', '--tumor_bam_file_prefix', required=True, help='prefix of the tumor bam file basename')
	parser.add_argument('-n', '--normal_bam_file_prefix', required=True, help='prefix of the normal bam file basename')
	parser.add_argument('-c', '--nCores', type=int, default=4, help="the max number of CPUs to use.")
	parser.add_argument('--call_snps',  action='store_true', dest='call_snps', help="If toggled, SNPs will be called.")
	parser.add_argument('-o', '--output_dir', type=str, required=True, help="output directory to hold the main output")
	parser.add_argument('-w', '--pyflow_dir', type=str, required=True, help="directory to hold pyflow output")
	parser.add_argument("--segment_stddev_divider", type=float, default=10.0, help="artificially reduce segment noise level")
	parser.add_argument('-r', '--ref_version', type=str, default="hs37d5", help="The version of the referece data")
	args = parser.parse_args()
	env_dist = os.environ
	wflow = TCGATestFlow(alignment_folder=args.alignment_folder, tumor_bam_file_prefix=args.tumor_bam_file_prefix,
	                     normal_bam_file_prefix=args.normal_bam_file_prefix, call_snps=args.call_snps, 
						 output_dir=args.output_dir,
	                     segment_stddev_divider=args.segment_stddev_divider,
	                     commit_id=env_dist.get('CI_BUILD_REF')[:8], nCores=args.nCores, ref_version=args.ref_version)
	retval = wflow.run(mode="local", nCores=args.nCores, dataDirRoot=args.pyflow_dir, isContinue='Auto',
	                   isForceContinue=True, retryMax=0)
	wflow.print_result_to_console()
	sys.exit(retval)
