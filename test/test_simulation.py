#!/usr/bin/env python
import argparse
import os
import shutil
import sys
from pyflow import WorkflowRunner
sys.path.insert(0, os.getcwd())
from software.main import MainFlow
binary_folder = "software"
class SampleTestFlow(WorkflowRunner):
    def __init__(self, simulation_type=None, purity=0.1,
        segment_stddev_divider=10, call_snps=False, 
        commit_id=None, debug=0, nCores=4, ref_version=None):
        self.simulation_type = simulation_type
        self.purity = purity
        self.segment_stddev_divider = segment_stddev_divider
        self.call_snps = call_snps
        self.commit_id = commit_id
        self.debug = debug
        self.nCores = nCores
        self.ref_version = ref_version

        self.simulation_data = "/y/hlab/AccurityTestData/simulation/%s" % \
            self.simulation_type
        self.sample_dir = self.simulation_type
        self.result_file = "%s/result_%s_%s_%s.tsv" % (self.sample_dir,
            self.commit_id, self.simulation_type, self.purity)
        self.configure_file_path = "configure"

        self.sample_name = "purity.%s"%self.purity
        self.main_output_dir = "%s/%s/infer_result_debug"%(self.sample_dir,
            self.sample_name)
        self.main_output_parent_dir = "%s/%s"%(self.sample_dir, self.sample_name)

    def workflow(self):
        ############################################
        # STEP 0: Generate Accurity configure file #
        ############################################
        sys.stderr.write("Generating the configure file\n")
        ref_folder_path = "/y/hlab/AccurityTestData/%s" % self.ref_version
        with open(self.configure_file_path,'w') as f:
            f.write("read_length\t101\n")
            f.write("window_size\t500\n")
            f.write("reference_folder_path\t%s\n" % ref_folder_path)
            f.write("samtools_path\t/y/program/bin/samtools\n")
            f.write("strelka_path\t/y/program/strelka-2.9.10.centos6_x86_64\n")
            f.write("binary_folder\t%s\n" % binary_folder)

        #####################################################
        # STEP 1: Generate sample directory and result file #
        #####################################################
        if os.path.isfile(self.sample_dir):
            os.remove(self.sample_dir)
            os.mkdir(self.sample_dir)
        elif os.path.isdir(self.sample_dir):
            shutil.rmtree(self.sample_dir)
            os.mkdir(self.sample_dir)
        else:
            os.mkdir(self.sample_dir)
        with open(self.result_file,'w') as f:
            f.write("commit_id\tsample_class\ttruth_purity\ttruth_ploidy\t"
                "Accurity_purity\tAccurity_ploidy\t"
                "delta_purity\tdelta_ploidy\tcp_recall\tcp_precision\t"
                "mcp_recall\tmcp_precision\t"
                "period\tfirst_peak_int\tlikelihood\n")

        ##############################
        # For each sub sample, do....#
        ##############################
        main_job_list = []

        sys.stderr.write("Working on %s ...\n" % self.main_output_parent_dir)
        ######################################################################
        # STEP 2: Generate directory for each sub sample and set symbol link #
        ######################################################################
        tumor_source_folder = "%s/%s" % (self.simulation_data, self.sample_name)
        normal_source_folder = "%s/normalSample" % (self.simulation_data)
        if not os.path.isdir(tumor_source_folder):
            sys.stderr.write("ERROR: no sample: %s %s.\n" % (
                self.simulation_data, self.sample_name))
            sys.exit(2)
        # link the data file
        if os.path.isfile(self.main_output_parent_dir):
            os.remove(self.main_output_parent_dir)
            os.mkdir(self.main_output_parent_dir)
        elif os.path.isdir(self.main_output_parent_dir):
            shutil.rmtree(self.main_output_parent_dir)
            os.mkdir(self.main_output_parent_dir)
        else:
            os.mkdir(self.main_output_parent_dir)

        normal_bam = self.main_output_parent_dir+"/normal.bam"
        tumor_bam = self.main_output_parent_dir+"/tumor.bam"
        os.symlink(os.path.join(normal_source_folder, 
            "normal.aligned.sort.bam"), normal_bam)
        os.symlink(os.path.join(normal_source_folder, 
            "normal.aligned.sort.bam.bai"),
            self.main_output_parent_dir+"/normal.bam.bai")
        os.symlink(tumor_source_folder+"/tumor.aligned.sort.bam", tumor_bam)
        os.symlink(tumor_source_folder+"/tumor.aligned.sort.bam.bai", 
            self.main_output_parent_dir+"/tumor.bam.bai")

		#default starting from step 1 (call SNPs)
        accurity_starting_step = 1
        snp_vcf_strelka_file_path = os.path.join(tumor_source_folder, 
            "strelka_snp.vcf.gz")
        if os.path.isfile(snp_vcf_strelka_file_path) and not self.call_snps:
            # if call_snps==True and the source SNP file exists.
            #create dir to put the symbol link to variants.vcf.gz
            os.makedirs(os.path.join(self.main_output_parent_dir, 
                "strelka_snp/results/variants"))
            os.symlink(snp_vcf_strelka_file_path, os.path.join(
                self.main_output_parent_dir,
                "strelka_snp/results/variants/variants.vcf.gz"))
            accurity_starting_step = 2

        ########################
        # STEP 3: Run the main job #
        ########################
        sys.stderr.write("##### Starting step %s.\n"%(accurity_starting_step))
        main_flow = MainFlow(self.configure_file_path, tumor_bam, normal_bam,
            output_dir = self.main_output_dir,
            segment_stddev_divider=self.segment_stddev_divider,
            snp_coverage_min=2,
            snp_coverage_var_vs_mean_ratio=10,
            clean=0, step=accurity_starting_step, debug=self.debug,
            nCores=self.nCores)
        main_flow.readConfigureFile(self.configure_file_path)
        main_flow.readDictFile()
        parent_job = None
        if len(main_job_list)>0:
            parent_job = main_job_list[-1]
        main_job = self.addWorkflowTask(("Main.%s"%self.sample_name).replace(".", "_"),
            main_flow, dependencies=parent_job)
        main_job_list.append(main_job)

        #####################################
        # STEP 4: Call recall and precision #
        #####################################
        sys.stderr.write("Calculate cnv_recall and cnv_precision\n")
        cmd = "%s/maestre recall_precision -t %s/truth_cnv.tab -p %s/cnv.output.tsv " \
            "-o %s/recall_precision.tsv || exit 1" % \
            (binary_folder, self.simulation_data, self.main_output_dir,
            self.main_output_dir)
        CallRecallAndPrecisionJob = self.addTask(
            ("CallRecallAndPrecision.%s"%self.sample_name).replace(".", "_"),
            cmd, dependencies=main_job)

        ###########################
        # STEP 5: Analysis Result #
        ###########################
        sys.stderr.write("Analyze result and write result to %s\n" % self.result_file)
        cmd = "cp_recall=$(cat %s/recall_precision.tsv|awk 'NR>1'|cut -f1);" \
            "cp_precision=$(cat %s/recall_precision.tsv|awk 'NR>1'|cut -f2);" \
            "mcp_recall=$(cat %s/recall_precision.tsv|awk 'NR>1'|cut -f3);" \
            "mcp_precision=$(cat %s/recall_precision.tsv|awk 'NR>1'|cut -f4);" \
            'Accurity_purity=$(sed -n "2p" %s/infer.out.tsv|cut -f1);' \
            'Accurity_ploidy=$(sed -n "2p" %s/infer.out.tsv|cut -f2);' \
            'truth_purity=%s;truth_ploidy=$(sed -n "2p" %s/ploidy.txt|cut -f3);' \
            "delta_purity=$(echo ${Accurity_purity} ${truth_purity}|awk '{print $1-$2}');" \
            "delta_ploidy=$(echo ${Accurity_ploidy} ${truth_ploidy}|awk '{print $1-$2}');" \
            "delta_purity=${delta_purity#-};delta_ploidy=${delta_ploidy#-};" \
            'period=$(sed -n "4p" %s/infer.out.tsv|cut -f2);' \
            'first_peak_int=$(sed -n "4p" %s/infer.out.tsv|cut -f4);' \
            'likelihood=$(sed -n "4p" %s/infer.out.tsv|cut -f1);' \
            'echo -e "%s\t%s\t${truth_purity}\t${truth_ploidy}\t${Accurity_purity}\t' \
            '${Accurity_ploidy}\t${delta_purity}\t${delta_ploidy}\t${cp_recall}\t'\
            '${cp_precision}\t${mcp_recall}\t${mcp_precision}\t${period}\t' \
            '${first_peak_int}\t${likelihood}" >> %s' % \
            (self.main_output_dir, self.main_output_dir, self.main_output_dir,
            self.main_output_dir, self.main_output_dir, self.main_output_dir,
            self.purity, self.simulation_data, self.main_output_dir,
            self.main_output_dir, self.main_output_dir, self.commit_id,
            self.simulation_type, self.result_file)
        AnalyzeResultJob = self.addTask(
            ("AnalyzeResult.%s"%self.sample_name).replace(".", "_"), cmd,
            dependencies=CallRecallAndPrecisionJob)

        ##################################################
        # STEP 6: Plot copy number and major copy number #
        ##################################################
        cmd = "%s/plotCPandMCP.py -t %s/truth_cnv.tab -i %s/cnv.output.tsv "\
            "-o %s/plot.cnv.vs.truth.purity%s.png -r %s || exit 2" % \
            (binary_folder, self.simulation_data, self.main_output_dir, 
            self.main_output_dir, self.purity, os.path.join(ref_folder_path, "genome.dict"))
        PlotCPAndMCPJob = self.addTask(
            ("PlotCPAndMCP.%s"%self.sample_name).replace(".", "_"), cmd,
            dependencies=main_job)

        ####################
        # STEP 7: Clean....#
        ####################
        sys.stderr.write("Clean ...\n")
        cmd = "mv %s/*ratio.w*.csv.gz %s/../"%(self.main_output_dir,
            self.main_output_dir)
        clean_job = self.addTask(
            ("Clean.%s"%self.sample_name).replace(".", "_"),
            cmd, dependencies=PlotCPAndMCPJob)

    def print_result_to_console(self):
        """
        print
            .../infer.out.tsv
            .../infer.out.details.tsv

        :return:
        """
        result_file = open(os.path.join(self.main_output_dir, "infer.out.tsv"))
        sys.stderr.write("### Content of infer.out.tsv:\n")
        for line in result_file:
            sys.stderr.write(line)
        result_file.close()

        result_file = open(os.path.join(self.main_output_dir, "infer.out.details.tsv"))
        sys.stderr.write("### First 15 lines of infer.out.details.tsv:\n")
        counter = 0
        for line in result_file:
            sys.stderr.write(line)
            counter += 1
            if counter>=15:
                break

        result_file = open(os.path.join(self.result_file))
        sys.stderr.write("### Content of %s:\n"%self.result_file)
        for line in result_file:
            sys.stderr.write(line)
        result_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-y', '--simulation_type', required=True,
        help='which type of simulated sample, i.e. singleclone')
    parser.add_argument('-p', '--purity', default=0.7,
        help='which purity for this type of simulated sample')
    parser.add_argument('-c', '--nCores', type=int, default=8,
        help="the max number of CPUs to use. Default is 8.")
    parser.add_argument("-d", "--debug", type=int, default=0,
        help="Set debug value. Default 0 means no debug info output."
		"Anything >0 enables more debug output.")
    parser.add_argument("--segment_stddev_divider", type=float, default=10.0,
        help="artificially reduce segment noise level")
    parser.add_argument('--call_snps',  action='store_true', dest='call_snps',
        help="If toggled, SNPs will be called.")
    parser.add_argument('-o', '--output_dir', type=str, required=True,
        help="pyflow output directory. Accurity output is in $simulation_type/purity*/")
    parser.add_argument('-r', '--ref_version', type=str, default="hs37d5",
        help="The version of referece data used in Accurity")
    args = parser.parse_args()
    simulation_type = args.simulation_type
    env_dist = os.environ
    commit_id = env_dist.get('CI_BUILD_REF')[:8]
    wflow = SampleTestFlow(simulation_type=simulation_type, purity=args.purity,
        segment_stddev_divider=args.segment_stddev_divider,
        call_snps=args.call_snps, commit_id=commit_id,
        debug=args.debug, nCores=args.nCores,
        ref_version=args.ref_version)
    retval = wflow.run(mode="local", nCores=args.nCores,
        dataDirRoot=args.output_dir,
        isContinue='Auto', isForceContinue=True, retryMax=0)
    wflow.print_result_to_console()
    sys.exit(retval)
