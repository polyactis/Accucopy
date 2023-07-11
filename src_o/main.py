#!/usr/bin/env python2
# -*- coding: future_fstrings -*-
"""
 Author:
 Yu S. Huang, polyactis@gmail.com
 Xinping Fan, 897488736@qq.com

"""
from argparse import ArgumentParser
import logging
import os
import sys
from subprocess import Popen
import shutil,re 
from datetime import datetime, timedelta
from pyflow import WorkflowRunner

class MainFlow(WorkflowRunner):
    def __init__(self, configure_filepath=None, tumor_bam=None,
        normal_bam=None, output_dir=None,
        snp_output_dir=None,
        segment_stddev_divider=20, snp_coverage_min=2,
        snp_coverage_var_vs_mean_ratio=10.0,
        no_of_autosomes=22,
        clean=False,
        step=0, debug=False, auto=1,
        max_no_of_peaks_for_logL=3, nCores=4, force_period_id=0, **keywords):
        self.configure_filepath = configure_filepath
        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam
        self.output_dir = output_dir
        self.snp_output_dir = snp_output_dir
        self.segment_stddev_divider = segment_stddev_divider
        self.snp_coverage_min = snp_coverage_min
        self.snp_coverage_var_vs_mean_ratio = snp_coverage_var_vs_mean_ratio
        self.no_of_autosomes = no_of_autosomes
        self.clean = clean
        self.step = step
        self.debug = debug
        self.auto = auto
        self.max_no_of_peaks_for_logL = max_no_of_peaks_for_logL
        self.nCores = nCores
        self.strelka_cores = max(1, self.nCores-2)
        self.force_period_id = force_period_id

        if not os.path.isdir(self.output_dir):
            os.mkdir(self.output_dir)
        #window_size will be read from the config file.
        self.samtools_path = ""
        self.strelka_path = ""
        self.read_len = None
        self.window_size = None
        self.ref_folder_path = ""
        self.binary_folder = None

        # normalization parameters
        self.max_coverage = 300
        self.smooth_window_half_size = 2

        # segmentation parameters
        self.min_segment_len = 50
        self.t_score_threshold = 30

        # self.chromosomeNames = ["chr1", "chr2", "chr3", "chr4", "chr5",
        #   "chr6", "chr7",
        #   "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
        #   "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
        #   "chr22", "chrX", "chrY"]
        # self.NUM_AUTO_CHR = 22

        self.chromosomeNames = None
        self.NUM_AUTO_CHR = None

        if self.snp_output_dir:
            self.strelka_output_dir = self.snp_output_dir
        else:
            self.strelka_output_dir = os.path.join(
                os.path.split(self.tumor_bam)[0], "strelka_snp")
        self.two_sample_snp_file = os.path.join(self.strelka_output_dir, 
            "results/variants/variants.vcf.gz")
        self.het_snp_filepath = os.path.join(self.output_dir,
            "het_snp.tsv.gz")
        self.het_snp_normal_filepath = os.path.join(self.output_dir,
            "het_snp.normal.tsv.gz")
        self.segment_data_filepath = os.path.join(self.output_dir,
            "all_segments.tsv.gz")
        self.infer_status_out_path = os.path.join(self.output_dir,
            "infer.status.txt")

        if os.path.isdir(self.output_dir):
            pyflowdir = os.path.join(self.output_dir, "pyflow.data")
            if os.path.isdir(pyflowdir):
                sys.stderr.write("Deleting pyflow folder ...")
                shutil.rmtree(pyflowdir)
                sys.stderr.write("Done.\n")

    def runShellCommand(self, cmdLine):
        sys.stderr.write("Running %s ...\n" % cmdLine)
        p = Popen(cmdLine, shell=True)
        p.wait()
        if p.returncode != 0:
            sys.stderr.write("Error encountered. Stop now.\n")
            self.final_log()
            sys.exit(p.returncode)
        return p.returncode

    def readConfigureFile(self, inputFname):
        # read in configure file
        if not os.path.isfile(inputFname):
            logging.error("Configure file %s does not exist!"%(inputFname))
            sys.exit(2)
        with open(inputFname, 'r') as f:
            lines = f.readlines()
            self.read_len = lines[0].strip().split("\t")[1]
            self.window_size = lines[1].strip().split("\t")[1]
            self.ref_folder_path = lines[2].strip().split("\t")[1]
            self.samtools_path = lines[3].strip().split("\t")[1]
            self.strelka_path = lines[4].strip().split("\t")[1]
            self.binary_folder = lines[5].strip().split("\t")[1]
    
    def readDictFile(self):
        ref_dict_filename = os.path.join(self.ref_folder_path, "genome.dict")
        chromosomeNames = []
        # dict file's chromosome id should be well sorted
        pattern = re.compile(r'^chr\d+$')
        with open(ref_dict_filename, 'r') as f:
            for line in f:
                if line.startswith('@SQ'):
                    records = line.split('\t')
                    chr_name = records[1].split(':')[1]
                    m = pattern.findall(chr_name)
                    if m:
                        chromosomeNames.append(m[0])
        self.chromosomeNames = chromosomeNames
        # chromosomeNames name did not contain sexual chromosome
        self.NUM_AUTO_CHR = len(chromosomeNames)

    def workflow(self):
        sys.stderr.write("Step=%s\n" % self.step)
        self.startTimeList = [datetime.now()]
        ########################################
        # STEP 0: preparation				   #
        ########################################
        # check whether index files exist
        status_string = "step 0: preparation (mkdir, index bam if bai is missing)."+\
            " start time: %s\n" % self.startTimeList[-1]
        sys.stderr.write(status_string)

        if os.path.isfile(self.output_dir):
            sys.stderr.write("Output dir %s is a file. Remove it.\n"%self.output_dir)
            os.remove(self.output_dir)
            os.mkdir(self.output_dir)
        elif os.path.isdir(self.output_dir):
            if self.clean:
                sys.stderr.write("Clean flag is on. Force remove %s and mkdir it.\n"%\
                    self.output_dir)
                shutil.rmtree(self.output_dir)
                os.mkdir(self.output_dir)
        else:
            os.mkdir(self.output_dir)

        model_select_output_dir = os.path.join(self.output_dir, "model_selection_log")
        if os.path.isfile(self.output_dir):
            sys.stderr.write("Output dir %s is a file. Remove it.\n" %\
                model_select_output_dir)
            os.remove(model_select_output_dir)
            os.mkdir(model_select_output_dir)
        elif os.path.isdir(model_select_output_dir):
            sys.stderr.write("%s exists, remove and mkdir it.\n" % \
                model_select_output_dir)
            shutil.rmtree(model_select_output_dir)
            os.mkdir(model_select_output_dir)
        else:
            os.mkdir(model_select_output_dir)

        #clean pyflow folder, otherwise it'll conflict with the next pyflow run.
        #pyflow_dir = os.path.join(self.output_dir, "pyflow.data")
        #if os.path.isdir(pyflow_dir):
        #	sys.stderr.write("Clean out pyflow folder, %s, to avoid conflicts.\n" % 
        #       pyflow_dir)
        #	shutil.rmtree(pyflow_dir)

        tumor_idx = self.tumor_bam + ".bai"
        normal_idx = self.normal_bam + ".bai"
        if not os.path.isfile(tumor_idx):
            cmd = self.samtools_path + " index " + self.tumor_bam
        else:
            cmd = None
        indexTumorBamJob = self.addTask("indexTumorBam", cmd)

        if not os.path.isfile(normal_idx):
            cmd = self.samtools_path + " index " + self.normal_bam
        else:
            cmd = None
        indexNormalBamJob = self.addTask("indexNormalBam", cmd)



        ############################################################
        # STEP 1: SNP calling                                      #
        ############################################################
        if self.step <= 1:
            self.startTimeList.append(datetime.now())
            status_string = "Last step time span: %s\n" % \
                (self.startTimeList[-1] - self.startTimeList[-2])
            status_string += "step 1: call SNPs.\n\tStart time: %s\n"%\
                self.startTimeList[-1]
            sys.stderr.write(status_string)
            oneThousandSNPFilepath = os.path.join(self.ref_folder_path, "snp_sites.gz")
            #input: tumor bam
            #output: self.vcf_tumor_file_path
            cmd = f"{self.strelka_path}/bin/configureStrelkaGermlineWorkflow.py "\
                f"--bam {self.normal_bam} "\
                f"--bam {self.tumor_bam} "\
                f"--ref {os.path.join(self.ref_folder_path, 'genome.fa')} "\
                f"--callRegions {oneThousandSNPFilepath} --runDir {self.strelka_output_dir}"
            strelka_prepare_job = self.addTask("strelka_prepare", cmd,
                dependencies=[indexNormalBamJob, indexTumorBamJob])
            cmd = f"{self.strelka_output_dir}/runWorkflow.py -m local -j {self.strelka_cores}"
            strelka_call_snp_job = self.addTask("strelka_call_snp", cmd,
                nCores=self.strelka_cores, dependencies=[strelka_prepare_job])
        else:
            strelka_prepare_job = self.addTask("strelka_prepare",
                dependencies=[indexNormalBamJob, indexTumorBamJob])
            strelka_call_snp_job = self.addTask("strelka_call_snp",
                dependencies=[strelka_prepare_job])



        ############################################################
        # STEP 2: GC normalization								 #
        ############################################################
        normalize_jobs = []
        normalize_output_file_ls = []
        if self.step <= 2:
            self.startTimeList.append(datetime.now())
            status_string = "Last step time span: %s\n" % \
                (self.startTimeList[-1] - self.startTimeList[-2])
            status_string += "step 2: GC normalization.\n\tstart time: %s\n" % \
                self.startTimeList[-1]
            sys.stderr.write(status_string)
            #input: tumor.bam and normal.bam
            #output: reg.in.txt, reg.out.txt in both tumor (tumor.reg.in.txt) and normal
            #output: tumor/"%s.ratio.w%s.csv.gz"%(chromosome, self.window_size)
            #output (GC-normalize adj factors): tumor/tumor.cov.adj.factor.txt,
            #   tumor/normal.cov.adj.factor.txt
            reg_input_base_filename = "reg.in.txt"
            reg_output_base_filename = "reg.out.txt"
            cmd = f'{os.path.join(self.binary_folder, "maestre")} normalize '\
                f'-t {self.tumor_bam} -n {self.normal_bam} '\
                f'--genome_dict_path {os.path.join(self.ref_folder_path, "genome.dict")} '\
                f'-w {self.window_size} -l {self.read_len} '\
                f'--smooth_window_half_size {self.smooth_window_half_size} '\
                f'--max_coverage {self.max_coverage} --debug {self.debug} '\
                f'--no_of_autosomes {self.no_of_autosomes} '\
                f'-o {self.output_dir} 2>&1 | tee -a {self.infer_status_out_path}'
            normalize_jobs.append(self.addTask("normalize", cmd,
                dependencies=[indexTumorBamJob, indexNormalBamJob]))
            #add a gzip job
            #cmd = "gzip %s/tumor.%s %s/normal.%s"%(self.output_dir, 
            #   reg_input_base_filename, self.output_dir, reg_input_base_filename)
            #self.addTask("gzip_regression_input", cmd, dependencies=normalize_jobs[-1])
        else:
            normalize_jobs.append(self.addTask("normalize", \
                dependencies=[indexTumorBamJob, indexTumorBamJob]))

        for chr_index in range(self.NUM_AUTO_CHR):
            chromosome = self.chromosomeNames[chr_index]
            normalized_output_file_path = os.path.join(self.output_dir, \
                "%s.ratio.w%s.csv.gz"%(chromosome, self.window_size))
            normalize_output_file_ls.append(normalized_output_file_path)

        if self.debug:
            #plot the coverage plot between tumor and normal by adjust
            cmd = f'{os.path.join(self.binary_folder, "plot_coverage_after_normalization.py")} '\
                f'-i {os.path.join(self.output_dir, "chr22.ratio.w%s.csv.gz"%self.window_size)} '\
                f'-o {os.path.join(self.output_dir, "plot.tumor_vs_normal.chr22.png")}'
            plot_coverage_job = self.addTask("plot_tumor_normal_coverage", cmd,
                dependencies=normalize_jobs)

            """
            #plot GC normalization png
            cmd = "%s -r %s -a %s -o %s" % (
                os.path.join(self.binary_folder, "plot_GC_normalization.py"),
                os.path.join(self.output_dir, "reg.in.txt"),
                os.path.join(self.output_dir, "cov.adj.factor.txt"),
                os.path.join(self.output_dir, "plot.gc.adj.png"))
            plot_gc_adjust_tumor_job = self.addTask("plot_gc_correction", 
                cmd, dependencies=normalize_jobs)
            """

        ############################################################
        # STEP 3: select heterozygous SNPs                 		#
        ############################################################

        if self.step <= 3:
            self.startTimeList.append(datetime.now())
            status_string = "Last step time span: %s\n" % \
                (self.startTimeList[-1] - self.startTimeList[-2])
            status_string += "step 3: select heterozygous SNPs.\n\t"\
                "Start time: %s\n" % \
                self.startTimeList[-1]
            sys.stderr.write(status_string)
            #input: self.vcf_tumor_file_path, self.vcf_normal_file_path
            #output: het_snp
            cmd = f"{os.path.join(self.binary_folder, 'maestre')} "\
                f"select_het_snp -s {self.two_sample_snp_file} -m 2 -x 200 "\
                f"--debug 0 -o {self.het_snp_filepath} 2>&1 "\
                f"| tee -a {self.infer_status_out_path}"
            call_het_snps_tumor_job = self.addTask("call_het_snps_tumor", cmd,
                dependencies=[strelka_call_snp_job])

            #cmd = "%s %s 15 | gzip > %s" % \
            #   (os.path.join(self.binary_folder, "snp_calling"), 
            #   self.vcf_normal_file_path, self.het_snp_normal_filepath)
            #call_het_snps_normal_job = self.addTask("call_het_snps_normal", 
            # cmd, dependencies=[call_snps_normal_job])
        else:
            call_het_snps_tumor_job = self.addTask("call_het_snps_tumor",
                        dependencies=[strelka_call_snp_job])
            #call_het_snps_normal_job = self.addTask("call_het_snps_normal",
            #   dependencies=[call_snps_normal_job])

        ############################################################
        # STEP 4: Segmentation									   #
        ############################################################
        if self.step <= 4:
            self.startTimeList.append(datetime.now())
            status_string = "Last step time span: %s\n" % \
                (self.startTimeList[-1] - self.startTimeList[-2])
            status_string += "step 4: Segmentation \n\tstart time: %s\n" % \
                self.startTimeList[-1]
            sys.stderr.write(status_string)

            segment_jobs = []
            segment_out_ls = []
            for chr_index in range(self.NUM_AUTO_CHR):
                chromosome = self.chromosomeNames[chr_index]
                segment_out_path = os.path.join(self.output_dir, \
                    f"{chromosome}.segments.M{self.min_segment_len}."\
                    f"T{self.t_score_threshold}.tsv")
                segment_out_ls.append(segment_out_path)
                cmd = f'{os.path.join(self.binary_folder, "GADA")} '\
                    f'--chromosome_id {chromosome} --window_size {self.window_size} '\
                    f'-M {self.min_segment_len} -T {self.t_score_threshold} '\
                    f'-i {normalize_output_file_ls[chr_index]} -o {segment_out_path} '\
                    f'2>&1 | tee -a { self.infer_status_out_path}'
                segment_jobs.append(self.addTask("segment_%s"%chromosome, cmd, 
                    dependencies=normalize_jobs))
            cmd = f"cat {' '.join(segment_out_ls)} | gzip > {self.segment_data_filepath}"
            reduce_all_segments_job = self.addTask("reduce_all_segments", cmd,
                dependencies=segment_jobs)
            #remove all individual segments files
            self.addTask("rm_individual_seg_files", "rm %s" % \
                " ".join(segment_out_ls), dependencies=reduce_all_segments_job)
        else:
            reduce_all_segments_job = self.addTask("reduce_all_segments")

        ############################################################
        # STEP 5: Infer purity, ploidy, etc.
        ############################################################
        if self.step <= 5:
            self.startTimeList.append(datetime.now())
            status_string = "Last step time span: %s\n" % \
                (self.startTimeList[-1] - self.startTimeList[-2])
            status_string += f"step 5: Infer tumor purity and ploidy.\n\t"\
                f"start time: {self.startTimeList[-1]}\n"
                
            sys.stderr.write(status_string)

            #input: self.segment_data_filepath (all_segments),
            #   self.het_snp_filepath (het_snp)
            #input: reg_coeff (to get depth of the of tumor bam)
            #output: infer.out.tsv, infer.out.details.tsv,
            #   rc_ratio_window_count_smoothed.tsv, peak_bounds.tsv
            #output: auto.tsv, cnv.output.tsv
            
            cmd = f"{os.path.join(self.binary_folder, 'infer')} "\
                f"{self.configure_filepath} {self.segment_data_filepath} "\
                f"{self.het_snp_filepath} {self.output_dir} "\
                f"{self.segment_stddev_divider} {self.snp_coverage_min} "\
                f"{self.snp_coverage_var_vs_mean_ratio} "\
                f"{self.max_no_of_peaks_for_logL} {self.debug} {self.auto} "\
                f"{os.path.join(self.ref_folder_path, 'genome.dict')} "\
                f"{self.force_period_id} "\
                f" 2>&1 | tee -a {self.infer_status_out_path}"
            infer_job = self.addTask("infer", cmd, 
                dependencies=[reduce_all_segments_job, call_het_snps_tumor_job])
            if self.debug:
                self.addTask("gzip_rc_ratio_no_of_windows_by_chr",
                    f"gzip {self.output_dir}/rc_ratio_no_of_windows_by_chr.tsv",
                    dependencies=infer_job)
        else:
            infer_job = self.addTask("infer")

        ############################################################
        # STEP 6: Make plots.
        ############################################################
        if self.step <= 6:
            self.startTimeList.append(datetime.now())
            status_string = "Last step time span: %s\n" % (
                self.startTimeList[-1] - self.startTimeList[-2])
            status_string += "step 6: Make plots.\n\tstart time: %s\n" % (
                self.startTimeList[-1])
            sys.stderr.write(status_string)
            # input: $(output_dir)/infer.out.tsv, infer.out.details.tsv,
            #   rc_ratio_window_count_smoothed.tsv, peak_bounds.tsv
            # output: plot.tre.jpg
            inferOutPath = os.path.join(self.output_dir, "infer.out.tsv")
            inferOutDetailsPath =os.path.join(self.output_dir,
                "infer.out.details.tsv")
            rcRatioSmoothedPath =os.path.join(self.output_dir,
                "rc_ratio_window_count_smoothed.tsv")
            peakBoundsPath =os.path.join(self.output_dir, "peak_bounds.tsv")

            # cmd = "%s %s %s %s %s"%(os.path.join(self.path, "plot.cnv.R"),
            #   self.configure_filepath, \
            #	tumor_samplename, self.output_dir,
            #   os.path.join(self.output_dir, "plot.cnv.jpg"))
            cmd = f"{os.path.join(self.binary_folder, 'plotCPandMCP.py')} "\
                f"-i {os.path.join(self.output_dir, 'cnv.output.tsv')} "\
                f"-r {os.path.join(self.ref_folder_path, 'genome.dict')} "\
                f"--no_of_autosomes {self.no_of_autosomes} "\
                f"-o {os.path.join(self.output_dir, 'plot.cnv.png')} "

            plot_cnv_job = self.addTask("plot_cnv", cmd, dependencies=infer_job)

            if self.debug:
                #plot the auto_cor diff program
                autocorPath =os.path.join(self.output_dir, "auto.tsv")
                cmd = f"{os.path.join(self.binary_folder, 'plot_autocor_diff.py')} "\
                    f"-i {os.path.join(self.output_dir, 'candidate.period.GADA.in.tsv')} "\
                    f"-o {os.path.join(self.output_dir, 'plot.tre.autocor.png')} "\
                    f"-s {os.path.join(self.output_dir, 'candidate.period.GADA.out.tsv')} "\
                    f"-a {autocorPath}"
                plot_autocor_diff_job = self.addTask("plot_autocor_diff", cmd, \
                    dependencies=infer_job)

                #tre job may fail. so run it after plot_autocor_diff_job
                #  (which will not fail if input files do not exist)
                cmd = f"{os.path.join(self.binary_folder, 'plot_tre.py')} "\
                    f"-i {rcRatioSmoothedPath} -p {peakBoundsPath} "\
                    f"-o {os.path.join(self.output_dir, 'plot.tre.png')}"
                plot_tre_job = self.addTask("plot_tre", cmd, dependencies=infer_job)

                #plot model selection result.
                hdf5_file = os.path.join(self.output_dir, \
                    "model_selection_log", "model_selection.h5")
                plot_output = os.path.join(self.output_dir, "model_selection_log")
                cmd = f"{os.path.join(self.binary_folder, 'plot_model_select_result.py')} "\
                    f"-f {hdf5_file} -o {plot_output}"
                plot_model_select_job = self.addTask("plot_model_select", cmd,
                    dependencies=infer_job)

                # input: $(output_dir)/infer.out.tsv, infer.out.details.tsv, auto.tsv
                # output: plot.tre.autocor.jpg
                #cmd = "%s %s %s %s %s" % (os.path.join(self.binary_folder, 
                #   "plot.tre.autocor.R"), inferOutPath,
                #	inferOutDetailsPath, autocorPath,
                #	os.path.join(self.output_dir, "plot.tre.autocor.jpg"))
                #plot_tre_autocor_job = self.addTask("plotTREAutocor", cmd,
                #   dependencies=infer_job)

        self.final_log()

    def final_log(self):
        self.startTimeList.append(datetime.now())
        statusString = "Last step time span: %s\n" % \
            (self.startTimeList[-1] - self.startTimeList[-2])
        statusString += "End time: %s\n" % self.startTimeList[-1]
        sys.stderr.write(statusString)


if __name__ == '__main__':
    ap = ArgumentParser(description='Please go to '
        'https://www.yfish.org/display/PUB/Accucopy for help or '
        'email polyactis@gmail.com.')
    ap.add_argument("-v", "--version", action="version", 
        version="32acfd1e-debug")
    ap.add_argument("-c", "--configure_filepath", type=str, required=True,
        help="the path to the configure file.")
    ap.add_argument("-t", "--tumor_bam", type=str, required=True,
        help="the path to the tumor bam file. "
        "If the bam is not indexed, an index file will be generated")
    ap.add_argument("-n", "--normal_bam", type=str, required=True,
        help="the path to the normal bam file. "
        "If the bam is not indexed, an index file will be generated")
    ap.add_argument("-o", "--output_dir", type=str, required=True,
        help="the output directory path.")
    ap.add_argument("--snp_output_dir", type=str, default=None,
        help="the directory to hold the SNP calling output. "
        "Default is the same folder as the bam file.")
    ap.add_argument("--clean", action='store_true',
        help="Toggle to remove the existing output folders and files.")
    ap.add_argument("--segment_stddev_divider", type=float, default=20.0,
        help="A factor that reduces the segment noise level. "
        "The default value (%(default)s) is recommended.")
    ap.add_argument("--no_of_autosomes", type=int, default=22,
        help="The number of autosome chromosomes for the species. "\
            "Sex chromosomes are excluded. "
        "Default is %(default)s.")
    ap.add_argument("--snp_coverage_min", type=int, default=2,
        help="The minimum SNP coverage in adjusting the expected SNP MAF. "
        "Default is %(default)s.")
    ap.add_argument("--snp_coverage_var_vs_mean_ratio", type=float, default=10.0,
        help="Instead of using the observed SNP coverage variance (not consistent), "
        "use coverage_mean X this-parameter as the variance for "
        "the negative binomial model "
        "which is used in adjusting the expected SNP MAF. "
        "Default is %(default)s.")
    ap.add_argument("--max_no_of_peaks_for_logL", type=int, default=3,
        help="the maximum number of peaks used in the log likelihood calculation. "
        "The final logL is average over the number of peaks used. "
        "Default is %(default)s")
    ap.add_argument("--nCores", type=int, default=8, 
        help="the max number of CPUs to use in parallel. "
            "Increase the number if you have many cores. Default is %(default)s.")
    ap.add_argument("-s", "--step", type=int, default=0,
        help='0: start from the very beginning (Default). '\
        '1: obtain the read positions and the major allele fractions. '\
        '2: normalization. '\
        '3: segmentation. '\
        '4: infer purity and ploidy only.')
    ap.add_argument("-l", "--lam", type=int, default=4,
        help="lambda for the segmentation algorithm. Default is %(default)s.")
    ap.add_argument("-d", "--debug", type=int, default=0,
        help="Set debug value. Default 0 means no debug info output."
            "Anything >0 enables more debug output.")
    ap.add_argument("--auto", type=int, default=1,
        help="The integer-valued argument that decides which method to use "
        "to detect the period in the read-count ratio histogram. "
        "0: the simple auto-correlation method. "
        "1: a GADA-based algorithm (recommended). Default is %(default)s.")
    ap.add_argument("--force_period_id", type=int, default=0,
        help="The non-negative integer-valued argument that force the program "
        "to use which period, instead of relying on max likelihood. "
        "0: use the one inferred automatically by the program; "
        "1: force-use the 1st period; "
        "2: force-use the 2nd period, etc. Max is 3. Default is %(default)s")
    args = ap.parse_args()
    if (args.force_period_id < 0):
        msg = (f"Argument `period` is non-negative integer-value, but you have "
               f"specified {args.force_period_id}. O will be used instead.")
        sys.stderr.write(msg)
    wflow = MainFlow(args.configure_filepath, args.tumor_bam, args.normal_bam,
        output_dir=args.output_dir,
        snp_output_dir=args.snp_output_dir,
        segment_stddev_divider=args.segment_stddev_divider,
        snp_coverage_min=args.snp_coverage_min,
        snp_coverage_var_vs_mean_ratio=args.snp_coverage_var_vs_mean_ratio,
        no_of_autosomes=args.no_of_autosomes,
        clean=args.clean, step=args.step, debug=args.debug, auto=args.auto,
        max_no_of_peaks_for_logL=args.max_no_of_peaks_for_logL,
        nCores=args.nCores, force_period_id=args.force_period_id)
    wflow.readConfigureFile(args.configure_filepath)
    wflow.readDictFile()
    retval = wflow.run(mode="local", nCores=args.nCores,
        dataDirRoot=args.output_dir, isContinue='Auto',
        isForceContinue=True, retryMax=0)
    sys.exit(retval)
