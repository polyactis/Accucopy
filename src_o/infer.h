/*
 Author:
 Yu S. Huang, polyactis@gmail.com
 */

#pragma once
#ifndef __INFER_H
#define __INFER_H

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <gsl/gsl_cdf.h>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>

#include "BaseGADA.h"
#include "read_para.h"
#include "prob.h"

using namespace std;

inline bool larger(double a, double b) {
    return a > b;
}

class OnePeak {
   public:
    OnePeak() : no_of_windows(0), no_of_snps(0) {
        peak_center_int=-1;
        peak_height = -1.0;
        peak_index = -1;
        lower_bound_int=-1;
        upper_bound_int=-1;
        half_width_int=-1;
    }
    OnePeak(int peak_center_int, double peak_height, int peak_index, int lower_bound_int, int upper_bound_int,
            int half_width_int)
        : peak_center_int(peak_center_int),
          peak_height(peak_height),
          peak_index(peak_index),
          lower_bound_int(lower_bound_int),
          upper_bound_int(upper_bound_int),
          half_width_int(half_width_int),
          no_of_windows(0),
          no_of_snps(0) {}
    ~OnePeak() {
        segment_rc_ratio_vector.clear();
        segment_obj_vector.clear();
    }
    int peak_center_int;
    double peak_height;
    int peak_index;
    int lower_bound_int;
    int upper_bound_int;
    int half_width_int;
    int no_of_periods_since_1st_peak;
    int no_of_windows;
    int no_of_snps;
    int no_of_logOR_peaks;
    float snp_coverage;
    vector<int> segment_rc_ratio_vector;
    // segments that belong to this peak
    vector<OneSegment> segment_obj_vector;

    int get_no_of_segments() {
        return segment_obj_vector.size();
    }
    void ResetCounters() {
        no_of_snps = 0;
        no_of_logOR_peaks = 0;
    }
    //sort by peak_center_int
    bool operator < (const OnePeak& other) const
    {
        return (peak_center_int < other.peak_center_int);
    }
};

class OnePeriod {
   public:
    OnePeriod() {
        no_of_snps = 0;
        no_of_logOR_peaks = 0;
        no_of_segments = 0;
        no_of_rc_peaks = 0;
        no_of_windows = 0;
        logL = 0.0;
        auto_cor_value = 0.0;
        best_purity = -1;
        best_ploidy = -1;
        rc_ratio_int_of_cp_2 = -1;
        rc_ratio_int_of_cp_2_corrected = -1;
        purity_corrected = -1.0;
        ploidy_corrected = -1.0;
        no_of_peaks_for_logL = 0;
    }

    OnePeriod(int period_int, int lower_bound_int, int upper_bound_int)
        : period_int(period_int),
          lower_bound_int(lower_bound_int),
          upper_bound_int(upper_bound_int) {
        no_of_snps = 0;
        no_of_logOR_peaks = 0;
        no_of_segments = 0;
        no_of_rc_peaks = 0;
        no_of_windows = 0;
        logL = 0;
        auto_cor_value = 0.0;
        best_purity = -1;
        best_ploidy = -1;
        rc_ratio_int_of_cp_2 = -1;
        rc_ratio_int_of_cp_2_corrected = -1;
        purity_corrected = -1.0;
        ploidy_corrected = -1.0;
        best_logL_snp = -1E-99;
        no_of_peaks_for_logL;
    }

    int getWidth() {
        return upper_bound_int - lower_bound_int + 1;
    }

    void ResetCounters() {
        no_of_snps = 0;
        no_of_logOR_peaks = 0;
        no_of_segments = 0;
        no_of_rc_peaks = 0;
        no_of_windows = 0;
        logL = 0.0;
    }

    void ResetSNPCounters() {
        no_of_snps = 0;
        no_of_logOR_peaks = 0;
    }

    int period_int;
    int lower_bound_int;
    int upper_bound_int;
    double auto_cor_value;
    float rc_ratio_int_of_cp_2;
    float rc_ratio_int_of_cp_2_corrected;
    int period_index;
    int no_of_peaks_for_logL;
    int no_of_segments;
    int no_of_rc_peaks;
    int no_of_snps;
    int no_of_logOR_peaks;
    double logL;
    double purity_corrected;
    double ploidy_corrected;
    double best_purity, best_ploidy, logL_rc, adj_logL_rc, logL_rc_penalty,
        best_logL_snp, best_logL_snp_penalty,
        best_logL_snp_no_of_parameters;
    int first_peak_int;
    int width;
    int best_no_of_copy_nos_bf_1st_peak;
    int no_of_windows;
    OnePeak first_peak_obj;
    vector<double> logL_snp_vector;
    vector<double> snp_penalty_vector;
    vector<double> snp_no_of_parameters_vector;
    vector<double> purity_vector;
    vector<double> ploidy_vector;
    vector<int> no_of_copy_nos_bf_1st_peak_vector;
    vector<OnePeak> peak_obj_vector;
    //sort in ascending order by auto_cor_value
    bool operator < (const OnePeriod& other_period) const {
        return (auto_cor_value < other_period.auto_cor_value);
    }
    //sort in descending order by auto_cor_value
    bool operator > (const OnePeriod& other_period) const {
        return (auto_cor_value > other_period.auto_cor_value);
    }

};



class Infer {
   public:
    Infer(string configFilepath, string segment_data_input_path,
          string snp_data_input_path, string output_dir,
          float segment_stddev_divider,
          int snp_coverage_min, float snp_coverage_var_vs_mean_ratio,
          int no_of_peaks_for_logL,
          int debug, int auto_, string refdictFilepath, int force_period_id);
    ~Infer();
    int run();

   private:
    int getSNPDataFromFile(string inputFname);
    int getSegmentDataFromFile(string inputFname);
    int output_segment_ratio(int **noOfWindowsByRatioAndChr);
    int output_snp_logOR_by_segment();
    int output_snp_logOR_by_peak(vector<OnePeak> &peak_obj_vector);
    int output_rc_ratio_of_peaks(vector<OnePeak> &peak_obj_vector);
    void kernel_smoothing(double mean_value, double stddev, int sample_size,
                          vector<double> &vec_to_hold_data);
    void calculate_autocor();
    int infer_candidate_period_by_autocor(OnePeriod &period_obj);
    void calc_autocor_shift_diff(double* all_diff, double &left_x, double &right_x);
    vector<OnePeriod> infer_candidate_period_by_GADA(double* all_diff, double left_x, double right_x, int run_type);
    OnePeak find_first_peak_ab_init(int candidate_period_int);
    OnePeak find_first_peak_given_bounds(int candidate_period_int,
                                         int first_peak_lower_bound_int,
                                         int first_peak_upper_bound_int);
    OnePeak refine_first_peak(int candidate_period_int,
                              OnePeak &first_peak_obj);
    int chrStr_to_index(string);
    int findSNPsWithinSegment(OneSegment &oneSegment);
    vector<OnePeak> find_peaks(OnePeriod &period_obj, OnePeak &first_peak_obj);
    int output_peak_bounds(vector<OnePeak> &peak_obj_vector);

    OnePeriod infer_best_period_by_logL(vector<OnePeriod> &candidate_period_vec);

    int output_logL(OnePeriod &best_period_obj,
                   vector<OnePeriod> &period_obj_vector);
    double calc_one_peak_logL_rc(float peak_center_float, OnePeak &peak_obj,
                                 OnePeriod &period_obj, double &adj_logL);
    void calc_purity_ploidy_from_period_and_cp_no_two(
        int cp_no_two_rc_ratio_int, int period_int, double &purity,
        double &ploidy);

    double infer_no_of_copy_nos_bf_1st_peak_for_one_period_by_logL_snp(
        OnePeriod &candidate_period);
    double getReadDepthFromRegCoeffFile(string inputFname);

    int refine_peak_center(OnePeak &peak_obj, vector<int> segment_rc_ratio_vector,
                                          int candidate_period_int,
                                          int first_peak_center_int);

    vector<double> call_subclone_peaks(double *, int);
    double output_copy_number_segments(OnePeriod &best_period_obj,
                                    vector<OnePeak> &peak_obj_vector);

    string _configFilepath;
    string _refdictFilepath;
    string _segment_data_input_path;
    string _snp_data_input_path;
    string _output_dir;
    float _segment_stddev_divider;
    int _no_of_peaks_for_logL;
    float _snp_maf_stddev_divider;
    //parameter used to adjust snp expected MAF, adjust_maf_expect(). Being same to 6, used in selecting heterozygous SNPs is no good.
    int _snp_coverage_min;
    float _snp_coverage_var_vs_mean_ratio;
    int _debug;
    int _auto;
    int _returnCode;
    int _force_period_id; // Force to select whichever period to use, not via logL
    int _min_no_of_snps_for_1_peak; // min #SNPs for one peak to be included in SNP logL calc
    int _period_min_segment_len; // minimal length of a GADA segment in searching for candidate period


    Config _config;
    RefDictInfo ref_dict_info;


    Prob _probInstance;
    vector<vector<OneSNP> > _SNPs;                   // indexed by chromosomes
    vector<vector<OneSegment> > _rc_ratio_segments;  // vector of segments at
                                                     // each
    // rc_ratio (high-resolution)
    vector<OnePeriod> _periodObjVector;
    // probability density function (in the number of windows) at any given ratio_int
    vector<double> _ratio_int_pdf_vec;
    double _cor_array[kPeriodMax + 1];

    unsigned long _genome_len_cnv_all;
    unsigned long _genome_len_clonal;
    int _period_discover_run_type;
    float _ploidy_cnv_all;
    float _ploidy_clonal;
    int _total_no_of_snps;
    int _total_no_of_snps_used;
    int _total_no_of_segments;
    int _total_no_of_segments_used;
    // maximum peak half width, consistent across different periods.
    int _max_peak_half_width;
    OnePeriod _period_obj_from_autocor;  // period_int, lower bound, upper
    // bound
    OnePeriod _period_obj_from_logL;
    OnePeak _first_peak_obj;  // best_first_peak,
                              // half_interval_width
                              // for 2
                              // deviations
    //  int _num_peak_less_one_half; 20161206
    //  const int _one_half = 1.5; 20161206
    ofstream _infer_outf;
    ofstream _infer_details_outf;

    ofstream _rc_logL_outf;
    ofstream _snp_logL_outf;
    ofstream _sub_outf;
    ofstream _sub_peak_outf;

    std::ofstream rc_ratio_by_chr_out_file;

    double _pool_hist[MAX_RATIO_HIGH_RES + 1];
    int _half_period_int;

    int _valley;
};
#endif
