/*
 Author:
 Yu S. Huang, polyactis@gmail.com
 */
#pragma once
#ifndef __READ_PARA_H
#define __READ_PARA_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <sys/stat.h>
#include <vector>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/copy.hpp>
#include "format.h"
#include <regex>

using namespace std;

#define TabMACRO << "\t" <<
#define FourSpaceMACRO << "    " <<
#define NewLineMACRO << "\n"

const int kMaxInputLineLength = 200;
// const int MAX_LEN = 250000000;
const int GC_RATIO_RESOLUTION = 1000;
// const int MAX_GC_TOTAL = 500;
const int MAX_FRAGMENT_LENGTH = 500;
const int MAX_READ_COUNT = 300;  // 20170124 was 200 before

//	const int
// LARGE=GC_RATIO_RESOLUTION*MAX_GC_TOTAL/5*MAX_GC_TOTAL/5*MAX_READ_COUNT;
//	const int COEF1=MAX_GC_TOTAL/5*MAX_GC_TOTAL/5*MAX_READ_COUNT;
//	const int COEF2=MAX_GC_TOTAL/5*MAX_READ_COUNT;
// const int COEF3 = MAX_READ_COUNT;
const int LARGE = GC_RATIO_RESOLUTION * MAX_READ_COUNT;

float const MIN_PLOIDY = 1.0;
float const MAX_PLOIDY = 4.0;
float const SHOULDER_RATIO = 0.8;
float const SHOULDER_RATIO_LEFT = 0.9;
float const kGaussianDensityFrontScalar =
    0.3989;  // 1/sqrt(2*pi), Gaussian density function constant
int const RESOLUTION = 1000;
float const FRESOLUTION = 1.0 * RESOLUTION;
int const CLONE_HALF_RESOLUTION = (int)(0.02 * RESOLUTION);  // 200
int const MAX_NUM_RECORD = 1000000;
int const kPeriodMin = (int)(0.1 * RESOLUTION);  // 100
int const kPeriodMax = 1.0 * RESOLUTION;         // 1000
int const kPeriodHalfWidthMax = 20;
// int const WIDTH=80;
int const MAX_RATIO = 3;
int const MAX_RATIO_HIGH_RES = (int)(MAX_RATIO * RESOLUTION);  // 3000
int const MAX_RATIO_RANGE = 100;  // maximum segment enrichment ratio allowed
int const MAX_RATIO_RANGE_HIGH_RES =
    (int)(MAX_RATIO_RANGE * RESOLUTION);  // 100K
int const kFirstPeakMin = 10;
int const kFirstPeakMax = 1.05 * RESOLUTION;            // 1050
int const kPeakHalfWidthMax = (int)(0.2 * RESOLUTION);  // 200
float const kPeakHeightMin = 1e2;
float const DEV1 =
    0.606;  // normal density with one standard deviation away from the mean
float const TAIL =
    0.2;  // normal density with one standard deviation away from the mean
float const DEV2 =
    0.135;  // normal density with two standard deviation away from the mean
// const string chromosomeNameArray[] = {
//     "chr1",  "chr2",  "chr3",  "chr4",  "chr5",  "chr6",  "chr7",  "chr8",
//     "chr9",  "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
//     "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX",  "chrY"};
// const string chromosomeNoArray[] = {
//     "1",  "2",  "3",  "4",  "5",  "6",  "7",  "8",  "9",  "10", "11", "12",
//     "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X",  "Y"};
// const int CHR_SIZE[] = {251, 244, 199, 193, 182, 172, 160, 147, 143,
//                         137, 136, 135, 116, 108, 104, 91,  83,  79,
//                         60,  64,  49,  52,  156, 57};  // in million bp, 24
//                                                        // chromosomes
// const int NUM_AUTO_CHR = 22;
// const int NUM_CHR = 23;

const int MAX_NUM_OF_COR_TO_SUM =
    (int)(0.5 * 0.1 * 10 *
          RESOLUTION);  // 0.5: 0.5 one interval; 0.1 peak half width; 3 peaks


class RefDictInfo
{
   public:
    RefDictInfo();
    RefDictInfo(int NUM_AUTO_CHR,vector <string> chromosomeNameArray,vector <long> CHR_SIZE);
    int NUM_AUTO_CHR;
    vector <string> chromosomeNameArray;
    vector <long> CHR_SIZE; 

};


class OneSegmentSNPs
{
   public:
    OneSegmentSNPs();
    OneSegmentSNPs(int no_of_snps, float coverage, vector<double>& logOR_list);
    int no_of_snps;
    float coverage;
    vector<double> logOR_list;
};

class OneSegment
{
   public:
    OneSegment();
    OneSegment(int chr_index, int start_pos, int end_pos, float rc_ratio,
               double stddev, int no_of_windows);
    int get_rc_ratio_high_res();
    int chr_index;
    int start_pos;
    int end_pos;
    float rc_ratio;
    double stddev;
    int no_of_windows;
    int no_of_snps;
    OneSegmentSNPs oneSegmentSNPs;
};

class OneSNP
{
   public:
    OneSNP(int chr_index, int position, double logOR, int coverage);
    int chr_index;
    int position;
    double logOR;
    int coverage;
};

class Config
{
   public:
    bool has_valid_para;
    string ref;
    int readlength;
    int window;
    string path_to_ref;
    string path;
};

RefDictInfo read_dict_from_file(string _refdictFilepath);

Config read_para(string file_conf);
inline bool isfile (const std::string& name) {
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
}


void calculate_median_mad(double *input_array, long start_index, long stop_index,
                                 float &median_value, float &mad_value);
void calculate_robust_mean_stddev(double *input_array, long start_index, long stop_index, int percent_to_exclude,
                                  float &mean_ref, float &stddev_ref);
void calculate_robust_mean_stddev(vector<float> float_vector, int percent_to_exclude,
                                  float &mean_ref, float &stddev_ref, double &squared_sum, int &sample_size);
vector<std::string> string_split(std::string str,std::string sep);
#endif
