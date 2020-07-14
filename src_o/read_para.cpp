#include "read_para.h"

using namespace std;

OneSegment::OneSegment()
    : chr_index(-1),
      start_pos(-1),
      end_pos(-1),
      rc_ratio(-1),
      stddev(-1),
      no_of_windows(-1)
{
}

OneSegment::OneSegment(int chr_index, int start_pos, int end_pos,
                       float rc_ratio, double stddev, int no_of_windows)
    : chr_index(chr_index),
      start_pos(start_pos),
      end_pos(end_pos),
      rc_ratio(rc_ratio),
      stddev(stddev),
      no_of_windows(no_of_windows)
{
}

int OneSegment::get_rc_ratio_high_res()
{
    return int(rc_ratio * RESOLUTION);
}

OneSNP::OneSNP(int chr_index, int position, double logOR, int coverage)
    : chr_index(chr_index), position(position), logOR(logOR), coverage(coverage)
{
}

OneSegmentSNPs::OneSegmentSNPs()
    : no_of_snps(-1),
      coverage(-1),
      logOR_list(){

}

//logOR_list will be moved
OneSegmentSNPs::OneSegmentSNPs(int no_of_snps, float coverage, vector<double>& logOR_list)
    : no_of_snps(no_of_snps),
      coverage(coverage),
      logOR_list(std::move(logOR_list)){
}

RefDictInfo::RefDictInfo()
    : NUM_AUTO_CHR(-1),
      chromosomeNameArray(),
      CHR_SIZE(){
      }

RefDictInfo::RefDictInfo(int NUM_AUTO_CHR,vector <string> chromosomeNameArray,
                         vector <long> CHR_SIZE)
    : NUM_AUTO_CHR(NUM_AUTO_CHR),
      chromosomeNameArray(chromosomeNameArray),
      CHR_SIZE(CHR_SIZE){
      }




RefDictInfo read_dict_from_file(string _refdictFilepath){

    cerr << "Reading in reference dict from " << _refdictFilepath << " ...\n";
    if (!isfile(_refdictFilepath)){
        cerr << _refdictFilepath << " does not exist. ERROR!" << endl;
        exit(1);
    }
    ifstream input_file;
    input_file.open(_refdictFilepath.c_str(), std::ios::in | std::ios::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> input_filter_stream_buffer;
    // input_filter_stream_buffer.push(boost::iostreams::gzip_decompressor());
    input_filter_stream_buffer.push(input_file);
    std::istream input_stream(&input_filter_stream_buffer);

    string line;
    int NUM_AUTO_CHR;
    vector <string> chromosomeNameArray;
    vector <long> CHR_SIZE;
    string first_string;

    std::getline(input_stream, line);
    while (!line.empty())
    {
        std::vector<std::string> element_vec = string_split(line, "\t");

        first_string = element_vec[0];
        std::getline(input_stream, line);

        if (first_string!="@SQ") {
            //ignore comments
            continue;
        }

        string second_string = element_vec[1];
        std::vector<std::string> element_vec2 = string_split(second_string, ":");
        string chr_string = element_vec2[1];
        regex pattern("^chr\\d+$");
        if (!regex_match(chr_string,pattern)) {
            continue;
        }
        string third_string = element_vec[2];
        std::vector<std::string> element_vec3 = string_split(third_string, ":");
        long chrosize = stol(element_vec3[1]);
        chromosomeNameArray.push_back(chr_string);
        CHR_SIZE.push_back(chrosize);
    }
    NUM_AUTO_CHR = CHR_SIZE.size();
    return RefDictInfo(NUM_AUTO_CHR=NUM_AUTO_CHR, chromosomeNameArray=chromosomeNameArray,
                         CHR_SIZE=CHR_SIZE);

//@SQ     SN:chr1 LN:248956422 

}
Config read_para(string file_conf)
{
    Config config;
    config.has_valid_para = false;

    string trash;
    try
    {
        if (!isfile(file_conf))
        {
            cout << fmt::format("ERROR: configure file {} does not exist.\n", file_conf);
            exit(1);
        }
        ifstream in(file_conf.c_str());
        in >> trash >> config.readlength;
        in >> trash >> config.window;
        in >> trash >> config.path_to_ref;
        in >> trash >> trash;  // samtools
        in >> trash >> trash;  // snp-caller
        in >> trash >> config.path;
        if (config.path_to_ref[config.path_to_ref.length() - 1] != '/')
            config.path_to_ref = config.path_to_ref + "/";
        if (config.path[config.path.length() - 1] != '/')
            config.path = config.path + "/";
        config.has_valid_para = true;
    }
    catch (...)
    {
        cout << fmt::format("ERROR: wrong format in configure file {}.\n", file_conf);
        exit(1);
    }
    cout << "Data from configure file:\n" TabMACRO
            "read length " TabMACRO config.readlength NewLineMACRO TabMACRO
            "window size " TabMACRO config.window NewLineMACRO TabMACRO
            "reference genome folder " TabMACRO config.path_to_ref NewLineMACRO TabMACRO
            "accurity path " TabMACRO config.path NewLineMACRO;
    return config;
}

void calculate_median_mad(double *input_array, long start_index, long stop_index,
                          float &median_value, float &mad_value) {
    vector<float> ratio_vector;
    for (int i = start_index; i < stop_index; i++){
        ratio_vector.push_back(input_array[i]);
    }
    sort(ratio_vector.begin(), ratio_vector.end(),
         [](float a, float b)->bool
         { return a > b; });
    int vector_size = ratio_vector.size();
    if (vector_size % 2 == 0)
        median_value = (ratio_vector[(int)vector_size / 2 - 1] +
                        ratio_vector[(int)vector_size / 2]) *
                       1.0 / 2;
    else
        median_value = ratio_vector[(int)(vector_size + 1) / 2 - 1];
    vector<float> median_diff_vector;
    for (uint i = 0; i < ratio_vector.size(); i++)
        median_diff_vector.push_back(
                fabs(ratio_vector[i] - median_value));
    sort(median_diff_vector.begin(), median_diff_vector.end(),
         [](float a, float b)->bool
         { return a > b; });
    if (vector_size % 2 == 0)
        mad_value = (median_diff_vector[(int)vector_size / 2 - 1] +
                     median_diff_vector[(int)vector_size / 2]) *
                    1.0 / 2;
    else
        mad_value = median_diff_vector[(int)(vector_size + 1) / 2 - 1];
}

/*
 * percent_to_exclude is a number from 0 to 100. usually 20.
 * For example, if percent_to_exclude=20, then the function will exclude 10% highest and 10% lowest values.
 */
void calculate_robust_mean_stddev(double *input_array, long start_index, long stop_index, int percent_to_exclude,
                                        float &mean_ref, float &stddev_ref) {
    vector<float> ratio_vector;
    for (int i = start_index; i < stop_index; i++){
        ratio_vector.push_back(input_array[i]);
    }
    sort(ratio_vector.begin(), ratio_vector.end(),
         [](float a, float b)->bool
         { return a > b; });
    int vector_size = ratio_vector.size();
    int lower_index = max(0, vector_size*percent_to_exclude/200);
    int upper_index = min(vector_size*(100-percent_to_exclude/2)/100+1, vector_size);
    int sample_size = upper_index-lower_index;
    float sum = 0;
    float sum_squared = 0;
    for (int i=lower_index; i<upper_index; i++){
        sum += ratio_vector[i];
        sum_squared += ratio_vector[i]*ratio_vector[i];
    }
    mean_ref = sum/(sample_size*1.0);
    float variance = sum_squared/(sample_size*1.0) - mean_ref*mean_ref;
    if (variance>=0) {
        stddev_ref = sqrt(variance);
    } else {
        stddev_ref = 0.0;
    }
}


void calculate_robust_mean_stddev(vector<float> float_vector, int percent_to_exclude,
                                  float &mean_ref, float &stddev_ref, double &squared_sum, int &sample_size) {
    sort(float_vector.begin(), float_vector.end(),
         [](float a, float b)->bool
         { return a > b; });
    int vector_size = float_vector.size();
    int lower_index = max(0, vector_size*percent_to_exclude/200);
    int upper_index = min(vector_size*(100-percent_to_exclude/2)/100+1, vector_size);
    float sum = 0;
    for (int i=lower_index; i<upper_index; i++){
        sample_size ++;
        sum += float_vector[i];
        squared_sum += float_vector[i]*float_vector[i];
    }
    mean_ref = sum/(sample_size*1.0);
    float variance = squared_sum/(sample_size*1.0) - mean_ref*mean_ref;
    if (variance>=0) {
        stddev_ref = sqrt(variance);
    } else {
        stddev_ref = 0.0;
    }
}

std::vector<std::string> string_split(std::string str,std::string sep){
    char* cstr=const_cast<char*>(str.c_str());
    char* current;
    std::vector<std::string> arr;
    current=strtok(cstr,sep.c_str());
    while(current!=NULL){
        arr.push_back(current);
        current=strtok(NULL,sep.c_str());
    }
    return arr;
}



// int main(int argc, char** argv) {
//	string f_conf="configure";
//	string hg;
//	int readlength, windowsize;
//	read_para(f_conf,hg,readlength,windowsize);
//	cout<< f_conf TabMACRO hg TabMACRO readlength TabMACRO windowsize
// NewLineMACRO;
//
//}
