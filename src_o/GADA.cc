/*
 eGADA enhanced Genome Alteration Detection Algorithm

 eGADA is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 any later version.

 eGADA is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with eGADA.  If not, see <http://www.gnu.org/licenses/>.

 Author:
 Yu S. Huang, polyactis@gmail.com

 */
#include <boost/program_options.hpp>  //for program options
#include <boost/tokenizer.hpp>
#include "BaseGADA.h"
#include "read_para.h"

using namespace std;
using namespace boost;
namespace po = boost::program_options;

typedef boost::tokenizer<boost::char_separator<char> > tokenizerCharType;	//to break a line into list by some delimiter

class GADA
{
   public:
    int argc;
    char **argv;  // 2013.08.20 somehow, can't use "char* argv[]" even though
                  // they are same.
    // otherwise "argv=_argv;" will not work for this error: incompatible types
    // in assignment of ‘char**’ to ‘char* [0]’
    string programName;
    boost::format usageDoc;
    boost::format examplesDoc;

    long input_array_len;
    double *input_array;
    string chromosome_id;
    std::vector<long> chr_start_pos_vector;

    int report;
    int reportIntervalDuringBE;  // how often to report progress during backward
                                 // elimination, default is 100K

    long debug;  // verbosity... set equal to 1 to see messages of SBLandBE(). 0
                 // to not see them
    double T;  // Backward elimination threshold
    // double T2=5.0; //Segment collapse to base Non-Alteration level threshold
    double BaseAmp;  // Base-level
    double a;  // SBL hyperprior parameter
    double sigma2;  // Variance observed, if negative value, it will be
                    // estimated by the mean of the differences
    // I would recommend to be estimated on all the chromosomes and as a trimmed
    // mean.
    long MinSegLen;  // Length in number of probes for a CNA segment to be
                     // called significan.
    long SelectClassifySegments;  // Classify segment into altered state (1),
                                  // otherwise 0
    long SelectEstimateBaseAmp;  // Estimate Neutral hybridization amplitude.
    double convergenceDelta;  // 1E-10 or 1E-8 seems to work well for this
                              // parameter. -- => ++ conv time
    // 1E8 better than 1E10 seems to work well for this parameter. -- => -- conv
    // time
    long maxNoOfIterations;      //=50000, //10000 is enough usually
    double convergenceMaxAlpha;  // 1E8 Maximum number of iterations to reach
                                 // convergence...
    double convergenceB;         // a number related to convergence = 1E-20
    int window_size;

    string input_file_path;
    string output_file_path;

    GADA(int _argc, char *_argv[]);  // 2013.08.28 commandline version

    virtual ~GADA()
    {
        chr_start_pos_vector.clear();
        free(input_array);
        // free(SegState);	//2013.08.30 SegState is not always allocated
        // with extra memory
    }

    // 2013.08.28 stream causes error " note: synthesized method ... required
    // here" because stream is noncopyable.
    po::options_description optionDescription;  //("Allowed options")
    po::positional_options_description positionOptionDescription;
    po::variables_map optionVariableMap;

    std::ifstream inputFile;
    std::ofstream outputFile;
    boost::iostreams::filtering_streambuf<boost::iostreams::input>
        inputFilterStreamBuffer;
    boost::iostreams::filtering_streambuf<boost::iostreams::output>
        outputFilterStreamBuffer;

    // 2013.08.28 stream causes error " note: synthesized method ... required
    // here" because stream is noncopyable.
    virtual void constructOptionDescriptionStructure();
    virtual void parseCommandlineOptions();
    virtual void openOneInputFile(
        string &input_file_path,
        boost::iostreams::filtering_streambuf<boost::iostreams::input> &
            inputFilterStreamBuffer);
    void readInputFile();

    virtual void openOutputFile();
    virtual void closeFiles();
    void run();

};

GADA::GADA(int _argc, char *_argv[]) : argc(_argc), argv(_argv)
{
    debug = 0;
    report = 0;
    programName = _argv[0];
    // strcpy(_argv[0], programName.c_str());	//2013.08.20 does no work

    cerr << "program name is " << programName << "." << endl;

    usageDoc = boost::format("%1% -i INPUTFNAME -o OUTPUTFNAME [OPTIONS]\n") %
               programName;
    examplesDoc = boost::format(
                      "%1% -i /tmp/input.tsv.gz -o /tmp/output.gz -M 1000 "
                      "--convergenceDelta 0.01 \n") %
                  programName;


    input_array_len = 1000;
    input_array = (double *)calloc(input_array_len, sizeof(double));

}



// 2013.08.28 stream causes error " note: synthesized method ... required here"
// because stream is noncopyable.
// boost smart pointer (or cxx11) solves the problem but too much work

void GADA::constructOptionDescriptionStructure()
{
    optionDescription.add_options()("help,h", "produce help message")
            ("TBackElim,T", po::value<double>(&T)->default_value(5.0),
             " is the backward elimination critical value for a breakpoint. i.e. "
                     "minimum (mean1-mean2)/stddev difference between two adjacent "
                     "segments.")
            ("aAlpha,a", po::value<double>(&a)->default_value(0.5),
             "is the SBL hyper prior parameter for a breakpoint. It is "
                     "the  shape parameter of the Gamma distribution. Higher "
                     "(lower) value means less (more) breakpoints expected a "
                     "priori.")
            ("MinSegLen,M", po::value<long>(&MinSegLen)->default_value(0),
             "is the minimum size in number of probes for a "
                     "segment to be deemed significant.")
            ("BaseAmp", po::value<double>(&BaseAmp)->default_value(0.0),
             "Mean amplitude associated to the Neutral state. If not "
                     "provided, and c option is used, then it is estimated as the "
                     "median value of all probes hybridization values after running "
                     "the algorithm. We recommend to estimate this on chromosomes "
                     "that are known to have a Neutral state on most areas. In some "
                     "cases this value may be known if we have been applied some "
                     "normalization, preprocessing or using another sample as ref.")
            ("sigma2,s", po::value<double>(&sigma2)->default_value(-1),
             "Variance observed, if negative value, it will be estimated by the "
                     "mean of the differences. "
                     "I would recommend to be estimated on all the chromosomes and as a "
                     "trimmed mean.")
            ("SelectClassifySegments,c", po::value<long>(&SelectClassifySegments)->default_value(0),
             "Classify segment into altered state (1), otherwise 0")
            ("SelectEstimateBaseAmp", po::value<long>(&SelectEstimateBaseAmp)->default_value(1),
             "toggle this to estimate BaseAmp from data, rather than "
                     "user-supplied.")
            ("convergenceDelta", po::value<double>(&convergenceDelta)->default_value(1E-8),
             "a delta number controlling convergence")
            ("maxNoOfIterations", po::value<long>(&maxNoOfIterations)->default_value(50000),
             "maximum number of iterations for EM convergence algorithm to run "
                     "before being stopped")
            ("convergenceMaxAlpha", po::value<double>(&convergenceMaxAlpha)->default_value(1E8),
             "one convergence related number, not sure what it does.")
            ("convergenceB", po::value<double>(&convergenceB)->default_value(1E-20),
             "one convergence related number, not sure what it does")
            ("debug,b", "toggle debug mode")
            ("report,r", "toggle report mode")
            ("reportIntervalDuringBE", po::value<int>(&reportIntervalDuringBE)->default_value(100000),
             "how often to report the break point to be removed during backward "
                     "elimination")
            ("chromosome_id", po::value<string>(&chromosome_id)->default_value("hello"), "chromosome ID for the input data")
            ("input_file_path,i", po::value<string>(&input_file_path),
             "input file path, csv file, gzipped or plain. Comment lines start with #."
                     " 4 columns with a header start,tumor_read_count,normal_read_count,read_count_ratio."
                     " This file can be an option or a positional argument.")
            ("output_file_path,o", po::value<string>(&output_file_path), "output filepath")
            ("window_size", po::value<int>(&window_size)->default_value(500), "the windows size used in GC normalization");
}

void GADA::parseCommandlineOptions()
{
    // all positional arguments are input files.
    positionOptionDescription.add("input_file_path", -1);

    po::store(po::command_line_parser(argc, argv)
                  .options(optionDescription)
                  .positional(positionOptionDescription)
                  .run(),
              optionVariableMap);

    // po::store(po::parse_command_line(argc, argv, optionDescription),
    // optionVariableMap);
    po::notify(optionVariableMap);
    if (optionVariableMap.count("help") || input_file_path.empty() ||
        output_file_path.empty())
    {
        cout << "Usage:" << endl << usageDoc << endl;

        cout << optionDescription << endl << endl;
        cout << "Examples:" << endl << examplesDoc << endl;
        exit(1);
    }
    if (optionVariableMap.count("debug"))
    {
        debug = 1;
    }
    else
    {
        debug = 0;
    }
    if (optionVariableMap.count("report"))
    {
        report = 1;
    }
    else
    {
        report = 0;
    }
}

void GADA::openOneInputFile(
    string &input_file_path, boost::iostreams::filtering_streambuf<
                            boost::iostreams::input> &inputFilterStreamBuffer)
{

    int inputFnameLength = input_file_path.length();
    if (input_file_path.substr(inputFnameLength - 3, 3) == ".gz")
    {
        inputFilterStreamBuffer.push(boost::iostreams::gzip_decompressor());
        inputFile.open(input_file_path.c_str(), std::ios::in | std::ios::binary);
    }
    else
    {
        inputFile.open(input_file_path.c_str(), std::ios::in);
    }
    inputFilterStreamBuffer.push(inputFile);
}

void GADA::readInputFile() {
    std::cerr << "Reading data from " << input_file_path << " ... ";
    openOneInputFile(input_file_path, inputFilterStreamBuffer);
    std::istream inputStream(&inputFilterStreamBuffer);

    boost::char_separator<char> sep(",");		//csv

    long array_moving_index = 0;
    int columnMovingIndex;
    string statInStr;
    std::string line;
    std::getline(inputStream, line);
    while (!line.empty())
    {
        if ( line.compare(0, 1, "#")!=0 && line.compare(0, 5, "start")!=0 ) {
            tokenizerCharType line_toks(line, sep);
            columnMovingIndex = 0;
            for(tokenizerCharType::iterator it = line_toks.begin(), ite = line_toks.end(); it!=ite; ++it) {
                statInStr = *it;
                if (columnMovingIndex==0){
                    chr_start_pos_vector.push_back(atol(statInStr.c_str()));
                }
                else if (columnMovingIndex==1) {
                    input_array[array_moving_index++] = atof(statInStr.c_str());
                    if (array_moving_index >= input_array_len) {
                        input_array_len = input_array_len + 10000;
                        input_array = (double *) realloc(input_array, input_array_len * sizeof(double));
                    }
                }
                columnMovingIndex++;
            }
        }
        std::getline(inputStream, line);
    }
    input_array_len = array_moving_index;
    input_array = (double *)realloc(input_array, input_array_len * sizeof(double));
    std::cerr << input_array_len << " data points for chromosome " << chromosome_id << "." << endl;
}

void GADA::openOutputFile()
{
    if (!output_file_path.empty())
    {
        std::cerr << "Open file " << output_file_path << " for writing ";
        int outputFnameLength = output_file_path.length();
        if (output_file_path.substr(outputFnameLength - 3, 3) == ".gz")
        {
            // boost::iostreams::gzip_compressor gzipCompressor;
            outputFilterStreamBuffer.push(boost::iostreams::gzip_compressor());
            // outputFilterStreamBuffer.push(boost::iostreams::base64_encoder());
            // outputFile.open(output_file_path.c_str(), std::ios::out |
            // std::ios::binary);
            outputFile.open(output_file_path.c_str(),
                            std::ios::out | std::ios::binary);
        }
        else
        {
            outputFile.open(output_file_path.c_str(), std::ios::out);
        }
        outputFilterStreamBuffer.push(outputFile);
        // outputStream.rdbuf(&outputFilterStreamBuffer);
        //(&outputFilterStreamBuffer);
        std::cerr << endl;
    }
    else
    {
        std::cerr << "Warning: Output file, " << output_file_path
                  << ", is an empty string." << endl;
    }
}

void GADA::closeFiles()
{
    //
    std::cerr << "Closing files ...";
    inputFile.close();
    if (!output_file_path.empty())
    {

        // delete outputFilterStreamBuffer;
        outputFile.flush();
        // outputFile.close();	//2013.08.21 if output is a .gz, closing it
        // here would result a truncated .gz file. don't know why.
        // maybe because the buffer or the gzip compressor filter needs to be
        // closed beforehand.
    }
    std::cerr << std::endl;
}




void GADA::run()
{
    constructOptionDescriptionStructure();
    parseCommandlineOptions();
    readInputFile();

    std::cerr << "Running SBLandBE ... " << endl;
    BaseGADA baseGADA =
        BaseGADA(input_array, input_array_len, sigma2, BaseAmp, a, T, MinSegLen, debug,
                 convergenceDelta, maxNoOfIterations, convergenceMaxAlpha,
                 convergenceB, reportIntervalDuringBE);
    baseGADA.SBLandBE();
    // K = SBLandBE(input_array, input_array_len, &sigma2, a, T, MinSegLen, &Iext, &Wext, debug ,
    // delta, numEMsteps, noOfBreakpointsAfterSBL, convergenceDelta,
    // maxNoOfIterations, convergenceMaxAlpha, convergenceB);
    std::cerr << boost::format(" %1% breakpoints after SBL, %2% breakpoints after BE.\n") %
                   baseGADA.noOfBreakpointsAfterSBL % baseGADA.K;
    // std::cerr<< boost::format("Backward elimination (T=%2%) and remove
    // segments that are shorter than %1% ... ") % MinSegLen % T;
    std::cerr << boost::format(
                         "Starting IextToSegLen() & IextWextToSegAmp() ... ");

    baseGADA.IextToSegLen();
    baseGADA.IextWextToSegAmp();

    std::cerr << " IextToSegLen() & IextWextToSegAmp() done." << endl;
    std::cerr << "Outputting final result ... ";
    openOutputFile();  // 2013.08.30 open this file here. Do not open it way
                       // before the main writing starts.
    // it will leave a long period of zero-writing-activity (due to
    // computation), which could hang the program sometimes on panfs system

    std::ostream outputStream(&outputFilterStreamBuffer);

    //outputStream << "# GADA Genome Alteration Detection Algorithm\n";
    //outputStream << "# Author: www.yfish.org polyactis@gmail.com. Originally from Roger Pique-Regi\n";
    outputStream << boost::format(
            "# Parameters: a=%1%,T=%2%,MinSegLen=%3%,sigma2=%4%,BaseAmp=%5%, convergenceDelta=%6%, maxNoOfIterations=%7%, "
                    "convergenceMaxAlpha=%8%, convergenceB=%9%.\n") %
            baseGADA.a % baseGADA.T % baseGADA.MinSegLen %
            baseGADA.sigma2 % baseGADA.BaseAmp %
            baseGADA.convergenceDelta % baseGADA.maxNoOfIterations %
            baseGADA.convergenceMaxAlpha % baseGADA.convergenceB;
    outputStream << boost::format("# %1% data points in input file\n") %
                        baseGADA._M_total_length;
    outputStream << boost::format("# Overall mean %1%\n") % baseGADA.Wext[0];
    outputStream << boost::format("# Sigma^2=%1%\n") % baseGADA.sigma2;
    outputStream << boost::format(
                        "# Convergence: delta=%1% after %2% EM iterations.\n") %
                        baseGADA.delta % baseGADA.numEMsteps;
    outputStream << boost::format("# Found %1% breakpoints after SBL\n") %
                        baseGADA.noOfBreakpointsAfterSBL;
    outputStream << boost::format("# Kept %1% breakpoints after BE\n") %
                        baseGADA.K;


    if (SelectClassifySegments == 0)
    {
        //outputStream << boost::format("Chromosome\tStart\tStop\tMean\tStddev\tNoOfValidWindows\n");
        for (int i = 0; i < baseGADA.K + 1; i++) {
            int chr_start_pos = chr_start_pos_vector[baseGADA.Iext[i]];
            int chr_stop_pos = chr_start_pos_vector[baseGADA.Iext[i+1]-1] + window_size - 1;
            float segment_mean;
            float segment_stddev;
            calculate_robust_mean_stddev(input_array, baseGADA.Iext[i], baseGADA.Iext[i+1], 40, segment_mean, segment_stddev);
            outputStream << chromosome_id << "\t" << chr_start_pos
                         << "\t" << chr_stop_pos
                         << "\t" << segment_mean
                         << "\t" << segment_stddev
                         << "\t" << baseGADA.SegLen[i]
                         << std::endl;
        }
    }
    else if (SelectClassifySegments == 1)
    {
        std::cerr << " Select and Classify Segments ...";
        if (SelectEstimateBaseAmp == 1)
        {
            outputStream << boost::format("# Estimating BaseAmp\n");
            baseGADA.CompBaseAmpMedianMethod();
        }
        outputStream << boost::format("# BaseAmp=%1% \n") % baseGADA.BaseAmp;
        // outputStream<< boost::format("# Classify Segments \n", BaseAmp);

        baseGADA.CollapseAmpTtest();

        std::cerr << " SelectClassifySegments done.\n";

        outputStream << boost::format("Start\tStop\tLength\tAmpl\tState\n");
        for (int i = 0; i < baseGADA.K + 1; i++)
        {
            outputStream << baseGADA.Iext[i] + 1 << "\t" << baseGADA.Iext[i + 1]
                         << "\t" << baseGADA.SegLen[i] << "\t"
                         << baseGADA.SegAmp[i];
            if (baseGADA.SegState[i] > baseGADA.BaseAmp)
                outputStream << "G";
            else if (baseGADA.SegState[i] < baseGADA.BaseAmp)
                outputStream << "L";
            else
                outputStream << "N";
            outputStream << endl;
        }
    }
    std::cerr << " output done." << endl;
    closeFiles();
}

int main(int argc, char *argv[])
{
    GADA instance(argc, argv);
    instance.run();
}
