#include <vector>
#include <iostream>
#include <tuple>
#include <string>
#include "H5Cpp.h"

struct Result {
    std::vector<std::pair<int,int> > model_list;
    std::vector<std::tuple<double,double,double,bool> > arg_list;
    std::vector<double> logL_list;
    std::pair<int,int> selection;
    double best_logL;
};

class HDF5_Log {
    private:
        H5::H5File file;
        size_t group_idx;
    public:
        HDF5_Log(std::string filename);
        void write(std::string &data_name, std::vector<double> &data, Result &result);
};

class Model_Selection {
    private:
        std::vector<double> &data;
        int cp;
        double purity;
        double tumor_depth;
        std::ostream &out;
    public:
        Result result;
        Model_Selection(std::vector<double> &data, int cp, double purity,
                        double tumor_depth, std::ostream &out);
        void run();
};