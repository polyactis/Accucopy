#include "model_selection.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/normal.hpp>
#include "H5Cpp.h"

struct Model {
    int minor;
    int major;
};

struct Arg {
    double alpha;
    double mu;
    double variance;
    int convergence;
};

Model_Selection::Model_Selection(std::vector<double> &data, int cp, 
    double purity, double tumor_depth, std::ostream &out):
        data(data), cp(cp), purity(purity), tumor_depth(tumor_depth), out(out) {}

void Model_Selection::run() {
    //generate model
    for(int minor = 0; minor < (cp/2 + 1); minor++) {
        result.model_list.push_back(std::make_pair(minor, cp - minor));
    }
    std::sort(data.begin(),data.end());
    double min_num = data[0];
    double max_num = data.back();
    int data_size = data.size();
    double mean = std::accumulate(data.begin(),data.end(),0.0)*1.0/data_size;

    double alpha1_new = 0;
    double mu = 0;
    double var_new = 0;
    bool convergence = false;

    for(auto &model: result.model_list) {
        int minor = model.first;
        int major = model.second;
        if(minor == major) { //for single Gaussian model
            alpha1_new = 1.0;
            mu = 0.0;
            var_new = std::accumulate(data.begin(), data.end(), 0.0,
                        [&mean](double x, double y){return x+std::pow(y - mean, 2);})/data_size;
            convergence = true;
            out << "\n\n\nFor model: " << minor << "-" << major << std::endl;
            out << "Single Gaussian model: var=" << var_new << std::endl;
        }
        else { // for gaussian mixture model
            //initialize parameters
            double alpha1 = (mean - min_num)/(max_num - min_num);
            double var = std::accumulate(data.begin(), data.end(), 0.0,
                           [&mean](double x, double y){return x+std::pow(y-mean, 2);})/data_size;
            //adjusted expected mean for logOR
            double lambda_a = (minor * purity + 1 - purity)/(cp * purity + 2 * (1 - purity))\
                               * tumor_depth;
            double lambda_b = (major * purity + 1 - purity)/(cp * purity + 2 * (1 - purity))\
                               * tumor_depth;
            boost::math::poisson_distribution<> poisson_a(lambda_a);
            boost::math::poisson_distribution<> poisson_b(lambda_b);
            double a_pmf_accumlate, b_pmf_accumlate;
            double pmf_accumlate = 0;
            double logOR_accumlate = 0;
            double a_pmf, b_pmf;
            int b, a;
            a = 1;
            a_pmf_accumlate = boost::math::pdf(poisson_a, 0);
            while(a_pmf_accumlate < 0.9999) {
                a_pmf = boost::math::pdf(poisson_a, a);
                a_pmf_accumlate += a_pmf;

                b = 1;
                b_pmf_accumlate = boost::math::pdf(poisson_b, 0);
                while(b_pmf_accumlate < 0.9999) {
                    b_pmf = boost::math::pdf(poisson_b, b);
                    b_pmf_accumlate += b_pmf;
                    pmf_accumlate += a_pmf * b_pmf;
                    logOR_accumlate += a_pmf * b_pmf * std::log(1.0 * a / b);
                    b += 1;
                }
                a += 1;
            }
            mu = logOR_accumlate / pmf_accumlate;

            std::vector<double> hidden_variable1(data_size, 0);
            convergence = true;
            out << "\n\n\nEM for model: " << minor << "-" << major << std::endl;
            out << "initialize parameters" << std::endl;
            out << "alpha1=" << alpha1 << ", alpha2=" << 1 - alpha1 << std::endl
                << "mu1=" << -mu << ", mu2=" << mu << std::endl
                << "var1=" << var << ", var2=" << var << std::endl;
            
            // define threshold
            double delta=1e-4;
            int iter_max = 10000;
            int iter_count = 0;

            while(true) {
                // E step
                boost::math::normal norm1(-mu, std::sqrt(var));
                boost::math::normal norm2(mu, std::sqrt(var));
                //data and hidden_variable1 has the same length
                for(auto data_iter = data.begin(), hidden_iter = hidden_variable1.begin();
                    data_iter != data.end() && hidden_iter != hidden_variable1.end();
                    data_iter++, hidden_iter++) {
                        *hidden_iter = alpha1 * boost::math::pdf(norm1, *data_iter)/\
                                       (alpha1 * boost::math::pdf(norm1, *data_iter)
                                        + (1 - alpha1) * boost::math::pdf(norm2, *data_iter));
                }
                // M step
                double sum1 = 0;
                double denominator = 0;
                for(auto data_iter = data.cbegin(), hidden_iter = hidden_variable1.cbegin();
                    data_iter != data.cend() && hidden_iter != hidden_variable1.cend();
                    data_iter++, hidden_iter++) {
                        sum1 += (*hidden_iter) * std::pow((*data_iter) - (-mu), 2)
                                 + (1 - (*hidden_iter)) * std::pow((*data_iter) - mu , 2);
                        denominator += (*hidden_iter);
                }
                var_new = sum1 / data_size;
                alpha1_new = denominator / data_size;
                //somethimes due to the precision, alpha1_new will greater than 1 when
                //alpha1_new is very close to 1
                if(alpha1_new > 1.0) {
                    alpha1_new = 1;
                }

                iter_count += 1;
                out << "Iteration " << iter_count << std::endl
                    << "alpha1=" << alpha1_new << ", alpha2=" << 1 - alpha1_new << std::endl
                    << "mu1=" << -mu << ", mu2=" << mu << std::endl
                    << "var1=" << var_new << ", var2=" << var_new << std::endl;
                // convergence ?
                if(std::abs((alpha1 - alpha1_new)/alpha1) < delta
                   && std::abs((var - var_new)/var) < delta) {
                       break;
                }
                //hit iter_max ?
                if(iter_count > iter_max) {
                    out << "Warning: the EM algorithm hit maximum iteration "
                        << "before convergence!!!" << std::endl;
                        convergence = false;
                        break;
                }
                //update parameters
                alpha1 = alpha1_new;
                var = var_new;
            }
            out << std::endl << "Optimal value:" << std::endl
                << "alpha1=" << alpha1_new << ", alpha2=" << 1 - alpha1_new << std::endl
                << "mu1=" << -mu << ", mu2=" << mu << std::endl
                << "var1=" << var_new << ", var2=" << var_new << std::endl;
        }
        result.arg_list.push_back(std::make_tuple(alpha1_new, -mu, var_new, convergence));
        // calculate logL
        double logL = 0;
        boost::math::normal norm1(-mu, std::sqrt(var_new));
        boost::math::normal norm2(mu, std::sqrt(var_new));
        for(auto &x: data) {
            logL += std::log(alpha1_new * boost::math::pdf(norm1, x) + (1 - alpha1_new)
                     * boost::math::pdf(norm2, x));
        }
        result.logL_list.push_back(logL);
    }
    // select model
    double max_logL = result.logL_list[0];
    int max_idx = 0;
    for(int i = 0; i < result.logL_list.size(); i++) {
        if(result.logL_list[i] > max_logL) {
            max_logL = result.logL_list[i];
            max_idx = i;
        }
    }
    result.best_logL = max_logL;
    result.selection = result.model_list[max_idx];
}

HDF5_Log::HDF5_Log(std::string filename):
    file(H5std_string(filename), H5F_ACC_TRUNC), group_idx(0) {}

void HDF5_Log::write(std::string &data_name, std::vector<double> &data, Result &result) {
    hsize_t dims[1] = {1};
    H5std_string groupname(std::string("/segment_") + std::to_string(group_idx++));
    H5::Group group(file.createGroup(groupname));

    H5std_string dataset_name("title");
    H5::DataSpace * dataspace = new H5::DataSpace(1, dims);
    H5::DataSet * dataset = new H5::DataSet(group.createDataSet(dataset_name,
        H5::StrType(0, H5T_VARIABLE), *dataspace));
    dataset->write(data_name, H5::StrType(0, H5T_VARIABLE));
    delete dataspace;
    delete dataset;

    dataset_name = "model_list";
    int model_size = result.model_list.size();
    dims[0] = model_size;
    Model * model_list = new Model[model_size];
    for(int i = 0; i < model_size; i++) {
        model_list[i].minor = result.model_list[i].first;
        model_list[i].major = result.model_list[i].second;
    }
    dataspace = new H5::DataSpace(1, dims);
    H5::CompType model_type(sizeof(Model));
    model_type.insertMember(std::string("minor"), HOFFSET(Model, minor), H5::PredType::NATIVE_INT);
    model_type.insertMember(std::string("major"), HOFFSET(Model, major), H5::PredType::NATIVE_INT);
    dataset = new H5::DataSet(group.createDataSet(dataset_name, model_type, *dataspace));
    dataset->write(model_list, model_type);
    delete [] model_list;
    delete dataset;
    delete dataspace;

    dataset_name = "arg_list";
    int arg_size = result.arg_list.size();
    dims[0] = arg_size;
    Arg * arg_list = new Arg[arg_size];
    for(int i = 0; i < arg_size; i++) {
        arg_list[i].alpha = std::get<0>(result.arg_list[i]);
        arg_list[i].mu = std::get<1>(result.arg_list[i]);
        arg_list[i].variance = std::get<2>(result.arg_list[i]);
        arg_list[i].convergence = std::get<3>(result.arg_list[i]) ? 1 : 0;
    }
    dataspace = new H5::DataSpace(1, dims);
    H5::CompType arg_type(sizeof(Arg));
    arg_type.insertMember(std::string("alpha"), HOFFSET(Arg, alpha), H5::PredType::NATIVE_DOUBLE);
    arg_type.insertMember(std::string("mu"), HOFFSET(Arg, mu), H5::PredType::NATIVE_DOUBLE);
    arg_type.insertMember(std::string("variance"), HOFFSET(Arg, variance), H5::PredType::NATIVE_DOUBLE);
    arg_type.insertMember(std::string("convergence"), HOFFSET(Arg, convergence), H5::PredType::NATIVE_INT);
    dataset = new H5::DataSet(group.createDataSet(dataset_name, arg_type, *dataspace));
    dataset->write(arg_list, arg_type);
    delete [] arg_list;
    delete dataspace;
    delete dataset;

    dataset_name = "logL_list";
    dims[0] = result.logL_list.size();
    dataspace = new H5::DataSpace(1, dims);
    dataset = new H5::DataSet(group.createDataSet(dataset_name, H5::PredType::NATIVE_DOUBLE, *dataspace));
    dataset->write(result.logL_list.data(), H5::PredType::NATIVE_DOUBLE);
    delete dataspace;
    delete dataset;

    dataset_name = "selection";
    int selection[2];
    selection[0] = result.selection.first;
    selection[1] = result.selection.second;
    dims[0] = 2;
    dataspace = new H5::DataSpace(1, dims);
    dataset = new H5::DataSet(group.createDataSet(dataset_name, H5::PredType::NATIVE_INT, *dataspace));
    dataset->write(selection, H5::PredType::NATIVE_INT);
    delete dataspace;
    delete dataset;

    dataset_name = "best_logL";
    dims[0] = 1;
    dataspace = new H5::DataSpace(1, dims);
    dataset = new H5::DataSet(group.createDataSet(dataset_name, H5::PredType::NATIVE_DOUBLE, *dataspace));
    dataset->write(&result.best_logL, H5::PredType::NATIVE_DOUBLE);
    delete dataspace;
    delete dataset;

    dataset_name = "data";
    dims[0] = data.size();
    dataspace = new H5::DataSpace(1, dims);
    dataset = new H5::DataSet(group.createDataSet(dataset_name, H5::PredType::NATIVE_DOUBLE, *dataspace));
    dataset->write(data.data(), H5::PredType::NATIVE_DOUBLE);
    delete dataspace;
    delete dataset;
    group.close();
}