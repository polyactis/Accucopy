#pragma once
#ifndef __PROB_H
#define __PROB_H

#include <algorithm>
#include <cmath>
#include <iostream>

using namespace std;

class Prob
{
   public:
    Prob();

    double PI;
    double E;
    double binomial_max(int n, double p);
    double binomial_max_log10(int n, double p);
    double factorial(double x);
    double gamma(double z);
    double log_factorial(double x);
    double nchoosek(double n, double k);
    double nchoosek2(double n, double k);
    double neg_bi(double p, double r, double i);
    void neg_bi_repara(double mn, double var, double &p, double &r);
    double moving_average(double *a, int size, double &mean, double &dev,
                          int idx, int moving_average_window,
                          int right_moving_average_window);
    void calc_window_average(double *a, double *smoothed, int sample_size,
                             int width);

   private:
    int counter;
};
#endif
