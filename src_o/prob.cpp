#include "prob.h"

Prob::Prob()
{
    PI = 3.141592653589793238462;
    E = 2.71828182845;
}

void Prob::neg_bi_repara(double mn, double var, double &p, double &r)
{
    p = 1 - mn / var;
    r = mn * mn / (var - mn);
    // cerr<<"neg_bi_repara " << mn SEP var SEP p SEP r NL;
}

double Prob::gamma(double z)
{
// 20170312 approximation of gamma function.
// A more accurate approximation can be obtained by using more terms from the asymptotic expansions of ln(Γ(z)) and Γ(z), which are based on Stirling's approximation.
    return sqrt(2 * PI / z) * pow( (z + 1 / (12 * z - 1 / 10.0 / z))/E, z);
}

double Prob::factorial(double x)
{
    return gamma(x + 1);
}

double Prob::log_factorial(double x)
{
// approximation of gamma
    x++;
    return log(sqrt(2 * PI / x)) +
           log(1 / E * (x + 1 / (12 * x - 1 / 10.0 / x))) * x;
}

double Prob::nchoosek(double n, double k)
{
    // cerr<<"nchoosek " SEP factorial(n)/factorial(k)/factorial(n-k) SEP
    // nchoosek2(n,k) NL;
    return factorial(n) / factorial(k) / factorial(n - k);
}

double Prob::nchoosek2(double n, double k)
{
    return exp(log_factorial(n) - log_factorial(k) - log_factorial(n - k));
}

double Prob::binomial_max(int n, double p)
{
    double freq = 0;
    if (n == 0) return 0.0;
    for (int i = 0; i <= n; i++)
        freq += (std::max(i, n - i) *1.0 / double(n) * nchoosek2(n, i) *
                pow(p, i) * pow(1 - p, n - i));
    return freq;
}

double Prob::binomial_max_log10(int n, double p)
{
    double freq = 0;
    if (n == 0) return 0.0;
    for (int i = 0; i <= n; i++) {
        //20171227 log10
        freq += (log10(std::max(i, n - i)*1.0 / double(n)) * nchoosek2(n, i) *
                pow(p, i) * pow(1 - p, n - i) );
    }
    return freq;
}

double Prob::neg_bi(double p, double r, double i)
{
    // cerr<<" neg_bi " SEP i SEP r SEP nchoosek(i+r-1,i) SEP pow(1-p,r) SEP
    // pow(p,i) NL;
    return nchoosek2(i + r - 1, i) * pow(1 - p, r) * pow(p, i);
}

double Prob::moving_average(double *a, int array_length, double &mean,
                            double &stddev, int idx, int moving_average_window,
                            int right_moving_average_window)
{
    //	int moving_average_window=(int)(array_length*0.30/2);
    //	int right_moving_average_window=(int)(array_length*0.20/2);
    int counter = 0;
    double sum = 0, square_sum = 0;
    for (int i = max(idx - moving_average_window, 0);
         i <= min(max(idx - moving_average_window, 0) +
                      right_moving_average_window,
                  array_length - 1);
         i++)
    {
        counter++;
        sum += a[i];
        square_sum += (a[i] * a[i]);
    }
    mean = sum / counter;
    stddev = sqrt((square_sum - counter * mean * mean) / (counter - 1));
    return 0.0;
}

void Prob::calc_window_average(double *a, double *smoothed, int sample_size,
                               int width)
{
    cerr << "Calculating window average (width=" << width << ")...";
    for (int i = 0; i < sample_size; i++) smoothed[i] = 0;

    for (int i = 0; i < sample_size; i++)
    {
        if (i < width)
        {
            int l = width + i + 1;
            for (int j = 0; j < l; j++) smoothed[j] += a[i] / l;
        }
        else if ((sample_size - i - 1) < width)
        {
            int l = width + 1 + (sample_size - 1 - i);
            for (int j = sample_size - 1; j > sample_size - 1 - l; j--)
                smoothed[j] += a[i] / l;
        }
        else
        {
            for (int j = i - width; j <= i + width; j++)
                smoothed[j] += a[i] / (2 * width + 1);
        }
    }
    cerr << "Done.\n";
}
