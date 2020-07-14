/*

 * recall_precision.cpp
 *
 *  Created on: 2017.2.17
 *      Author: fanxp
 */

#include "singly_linked_list.h"
#include <string>
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

#define max(a, b) (((a) > (b)) ? (a) : (b))
#define min(a, b) (((a) < (b)) ? (a) : (b))
#define LEN 22

using namespace std;

long chr_len[22] = {249250621, 243199373, 198022430, 191154276, 180915260,
                    171115067, 159138663, 146364022, 141213431, 135534747,
                    135006516, 133851895, 115169878, 107349540, 102531392,
                    90354753,  81195210,  78077248,  59128983,  63025520,
                    48129895,  51304566};

class Segment
{
   public:
    int chr_id;
    long start;
    long end;
    double copynumber;

    Segment() : chr_id(0), start(0), end(0), copynumber(0.0) {}
    Segment(int chr, long start, long end, double copynumber)
        : chr_id(chr), start(start), end(end), copynumber(copynumber)
    {
    }
    bool operator<(const Segment& seg) const;
    friend std::ostream& operator<<(std::ostream& os, const Segment& seg);
};

class GenomicSegments
{
   public:
    SLinkedList<Segment> segments[LEN];

    GenomicSegments(){};
    ~GenomicSegments(){};
    void add(Segment& seg);
    friend std::ostream& operator<<(std::ostream& os, GenomicSegments& gs);
    SLinkedList<Segment>& operator[](int chr_idx);
};

void ReadFromActualResult(string& filename, GenomicSegments& gs);
void ReadFromAccurityResult(string& filename, GenomicSegments& gs);
double SumOfArea(GenomicSegments& gs, int length);
double SimilarityOfXtoY(GenomicSegments& XonY, int length);
void CallRecallAndPrecision(GenomicSegments& actual, GenomicSegments& accurity,
                            double& recall, double& precision);
void YSegmentsMapOnX(GenomicSegments& X, GenomicSegments& Y,
                     GenomicSegments& XandY, int length);

int main(int argc, char* argv[])
{
    GenomicSegments actual;
    GenomicSegments accurity;
    double recall, precision;

    string resultfile = argv[1];
    string accurityfile = argv[2];
    string outputfile = argv[3];

    ReadFromActualResult(resultfile, actual);
    ReadFromAccurityResult(accurityfile, accurity);
    CallRecallAndPrecision(actual, accurity, recall, precision);
    cout << "recall: " << recall << " precision: " << precision << endl;

    ofstream outfile(outputfile.c_str());
    outfile << "recall"
            << "\t"
            << "precision" << endl;
    outfile << recall << "\t" << precision << endl;
    outfile.close();
    return 0;
}

void GenomicSegments::add(Segment& seg)
{
    int chr_idx = seg.chr_id;
    segments[chr_idx - 1].Add(seg);
}

std::ostream& operator<<(std::ostream& os, GenomicSegments& gs)
{
    int i = 0;
    for (; i < LEN; i++) gs[i].PrintSLList();
    return os;
}

SLinkedList<Segment>& GenomicSegments::operator[](int chr_idx)
{
    return segments[chr_idx];
}

bool Segment::operator<(const Segment& seg) const
{
    if (start < seg.start)
        return true;
    else
        return false;
}

std::ostream& operator<<(std::ostream& os, const Segment& seg)
{
    os << "chr_id"
       << "\t" << seg.chr_id << "\t"
       << "start"
       << "\t" << seg.start << "\t"
       << "end"
       << "\t" << seg.end << "\t"
       << "copy number"
       << "\t" << seg.copynumber;
    return os;
}

void ReadFromActualResult(string& filename, GenomicSegments& gs)
{
    ifstream infile(filename.c_str());
    string trash, chr;
    double copynumber;
    long start, end;
    int chr_idx;
    stringstream ss;
    if (!infile)
    {
        cout << "no such file: " << filename << endl;
        exit(0);
    }
    getline(infile, trash);  // skip header
    while (infile >> chr >> start >> end >> copynumber >> trash)
    {
        ss << chr.substr(3);
        ss >> chr_idx;
        ss.clear();
        if (chr_idx < 1 || chr_idx > 22) continue;  // skip non autosome
        Segment seg(chr_idx, start, end, copynumber);
        gs.add(seg);
    }
    infile.close();
}

void ReadFromAccurityResult(string& filename, GenomicSegments& gs)
{
    ifstream infile(filename.c_str());
    string trash;
    double copynumber;
    long start, end;
    int chr_idx;
    if (!infile)
    {
        cout << "no such file: " << filename << endl;
        exit(0);
    }
    getline(infile, trash);  // skip header
    while (infile >> chr_idx >> trash >> trash >> copynumber >> trash >>
           trash >> trash >> trash >> trash >> trash >> start >> end)
    {
        if (copynumber == 2) continue;
        Segment seg(chr_idx, start, end, copynumber);
        gs.add(seg);
    }
    infile.close();
}

double SumOfArea(GenomicSegments& gs, int length)
{
    int i;
    double sum = 0;
    for (i = 0; i < length; i++)
    {
        if (gs[i].IsEmpty()) continue;
        SLinkedList<Segment>::Iterator it = gs[i].begin();
        for (; it != gs[i].end(); ++it) sum += (*it).end - (*it).start + 1;
    }
    return sum;
}

double SimilarityOfXtoY(GenomicSegments& XonY, int length)
{
    int i;
    double similarity = 0;
    for (i = 0; i < length; i++)
    {
        if (XonY[i].IsEmpty()) continue;
        SLinkedList<Segment>::Iterator it = XonY[i].begin();
        for (; it != XonY[i].end(); ++it)
            similarity +=
                ((*it).end - (*it).start + 1) * exp(-(*it).copynumber);
    }
    return similarity;
}

void YSegmentsMapOnX(GenomicSegments& X, GenomicSegments& Y,
                     GenomicSegments& XandY, int length)
{
    int i;
    long start, end;

    for (i = 0; i < length; i++)
    {
        if (X[i].IsEmpty() || Y[i].IsEmpty()) continue;
        SLinkedList<Segment>::Iterator Yi = Y[i].begin();
        for (; Yi != Y[i].end(); ++Yi)
        {
            SLinkedList<Segment>::Iterator Xi = X[i].begin();
            for (; Xi != X[i].end(); ++Xi)
            {
                start = max((*Yi).start, (*Xi).start);
                end = min((*Yi).end, (*Xi).end);
                if (start < end + 1)
                {
                    Segment seg(i + 1, start, end,
                                fabs((*Yi).copynumber - (*Xi).copynumber));
                    XandY.add(seg);
                }
            }
        }
    }
}

void CallRecallAndPrecision(GenomicSegments& actual, GenomicSegments& accurity,
                            double& recall, double& precision)
{
    double similarity;
    GenomicSegments XonY;
    YSegmentsMapOnX(actual, accurity, XonY, LEN);
    /*cout << "actual\n" << actual << "accurity\n" << accurity
             << "XonY\n" << XonY << endl;*/
    similarity = SimilarityOfXtoY(XonY, LEN);
    recall = similarity / SumOfArea(actual, LEN);
    precision = similarity / SumOfArea(accurity, LEN);
}
