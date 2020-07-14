/*=================================================================
 * BaseGenomeBreaks.h
 *
 *=================================================================*/
/*
 This File is part of GADA

 GADA v1.0 Genome Alteration Detection Algorithm
 Copyright (C) 2008  Childrens Hospital of Los Angeles

 GADA is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 any later version.

 GADA is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with GADA.  If not, see <http://www.gnu.org/licenses/>.

 Author:
 Roger Pique-Regi    piquereg@usc.edu

 */

#ifndef _BaseGADA_H_
#define _BaseGADA_H_

//#include "matlabdefines.h"
#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <map>	//for hash_map
#include <set>	//for set
#include <functional>	//2013.09.11 for customize hash
#include <boost/functional/hash.hpp>	//2013.09.10 yh: for customize boost::hash
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include "RedBlackTree.h"	//2013.09.19 red-black tree to store segment breakpoint, score, etc.

#define log2(x) log(x)/log(2)
//#define min(x,y) x<y?x:y


using namespace std;


class BreakPointKey{
	/*
	 * 2013.09.22 class defines the key of a breakpoint in a red-black tree.
	 * 	content is similar to BreakPoint but different ordering function.
	 */
	// for output
	friend ostream& operator<<(ostream& out, BreakPointKey& bpKey){
		/*
		 * 2013.09.22 something wrong here, it can't be streamed to an ostream
		 */
		out << boost::format("position=%1%, tscore=%2%, weight=%3%, length=%4%, MinSegLen=%5%, totalLength=%6%")%
				bpKey.position % bpKey.tscore % bpKey.weight %
				bpKey.segmentLength % bpKey.MinSegLen % bpKey.totalLength;
		return out;
	}
	// define comparison operator for the RBTree, how to order them
	// For all break points below MinSegLen, order by tscore
	// For all break points above MinSegLen, order by tscore
	// for break points crossing MinSegLen, order by segmentLength
	friend bool operator<(const BreakPointKey& a, const BreakPointKey& b){
		/*
		 * make sure no duplicate (or close keys, for floating values) keys exist in red-black trees.
		 */
		//return a.tscore < b.tscore;

		if(a.tscore<a.minTScore && b.tscore<b.minTScore){	//both below minTScore, then order them by score
			if (a.tscore==b.tscore){
				return a.segmentLength<b.segmentLength;
			}
			else{
				return a.tscore<b.tscore;
			}
		}
		else if ((a.tscore<a.minTScore && b.tscore>=b.minTScore) || (a.tscore>=a.minTScore && b.tscore<b.minTScore)){	//one below, one above
			if (a.tscore==b.tscore){
				return a.segmentLength<b.segmentLength;
			}
			else{
				//order by tscore
				return a.tscore<b.tscore;
			}
		}
		else{	//both above minTScore (or one equal, one above or both equal), then order them by segment length, unless they are identical
			if (a.segmentLength==b.segmentLength){
				return a.tscore<b.tscore;
			}
			else{
				return a.segmentLength<b.segmentLength;
			}
		}

	}

	friend bool operator>(const BreakPointKey& a, const BreakPointKey& b){
		/*
		 * make sure no duplicate (or close keys, for floating values) keys exist in red-black trees.
		 */
		//return a.tscore > b.tscore;
		if(a.tscore<a.minTScore && b.tscore<b.minTScore){	//both below minTScore, then order them by score, unless identical
			if (a.tscore==b.tscore){
				return a.segmentLength>b.segmentLength;
			}
			else{
				return a.tscore>b.tscore;
			}
		}
		else if ((a.tscore<a.minTScore && b.tscore>=b.minTScore) || (a.tscore>=a.minTScore && b.tscore<b.minTScore)){	//one below, one above
			if (a.tscore==b.tscore){
				return a.segmentLength>b.segmentLength;
			}
			else{
				//order by tscore
				return a.tscore>b.tscore;
			}
		}
		else{	//both above minTScore (or one equal, one above or both equal), then order them by segment length, unless they are identical
			if (a.segmentLength==b.segmentLength){
				return a.tscore>b.tscore;
			}
			else{
				return a.segmentLength>b.segmentLength;
			}
		}
	}

	friend bool operator==(const BreakPointKey& a, const BreakPointKey& b) {
		return a.tscore == b.tscore;
	}

public:
	long position;
	double weight;
	double tscore;
	long segmentLength;
	long MinSegLen;
	double minTScore;
	long totalLength;
	void* nodePtr;
	BreakPointKey(){
		position=1;
		weight=0;
		tscore=0;
		segmentLength=0;
		MinSegLen=0;
		minTScore = 0;
		totalLength=0;
		nodePtr = NULL;
	}
	BreakPointKey(long _position, double _weight, double _tscore, long _segmentLength, long _MinSegLen, double _minTScore, long _totalLength):
		position(_position),weight(_weight), tscore(_tscore), segmentLength(_segmentLength), MinSegLen(_MinSegLen), minTScore(_minTScore), totalLength(_totalLength){
		nodePtr = NULL;
	};
	~BreakPointKey(){
	}
};


class BreakPoint{
	// for output
	friend ostream& operator<<(ostream& out, BreakPoint& breakPoint){
		out << boost::format("position=%1%, tscore=%2%, weight=%3%, length=%4%, MinSegLen=%5%, totalLength=%6%")%
				breakPoint.position % breakPoint.tscore % breakPoint.weight %
				breakPoint.segmentLength % breakPoint.MinSegLen % breakPoint.totalLength;
		return out;
	}
	friend bool operator<(const BreakPoint& a, const BreakPoint& b){
		/*
		 * to be used in std::sort
		 */
		return a.position<b.position;
	}

	friend bool operator==(const BreakPoint& a, const BreakPoint& b) {
		return a.position==b.position;
	}

public:
	long position;
	double weight;
	double tscore;
	long segmentLength;
	long MinSegLen;
	double minTScore;
	long totalLength;
	BreakPoint* leftBreakPointPtr;
	BreakPoint* rightBreakPointPtr;
	void* nodePtr;
	BreakPoint(){
		position=1;
		weight=0;
		tscore=0;
		segmentLength=0;
		MinSegLen=0;
		totalLength=0;
		leftBreakPointPtr=NULL;
		rightBreakPointPtr=NULL;
		nodePtr = NULL;
	}
	BreakPoint(long _position, double _weight, double _tscore, long _segmentLength, long _MinSegLen, double _minTScore, long _totalLength):
		position(_position),weight(_weight), tscore(_tscore), segmentLength(_segmentLength), MinSegLen(_MinSegLen), minTScore(_minTScore), totalLength(_totalLength){
		leftBreakPointPtr=NULL;
		rightBreakPointPtr=NULL;
		nodePtr = NULL;
	};
	~BreakPoint(){};
	void setLeftBreakPoint(BreakPoint* bpPtr){
		this->leftBreakPointPtr=bpPtr;
	}
	BreakPoint* getLeftBreakPoint(){
		return this->leftBreakPointPtr;
	}
	void setRightBreakPoint(BreakPoint* bpPtr){
		this->rightBreakPointPtr=bpPtr;
	}
	BreakPoint* getRightBreakPoint(){
		return this->rightBreakPointPtr;
	}
	void setRBTreeNodePtr(void* ndPtr){
		this->nodePtr = ndPtr;
	}
	void* getRBTreeNodePtr(){
		return this->nodePtr;
	}
	BreakPointKey* getKeyPointer(){
		BreakPointKey* bpKeyPtr = new BreakPointKey(position, weight, tscore, segmentLength, MinSegLen, minTScore, totalLength);
		return bpKeyPtr;
	}
	BreakPointKey getKey(){
		BreakPointKey bpKey = BreakPointKey(position, weight, tscore, segmentLength, MinSegLen, minTScore, totalLength);
		return bpKey;
	}
	void removeItself(){
		/*
		 * updating left and right breakpoint
		 */
		double iC, iL, iR, totalLength_double;
		double h0;

		totalLength_double = (double) totalLength;
		//initialize
		BreakPoint* leftLeftBreakPointPtr = NULL;
		BreakPoint* rightRightBreakPointPtr = NULL;

		if (leftBreakPointPtr!=NULL && rightBreakPointPtr!=NULL){
			iL = (double) leftBreakPointPtr->position;
			iC = (double) this->position;
			iR = (double) rightBreakPointPtr->position;

			leftLeftBreakPointPtr = leftBreakPointPtr->leftBreakPointPtr;
			rightRightBreakPointPtr = rightBreakPointPtr->rightBreakPointPtr;
			//update the weights first
			if (leftLeftBreakPointPtr!=NULL){	//leftBreakPointPtr is NOT the left most.
				leftBreakPointPtr->weight = leftBreakPointPtr->weight
						+ sqrt((totalLength_double - iL) / (totalLength_double - iC) * iL / iC) * (iR - iC) / (iR - iL) *weight;
			}
			if (rightRightBreakPointPtr!=NULL){	//rightBreakPointPtr is NOT the right most break point
				rightBreakPointPtr->weight = rightBreakPointPtr->weight
						+ sqrt((totalLength_double - iR) / (totalLength_double - iC) * iR / iC) * (iC - iL) / (iR - iL) * weight;
			}
			//update the tscoe and segmentLength, which needs the updated weights
			if (leftLeftBreakPointPtr!=NULL){
				h0 = (double) (totalLength_double - leftBreakPointPtr->position) * (double) leftBreakPointPtr->position / totalLength_double * (double) (rightBreakPointPtr->position - leftLeftBreakPointPtr->position)
							/ (double) (rightBreakPointPtr->position - leftBreakPointPtr->position) / (double) (leftBreakPointPtr->position - leftLeftBreakPointPtr->position);
				leftBreakPointPtr->tscore = fabs(leftBreakPointPtr->weight) / sqrt(h0);
				leftBreakPointPtr->segmentLength = min(leftBreakPointPtr->position-leftLeftBreakPointPtr->position ,
						rightBreakPointPtr->position-leftBreakPointPtr->position);
			}
			if (rightRightBreakPointPtr!=NULL){
				h0 = (double) (totalLength_double - rightBreakPointPtr->position) * (double) rightBreakPointPtr->position / totalLength_double * (double) (rightRightBreakPointPtr->position - leftBreakPointPtr->position)
							/ (double) (rightRightBreakPointPtr->position - rightBreakPointPtr->position) / (double) (rightBreakPointPtr->position - leftBreakPointPtr->position);
				rightBreakPointPtr->tscore = fabs(rightBreakPointPtr->weight) / sqrt(h0);
				rightBreakPointPtr->segmentLength= min(rightRightBreakPointPtr->position-rightBreakPointPtr->position,
						rightBreakPointPtr->position-leftBreakPointPtr->position);
			}
			/*
			for (j = start; j <= end; j++) {
				h0 = (double) (totalLength_double - Iext[j]) * (double) Iext[j] / totalLength_double * (double) (Iext[j + 1] - Iext[j - 1])
						/ (double) (Iext[j + 1] - Iext[j]) / (double) (Iext[j] - Iext[j - 1]);
				Scores[j] = fabs(Wext[j]) / sqrt(h0);
			}
			*/
		}
		//update the left & right of the left & right
		if (leftBreakPointPtr!=NULL){
			leftBreakPointPtr->rightBreakPointPtr = rightBreakPointPtr ;
		}
		if (rightBreakPointPtr!=NULL){
			rightBreakPointPtr->leftBreakPointPtr = leftBreakPointPtr;
		}
		//release the memory
		//delete leftBreakPointPtr;
		//delete rightBreakPointPtr;
		//delete nodePtr;
	}
	//define methods to compute /update score/length/weight/position based on neighboring breakpoints
};



namespace std {
	template<> struct hash<BreakPoint> {
		/*
		 * 2013.09.22 hash function for BreakPoint, this requires g++ flag "-std=c++0x"
		 */
		size_t operator()(const BreakPoint& bp) const {

			// Start with a hash value of 0    .
			size_t seed = 0;

			// Modify 'seed' by XORing and bit-shifting in
			// one member of 'Key' after the other:
			boost::hash_combine(seed, boost::hash_value(bp.position));
			//boost::hash_combine(seed, boost::hash_value(bp.segmentLength));
			return seed;
		}
	};

};

typedef set<BreakPoint*> rbNodeDataType;
typedef RedBlackTree<BreakPointKey, rbNodeDataType > treeType;
typedef RedBlackTreeNode<BreakPointKey, rbNodeDataType > rbNodeType;

class BaseGADA{


public:
	double *Wext; //IO Breakpoint weights in extended notation...
	long *Iext; //IO Breakpoint positions in extended notation...
	double *_tscore_array;
	long *pK; //IO Number breakpoint positions remaining.
	double T; //IP  Threshold to prune
	long MinSegLen;	//minimum segment length

	long debug; //verbosity... set equal to 1 to see messages of SBLandBE(). 0 to not see them
	int report;
	long _M_total_length;	// total/max length

	long i;
	long K;
	double *normalized_data_array;	//would store normalized array of input data (raw-mean)
	double *inputDataArray;
	long *SegLen;
	double *SegAmp;
	double *SegState;
	double delta;
	long numEMsteps;
	long noOfBreakpointsAfterSBL;

	double *_alpha_array, *_aux_array;

	//double T2=5.0; //Segment collapse to base Non-Alteration level threshold
	double BaseAmp; //Base-level
	double a; //SBL hyperprior parameter
	double sigma2; //Variance observed, if negative value, it will be estimated by the mean of the differences
				   // I would recommend to be estimated on all the chromosomes and as a trimmed mean.
	long SelectClassifySegments; //Classify segment into altered state (1), otherwise 0
	long SelectEstimateBaseAmp; //Estimate Neutral hybridization amplitude.

	long maxNoOfIterations;	//=50000, //10000 is enough usually
	double convergenceDelta;	//1E-10 or 1E-8 seems to work well for this parameter. -- => ++ conv time
			//1E8 better than 1E10 seems to work well for this parameter. -- => -- conv time
	double convergenceMaxAlpha;	// 1E8 Maximum number of iterations to reach convergence...
	double convergenceB;	// a number related to convergence = 1E-20
	double ymean;	//mean of inputDataArray
	int reportIntervalDuringBE;	// how often to report progress during backward elimination, default is 100K

	BaseGADA(double* _inputDataArray, long _M, double _sigma2, double _BaseAmp, double _a, double _T, long _MinSegLen,
			long _debug , double _convergenceDelta,
			long _maxNoOfIterations, double _convergenceMaxAlpha, double _convergenceB, int _reportIntervalDuringBE):
			inputDataArray(_inputDataArray), _M_total_length(_M),  sigma2(_sigma2), BaseAmp(_BaseAmp), a(_a), T(_T), MinSegLen(_MinSegLen),
			debug(_debug),
			convergenceDelta(_convergenceDelta),
			maxNoOfIterations (_maxNoOfIterations), convergenceMaxAlpha(_convergenceMaxAlpha),
			convergenceB(_convergenceB), reportIntervalDuringBE(_reportIntervalDuringBE){
		noOfBreakpointsAfterSBL = 0;
	}
	~BaseGADA(){
		//free(SegLen);
		//free(SegAmp);
		//free(SegState);
	}
	void reconstruct(double *wr, long M_total_length, double *aux_vec);
	void BubbleSort(long *I, long L);
	void doubleBubbleSort(double *D, long *I, long L);
	void TrisolveREG(double *t0, double *tu, double *tl, double *coef, double *sol,
			long sizeh0);
	void DiagOfTriXTri(double *ll, double *l0, double *lu, double *rl, double *r0,
			double *ru, double *d, long N);
	void tridiagofinverse(double *t0, double *tl, double *itl, double *it0,
			double *itu, long N, double *d, double *e);
	void ForwardElimination(double *A, long N);
	void BackSubstitution(double *A, long N);
	void BackwardElimination(double *A, long N);
	void TriSolveINV(double *AA, long N, double *x, double *d, double *e);
	void ComputeH(double *h0, double *h1, long M_total_length);
	void ComputeFdualXb(long M_total_length, double *convergenceB);
	// 20080119 REMOVED void ComputeHs(long *s,double *a,long M_total_length,long Ms,double *h0,double *h1);
	void ComputeHs(long *s, long M_total_length, long Ms, double *h0, double *h1);
	void TriSymGaxpy(double *t0, double *t1, double *x, long M_total_length, double *y);
	void ComputeT(double *h0, double *h1, long M_total_length, double *alfa, double sigma,
			double *t0, double *tl, double *tu);
	long findminus(double *alpha_array, long Ms, double convergenceMaxAlpha, long *sel);
	long simpletresholding(double *inputvector, long N, double thres, double *disc);
	void computesegmentmeans(double *inputvector, long N, double *disc,
			long numdisc, double *amp);
	void reconstructoutput(double *rec, long N, double *disc, long numdisc,
			double *amp);
	long SBL(double *y, //I -- 1D array with the input signal
			long *I, //IO -- 1D array with the initial (final) candidate breakpoints
			double *alpha_array, //I -- 1D array with the initial (final) hyperparameter inv. varainces.
			double *w, //O -- 1D array containing the breakpoint weigths or posterior mean.
			double *sigw, //O -- Posterior variances, I would not really need them since they can be computed from alpha and H
			long M_total_length, //Initial size of the array in y
			long *K, //Size of the I alpha w

			//Algorithm parameters:
			double sigma2, //Noise estimated
			double a, //
			double convergenceB, double convergenceMaxAlpha, //Basis reduction parameter
			long maxNoOfIterations, //Max number of iterations
			double convergenceDelta, //Tolerance for convergence
			long debug //verbosity... set equal to 1 to see messages  0 to not see them
			);

	long BEthresh(
			//To eliminate...
			double *Scores, long Nscores, double *wr, long *indsel,
			long *pointNumRem, double *pointTau);
	//Returns breakpoint list lenght.
	long SBLandBE();

	void Project(double *y, long M_total_length, long *I, long L, double *xI, double *wI);
	void IextToSegLen();
	void IextWextToSegAmp();
	void CompZ( // computes z=F'y for entire possilbe breakpoint positions (normalized PWC)
			double *y, double *z, long M_total_length);
	void ComputeHsIext(
			long *Iext, // Indices of selection,
			long K, // Length of the indices,
			double *h0, // Returning diagonal of H,
			double *h1 // Returning upper diagonal of H
			);
	void ProjectCoeff( //IextYobs2Wext
			double *y, long M_total_length, long *Iext, long K, double *Wext);
	void CollapseAmpTtest( //Uses a T test to decide which segments collapse to neutral
			);

	double // Returns BaseAmp corresponding to the base level.
	CompBaseAmpMedianMethod();

	void ClassifySegments(double *SegAmp, long *SegLen, double *SegState, long K,
			double BaseAmp, double ploidy, double sigma2, //Reference noise
			double T //Critical value that decides when to colapse
			);

	void ComputeTScores(const double *Wext, const long *Iext, double *Scores, long K,
			long start, long end);

	long BEwTscore(double *Wext, //IO Breakpoint weights extended notation...
			long *Iext, //IO Breakpoint positions in extended notation...
			double *tscore_array, long *pK, //IO Number breakpoint positions remaining.
			double T, //IP  Threshold to prune
			long MinSegLen=0,	//minimum segment length
			long debug=0
			);

	long BEwTandMinLen( //Returns breakpoint list lenght. with T and MinSegLen
			double *Wext, //IO Breakpoint weights extended notation...
			long *Iext, //IO Breakpoint positions in extended notation...
			long *pK, //IO Number breakpoint positions remaining.
			double sigma2, //IP If sigma2
			double T, //IP  Threshold to prune,  T=T*sqrt(sigma2);
			long MinSegLen, //IP Minimum length of the segment.
			long debug
			);
	long RemoveBreakpoint(double *Wext, long *Iext, double *tscore_array, long K, long indexOfSegmentToRemove);

};

#endif //_BaseGADA_H_
