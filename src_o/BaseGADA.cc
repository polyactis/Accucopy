/*=================================================================
 * BaseGADA.c
 * All basic functions for SBL, BE, etc.
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
 Jordi Monso-Varona

 */
#include "BaseGADA.h"


//define how to output rbNodeDataType
std::ostream &operator<<(std::ostream &stream, const rbNodeDataType &b) {
    return stream << boost::format("%1%-elements-in-this-set")% b.size();
}


/* General Notation definition
 K		= number of breakpoints, segments is K+1
 M_total_length		= number of probes in chromosome / unit of analysis (can be a chr arm)
 I Ir Ired	Breakpoint set, double array I[0] contains first breakpoint I[K-1] contains Kth breakpoint length 0 if no breakpoints
 Ie Iext		Breakpoint set, extended notation I[0] contains 0 I[K] contains Kth breakpoint I[K+1]=M_total_length
 We		Aligned with Ie We[0] contains the overall mean and We[K+1] nothing
 SegAmp  Segment amplitudes length K+1
 SegLen  Segment longitudes length K+1 and should sum up to M_total_length?
 */

void BaseGADA::IextToSegLen() {
	SegLen = (long*) calloc(K + 1, sizeof(long));
	SegAmp = (double *) calloc(K + 1, sizeof(double));
	long k;
	for (k = 1; k <= K + 1; k++)
		SegLen[k - 1] = Iext[k] - Iext[k - 1];
}

void BaseGADA::IextWextToSegAmp() {
	long k;
	double M_total_length, TotalMean, AuxMean;
	M_total_length = Iext[K + 1];

	TotalMean = Wext[0];
	SegAmp[0] = 0;
	for (k = 1; k <= K; k++)
		SegAmp[k] = Wext[k]
				/ sqrt((double) (M_total_length - Iext[k]) * (double) Iext[k] / M_total_length)
				+ SegAmp[k - 1];
	AuxMean = 0;
	for (k = 0; k <= K; k++)
		AuxMean = AuxMean + SegAmp[k] * (double) (Iext[k + 1] - Iext[k]);
	AuxMean = AuxMean / M_total_length;
	for (k = 0; k <= K; k++)
		SegAmp[k] = SegAmp[k] - AuxMean + TotalMean;

}

// Implements: Reconstruct.M_total_length -- Recontructs x from w with zero mean
void BaseGADA::reconstruct(double *wr, long M_total_length, double *aux_vec) {
	long i = 0;
	double aux_double = 0;

	for (i = 0; i < M_total_length; i++) {
		aux_vec[i] = 0;
	}
	for (i = 0; i < M_total_length; i++) {
		aux_vec[i] = (-1)
				* sqrt(
						((double) (M_total_length - i) * (double) (i + 1))
								/ (double) (M_total_length + 1));
	}
	for (i = 0; i < M_total_length; i++) {
		wr[i] = (wr[i] / aux_vec[i]);
	}
	for (i = 0; i < M_total_length - 1; i++) {
		wr[M_total_length - (i + 2)] = wr[M_total_length - (i + 2)] + wr[M_total_length - (i + 1)];
	}
	aux_double = 0;
	for (i = 0; i < M_total_length; i++) {
		aux_double = aux_double + wr[i];
	}
	aux_double = (aux_double / M_total_length);
	for (i = 0; i < M_total_length; i++) {
		wr[i] = wr[i] - aux_double; // I have y mean in aux_double
	}
}

/* Bubble Sort */
void BaseGADA::BubbleSort(long *I, long L) {
	long i = 0;
	long aux_local = 0;

	for (i = 0; i < L - 1; i++) {
		if (I[i] > I[i + 1]) {
			aux_local = I[i + 1];
			I[i + 1] = I[i];
			I[i] = aux_local;
			i = -1;
		}
	}
}
void BaseGADA::doubleBubbleSort(double *D, long *I, long L) {
	long i = 0;
	double Daux = 0;
	long Iaux = 0;

	for (i = 0; i < L - 1; i++) {
		if (D[i] > D[i + 1]) {
			Daux = D[i + 1];
			D[i + 1] = D[i];
			D[i] = Daux;

			if (I != NULL) {
				Iaux = I[i + 1];
				I[i + 1] = I[i];
				I[i] = Iaux;
			}

			i = -1;
		}
	}
}

/////////////  TRIDIAGONAL MATRIX OPERATION FUNCTIONS ////////////////////
void BaseGADA::TrisolveREG(
		//input variables:
		double *t0, double *tu, double *tl, double *coef, double *sol,
		long sizeh0) {
	double M_total_length = 0;
	long i = 0;

	for (i = 0; i < sizeh0 - 1; i++) {
		//forward elimination
		M_total_length = tl[i] / t0[i];
		t0[i + 1] = t0[i + 1] - M_total_length * tu[i];
		coef[i + 1] = coef[i + 1] - coef[i] * M_total_length;
	}
	sol[sizeh0 - 1] = coef[sizeh0 - 1] / t0[sizeh0 - 1];
	for (i = 0; i < sizeh0 - 1; i++) {
		sol[sizeh0 - 2 - i] = (coef[sizeh0 - 2 - i]
				- (tu[sizeh0 - 2 - i] * sol[sizeh0 - 1 - i]))
				/ t0[sizeh0 - 2 - i];
	}
}

// DiagOfTriXTri  Diagonal of L*R where L,R tridiagonal (MexDiagOfTriXTri)
void BaseGADA::DiagOfTriXTri(
//Input variables:
		double *ll, // ll Lower diagonal
		double *l0, // l0 Central diagonal
		double *lu, // lu Upper diagonal
		double *rl, // rl Lower diagonal
		double *r0, // r0 Central diagonal
		double *ru, // ru Upper diagonal
		double *d, // d=DiagOfTriXTri(ll,l0,lu,rl,r0,ru)
		long N //Number of variables or length of l0,
		) {
	long i;

	d[0] = l0[0] * r0[0] + lu[0] * rl[0];
	for (i = 1; i < N - 1; i++) {
		d[i] = ll[i - 1] * ru[i - 1] + l0[i] * r0[i] + lu[i] * rl[i];
	}
	d[N - 1] = ll[N - 2] * ru[N - 2] + l0[N - 1] * r0[N - 1];

}

/* Tridiagonal de la inversa imput parameters: (tl,t0,d,e)
 *
 */
void BaseGADA::tridiagofinverse(
//Input variables:
		double *t0, double *tl, double *itl, double *it0, double *itu,

		long N, //Number of variables
		double *d, //Array cointaining the upper diagonal of pseudo-inverse
		double *e //Array cointaining the lower diagonal of pseudo-inverse
		) {
	long i;
	it0[0] = 1 / (t0[0] * (1 - e[0] * d[0]));

	for (i = 1; i < (N - 1); i++) {
		it0[i] = 1 / ((t0[i] - (tl[i - 1] * d[i - 1])) * (1 - e[i] * d[i]));
	}
	it0[N - 1] = 1 / (t0[N - 1] - tl[N - 2] * d[N - 2]);

	for (i = 0; i < (N - 1); i++) {
		itl[i] = -e[i] * it0[i];
		itu[i] = -d[i] * it0[i + 1];
	}
}

//////////               MEXTRISOLVE            ////////////////////
void BaseGADA::ForwardElimination(
//Input variables:
		double *A, //2D Array containing the [tu';t0';tl';b']
				   // Retuns [d';t0';tl';z']
		long N //Number of variables
		) {
	long i = 0;
	double *d, *tl, *t0, *z, *b;
	double M_total_length, m2;

	d = &A[0];
	t0 = &A[1];
	tl = &A[2];
	b = &A[3];
	z = b;

	//First iteration
	//%normalization of the first row
	m2 = *t0;
	*d /= m2; //d=d/M_total_length
	*b /= m2; //b=b/M_total_length
	//Pointer relocation
	b = b + 4;
	t0 = t0 + 4;

	for (i = 1; i < (N - 1); i++) {
		M_total_length = *tl;
		tl += 4;
		m2 = (*t0) - (*d) * M_total_length;
		t0 += 4;
		d += 4;
		*d /= m2; //tu=tu/m2
		*b = (*b - M_total_length * (*z)) / m2;
		b += 4;
		z += 4;
	}

	//Last step
	M_total_length = *tl;
	m2 = (*t0) - (*d) * M_total_length;
	*b = (*b - M_total_length * (*z)) / m2;
}

void BaseGADA::BackSubstitution(
		double *A, //2D Array containing the [d';!t0';!tl';z'] from Forward Elimination
				   // Retuns [d';t0';tl';x']

		long N //Number of variables
		) {
	long i = 0;
	double *z, *x;

	// +4 jumps pointer to the element in next column of A
	// Pointer initialization
	z = &A[N * 4 - 1];
	x = z - 4;

	for (i = 0; i < (N - 1); i++) {
		x[0] -= x[4] * x[-3];
		x -= 4;
	}

}
void BaseGADA::BackwardElimination(
		double *A, //2D Array containing the [tu';t0';tl';!b']
				   // Retuns [tu';t0';e';b']
		long N //Number of variables
		) {
	long i = 0;
	double *e;
	double *t0;
	e = &A[4 * N - 6];
	t0 = &A[4 * N - 3];
	*e /= *t0;
	t0 -= 4;

	for (i = 1; i < (N - 1); i++) {
		t0[-3] /= t0[0] - (t0[-1]) * (t0[+1]);
		t0 -= 4;
	}

}

void BaseGADA::TriSolveINV(double *AA, long N/*columns*/, double *x,
		double *d, double *e) {
	long i;

	double aux_local;

	/* Save backup copy of tu d*/
	for (i = 0; i < N - 1; i++) {
		d[i] = AA[(i << 2)];
//        x[i]=A[(i<<2)+3];
	}

	/*Run algorithm*/
	ForwardElimination(AA, N);
	BackSubstitution(AA, N);

	/*Store solution and Recompose A*/
	for (i = 0; i < N - 1; i++) {
		aux_local = AA[(i << 2)];
		AA[(i << 2)] = d[i];
		d[i] = aux_local;

		x[i] = AA[(i << 2) + 3];
	}
	x[i] = AA[(i << 2) + 3]; //x[N-1]=A[N<<2-1];

	BackwardElimination(AA, N);

	/*Store solution on the output vectors*/
	for (i = 0; i < N - 1; i++) {
		e[i] = AA[(i << 2) + 2];
	}
}

////////////////////////////////////////////////////////////////////////

void BaseGADA::ComputeFdualXb(
//imput variables:
		long M_total_length, double *b) {
	long i = 0;
	double a = 0;

	/* diff */
	for (i = 0; i < M_total_length - 1; i++) {
		b[i] = b[i + 1] - b[i];

	}
	b[M_total_length - 1] = 0;
	for (i = 0; i < M_total_length - 1; i++) {
		a = (double) (M_total_length - 1 - i) * (double) (i + 1) / (double) M_total_length;
		b[i] = b[i] * sqrt(a);
	}
}

void BaseGADA::CompZ( // computes z=F'y for entire possilbe breakpoint positions (normalized PWC)
		double *y, double *z, long M_total_length) {
	double LeftSum;
	double RightSum;
	double *aux_local = NULL;

	aux_local = (double*) calloc(M_total_length-1,sizeof(double));

	LeftSum = 0;
	RightSum = 0;
	for (M_total_length = 1; M_total_length < M_total_length; M_total_length++) {
		LeftSum = LeftSum + y[M_total_length - 1];
		aux_local[M_total_length - 1] = (-1) * sqrt((double) (M_total_length - M_total_length) / ((double) M_total_length * (double) M_total_length))
				* LeftSum;
	}
	for (M_total_length = M_total_length - 1; M_total_length >= 1; M_total_length--) {
		RightSum = RightSum + y[M_total_length];
		z[M_total_length - 1] = aux_local[M_total_length - 1]
				+ sqrt((double) M_total_length / ((double) (M_total_length - M_total_length) * (double) (M_total_length)))
						* RightSum;
	}
	free(aux_local);
}

/* Compute H   function [h0,h1]=CompH(dim) */
void BaseGADA::ComputeH(
//Input variables:
		double *h0, double *h1, long M_total_length //Number of variables
		) {
	long i;

	for (i = 0; i < M_total_length - 1; i++) {
		h0[i] = ((double) (M_total_length - 1 - i) * (double) (i + 1)) / (double) M_total_length;
	}
	for (i = 0; i < M_total_length - 2; i++) {
		h1[i] = -sqrt((h0[i + 1] * h0[i]));
	}
	for (i = 0; i < M_total_length - 1; i++) {
		h0[i] = 2 * h0[i];
	}
}

// Computes the H at the vector of indices...
void BaseGADA::ComputeHs(
//input variables:
		long *s, // Indices of selection,
		long M_total_length, // Length of the chormosome,
		long K, // Length of the indices,
		double *h0, // Returning diagonal of H,
		double *h1 // Returning upper diagonal of H
		) {
	long i;
	double iC, iL, iR;
	//double M_total_length;
	//M_total_length=(double)MM;

	i = 0;

	if (K == 1) {
		iL = 0;
		iC = (double) (s[i] + 1);
		iR = M_total_length;
		h0[i] = ((M_total_length - iC) * iC / M_total_length) * (iR - iL) / ((iR - iC) * ((iC - iL)));
		//h1[i]=0;
	} else {
		iL = 0;
		iC = (double) (s[i] + 1);
		iR = (double) (s[i + 1] + 1);
		h0[i] = ((M_total_length - iC) * iC / M_total_length) * (iR - iL) / ((iR - iC) * ((iC - iL)));
		h1[i] = -sqrt((M_total_length - iC) * iC * (M_total_length - iR) * iR) / (M_total_length * (iR - iC));
		for (i = 1; i < K - 1; i++) {
			//        h0[i]=((double)((M_total_length-s[i]-1)*(s[i]+1))/(double)M_total_length)*(double)(s[i+1]-s[i-1])/(double)((s[i+1]-s[i])*((s[i]-s[i-1])));
			//        h1[i]=sqrt((double)(((M_total_length-s[i]-1)*(s[i]+1))*(M_total_length-s[i+1]-1)*(s[i+1]+1)))/(double)(M_total_length*(s[i+1]-s[i]));
			iL = iC;
			iC = iR;
			iR = (double) (s[i + 1] + 1);
			h0[i] = ((M_total_length - iC) * iC / M_total_length) * (iR - iL) / ((iR - iC) * ((iC - iL)));
			h1[i] = -sqrt((M_total_length - iC) * iC * (M_total_length - iR) * iR) / (M_total_length * (iR - iC));
		}
		//i=K-1
		iL = iC;
		iC = iR;
		iR = M_total_length;
		h0[i] = ((M_total_length - iC) * iC / M_total_length) * (iR - iL) / ((iR - iC) * ((iC - iL)));
		//    h0[i]=((double)((M_total_length-s[i]-1)*(s[i]+1))/(double)M_total_length)*(double)(M_total_length-s[i-1]-1)/(double)((M_total_length-s[i]-1)*((s[i]-s[i-1]))); //s[K]=M_total_length;
	}

}

void BaseGADA::ComputeHsIext(
//input variables:
		long *Iext, // Indices of selection,
		long K, // Length of the indices,
		double *h0, // Returning diagonal of H,
		double *h1 // Returning upper diagonal of H
		) {
	long i, M_total_length;
	double iC, iL, iR;
	//double M_total_length;
	//M_total_length=(double)MM;

	M_total_length = Iext[K + 1];

	//iL=0;iC=(double)Iext[1];iR=(double)s[2];
	for (i = 1; i < K; i++) {
		//iL=iC;iC=iR;iR=(double)(Iext[i+1]);
		iL = Iext[i - 1];
		iC = Iext[i];
		iR = (double) Iext[i + 1];
		h0[i - 1] = ((M_total_length - iC) * iC / M_total_length) * (iR - iL) / ((iR - iC) * ((iC - iL)));
		h1[i - 1] = -sqrt((M_total_length - iC) * iC * (M_total_length - iR) * iR) / (M_total_length * (iR - iC));
	}
	//i=K
	//iL=iC;iC=iR;iR=M_total_length;
	iL = Iext[i - 1];
	iC = Iext[i];
	iR = (double) Iext[i + 1];
	h0[i - 1] = ((M_total_length - iC) * iC / M_total_length) * (iR - iL) / ((iR - iC) * ((iC - iL)));
}

void BaseGADA::TriSymGaxpy(
//input variables:
		double *t0, double *t1, double *x, long M_total_length, double *y

		) {
	long i;

	if (M_total_length == 1) {
		y[0] = t0[0] * x[0];
	} else {
		y[0] = t0[0] * x[0] + t1[0] * x[1];
		for (i = 1; i < M_total_length - 1; i++) {
			y[i] = t0[i] * x[i] + t1[i] * x[i + 1] + t1[i - 1] * x[i - 1];
		}
		y[M_total_length - 1] = t0[M_total_length - 1] * x[M_total_length - 1] + t1[M_total_length - 2] * x[M_total_length - 2];
	}
}

void BaseGADA::ComputeT(double *h0, double *h1, long M_total_length, double *alfa, double sigma, /*pass 1 if scale not available*/
double *t0, double *tl, double *tu
/* in theory there is a parameter called 'scale' but we don't use it */
) {
	long i = 0;

	for (i = 0; i < M_total_length; i++) {
		t0[i] = (h0[i] * alfa[i] * sigma) + 1;
		if (i < M_total_length - 1) {
			tu[i] = h1[i] * alfa[i + 1] * sigma;
			tl[i] = h1[i] * alfa[i] * sigma;
		}
	}
}

/*allocates all those lower-than-convergenceMaxAlpha alpha_array's INDEXES in sel vector  */
long BaseGADA::findminus(
		double *alpha_array, long K, double convergenceMaxAlpha, long *sel //index vector
		) {
	long i, j = 0;
	for (i = 0; i < K; i++) {
		if (alpha_array[i] < convergenceMaxAlpha) {
			sel[j] = i;
			j++;
		}
	}
	return (j);
}

/* simpletresholding
 * Applies a simple tresholding algorithm to find the discontinuities
 */
long //Returns the number of discontinuities
BaseGADA::simpletresholding(
//Input variables:
		double *inputvector, //1D Array containing the input vector
		long N, //Vector length
		double thres, //Treshold value
		//Output variables:
		double *disc //1D empty array, with memory already allocated for finding up to N discontinuities positions
		) {
	long i = 0;
	long numdisc = 0; //Number of discontinuities
	long state = 0; //To keep up the state 0(<tresh),1(>tresh)

	//First step: Find discontinuities
	if (inputvector[0] >= thres)
		state = 1;
	for (i = 0; i < N; i++) {
		if (state == 0) {
			if (inputvector[i] >= thres) {
				disc[numdisc] = (double) (i + 1); //Matlab uses 1..N indices instead of 0..(N-1)
				numdisc++;
				state = 1;
//                printf("i%ld numdisc%ld State change state%ld.\n",i,numdisc,state);
			}
		} else {
			if (inputvector[i] < thres) {
				disc[numdisc] = (double) (i + 1); //Matlab uses 1..N indices instead of 0..(N-1)
				numdisc++;
				state = 0;
				//              printf("i%ld numdisc%ld State change state%ld.\n",i,numdisc,state);
			}
		}
	}

	return numdisc;
}

/* computesegmentmeans
 * Computed segment means.
 */
void BaseGADA::computesegmentmeans(
//Input variables:
		double *inputvector, //1D Array containing the input vector
		long N, //Vector length
		double *disc, //1D empty array, with memory already allocated for finding up to N discontinuities positions
		long numdisc, //Length of the discontinuity vector
		//Output variables:
		double *amp //1D empty array, containing the amplitudes of the segments between discotinuities,
		) {
	long i, j;
	long prevdisc = 0;

	for (i = 0; i < numdisc; i++) {
		amp[i] = 0;
		for (j = prevdisc; j < (long) (disc[i] - 1); j++)
			amp[i] += inputvector[j];
		amp[i] = amp[i] / ((double) (disc[i] - 1 - prevdisc));
		prevdisc = (long) (disc[i] - 1);
	}
	amp[i] = 0;
	for (j = prevdisc; j < N; j++)
		amp[i] += inputvector[j];
	amp[i] = amp[i] / (double) (N - prevdisc);
}

/* reconstructoutput
 * Reconstruct signal from the discontinuity location, and segment
 * amplitudes
 */
void BaseGADA::reconstructoutput(
//Output variables:
		double *rec, //1D Array returning the reconstructed vector
		//Input variable:
		long N, //Vector length
		double *disc, //1D empty array, with memory already allocated for finding up to N discontinuities positions
		long numdisc, //Length of the discontinuity vector
		//Output variables:
		double *amp //1D empty array, containing the amplitudes of the segments between discotinuities,
		) {
	long i, j;

	j = 0;
	for (i = 0; i < numdisc; i++)
		for (; j < (long) (disc[i] - 1); j++)
			rec[j] = amp[i];
	//If no discontinuities or from the last to the end
	for (; j < N; j++)
		rec[j] = amp[i];
}

/* +++++++++++============================================++++++++++++++++++++++++++++++************************************** */

/* SBL function
 * returns number of EM iterations
 */
long BaseGADA::SBL(double *y, //I -- 1D array with the input signal
		long *I, //IO -- 1D array with the initial (final) candidate breakpoints
		double *alpha_array, //I -- 1D array with the initial (final) hyperparameter inv. varainces.
		double *w, //O -- 1D array containing the breakpoint weigths or posterior mean.
		double *sigw, //O -- Posterior variances, I would not really need them since they can be computed from alpha_array and H
		long M_total_length, //Initial size of the array in y
		long *pK, //Size of the I alpha_array w

		//Algorithm parameters:
		double sigma2, //Noise estimated
		double a, //
		double convergenceB, double convergenceMaxAlpha, //Basis reduction parameter
		long maxNoOfIterations, //Max number of iterations
		double convergenceDelta, //Tolerance for convergence
		long debug //verbosity... set equal to 1 to see messages  0 to not see them
		) {
	long n, i, K;
	long M0, sizesel;
	double *w0;
	double myaux;

	//Extra memory necessary to store variables during the algorithm.... (require memory initialization)
	double *h0 = NULL;
	double *h1 = NULL;
	double *xx = NULL;
	double *z = NULL;
	double *t0 = NULL;
	double *tl = NULL;
	double *tu = NULL;
	double *wpred = NULL;
	double *d = NULL;
	double *e = NULL;
	long *sel = NULL;
	double *AA = NULL;
	double *yy = NULL;

	K = *pK;
	M0 = M_total_length - 1;

	//non discontinuities case
	if (K == 0) {
		w[0] = -1;
		*pK = K;
		sigw[0] = -1;
		return 0;
	}
	//Memory initialization of the outputs. (Already with mem assigned)
//    w=calloc(K,sizeof(double));
//    sigw=calloc(K,sizeof(double));

	//Memory initialization (internal to be freed)
	yy = (double*) calloc(M_total_length,sizeof(double));
	t0 = (double*) calloc(M_total_length,sizeof(double));
	tl = (double*) calloc(M_total_length-1,sizeof(double));
	tu = (double*) calloc(M_total_length-1,sizeof(double));
	AA = (double*) calloc(4*M_total_length,sizeof(double));
	d = (double*) calloc(M0-1,sizeof(double));
	e = (double*) calloc(M0-1,sizeof(double));
	h0 = (double*) calloc(M0,sizeof(double));
	h1 = (double*) calloc(M0-1,sizeof(double));
	wpred = (double*) calloc(K,sizeof(double));
	sel = (long*) calloc(K,sizeof(long));
	z = (double*) calloc(M0,sizeof(double));
	xx = (double*) calloc(K,sizeof(double)); //myDoubleMAlloc(K);

	//Create a copy of the input
	for (i = 0; i < M_total_length; i++)
		yy[i] = y[i];

	//printf("\n\nOPERATIONS BEFORE LOOP:\n");

	// printf("\nCOMPUTE H\n");

	ComputeH(h0, h1, M_total_length);
	//printf("H0\n\nh0[0]:%g\nh0[1]:%g\nh0[2]:%g\nh0[3]:%g\nh0[%ld]:%g\n",h0[0],h0[1],h0[2],h0[3],M0-1,h0[M0-1]);
	//printf("H1\n\nh1[0]:%g\nh1[1]:%g\nh1[2]:%g\nh1[3]:%g\nh1[%ld]:%g\n",h1[0],h1[1],h1[2],h1[3],M0-2,h1[M0-2]);
	//printf("\nCOMPUTE F DUAL\n");

	ComputeFdualXb(M_total_length, yy); //checked

	w0 = yy; //w0 now has M_total_length-1 dimmension
	//printf("W0\n\nw0[0]:%g\nw0[1]:%g\nw0[2]:%g\nw0[3]:%g\nw0[%ld]:%g\n",w0[0],w0[1],w0[2],w0[3],M0-1,w0[M0-1]);
	for (i = 0; i < K; i++) {
		xx[i] = w0[i];
		t0[i] = h0[i];
	}

	//printf("H0\n\nh0[0]:%g\nh0[1]:%g\nh0[2]:%g\nh0[3]:%g\nh0[%ld]:%g\n",t0[0],t0[1],t0[2],t0[3],M0-1,t0[M0-1]);

	//printf("\nTRISOLVE REGULAR\n");
	TrisolveREG(t0, h1, h1, xx, z, M0); //checked
	//printf("Z\n\nz[0]:%g\nz[1]:%g\nz[2]:%g\nz[3]:%g\nz[%ld]:%g\n",z[0],z[1],z[2],z[3],M0-1,z[M0-1]);
	//printf("Z\n\nz[0]:%g\nz[1]:%g\nz[2]:%g\nz[3]:%g\nz[%ld]:%g\n",xx[0],xx[1],xx[2],xx[M0-2],M0-1,xx[M0-1]);

	//printf("\n\nSUBSELECTION STARTS\n");
	/* initialize if subselection */
	if (K < M0) {
		//printf("\nCOMPUTE HS SUBSELECTION\n");
		ComputeHs(I, M_total_length, K, h0, h1); //checked

		//printf("h0s\n\nh0[0]:%g\nh0[1]:%g\nh0[2]:%g\nh0[3]:%g\nh0[%ld]:%g\n",h0[0],h0[1],h0[2],h0[3],K-1,h0[K-1]);
		//printf("h1s\n\nh1[0]:%g\nh1[1]:%g\nh1[2]:%g\nh1[3]:%g\nh1[%ld]:%g\n",h1[0],h1[1],h1[2],h1[3],K-2,h1[K-2]);
		//z(I)
		//printf("I[0]:%ld\n",I[0]);

		for (i = 0; i < K; i++){
			t0[i] = z[I[i]];
			xx[i] = 0;
		}

		//for (i = 0; i < K; i++)

		//printf("\nTRISYMGAXPY SUBSELECTION\n");
		TriSymGaxpy(h0, h1, t0, K, xx); //checked
		w0 = xx; // error??? nono I think it's perfecly fine!!!
		//printf("W0 TRISYMGAXPY\n\nw0[0]:%g\nw0[1]:%g\nw0[2]:%g\nw0[3]:%g\nd[%ld]:%g\n",w0[0],w0[1],w0[2],w0[3],K-1,w0[K-1]);
	}
	//printf("\n\n*****************\nITERATION STARTS\n*****************\n\n");

	/********************************/
	/******* EM loop         ********/
	/********************************/

	for (n = 0; n < maxNoOfIterations; n++) {

		for (i = 0; i < K; i++) {
			wpred[i] = w[i];
		}
		ComputeT(h0, h1, K, alpha_array, sigma2, t0, tl, tu); //checked
//        if (n==0){
//        //printf("COMPUTE T FIRST ITERATION:\n");
//        //printf("T0\n\nt0[0]:%g\nt0[1]:%g\nt0[2]:%g\nt0[3]:%g\nt0[%ld]:%g\n",t0[0],t0[1],t0[2],t0[3],K-1,t0[K-1]);
//        //printf("TL\n\ntl[0]:%g\ntl[1]:%g\ntl[2]:%g\ntl[3]:%g\ntl[%ld]:%g\n",tl[0],tl[1],tl[2],tl[3],K-2,tl[K-2]);
//        //printf("TU\n\ntu[0]:%g\ntu[1]:%g\ntu[2]:%g\ntu[3]:%g\ntu[%ld]:%g\n",tu[0],tu[1],tu[2],tu[3],K-2,tu[K-2]);
//        }

		if (K == 1) {
			w[0] = w0[0] / t0[0];
			sigw[0] = 1 / t0[0] * h0[0] * sigma2;
		} else {
			for (i = 0; i < K; i++) {
				if (i != K - 1) {
					AA[(i << 2)] = tu[i];
					AA[(i << 2) + 2] = tl[i];
				}
				AA[(i << 2) + 1] = t0[i];
				AA[(i << 2) + 3] = w0[i];
			}

			TriSolveINV(AA, K, w, d, e); //checked
//            if (n==0){
//            printf("TRISOLVE INVERSE FIRST ITERATION:\n");
//            printf("W TRISOLVE\n\nw[0]:%g\nw[1]:%g\nw[2]:%g\nw[3]:%g\nw[%ld]:%g\n",w[0],w[1],w[2],w[3],K-1,w[K-1]);
//            printf("D TRISOLVE\n\nd[0]:%g\nd[1]:%g\nd[2]:%g\nd[3]:%g\nd[%ld]:%g\n",d[0],d[1],d[2],d[3],K-2,d[K-2]);
//            printf("E TRISOLVE\n\ne[0]:%g\ne[1]:%g\ne[2]:%g\ne[3]:%g\ne[%ld]:%g\n",e[0],e[1],e[2],e[3],K-2,e[K-2]);
//            }

			tridiagofinverse(t0, tl, tl, t0, tu, K, d, e); // itl it0 itu on the sames without i to save memory can be recycled.
			// now itl itu it0 cointain the tridiagonal of the inverse....
//            if (n==0){
//            printf("TRISOLVE FIRST ITERATION:\n");
//            printf("IT0\n\nit0[0]:%g\nit0[1]:%g\nit0[2]:%g\nit0[3]:%g\nit0[%ld]:%g\n",it0[0],it0[1],it0[2],it0[3],K-1,it0[K-1]);
//            printf("ITL\n\nitl[0]:%g\nitl[1]:%g\nitl[2]:%g\nitl[3]:%g\nitl[%ld]:%g\n",itl[0],itl[1],itl[2],itl[3],K-2,itl[K-2]);
//            printf("ITU\n\nitu[0]:%g\nitu[1]:%g\nitu[2]:%g\nitu[3]:%g\nitu[%ld]:%g\n",itu[0],itu[1],itu[2],itu[3],K-2,itu[K-2]);
//            }

			DiagOfTriXTri(tl, t0, tu, h1, h0, h1, sigw, K); //checked
//            if (n==0){
//            printf("TRIDIAGOFINVERSE FIRST ITERATION:\n");
//            printf("DIAG\n\ndiag[0]:%g\ndiag[1]:%g\ndiag[2]:%g\ndiag[3]:%g\ndiag[%ld]:%g\n",diag[0],diag[1],diag[2],diag[3],K-1,diag[K-1]);
//            printf("DIAGOFTRIXTRI ITERATION:%ld\n",n);
//            }

			for (i = 0; i < K; i++) {
				sigw[i] = sigma2 * sigw[i];
			}

		}
		// 2013.09.01 merged to below
		//for (i = 0; i < K; i++) {
		//	alpha_array[i] = (1 + 2 * a) / (w[i] * w[i] + sigw[i] + 2 * convergenceB);
		//}

//        if (n==0){
//            printf("ALPHA & DIAGSIGMA FIRST ITERATION:\n");
//            printf("\n\nMS:%ld\n\n",K);
//            printf("ALPHA\n\nalpha_array[0]:%g\nalpha_array[1]:%g\nalpha_array[2]:%g\nalpha_array[3]:%g\nalpha_array[%ld]:%g\n",alpha_array[0],alpha_array[1],alpha_array[2],alpha_array[3],K-1,alpha_array[K-1]);
//            printf("DIAGSIGMA\n\nd[0]:%g\ndiagsigma[1]:%g\ndiagdigma[2]:%g\ndiagsigma[3]:%g\ndiagsigma[%ld]:%g\n",sigw[0],sigw[1],sigw[2],sigw[3],K-1,sigw[K-1]);
//        }

		//euclidean norm of wpred-w  CRITERIUM of CONVERGENCE
		/*
		delta=0;
		for (i=0;i<K;i++){
		delta=delta+(wpred[i]-w[i])*(wpred[i]-w[i]);
		}
		delta=sqrt(delta);
		*/

		//max diff of wpred-w  CRITERIUM of CONVERGENCE
		delta = 0;
		for (i = 0; i < K; i++) {
			alpha_array[i] = (1 + 2 * a) / (w[i] * w[i] + sigw[i] + 2 * convergenceB);
			//myaux=abs(wpred[i]-w[i]);
			myaux = (wpred[i] - w[i]);
			if (myaux < 0)
				myaux = -myaux;
			if (myaux > delta)
				delta = myaux;
		}

//        if (n==0){
//            printf("EUCLIDEAN NORM FIRST ITERATION:\n");
//            printf("\n\neuclidean_norm:%g\n\n",delta);
//        }
		if (delta < convergenceDelta) {
			if (debug > 0) {
				std::cerr<< boost::format("# \t SBL: Converged after %1% iterations, delta=%4%, within tolerance %2%, M_total_length=%3% \n") %
						n % convergenceDelta % K % delta;
			}
			break;
		}

		sizesel = findminus(alpha_array, K, convergenceMaxAlpha, sel);
//        if (n==0){
//            printf("\n\nSIZESEL FIRST ITERATION:\n");
//            printf("SEL\n\nsel[0]:%ld\nsel[1]:%ld\nsel[2]:%ld\nsel[3]:%ld\nsel[%ld]:%ld\n",sel[0],sel[1],sel[2],sel[3],K-1,sel[K-1]);
//        }
		if (sizesel == 0) {
			K = 0;
			w[0] = -1;
			I[0] = -1;
			sigw[0] = -1;
			alpha_array[0] = -1;
			std::cerr << boost::format("#      SBL: After %1% iterations, No disconinuities found M_total_length=%2% \n") %
					n % K;
			break;
		}
		if (sizesel < K) {
			K = sizesel;

			//I(sel)
			for (i = 0; i < sizesel; i++){
				I[i] = I[sel[i]];
			}
			//printf("I[%ld]=%ld\n",i,I[i]);

			ComputeHs(I, M_total_length, K, h0, h1);

//            if (k==0){
//            printf("\n\nCOMPUTE Hs FIRST TIME REDUCTION:\n");
//            printf("h0s sizesel<K\n\nh0[0]:%g\nh0[1]:%g\nh0[2]:%g\nh0[3]:%g\nh0[%ld]:%g\n",h0[0],h0[1],h0[2],h0[3],K-1,h0[K-1]);
//            printf("h1s sizesel<K\n\nh1[0]:%g\nh1[1]:%g\nh1[2]:%g\nh1[3]:%g\nh1[%ld]:%g\n",h1[0],h1[1],h1[2],h1[3],K-2,h1[K-2]);
//            }

			//alpha_array(sel)

			//printf("I[0]:%ld\n",I[0]);
			for (i = 0; i < K; i++){
				alpha_array[i] = alpha_array[sel[i]];
				t0[i] = z[I[i]];
				xx[i] = 0;
			}


			TriSymGaxpy(h0, h1, t0, K, xx); //checked until here
			w0 = xx;
//            if (k==0){
//                printf("\n\nTRISYMGAXPY FIRST TIME REDUCTION:\n");
//                printf("W0 trisymgaxpy\n\nw0[0]:%g\nw0[1]:%g\nw0[2]:%g\nw0[3]:%g\nw0[%ld]:%g\n",w0[0],w0[1],w0[2],w0[3],K-1,w0[K-1]);
//            }

//            if (k==0){
//                printf("W FIRSt TIME REDUCTION\n\nw[0]:%g\nw[1]:%g\nw[2]:%g\nw[3]:%g\nw[%ld]:%g\n",w[0],w[1],w[2],w[3],K-1,w[K-1]);
//                printf("SEL FIRSt TIME REDUCTION\n\nsel[0]:%ld\nsel[1]:%ld\nsel[2]:%ld\nsel[3]:%ld\nsel[%ld]:%ld\n",sel[0],sel[1],sel[2],sel[3],K-1,sel[K-1]);
//            }
			for (i = 0; i < K; i++) {
				w[i] = w[sel[i]];
				wpred[i] = wpred[sel[i]];
			}
//             if (k==0){
//                printf("W FIRST TIME REDUCTION AFTER W[SEL]\n\nw[0]:%g\nw[1]:%g\nw[2]:%g\nw[3]:%g\nw[%ld]:%g\n",w[0],w[1],w[2],w[3],K-1,w[K-1]);
//                k=1;
//             }

		} //End if of column elimination
	} //End loop MAXIT
	/*********************/

	if (debug > 0) {
		if (i >= maxNoOfIterations) {
			std::cerr << boost::format("# \t SBL: Converged??? Stopped after %1% iterations with delta=%2%, M_total_length=%3%, convergenceDelta=%4% \n") %
					n % delta % K % convergenceDelta;
		}
	}

	if (debug){
		std::cerr <<  boost::format("_SBL_ Reffiting K=%1%\n")%K;
	}

	// Projection onto Breakpoint set... (REFITTING)
	if (K != 0) {

		//printf("I[0]:%ld\n",I[0]);
		for (i = 0; i < K; i++){
			t0[i] = z[I[i]];
			sigw[i] = h0[i];
			w[i] = 0;
		}
		TriSymGaxpy(h0, h1, t0, K, w);
	} else {
		h0[0] = -1;
	}

	//Memory freeing....
	free(yy);
	free(h0);
	free(h1);
	free(xx);
	free(z);
	free(t0);
	free(tl);
	free(tu);
	free(wpred);
	free(d);
	free(e);
	free(sel);
	free(AA);

	*pK = K;
	return n;
}

/**************************************************************************************************************************/

/* Depending on the Tau value, this function will behave in two different
 ways. First, if Tau is lower than zero (-1 will be our agreement) the
 algorithm will remove the less-important break-point left. On the other
 hand, if it is greater than zero, it will work as a regular threshold
 algorithm i.e.,it will remove all those break-point who don't reach the
 selected value TAU.*/
long BaseGADA::BEthresh(double *Scores, long Nscores, double *wr, long *indsel,
		long *pointNumRem, double *pointTau) {
	long M_total_length = Nscores + 1; //genome length
	double vmin = -1e100; //Forces to enter at least once for tau<0
	long imin = -1;
	long inrem = -1;
	long i = 0;
	long aux_local = 0;
	long iC = 0;
	long iR = 0;
	long iL = 0;
	double h0R = 0;
	double h0L = 0;
	double Tau = *pointTau;
	long NumRem = *pointNumRem;

#ifdef _DEBUG_
	//checked-> inputs are OK
	long n=0;
	printf("INPUT CHECKING:\n\nNscores=%ld\n\nSCORES[indsel[0]]=%g\nSCORES[indsel[1]]=%g\nSCORES[indsel[2]]=%g\nSCORES[indsel[%ld]]=%g\n",Nscores,Scores[indsel[0]],Scores[indsel[1]],Scores[indsel[2]],NumRem-1,Scores[indsel[NumRem-1]]);
	printf("\nwr[indsel[0]]=%g\nwr[indsel[1]]=%g\nwr[indsel[2]]=%g\nwr[indsel[%ld]]=%g\n\nindsel[0]=%ld\nindsel[1]=%ld\nindsel[2]=%ld\nindsel[%ld]=%ld\n",wr[indsel[0]],wr[indsel[1]],wr[indsel[2]],NumRem,wr[indsel[NumRem-1]],indsel[0],indsel[1],indsel[2],NumRem-1,indsel[NumRem-1]);
	printf("\nNumRem=%ld\n\nTau=%g\n\nNscores=%ld\n\n",NumRem,Tau,Nscores);

#endif

	while (NumRem > 0) {
#ifdef _DEBUG_
		n++;
		printf("ITERATION NUMBER:\n\nn=%ld\n\n",n);
#endif

		//search for the minimum score>0
		imin = 0;
		vmin = Scores[indsel[imin]];
		for (i = 1; i < NumRem; i++) //search of the lowest Score (less important)
				{
			aux_local = indsel[i];
			if (Scores[aux_local] < vmin) {
				imin = i;
				vmin = Scores[aux_local];
			}
		}
		inrem = indsel[imin]; // Returns the absolute position of the disc to remove.

#ifdef _DEBUG_
		printf("NEXT COMPONENT OUT:\n\nimin=%ld\nvmin=%g\n\n",imin,vmin);
#endif

		if ((vmin > Tau) && (Tau >= 0)) //if the lowest score is greater than TAU, we don't do anything
				{

#ifdef _DEBUG_
			printf("Threshold is not enough --> THERE IS NO MORE TO EXTRACT\n\n");
#endif

			break;
		}
		if (imin == -1) {

#ifdef _DEBUG_
			printf("No discontinuities left --> THERE IS NO MORE TO EXTRACT\n\n");
#endif

			//there is nothing left
			NumRem = 0;
			vmin = 0;
			break;
		}
		/*If I'M_total_length here then imin is what I want to remove... because there is still
		 *something to remove (non zero scores) and they are lower than my threshold*/

		iC = indsel[imin];

#ifdef _DEBUG_
		printf("IC Component out:\n\niC=%ld\n\n",iC);
#endif

		//only one left
		if (NumRem == 1) {
			aux_local = indsel[imin];
			Scores[aux_local] = 0;
			wr[aux_local] = 0;
			NumRem = 0;

#ifdef _DEBUG_
			printf("There is one left to extract:\n\nindsel[imin]=%ld\n\n", aux_local);
#endif

			break;
		}
		/*is it efficient to look for them, or will it be better to mantain an
		 indsel like in the matlab version??? Mantain indsel
		 Update left weight if it exist, and right weight if it exist.*/

		if (imin == 0) //if I should remove the leftmost w
				{
#ifdef _DEBUG_
			printf("LEFTMOST REMOVING:\n\n",imin,vmin);
#endif

			iR = indsel[imin + 1]; //indsel[1]
			iL = 0;
			//recompute weights
			wr[iR] = wr[iR]
					+ (double) sqrt(
							(double) (M_total_length - iR - 1) / (M_total_length - iC - 1)
									* (double) (iR + 1) / (iC + 1))
							* (double) (iC + 1 - iL) / (iR + 1 - iL)
							* (double) wr[iC];
			wr[iC] = 0;

#ifdef _DEBUG_
			printf("wr[iR]=%g:\niR=%ld\n",wr[iR],iR);
#endif

			//I don't modify wr[iL] coz there is no left breakpoint
			//recompute Scores
			Scores[iC] = 0;
			iC = indsel[imin + 1];
			if (imin + 1 == NumRem - 1) ////imin+1 or imin+2???????
					{
				/*if there are just two discontinuitues so now we will
				 work with the rightmost w (or Score)*/
				iR = M_total_length - 1;
			} else {
				iR = indsel[imin + 2];
			}
			h0R = (double) (M_total_length - iC - 1) * (iC + 1) / M_total_length * (iR + 1 - iL)
					/ (iR - iC) / (iC + 1 - iL);
			if (wr[iC] >= 0) {
				Scores[iC] = (double) (wr[iC] / sqrt(h0R));
			} else {
				Scores[iC] = (double) (-wr[iC] / sqrt(h0R));
			}

#ifdef _DEBUG_
			printf("h0R=%g\niC=%ld\nScores[iC]=%g\n\n",h0R,iC,Scores[iC]);
#endif

		} else if (imin == NumRem - 1) {
#ifdef _DEBUG_
			printf("RIGHTMOST REMOVING:\n\n",imin,vmin);
#endif

			iL = indsel[imin - 1];
			iR = M_total_length;
			//recompute weights
			wr[iL] = wr[iL]
					+ sqrt(
							(double) (M_total_length - iL - 1) / (M_total_length - iC - 1)
									* (double) (iL + 1) / (iC + 1))
							* (double) (iR - iC - 1) / (iR - iL - 1) * wr[iC];
			wr[iC] = 0;
#ifdef _DEBUG_
			printf("wr[iL]=%g:\niL=%ld\n",wr[iL],iL);
#endif

			//Scores
			Scores[iC] = 0;
			/*new iC is gonna be placed in imin-1 so look if is the most
			 left breakpoint*/
			if (imin - 1 == 0) {
				iL = -1;
			} else {
				iL = indsel[imin - 2];
			}
			iC = indsel[imin - 1];
			h0L = (double) (M_total_length - iC - 1) * (iC + 1) / M_total_length * (iR - iL - 1)
					/ (iR - iC - 1) / (iC - iL);
			if (wr[iC] >= 0) {
				Scores[iC] = (double) (wr[iC] / sqrt(h0L));
			} else {
				Scores[iC] = (double) (-wr[iC] / sqrt(h0L));
			}
#ifdef _DEBUG_
			printf("h0L=%g\niC=%ld\nScores[iC]=%g\n\n",h0L,iC,Scores[iC]);
#endif
		} else {
#ifdef _DEBUG_
			printf("LEFT&RIGHT REMOVING:\n\n",imin,vmin);
#endif
			//removing left and right
			iR = indsel[imin + 1];
			iL = indsel[imin - 1];
			wr[iR] = wr[iR]
					+ sqrt(
							(double) (M_total_length - iR - 1) / (M_total_length - iC - 1)
									* (double) (iR + 1) / (iC + 1))
							* (double) (iC - iL) / (iR - iL) * wr[iC];
			wr[iL] = wr[iL]
					+ sqrt(
							(double) (M_total_length - iL - 1) / (M_total_length - iC - 1)
									* (double) (iL + 1) / (iC + 1))
							* (double) (iR - iC) / (iR - iL) * wr[iC];
			wr[iC] = 0;
#ifdef _DEBUG_
			printf("M_total_length-iR-1=%ld:\nM_total_length-iC-1=%ld\niR+1/iC+1=%g:\nsqrt(all)=%g\n",M_total_length-iR-1,M_total_length-iC-1,(iR+1)/(iC+1),sqrt(((M_total_length-iR-1)/(M_total_length-iC-1))*((iR+1)/(iC+1))));
			printf("wr[iL]=%g:\niL=%ld\nwr[iR]=%g:\niR=%ld\n",wr[iL],iL,wr[iR],iR);
#endif

			//Scores
			Scores[iC] = 0;
			iL = indsel[imin - 1];
			iC = indsel[imin + 1];
			if (imin + 1 == NumRem - 1) //iC is now the most right point
					{
				iR = M_total_length - 1;
			} else {
				iR = indsel[imin + 2];
			}
			h0R = (double) (M_total_length - iC - 1) * (iC + 1) / M_total_length * (iR - iL) / (iR - iC)
					/ (iC - iL);
			if (wr[iC] >= 0) {
				Scores[iC] = (double) (wr[iC] / sqrt(h0R));
			} else {
				Scores[iC] = (double) (-wr[iC] / sqrt(h0R));
			}
#ifdef _DEBUG_
			printf("h0R=%g\niCfirst=%ld\nScores[iC]=%g\n\n",h0R,iC,Scores[iC]);
#endif
			//scores
			if (imin - 1 == 0) //iC is now the most left point
					{
				iL = -1;
			} else {
				iL = indsel[imin - 2];
			}
			iC = indsel[imin - 1];
			iR = indsel[imin + 1];
			h0L = (double) (M_total_length - iC - 1) * (iC + 1) / M_total_length * (iR - iL) / (iR - iC)
					/ (iC - iL);
			if (wr[iC] >= 0) {
				Scores[iC] = (double) (wr[iC] / sqrt(h0L));
			} else {
				Scores[iC] = (double) (-wr[iC] / sqrt(h0L));
			}
#ifdef _DEBUG_
			printf("h0L=%g\niCsecond=%ld\nScores[iC]=%g\n\n",h0L,iC,Scores[iC]);
#endif
		}
		NumRem = NumRem - 1;
		for (i = imin; i < NumRem; i++) {
			indsel[i] = indsel[i + 1];

		}
#ifdef _DEBUG_
		printf("INDSEL RECOMPUTE:\n\nindsel[0]=%ld\nindsel[1]=%ld\nindsel[2]=%ld\nindsel[%ld]=%ld\n\n",indsel[0],indsel[1],indsel[2],NumRem-1,indsel[NumRem-1]);
#endif
		if (Tau < 0) {
#ifdef _DEBUG_
			printf("EXIT BECAUSE TAU IS LOWER THAN ZERO\n\n");
#endif
			break;
		}
	}

#ifdef _DEBUG_
	printf("OUTPUT Parameters:\n\nvmin=%g\nNumRem=%ld\nimin=%ld",vmin,NumRem,imin);
#endif

	*pointTau = vmin;
	*pointNumRem = NumRem;
	return (inrem);
}

/*************************************************************************/
/*
 % CollapseAmp -- Collapses the amplitudes of the non significantly altered segments into a base level.
 % [OutAmp,State]=CollapseAmp(SegAmp,SegLen,BaseAmp,sigma2,T);
 */

void BaseGADA::CollapseAmpTtest() {
	//Uses a T test to decide which segments collapse to neutral
	long k;

//	if(L==1)
//	{
//		if( fabs(SegAmp[0]-BaseAmp)/sqrt(sigma2/(double)SegLen[0]) < T)
//			SegAmp[0]=BaseAmp;
//	}
	SegState = (double *) calloc(K + 1, sizeof(double));
	for (i = 0; i <= K; i++)
		SegState[i] = SegAmp[i];
	for (k = 0; k < K; k++)
		if (fabs(SegAmp[k] - BaseAmp) / sqrt(sigma2 / (double) SegLen[k]) < T) {
			//but do it only if one of the neigboring ones have been collapsed
			if ((k > 0)
					&& (fabs(SegAmp[k - 1] - BaseAmp)
							/ sqrt(sigma2 / (double) SegLen[k - 1]) < T))
				SegAmp[k] = BaseAmp;
			if ((k < K - 1)
					&& (fabs(SegAmp[k + 1] - BaseAmp)
							/ sqrt(sigma2 / (double) SegLen[k + 1]) < T))
				SegAmp[k] = BaseAmp;
			//or it is the initial segment or the final segment of the unit, since we assume that the unseen neighbors where in collapsed state
			if ((k == 0) || (k == K - 1))
				SegAmp[k] = BaseAmp;
			//or it has larger size than the neighboring segments
			if ((k > 0) && (SegLen[k] > SegLen[k - 1]))
				SegAmp[k] = BaseAmp;
			if ((k < K - 1) && (SegLen[k] > SegLen[k + 1]))
				SegAmp[k] = BaseAmp;

		}

	/*
	 //Classify amplitudes...
	 for(k=0;k++;k<K)
	 {
	 if(SegAmp[k]>BaseAmp)
	 SegAmp=+1; //Gain
	 else if(SegAmp[k]<BaseAmp)
	 SegAmp=-1;
	 else
	 SegAmp=0;
	 }
	 */
}
/*************************************************************************/
/*
 % ClassifySegments -- Classifies/Collapses reconstructed segments into Altered/Gain/Loss
 % Gain uses positive numbers log2(3)-log2(2)
 */

void BaseGADA::ClassifySegments(double *SegAmp, long *SegLen, double *SegState, long K,
		double BaseAmp, double ploidy, double sigma2, //Reference noise
		double T //Critical value that decides when to colapse
		) {
	long k;
	double c;
	double aux_local;

	for (k = 0; k <= K; k++) {
		if (fabs(SegAmp[k] - BaseAmp) / sqrt(sigma2 / (double) SegLen[k]) < T) {
			SegState[k] = ploidy;
		} else if ((SegAmp[k] - BaseAmp) > 0) {
			for (c = ploidy; c < 100; c++) {
				aux_local = (SegAmp[k] - BaseAmp - log2(c) + log2(ploidy))
						/ sqrt(sigma2 / (double) SegLen[k]);
				if (aux_local < T) {
					break;
				}
			}
			SegState[k] = c;
		} else if ((SegAmp[k] - BaseAmp) < 0) {
			for (c = ploidy; c > 0; c--) {
				if ((SegAmp[k] - BaseAmp - log2(c) + log2(ploidy))
						/ sqrt(sigma2 / (double) SegLen[k]) > -T) {
					break;
				}
			}
			SegState[k] = c;
		}
	}
}


double // Returns BaseAmp corresponding to the base level.
BaseGADA::CompBaseAmpMedianMethod() {
	//Computes the median recontruction level, as baseline level.
	long M_total_length, k, RunLen;
	BaseAmp = 0;

	//If they need to be sorted use the following...
	double *D;
	long *I;
	D = (double*) calloc(K+1,sizeof(double));
	I = (long*) calloc(K+1,sizeof(long));
	for (k = 0; k < K + 1; k++) {
		D[k] = SegAmp[k];
		I[k] = SegLen[k];
	}
	doubleBubbleSort(D, I, K + 1); //I need indexes of the sort
	SegAmp = D;
	SegLen = I;

	M_total_length = 0;
	for (k = 0; k <= K; k++)
		M_total_length = M_total_length + SegLen[k];
#ifdef _DebugCompBaseAmpMedianMethod_
	printf("_DebugCompBaseAmpMedianMethod_: M_total_length%ld K%ld M_total_length/2%ld\n",M_total_length,K,M_total_length/2);
#endif

	RunLen = 0;
	k = 0;
	while (RunLen < M_total_length / 2)
		RunLen = RunLen + SegLen[k++];
#ifdef _DebugCompBaseAmpMedianMethod_
	printf("_DebugCompBaseAmpMedianMethod_: k%ld RunLen=%ld SegAmp[k-1]%g SegAmp[k]%g\n",k,RunLen,SegAmp[k-1],SegAmp[k]);
#endif

	BaseAmp = SegAmp[k - 1];

#ifdef _DebugCompBaseAmpMedianMethod_
	printf("_DebugCompBaseAmpMedianMethod_: BaseAmp=%g\n",BaseAmp);
#endif

	return BaseAmp;
}

/**************************************************************************************************************************/


long BaseGADA::SBLandBE() {
	//double convergenceDelta, convergenceMaxAlpha, convergenceB;
	//long maxNoOfIterations;
	// SBL optimization parameters
	//convergenceDelta = 1E-10; //1E-10 or 1E-8 seems to work well for this parameter. -- => ++ conv time
	//convergenceMaxAlpha = 1E8; //1E8 better than 1E10 seems to work well for this parameter. -- => -- conv time
	//maxNoOfIterations = 100000; //Maximum number of iterations to reach convergence...
	//convergenceB = 1E-20; //
	//sigma2 = *Psigma2;

    //2013.08.28 no more copying of input data. to reduce memory usage.
	//normalized_data_array = inputDataArray;

    normalized_data_array = (double* )calloc(_M_total_length, sizeof(double));
    for(i=0;i<_M_total_length;i++)
        normalized_data_array[i]=inputDataArray[i];

    long i;
    double delta;
    if (sigma2 < 0) {
        //If sigma2 < 0, estimate sigma2 from data
		sigma2 = 0;
		for (i = 1; i < _M_total_length; i++) {
			delta = normalized_data_array[i] - normalized_data_array[i - 1];
			sigma2 += (0.5 * delta * delta);
		}
		sigma2 = sigma2 / (_M_total_length - 1);
	}
	if (debug){
		std::cerr << boost::format("sigma^2 (of difference of adjacent values) = %1% .\n")% sigma2;
	}
	//Mean removal
	ymean = 0;
	for (i = 0; i < _M_total_length; ++i)
		ymean += normalized_data_array[i];
	ymean = ymean / _M_total_length;
	for (i = 0; i < _M_total_length; ++i)
		normalized_data_array[i] = normalized_data_array[i] - ymean;



	//Call to SBL
	if (debug){
		std::cerr << "_SBLBE_ Memory initialization\n";
	}

	K = _M_total_length - 1;
	Iext = (long*) calloc(_M_total_length,sizeof(long)); //R
	Wext = (double*) calloc(_M_total_length,sizeof(double)); //R
	_alpha_array = (double*) calloc(_M_total_length,sizeof(double));
	_aux_array = (double*) calloc(_M_total_length,sizeof(double));

	if (debug){
		std::cerr << "_SBLBE_ Breakpoint Initialization\n";
	}
	//Initialize breakpoints
	for (i = 0; i < _M_total_length; i++)
		_alpha_array[i] = 0.0;
	for (i = 0; i < _M_total_length; i++)
		Wext[i] = 0.0;
	for (i = 0; i < _M_total_length; i++)
		Iext[i] = i;

	if (debug){
		std::cerr << "_SBLBE_ SBL starts\n";
	}
	numEMsteps = SBL(normalized_data_array, Iext, _alpha_array, Wext + 1, _aux_array, _M_total_length, &K, sigma2, a, convergenceB,
			convergenceMaxAlpha, maxNoOfIterations, convergenceDelta, debug);
	noOfBreakpointsAfterSBL = K;	//2013.08.31 K would be changed later on.

	//Convert Iext and Wext to the extended notation.
	Iext = (long*) realloc(Iext, (K + 2) * sizeof(long));
	for (i = (K + 1); i > 0; i--)
		Iext[i] = Iext[i - 1] + 1;
	Iext[0] = 0;
	Iext[K + 1] = _M_total_length;

	Wext = (double *) realloc(Wext, (K + 1) * sizeof(double));
	//for(i=K;i>0;i++)
	//	Wext[i]=Wext[i-1];
	Wext[0] = ymean;

	if (debug){
		std::cerr << "_SBLBE_ Backward Elimination T=" << T << " MinSegLen=" << MinSegLen << std::endl;
	}
	BEwTandMinLen(Wext, Iext, &K, sigma2, T, MinSegLen, debug);

	Iext = (long*) realloc(Iext, (K + 2) * sizeof(long));
	Wext = (double *) realloc(Wext, (K + 1) * sizeof(double));

	if (debug){
		std::cerr << "_SBLBE_ After BE K=" <<K << endl;
	}

	//Freeing memory
	//free(_alpha_array);
	//free(normalized_data_array);	#2013.08.30
	//free(_aux_array);

	//*pIext = Iext;
	//*pWext = Wext;
	//*Psigma2 = sigma2;
	return K;
}

/**************************************************************************************************************************/
long BaseGADA::BEwTandMinLen( //Returns breakpoint list length. with T and MinSegLen
		double *Wext, //IO Breakpoint weights extended notation...
		long *Iext, //IO Breakpoint positions in extended notation...
		long *pK, //IO Number breakpoint positions remaining.
		double sigma2, //IP If sigma2
		double T, //IP  Threshold to prune,  T=T*sqrt(sigma2);
		long MinSegLen, //IP Minimum length of the segment.
		long debug
		) {
	long i, K;
	double *tscore_local; //Statistical Scores,

	K = *pK; //Number of breakpoints
	_M_total_length = Iext[K + 1]; //Total length
	T = T * sqrt(sigma2); //Adjusting T to the noise power

	tscore_local = (double*) calloc(K+1,sizeof(double)); //

	//Computing scores
	if (debug){
			std::cerr << boost::format("_BEwTandMinLen_() Computing scores\n");
	}
	for (i = 0; i < K + 1; i++)
		tscore_local[i] = 0;

	ComputeTScores(Wext, Iext, tscore_local, K, 1, K);

	if (debug){
		std::cerr << boost::format("_BEwTandMinLen_:K%1% w[0]=%2%,w[1]=%3%,w[2]=%4%,w[K-1]=%5%,w[K]=%6%\n") % K % Wext[0] % Wext[1] % Wext[2] % Wext[K-1] % Wext[K];
		std::cerr << boost::format("_BEwTandMinLen_:K%1% t[0]=%2%,t[1]=%3%,t[2]=%4%,t[K-1]=%5%,t[K]=%6%\n") % K % tscore_local[0] % tscore_local[1] % tscore_local[2] % tscore_local[K-1] % tscore_local[K];
		std::cerr << boost::format("_BEwTandMinLen_() Starting Backward Elimination: sigma2=%1% T=%2% MinSegLen %3% ... \n") % sigma2 % T % MinSegLen;
	}

	BEwTscore(Wext, Iext, tscore_local, &K, T, MinSegLen, debug);
	if (debug){
		std::cerr << boost::format("_BEwTandMinLen_:K%1% w[0]=%2%,w[1]=%3%,w[2]=%4%,w[K-1]=%5%,w[K]=%6%\n") % K % Wext[0] % Wext[1]%Wext[2]%Wext[K-1]%Wext[K];
		std::cerr << boost::format("_BEwTandMinLen_:K%1% t[0]=%2%,t[1]=%3%,t[2]=%4%,t[K-1]=%5%,t[K]=%6%\n") % K%tscore_local[0]%tscore_local[1]%tscore_local[2]%tscore_local[K-1]%tscore_local[K];
	}

	if (debug>0){
		std::cerr << boost::format("_BEwTandMinLen_() finished. number of breakpoints=%1% \n") % K;
	}

	//Do not free tscore_local as it will trigger memory free error during program exit. 
	//  tscore_local is realloc-ed in BEwTscore().
	//free(tscore_local); 

	*pK = K;
	return K;
}

/******************************************************/
//BEwTscore(Iext,Wext,h0,h1,tscore_array,&K,T);  //Need to update BEthres to operate on the Iext Wext notation...
long BaseGADA::BEwTscore(double *Wext, //IO Breakpoint weights extended notation...
		long *Iext, //IO Breakpoint positions in extended notation...
		double *tscore_array, long *pK, //IO Number breakpoint positions remaining.
		double T, //IP  Threshold to prune
		long MinSegLen,	//minimum segment length
		long debug
		) {
	long i, K, M_total_length;

	K = *pK; //Number of breakpoints
	M_total_length = Iext[K + 1]; //Total length

	if (debug>0){
		std::cerr << boost::format("BEwTscore(): BE starts K=%1% M_total_length=%2% T=%3% MinSegLen=%4% ... \n")% K % M_total_length % T % MinSegLen;
	}
	if (debug>0){
		std::cerr << boost::format("Adding %1% breakpoints into a red-black tree ... ")% K ;
	}
	treeType rbTree = treeType();

	BreakPoint* leftBreakPointPtr=NULL;
	BreakPoint* rightBreakPointPtr=NULL;

	rbNodeType* currentNodePtr=rbTree.nil;
	rbNodeType* genomeLeftNodePtr=rbTree.nil;
	rbNodeType* genomeRightNodePtr=rbTree.nil;
	BreakPoint *leftMostBreakPointPtr = new BreakPoint(Iext[0], Wext[0], tscore_array[0], 0, MinSegLen, T, Iext[K+1]);
	leftMostBreakPointPtr->nodePtr = rbTree.nil;
	BreakPoint *rightMostBreakPointPtr = new BreakPoint(Iext[K+1], 0, 0, 0, MinSegLen, T, Iext[K+1]);
	rightMostBreakPointPtr->nodePtr = rbTree.nil;
	int maxBPSetSize =0;
	for (i = 1; i < K + 1; i++){
		long segLength = min(Iext[i]-Iext[i-1],Iext[i+1]-Iext[i]);	//shorter of two neighboring segments as length for the breakpoint
		BreakPoint* bpPtr = new BreakPoint(Iext[i], Wext[i], tscore_array[i], segLength, MinSegLen, T, Iext[K+1]);
		BreakPointKey bpKey = bpPtr->getKey();
		//cerr<< *bpPtr << endl;
		//cerr << boost::format("i=%1%, tree size=%2%, tree valid=%3%")% i % rbTree.size() % rbTree.isValidRedBlackTree() << endl;
		if (i==1){
			bpPtr->setLeftBreakPoint(leftMostBreakPointPtr);
			leftMostBreakPointPtr->setRightBreakPoint(bpPtr);
		}
		else if (leftBreakPointPtr!=NULL){
			bpPtr->setLeftBreakPoint(leftBreakPointPtr);
			leftBreakPointPtr->setRightBreakPoint(bpPtr);
			if (i==K){	//last point
				bpPtr->setRightBreakPoint(rightMostBreakPointPtr);
				rightMostBreakPointPtr->setLeftBreakPoint(bpPtr);
			}
		}
		currentNodePtr = rbTree.queryTree(bpKey);
		if (rbTree.isNULLNode(currentNodePtr)){
			rbNodeDataType* dataPtr = new rbNodeDataType();
			currentNodePtr = rbTree.insertNode(bpKey, dataPtr);
		}
		currentNodePtr->getDataPtr()->insert(bpPtr);
		if (currentNodePtr->getDataPtr()->size()>maxBPSetSize){
			maxBPSetSize = currentNodePtr->getDataPtr()->size();
		}
		bpPtr->setRBTreeNodePtr(currentNodePtr);
		//reset the left break point pointer
		leftBreakPointPtr = bpPtr;
	}
	if (debug>0){
		//rbTree.printTree();
		std::cerr << boost::format(" noOfNodesInTree=%1%, tree max depth =%2%, tree valid=%3%, maxBPSetSize=%4%.\n") %
				rbTree.noOfNodes() % rbTree.maxDepth() % rbTree.isValidRedBlackTree() % maxBPSetSize;
	}
	long toRemoveSegmentLength=-1;	//2013.08.31 initial value =-1, so that it is <MinSegLen
	long previousToRemoveSegmentLength = -1;
	double currentMinScore;
	double previousRoundMinScore=-1;
	long counter = 0;

	rbNodeType*  minNodePtr = NULL;
	minNodePtr = rbTree.getMinimum();
	BreakPointKey minBPKey=minNodePtr->getKey();
	rbNodeDataType* setOfBPPtr = minNodePtr->getDataPtr();
	rbNodeDataType::iterator setOfBPIterator=(*setOfBPPtr).begin();
	//reset
	leftBreakPointPtr=NULL;
	rightBreakPointPtr=NULL;
	genomeLeftNodePtr=rbTree.nil;
	genomeRightNodePtr=rbTree.nil;

	currentMinScore = minBPKey.tscore;
	toRemoveSegmentLength = minBPKey.segmentLength;
	while (rbTree.noOfNodes()>0 && (currentMinScore<T || toRemoveSegmentLength<MinSegLen)){
		minBPKey = minNodePtr->getKey();
		setOfBPPtr = minNodePtr->getDataPtr();
		if (debug>0 && counter%reportIntervalDuringBE==0){
			std::cerr << boost::format("BEwTscore(): iteration no=%1% T=%2% MinSegLen=%3%: minimum break point key: %4% previousRoundMinScore=%5% previousToRemoveSegmentLength=%6% noOfSegments=%7% setOfBPPtr.size=%8% \n") %
					counter % T % MinSegLen % minBPKey %  previousRoundMinScore % previousToRemoveSegmentLength %
					rbTree.noOfNodes() % (*setOfBPPtr).size();
			cerr << boost::format("\t currentMinScore=%1%, toRemoveSegmentLength=%2% \n")%
							currentMinScore % toRemoveSegmentLength;
		}
		for (setOfBPIterator =(*setOfBPPtr).begin(); setOfBPIterator!=(*setOfBPPtr).end(); setOfBPIterator++){
			//remove all breakpoints in this node's data (they have same tscore and length)
			BreakPoint* minBPPtr = *setOfBPIterator;	//get address of BreakPoint
			leftBreakPointPtr = minBPPtr->leftBreakPointPtr;
			rightBreakPointPtr = minBPPtr->rightBreakPointPtr;

			if ((debug>0) && counter%reportIntervalDuringBE==0){
				std::cerr << boost::format("\t BEwTscore(): iteration no=%1% T=%2% MinSegLen=%3%: break point to be removed: %4% previousRoundMinScore=%5% previousToRemoveSegmentLength=%6% noOfSegments=%7% \n") %
						counter % T % MinSegLen % *minBPPtr %  previousRoundMinScore % previousToRemoveSegmentLength % rbTree.noOfNodes();
			}
			//update two neighboring break points.
			minBPPtr->removeItself();

			//modify genome left & right key, delete their nodes from tree and re-add them with new key and updated break point info
			if (leftBreakPointPtr!=NULL &&  leftBreakPointPtr->nodePtr!=rbTree.nil && leftBreakPointPtr->nodePtr!=NULL){
				//delete the outdated left node
				genomeLeftNodePtr = (rbNodeType*)leftBreakPointPtr->nodePtr;
				genomeLeftNodePtr->getDataPtr()->erase(leftBreakPointPtr);
				if (genomeLeftNodePtr->getDataPtr()->size()==0){
					//delete this node altogether if its vector is empty
					rbTree.deleteNode(genomeLeftNodePtr);
				}
				//new genomeLeftNodePtr that matches the new key
				genomeLeftNodePtr = rbTree.queryTree(leftBreakPointPtr->getKey());
				if (rbTree.isNULLNode(genomeLeftNodePtr)){
					//create an new node
					genomeLeftNodePtr = rbTree.insertNode(leftBreakPointPtr->getKey(), new rbNodeDataType() );
				}
				genomeLeftNodePtr->getDataPtr()->insert(leftBreakPointPtr);
				leftBreakPointPtr->nodePtr = genomeLeftNodePtr;
			}
			if (rightBreakPointPtr!=NULL && !rbTree.isNULLNode((rbNodeType*)rightBreakPointPtr->nodePtr) && rightBreakPointPtr->nodePtr!=NULL){
				//delete the outdated right node
				genomeRightNodePtr = (rbNodeType*)rightBreakPointPtr->nodePtr;
				genomeRightNodePtr->getDataPtr()->erase(rightBreakPointPtr);
				if (genomeRightNodePtr->getDataPtr()->size()==0){
					//delete this node altogether if its vector is empty
					rbTree.deleteNode(genomeRightNodePtr);
				}
				//new genomeRightNodePtr that matches the new key
				genomeRightNodePtr = rbTree.queryTree(rightBreakPointPtr->getKey());
				if (rbTree.isNULLNode(genomeRightNodePtr)){
					//create an new node
					genomeRightNodePtr = rbTree.insertNode(rightBreakPointPtr->getKey(), new rbNodeDataType() );
				}
				genomeRightNodePtr->getDataPtr()->insert(rightBreakPointPtr);
				rightBreakPointPtr->nodePtr = genomeRightNodePtr;
			}
		}
		(*setOfBPPtr).clear();
		//delete this minimum node after its data is all tossed out
		rbTree.deleteNode(minNodePtr);

		counter ++;
		previousRoundMinScore = minBPKey.tscore;
		previousToRemoveSegmentLength = minBPKey.segmentLength;
		if (rbTree.noOfNodes()>0){
			//get a new minimum
			minNodePtr = rbTree.getMinimum();
			minBPKey = minNodePtr->getKey();
			currentMinScore = minBPKey.tscore;
			toRemoveSegmentLength = minBPKey.segmentLength;
		}
		else{
			break;
		}

	}
	if (debug>0){
		std::cerr << boost::format("BEwTscore(): last iteration no=%1% T=%2% MinSegLen=%3%: minimum break point key: %4%, previousRoundMinScore=%5% previousToRemoveSegmentLength=%6% tree size=%7%, noOfNodesInTree=%8% \n") %
				counter % T % MinSegLen % minBPKey % previousRoundMinScore %
				previousToRemoveSegmentLength % rbTree.size() % rbTree.noOfNodes();
		cerr << boost::format("\t currentMinScore=%1%, toRemoveSegmentLength=%2% \n")%
				currentMinScore % toRemoveSegmentLength;
	}

	// convert data back to old data structures
	K = rbTree.noOfNodes();	//update the number of break points
	vector<BreakPoint> breakPointVector;
	breakPointVector.push_back(*leftMostBreakPointPtr);	//add this first
	while (rbTree.noOfNodes()>0){
		minNodePtr = rbTree.getMinimum();
		setOfBPPtr = minNodePtr->getDataPtr();
		for (setOfBPIterator =(*setOfBPPtr).begin(); setOfBPIterator!=(*setOfBPPtr).end(); setOfBPIterator++){
			breakPointVector.push_back(**setOfBPIterator);
		}
		rbTree.deleteNode(minNodePtr);
	}
	breakPointVector.push_back(*rightMostBreakPointPtr);	//add the right most.
	if (debug){
		cerr << boost::format("breakPointVector size=%1%, noOfNodesInTree=%2%, tree size=%3%, tree valid=%4%")%
				breakPointVector.size() % rbTree.noOfNodes() % rbTree.size() % rbTree.isValidRedBlackTree() << endl;
	}
	//sort all the breakpoints in the tree by chromosomal position, reconstruct Iext, Wext, tscore, pK
	sort(breakPointVector.begin(),breakPointVector.end());
		/*
		 * a c++0x sorting lambda function, not universally accepted
		 	[](const BreakPoint& a, const BreakPoint& convergenceB) -> bool
			{return a.position < b.position;});
		 *
		 */
	for (i=1; i<K+1; i++){
		Wext[i]=breakPointVector[i].weight;
		Iext[i]=breakPointVector[i].position;
		tscore_array[i]=breakPointVector[i].tscore;
	}
	Iext[K+1] = leftMostBreakPointPtr->totalLength;
	Iext = (long*) realloc(Iext, (K + 2) * sizeof(long));
	Wext = (double *) realloc(Wext, (K + 1) * sizeof(double));
	tscore_array = (double *) realloc(tscore_array, (K + 1) * sizeof(double));

	//free(leftMostBreakPointPtr);
	//free(rightMostBreakPointPtr);
	leftBreakPointPtr = NULL;
	rightMostBreakPointPtr = NULL;
	*pK = K;
	return K;
}

long BaseGADA::RemoveBreakpoint(double *Wext, long *Iext, double *tscore_array, long K, long indexOfSegmentToRemove) {
	long j;
	double iC, iL, iR, M_total_length;

	M_total_length = (double) Iext[K + 1];
	iL = (double) Iext[indexOfSegmentToRemove - 1];
	iC = (double) Iext[indexOfSegmentToRemove];
	iR = (double) Iext[indexOfSegmentToRemove + 1];

	//Change coefficients
	if (indexOfSegmentToRemove > 1){
		Wext[indexOfSegmentToRemove - 1] = Wext[indexOfSegmentToRemove - 1]
				+ sqrt((M_total_length - iL) / (M_total_length - iC) * iL / iC) * (iR - iC) / (iR - iL) * Wext[indexOfSegmentToRemove];
	}
	if (indexOfSegmentToRemove < K){
		Wext[indexOfSegmentToRemove + 1] = Wext[indexOfSegmentToRemove + 1]
				+ sqrt((M_total_length - iR) / (M_total_length - iC) * iR / iC) * (iC - iL) / (iR - iL) * Wext[indexOfSegmentToRemove];
	}
	Wext[indexOfSegmentToRemove] = 0; //
	//Shorten list
	for (j = indexOfSegmentToRemove; j < K; j++){
		tscore_array[j] = tscore_array[j + 1];
		Wext[j] = Wext[j + 1];
		Iext[j] = Iext[j + 1];
	}
	//2013.09.01 last one for Iext
	if (indexOfSegmentToRemove<=K){
		Iext[K] = Iext[K + 1];
	}
	//for (j = indexOfSegmentToRemove; j < K + 1; j++)
	//	Iext[j] = Iext[j + 1];
	return K - 1;
}

void BaseGADA::ComputeTScores(const double *Wext, const long *Iext, double *Scores,
		long K, long start, long end) {
	long j;
	double h0, M_total_length;

	M_total_length = (double) Iext[K + 1];

	for (j = start; j <= end; j++) {
		h0 = (double) (M_total_length - Iext[j]) * (double) Iext[j] / M_total_length * (double) (Iext[j + 1] - Iext[j - 1])
				/ (double) (Iext[j + 1] - Iext[j]) / (double) (Iext[j] - Iext[j - 1]);
		Scores[j] = fabs(Wext[j]) / sqrt(h0);
	}
}

void BaseGADA::Project(double *y, long M_total_length, long *I, long L, double *xI, double *wI) {
	// Intern variables declaration
	double aux_double = 0;
	double ymean = 0;
	long i = 0;
	double *h0 = NULL;
	double *h1 = NULL;
	double *w0 = NULL;
	double *z = NULL;
	double *temp = NULL;
	double *wr = NULL;
	double *aux_vec = NULL;

	// Variables inizialitation
	h0 = (double*) calloc(M_total_length-1,sizeof(double));
	h1 = (double*) calloc(M_total_length-2,sizeof(double));
	z = (double*) calloc(M_total_length-1,sizeof(double));
	temp = (double*) calloc(L,sizeof(double));
	wr = (double*) calloc(M_total_length,sizeof(double));
	aux_vec = (double*) calloc(M_total_length,sizeof(double));

//     h0=mxGetPr(mxCreateDoubleMatrix(1,M_total_length-1,mxREAL));
//     h1=mxGetPr(mxCreateDoubleMatrix(1,M_total_length-2,mxREAL));
//     z=mxGetPr(mxCreateDoubleMatrix(1,M_total_length-1,mxREAL));
//     temp=mxGetPr(mxCreateDoubleMatrix(1,L,mxREAL));
//     wr=mxGetPr(mxCreateDoubleMatrix(1,M_total_length,mxREAL));
//     aux_vec=mxGetPr(mxCreateDoubleMatrix(1,M_total_length,mxREAL));

	// Remove y's mean
	aux_double = 0;

	for (i = 0; i < M_total_length; i++) {
		aux_double = aux_double + y[i];
	}
	ymean = (double) (aux_double / M_total_length);
	for (i = 0; i < M_total_length; i++) {
		y[i] = y[i] - ymean; // I have y mean in aux_double
	}
#ifdef _DEBUG_
	//checked-> inputs are OK
	printf("MEAN CHECKING:\n\nymean=%g\n",ymean);
	printf("\ny[0]=%g\ny[1]=%g\ny[2]=%g\ny[3]=%g\ny[%ld]=%g\nM_total_length=%ld\n",y[0],y[1],y[2],y[3],M_total_length-1,y[M_total_length-1],M_total_length);

#endif

	//start reconstruction
	if (L > 0) {
		//I sort -> maybe desorded
		BubbleSort(I, L);
		//
#ifdef _DEBUG_
		//checked-> inputs are OK
		printf("BubbleSort CHECKING:\n");
		printf("\nI[0]=%ld\nI[1]=%ld\nI[2]=%ld\nI[3]=%ld\nI[%ld]=%ld\nL=%ld\n",I[0],I[1],I[2],I[3],L-1,I[L-1],L);
#endif
		//

		ComputeH(h0, h1, M_total_length);
#ifdef _DEBUG_
		//checked-> inputs are OK
		printf("COMPUTE H CHECKING:\n\nh0[0]=%g\nh0[1]=%g\nh0[2]=%g\nh0[3]=%g\nh0[%ld]=%g\nsizeh0=%ld\n",h0[0],h0[1],h0[2],h0[3],M_total_length-2,h0[M_total_length-2],M_total_length-1);
		printf("\nh1[0]=%g\nh1[1]=%g\nh1[2]=%g\nh1[3]=%g\nh1[%ld]=%g\nsizeh0-1=%ld\n",h1[0],h1[1],h1[2],h1[3],M_total_length-3,h1[M_total_length-3],M_total_length-2);

#endif

		ComputeFdualXb(M_total_length, y);
		w0 = y; //careful I just erased y's value
#ifdef _DEBUG_
		//checked-> inputs are OK
		printf("Compute F dual CHECKING:\n\nw0[0]=%g\nw0[1]=%g\nw0[2]=%g\nw0[3]=%g\nw0[%ld]=%g\nM_total_length=%ld\n",w0[0],w0[1],w0[2],w0[3],M_total_length-1,w0[M_total_length-2],M_total_length-1);

#endif

		TrisolveREG(h0, h1, h1, w0, z, M_total_length - 1);

#ifdef _DEBUG_
		//checked-> inputs are OK
		printf("Trisolve Checking:\n\nz[0]=%g\nz[1]=%g\nz[2]=%g\nz[3]=%g\nz[%ld]=%g\nM_total_length=%ld\n",z[0],z[1],z[2],z[3],M_total_length-1,z[M_total_length-2],M_total_length-1);

#endif
		ComputeHs(I, M_total_length, L, h0, h1);
#ifdef _DEBUG_
		//checked-> inputs are OK
		printf("Compute Hs2 CHECKING:\n\nh0[0]=%g\nh0[1]=%g\nh0[2]=%g\nh0[3]=%g\nh0[%ld]=%g\nL=%ld\n",h0[0],h0[1],h0[2],h0[3],L-1,h0[L-1],L);
		printf("\nh1[0]=%g\nh1[1]=%g\nh1[2]=%g\nh1[3]=%g\nh1[%ld]=%g\nL-1=%ld\n",h1[0],h1[1],h1[2],h1[3],L-2,h1[L-2],L-1);

#endif
		for (i = 0; i < L; i++) {
			temp[i] = z[I[i]];
		}

		for (i = 0; i < L; i++) {
			y[i] = 0;
		}
#ifdef _DEBUG_
		//checked-> inputs are OK
		printf("Z(I) CHECKING:\n\ntemp[0]=%g\ntemp[1]=%g\ntemp[2]=%g\ntemp[3]=%g\ntemp[%ld]=%g\nL=%ld\n",temp[0],temp[1],temp[2],temp[3],L-1,temp[L-1],L);

#endif
		TriSymGaxpy(h0, h1, temp, L, wI);
#ifdef _DEBUG_
		//checked-> inputs are OK
		printf("TriSymGaxpy CHECKING:\n\nwI[0]=%g\nwI[1]=%g\nwI[2]=%g\nwI[3]=%g\nwI[%ld]=%g\nL=%ld\n",wI[0],wI[1],wI[2],wI[3],L-1,wI[L-1],L);

#endif
		for (i = 0; i < L; i++) {
			wr[I[i]] = wI[i];
#ifdef _DEBUG_
			//checked-> inputs are OK
			printf("wr CHECKING:\n\nwr[I[%ld]]=%g\nI[i]=%ld\n",i,wr[I[i]],I[i]);

#endif
		}
#ifdef _DEBUG_
		//checked-> inputs are OK
		printf("wr CHECKING:\n\nwr[I[0]]=%g\nwr[I[1]]=%g\nwr[%ld]=%g\nM_total_length=%ld\n",wr[I[0]],wr[I[1]],I[L-1],wr[I[L-1]],M_total_length);
#endif

		reconstruct(wr, M_total_length, aux_vec);
#ifdef _DEBUG_
		//checked-> inputs are OK
		printf("RECONSTRUCT CHECKING:\n\nwr[0]=%g\nwr[1]=%g\nwr[2]=%g\nwr[3]=%g\nwr[%ld]=%g\nM_total_length=%ld\naux_double=%g\nymean=%g\n",wr[0],wr[1],wr[2],wr[3],M_total_length-1,wr[M_total_length-1],M_total_length,aux_double,ymean);

#endif
		for (i = 0; i < M_total_length; i++) {
			xI[i] = ymean + wr[i];
		}
	} else {
		L = 0;
		for (i = 0; i < M_total_length; i++) {
			xI[i] = ymean;
		}
	}

	free(h0);
	free(h1);
	free(z);
	free(temp);
	free(wr);
	free(aux_vec);

}

void BaseGADA::ProjectCoeff(double *y, long M_total_length, long *Iext, long K, double *Wext) {
	//IextYobs2Wext
	// Intern variables declaration
	double ymean = 0;
	long i = 0;
	double *h0;
	double *h1;
	double *z;

	// Variables inizialitation
	h0 = (double*) calloc(K,sizeof(double));
	h1 = (double*) calloc(K-1,sizeof(double));
	z = (double*) calloc(M_total_length-1,sizeof(double));

	ymean = 0;
	for (i = 0; i < M_total_length; i++)
		ymean = ymean + y[i];
	ymean = ymean / M_total_length;

	// Remove y's mean, Not necessary
//    for (i=0;i<M_total_length;i++)
//        y[i]=y[i]-ymean;

	if (K > 0) {
		//I sort -> assumed that I is already ordered
		// BubbleSort(Iext,K+2);

		CompZ(y, z, M_total_length);
//		printf("\n CHECKING: ymean=%g,M_total_length=%ld,K=%ld\n",ymean,M_total_length,K);//
//		printf("\n z CHECKING: z[0]=%g,z[1]=%g,z[2]=%g,z[M_total_length-3]=%g,z[M_total_length-2]=%g\n",z[0],z[1],z[2],z[M_total_length-3],z[M_total_length-2]);

		for (i = 1; i <= K; i++)
			z[i - 1] = z[Iext[i] - 1];
//    	printf("\n z CHECKING: z[0]=%g,z[1]=%g,z[2]=%g,z[K-2]=%g,z[K-1]=%g\n",z[0],z[1],z[2],z[K-2],z[K-1]);

		ComputeHsIext(Iext, K, h0, h1);
		//printf("\n h0 CHECKING: h0[0]=%g,h0[1]=%g,h0[2]=%g,h0[K-2]=%g,h0[K-1]=%g,h0[K]=%g\n",h0[0],h0[1],h0[2],h0[K-2],h0[K-1]);
		//printf("\n h1 CHECKING: h1[0]=%g,h1[1]=%g,h1[2]=%g,h1[K-2]=%g,h1[K-1]=%g\n",h1[0],h1[1],h1[2],h1[K-2]);

		for (i = 0; i < K + 1; i++)
			Wext[i] = 0.0;

		TriSymGaxpy(h0, h1, z, K, Wext + 1);
		//   	printf("\n w CHECKING: w[0]=%g,w[1]=%g,w[2]=%g,w[K-1]=%g,w[K]=%g\n",Wext[0],Wext[1],Wext[2],Wext[K-1],Wext[K]);
	}
	Wext[0] = ymean;


	free(h0);
	free(h1);
	free(z);


}

