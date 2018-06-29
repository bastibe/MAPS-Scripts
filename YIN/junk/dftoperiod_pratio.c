/* 
* dftoperiods.c
* estimate period from df
* 
* Alain de Cheveigné, CNRS/Ircam	Dec 2001
* (c) 2001 CNRS
*/

#include <math.h>
#include "mex.h"

//#define MWRAP
#include "mwrap_check.h"

/* Input Arguments */
#define D_IN	prhs[0]
#define B_IN	prhs[1]
#define T_IN	prhs[2]
#define R_IN	prhs[3]

/* output arguments */
#define PRD_OUT	plhs[0]
#define GOOD_OUT	plhs[1]

static void ddftoperiods(
				double *dp,			/* input df */
				double *bp,			/* bounds */
				double *prdp, 		// periods
				double *goodp, 		// flag indicating that a good estimate was found
				double thresh,		// threshold
				double pratio,		// ratio of min to global min
				int m,				// nrows of df
				int n				// ncols of df
			) 
{
	int j,k, p, lo, hi, maxj, goodflag;
	double z, min, globalmin, thresh2;
	
	// bounds
	lo= GET(bp[0]);	
	hi= GET(bp[1]);
	if (lo<0 || hi>m ) {
		mexPrintf("%d %d\n", lo, hi);
		mexErrMsgTxt("DFTOPERIOD: bounds are out of bounds");
	}

	#define MARGIN 1.5

	for (k=0;k<n;k++) {											// for each col
	
		goodflag=0;
	
		// estimate global min and ajust threshold
		globalmin=GET(dp[k*m+lo]);
		for (j=lo;j<hi; j++) {									// for each row					
			z = GET(dp[k*m+j]);
			if (z<globalmin) {									// update min & period estimate
				globalmin = z;
			}
		}
		thresh2=globalmin*pratio;
		
		// estimate period
		p=lo;
		min=GET(dp[k*m+lo]);
		maxj=hi;
		for (j=lo;j<maxj; j++) {								// for each row					
		
			z = GET(dp[k*m+j]);
			
			if (z<min) {										// update min & period estimate
				min = z;
				p=j;
			}
			
			if (globalmin < thresh) {
				if (z<thresh) {									// good candidate, restrict search range
					maxj=j*MARGIN; 
					if (maxj> hi) { maxj = hi; }
					if (z<0.1) { goodflag=1; }
				}
			} else {
				if (z<thresh2) {								// good candidate, restrict search range
					maxj=j*MARGIN; 
					if (maxj> hi) { maxj = hi; }
					if (z<0.1) { goodflag=1; }
				}
			}
		}
		// mexPrintf("%d\n", p);

		SET(prdp[k])=p;
		SET(goodp[k])=goodflag;
	}
	
	return;						
}

void mexFunction(
				int nlhs, mxArray *plhs[],
				int nrhs, const mxArray *prhs[]
				)
{
	double	*dp, *bp, *tp, *rp, *prdp, *goodp, thresh, pratio;
	int m, n;
	
	/* Check for proper number of arguments */
	if (nrhs != 4){
		mexErrMsgTxt("DFTOPERIOD requires 4 input arguments");
	}
	
	/* Check type of input */
	if (!mxIsNumeric(D_IN) || mxIsComplex(D_IN) || 
		mxIsSparse(D_IN)	|| !mxIsDouble(D_IN) ) {
		mexErrMsgTxt("DFTOPERIOD: D should be doubles");
	}
	if (!mxIsNumeric(B_IN) || mxIsComplex(B_IN) || 
		mxIsSparse(B_IN)	|| !mxIsDouble(B_IN) ) {
		mexErrMsgTxt("DFTOPERIOD: B should be doubles");
	}
	if (!mxIsNumeric(T_IN) || mxIsComplex(T_IN) || 
		mxIsSparse(T_IN)	|| !mxIsDouble(T_IN) ) {
		mexErrMsgTxt("DFTOPERIOD: T should be double");
	}
	if (!mxIsNumeric(R_IN) || mxIsComplex(R_IN) || 
		mxIsSparse(R_IN)	|| !mxIsDouble(R_IN) ) {
		mexErrMsgTxt("DFTOPERIOD: R should be double");
	}
	m=mxGetM(D_IN);		// rows
	n=mxGetN(D_IN);		// columns

	if (mxGetM(B_IN) * mxGetN(B_IN) != 2) {
		mexErrMsgTxt("DFTOPERIOD: B should be a size 2 vector");
	}
	if (mxGetM(T_IN)*mxGetN(T_IN) != 1) {
		mexErrMsgTxt("DFTOPERIOD: T should be scalar");
	}
	if (mxGetM(R_IN)*mxGetN(R_IN) != 1) {
		mexErrMsgTxt("DFTOPERIOD: R should be scalar");
	}

	/* Create matrix to return */
	PRD_OUT = mxCreateDoubleMatrix(1,n, mxREAL);	 
	GOOD_OUT = mxCreateDoubleMatrix(1,n, mxREAL);	 

	dp = mxGetPr(D_IN);
	bp = mxGetPr(B_IN);
	tp = mxGetPr(T_IN); thresh=tp[0];
	rp = mxGetPr(R_IN); pratio=rp[0];
	prdp = mxGetPr(PRD_OUT);
	goodp = mxGetPr(GOOD_OUT);
	
	checkin_matrix((mxArray *) D_IN);
	checkin_matrix((mxArray *) B_IN);
	checkin_matrix((mxArray *) PRD_OUT);
	checkin_matrix((mxArray *) GOOD_OUT);
	
	/* Do the actual computations in a subroutine */
	ddftoperiods(dp,bp,prdp,goodp,thresh,pratio,m,n);

	checkout_matrix((mxArray *) D_IN);
	checkout_matrix((mxArray *) B_IN);
	checkout_matrix((mxArray *) PRD_OUT);
	checkout_matrix((mxArray *) GOOD_OUT);
	return;
}

