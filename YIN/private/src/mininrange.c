/* 
* minInRange.c
* for each sample, return index of the smallest value within an interval 
* of that sample
* 
* Alain de Cheveigné, CNRS/Ircam	Dec 2001
* (c) 2001 CNRS
*/

#include <math.h>
#include "mex.h"

/* #define MWRAP */
#include "mwrap_check.h"

/* Input Arguments */

#define	X_IN	       prhs[0]
#define	INTERVAL_IN	   prhs[1]

/* Output Arguments */

#define	IDX_OUT	plhs[0]

static void mininrange(
            double *xp,
            double *intp,
            double *idxp,
            unsigned int n
            )
{
    int h,k, idx, left, right, interval;
    double min;
	double max;
    
	
    for (k=0; k<n; k++) {
    	interval=(int) GET(intp[k]);

        left = k - interval/2;
        right = left+interval;  
        
		if (left<0) left=0; 	/* clip */
		if (right>n) right=n;
			
        min = GET(xp[k]);
        idx = k;
        for (h=left;h<right; h++) {
            if (GET(xp[h]) < min) {  /* update min */
                min = GET(xp[h]);
                idx = h;
            }
        }
        SET(idxp[k])=idx+1;
     }
     return;           
}

void mexFunction(
                 int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]
		 )
{
  double	*xp, *idxp, *intp;
  unsigned int	m,n;
  
  /* Check for proper number of arguments */
  
  if (nrhs != 2) {
    mexErrMsgTxt("MININRANGE requires two input arguments");
  } else if (nlhs != 1) {
    mexErrMsgTxt("MININRANGE requires one output argument");
  }
  
  /* Check type of input */
  
  if (!mxIsNumeric(X_IN) || mxIsComplex(X_IN) || 
      mxIsSparse(X_IN)  || !mxIsDouble(X_IN)) {
    mexErrMsgTxt("MININRANGE: X should be matrix of doubles");
  }
  if (!mxIsNumeric(INTERVAL_IN) || mxIsComplex(INTERVAL_IN) || 
      mxIsSparse(INTERVAL_IN)  || !mxIsDouble(INTERVAL_IN)) {
    mexErrMsgTxt("MININRANGE X should be matrix of doubles");
  }

  m = mxGetM(X_IN); /* rows */
  n = mxGetN(X_IN); /* columns */ 
  if (m>1 || n <=1) {
      mexErrMsgTxt("MININRANGE X should be row vector");
  }
  if (m != mxGetM(INTERVAL_IN) || n != mxGetN(INTERVAL_IN)) {
	  mexErrMsgTxt("MININRANGE: INTERVAL should be of same size as X");
  }
  
  /* Create matrix to return */
  
  IDX_OUT = mxCreateDoubleMatrix(1, n, mxREAL);  
  
  /* Assign pointers to the various parameters */
  
  xp = mxGetPr(X_IN);
  intp = mxGetPr(INTERVAL_IN);
  idxp = mxGetPr(IDX_OUT);
  
  checkin_matrix((mxArray *) X_IN);
  checkin_matrix((mxArray *) INTERVAL_IN);
  checkin_matrix(IDX_OUT);
  
  /* Do the actual computations in a subroutine */
  
  mininrange(xp,intp,idxp,n);

  checkout_matrix((mxArray *) X_IN);
  checkout_matrix((mxArray *) INTERVAL_IN);
  checkout_matrix(IDX_OUT);
  return;
}


