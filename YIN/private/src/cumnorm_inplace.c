/* 
* cumnorm_inplace.c
* cumulative mean-normalization of diff function
* 
* Alain de Cheveigné, CNRS/Ircam	Dec 2001
* (c) 2001 CNRS
*/

#include <math.h>
#include "mex.h"

/* #define MWRAP */
#include "mwrap_check.h"

/* Input Arguments */
#define X_IN	prhs[0]

/* Output Arguments */

static void cumsum_inplace(
				double *xp,		/* matrix to cumulative-mean-normalize (along columns) */
				int m, 			/* rows */
				int n	 		/* columns */ 
			) 
{
	int j,k;
	double z, sum, mean;
		
	
	for (k=0; k<n; k++) {							/* for each column */
		SET(xp[k*m+0])=1;
		sum = 0.0;
		for (j=1;j<m;j++) {							/* for each row */
			z = GET(xp[k*m+j]);
			z = (z > 0) ? z:0;						/* clip to remove numerical artifacts */
			sum += z;
			z = (sum>0) ? (z / (sum/j)) : 1;		/* cumulative-mean-normalize */
			SET(xp[k*m+j])=z;
		}
	}
/*
	for (j=0;j<m;j++) {									// for each row (time) 
		sum=0;
		SET(xp[j])=1;
		for (k=1; k<n; k++) {							// for each column (lag)
			z=GET(xp[m*k+j]);
			sum+=z;
			mean=sum/k;
			SET(xp[m*k+j]) = (mean > 0) ? z/mean : 1;
		}
	}
*/
		
	
	return;						
}

void mexFunction(
				int nlhs, mxArray *plhs[],
				int nrhs, const mxArray *prhs[]
				)
{
	double	*xp; 
	int nx, mx;
	
	/* Check for proper number of arguments */
	if (nrhs != 1) {
		mexErrMsgTxt("CUMSUM_INPLACE takes 1 input argument");
	}
	
	/* Check type of input */
	if (!mxIsNumeric(X_IN) || mxIsComplex(X_IN) || 
		mxIsSparse(X_IN)	|| !mxIsDouble(X_IN) ) {
		mexErrMsgTxt("CUMSUM_INPLACE: X should be doubles");
	}
	mx=mxGetM(X_IN);		/* rows */
	nx=mxGetN(X_IN);		/* columns */
	if (nx*mx == 0) {
		mexErrMsgTxt("CUMSUM_INPLACE: input matrix is empty");
	}	
		
	xp = mxGetPr(X_IN);	
	checkin_matrix((mxArray *) X_IN);
			
	
	/* Do the actual computations in a subroutine */
	cumsum_inplace(xp,mx,nx);

	checkout_matrix((mxArray *) X_IN);
	return;
}

