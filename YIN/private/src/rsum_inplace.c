/* 
* rsum_inplace.c
* in-place running sum
* 
* Alain de Cheveigné, CNRS/Ircam	Jun 2002
* (c) 2002 CNRS
*/

/* Replaces each sample of the input matrix by the sum of itself and its N-1 
* neighbors. The operation is done in-place on the input argument.  The last
* N-1 samples of each row are invalid. */

#include <math.h>
#include "mex.h"

/* #define MWRAP */
#include "mwrap_check.h"

/* Input Arguments */
#define X_IN	prhs[0]
#define N_IN	prhs[1]

/* Output Arguments */

static void rsmth(
				double *xp,		/* matrix to sum */
				double N,		/* window size */
				int m, 			/* rows */ 
				int n	 		/* columns */
			) 
{
	int j,k, Ni;
	double Nr, tmp, sum, *x, *xx, *bump;
	
	Ni = (int) floor(N);
	Nr = N - (double) Ni;
	
	
	if (Nr == 0.0) {
		/* N is integer: simple summation */
		for (j=0;j<n;j++) {								/* for each column */
			/* prefill running sum */
			sum=0.0;
			x=&xp[j*m];
			bump=x+Ni-1;
			while (x<bump) {				/* for each row */
				sum += GET(*x);
				x++;
			}
			/* N */
			x=&xp[j*m];
			bump=x+m-(Ni-1);
			xx=x+Ni-1;
			while(x<bump) {
				tmp = GET(*x);
				sum += GET(*xx);
				SET(*x) = sum;
				sum -= tmp;
				x++; xx++;
			}
		}
	} else {
		/* N is fractionary: interpolate */
		for (j=0;j<n;j++) {								/* for each column */
			/* prefill running sum */
			sum=0.0;
			x=&xp[j*m];
			bump=x+Ni-1;
			while (x<bump) {				/* for each row */
				sum += GET(*x);
				x++;
			}
			/* N */
			x=&xp[j*m];
			bump=x+m-(Ni-1);
			xx=x+Ni-1;
			while(x<bump) {
				tmp = GET(*x);
				sum += GET(*xx);
				SET(*x) = sum + Nr * GET(*xx);
				sum -= tmp;
				x++; xx++;
			}
		}
	}	
	return;						
}

void mexFunction(
				int nlhs, mxArray *plhs[],
				int nrhs, const mxArray *prhs[]
				)
{
	double	*xp, *Np, N; 
	int n, m;
	
	/* Check for proper number of arguments */
	if (nrhs != 2) {
		mexErrMsgTxt("RSUM_INPLACE takes 2 input arguments");
	}
	
	/* Check type of input */
	if (!mxIsNumeric(X_IN) || mxIsComplex(X_IN) || 
		mxIsSparse(X_IN)	|| !mxIsDouble(X_IN) ) {
		mexErrMsgTxt("RSUM_INPLACE: X should be doubles");
	}
	if (!mxIsNumeric(N_IN) || mxIsComplex(N_IN) || 
		mxIsSparse(N_IN)	|| !mxIsDouble(N_IN) ) {
		mexErrMsgTxt("RSUM_INPLACE: N should be doubles");
	}
	m=mxGetM(X_IN);		/* rows */
	n=mxGetN(X_IN);		/* columns */

	if (mxGetM(N_IN)*mxGetN(N_IN) != 1) {
		mexErrMsgTxt("RSUM_INPLACE: N should be a scalar");
	}
			
	xp = mxGetPr(X_IN);
	Np = mxGetPr(N_IN);
	N = *Np;

	if (m < ceil(N)) {
		mexErrMsgTxt("RSUM_INPLACE: X should have at least N rows");
	}
		
	checkin_matrix((mxArray *) X_IN);
	
	/* Do the actual computations in a subroutine */
	rsmth(xp, N, m, n);

	checkout_matrix((mxArray *) X_IN);
	return;
}

