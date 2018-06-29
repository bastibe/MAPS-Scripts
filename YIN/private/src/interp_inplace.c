/* 
* interp_inplace.c
* linear interpolation
* 
* Alain de Cheveigné, CNRS/Ircam	
* (c) 2003 CNRS
*/

#include <math.h>
#include "mex.h"

/* #define MWRAP */
#include "mwrap_check.h"

/* Input Arguments */
#define X_IN	prhs[0]
#define IDX_IN	prhs[1]

/* Output Arguments */

static void interp_inplace(
				double *xp,		/* input vector */
				double *idxp,	/* index vector */
				int m, 			/* rows input vector */
				int midx		/* rows index vector */
			) 
{
	int j,idxi;
	double a,b,idx,idxf;
		
	
	for (j=0;j<midx;j++) {								/* for each index */
		idx=GET(idxp[j]);
		idxi=floor(idx);
		idxf=idx-idxi;
		if (idxi<0) {
			a=GET(xp[0]);
			b=GET(xp[1]);
			SET(idxp[j]) = a + idx*(b-a);
		} else if (idxi>=m) {
			a=GET(xp[m-2]);
			b=GET(xp[m-1]);
			SET(idxp[j]) = b + (idx-m+1)*(b-a);
		} else {
			a=GET(xp[idxi]);
			b=GET(xp[idxi+1]);
			SET(idxp[j]) = a + idxf*(b-a);
		}
	}
	
	return;						
}

void mexFunction(
				int nlhs, mxArray *plhs[],
				int nrhs, const mxArray *prhs[]
				)
{
	double	*xp, *idxp; 
	int nx, mx, nidx, midx;
	
	/* Check for proper number of arguments */
	if (nrhs !=2 ) {
		mexErrMsgTxt("INTERP_INPLACE takes 2 input arguments");
	}
	
	/* Check type of input */
	if (!mxIsNumeric(X_IN) || mxIsComplex(X_IN) || 
		mxIsSparse(X_IN)	|| !mxIsDouble(X_IN) ) {
		mexErrMsgTxt("INTERP_INPLACE: X should be doubles");
	}
	if (!mxIsNumeric(IDX_IN) || mxIsComplex(IDX_IN) || 
		mxIsSparse(IDX_IN)	|| !mxIsDouble(IDX_IN) ) {
		mexErrMsgTxt("INTERP_INPLACE: Y should be doubles");
	}
	mx=mxGetM(X_IN);		/* rows */
	nx=mxGetN(X_IN);		/* columns */
	midx=mxGetM(IDX_IN);		
	nidx=mxGetN(IDX_IN);		

	if (nx>1 || nidx>1) {
		mexErrMsgTxt("INTERP_INPLACE: X and IDX should be column vectors");
	}
	if (mx<1) {
		mexErrMsgTxt("INTERP_INPLACE: X should have at least two samples");
	}
		
	xp = mxGetPr(X_IN);
	idxp = mxGetPr(IDX_IN);
	
	checkin_matrix((mxArray *) IDX_IN);
	checkin_matrix((mxArray *) X_IN);

	/* Do the actual computations in a subroutine */
	interp_inplace(xp,idxp,mx,midx);

	checkout_matrix((mxArray *) X_IN);
	checkout_matrix((mxArray *) IDX_IN);
	return;
}

