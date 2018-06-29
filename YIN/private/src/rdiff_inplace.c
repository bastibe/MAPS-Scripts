/* 
* rdiff_inplace.c
* running difference function
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
#define Y_IN	prhs[1]
#define R_IN	prhs[2]
#define L_IN 	prhs[3]
#define N_IN	prhs[4]

/* Output Arguments */

static void rdiff_inplace(
				double *xp,		/* input vector 1 */
				double *yp,		/* input vector 2 */
				double *rp,		/* df matrix (time X lag) */
				double *lp,		/* lag matrix */
				int N,			/* size of summation window */
				int m, 			/* rows in corr matrix */
				int n	 		/* columns in corr matrix */
			) 
{
	int j,k,i, xlag, ylag;
	double z, d, *r, *x, *y, *bump;
		
	
	for (j=0;j<n;j++) {									/* for each column (lag pair) */
		xlag = GET(lp[2*j]);
		ylag = GET(lp[2*j+1]);
		if (N==1) {
			r=&rp[j*m];
			x=&xp[xlag];
			y=&yp[ylag];
			bump=r+m;
			while(r<bump) {								/* for each row */
				d = GET(*x)-GET(*y);
				SET(*r) = d*d;	
				r++; x++; y++;
			}
		} else {
			for (k=0; k<m; k++) {						/* for each row */
				z=0;
				x=&xp[xlag+k*N];
				y=&yp[ylag+k*N];
				bump=x+N;
				while(x<bump) {
					d = GET(*x)-GET(*y);
					z+=d*d;
					x++; y++;
				}
				SET(rp[j*m+k]) = z;
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
	double	*xp, *yp, *rp, *lp, *up, *np; 
	int nx, mx, ny, my, nr, mr, nl, ml, xmax, ymax, j, N, xlag, ylag;
	
	/* Check for proper number of arguments */
	if (nrhs < 4 || nrhs > 5) {
		mexErrMsgTxt("RDIFF_INPLACE takes 4 or 5 input arguments");
	}
	
	/* Check type of input */
	if (!mxIsNumeric(X_IN) || mxIsComplex(X_IN) || 
		mxIsSparse(X_IN)	|| !mxIsDouble(X_IN) ) {
		mexErrMsgTxt("RDIFF_INPLACE: X should be doubles");
	}
	if (!mxIsNumeric(Y_IN) || mxIsComplex(Y_IN) || 
		mxIsSparse(Y_IN)	|| !mxIsDouble(Y_IN) ) {
		mexErrMsgTxt("RDIFF_INPLACE: Y should be doubles");
	}
	if (!mxIsNumeric(R_IN) || mxIsComplex(R_IN) || 
		mxIsSparse(R_IN)	|| !mxIsDouble(R_IN) ) {
		mexErrMsgTxt("RDIFF_INPLACE: R should be doubles");
	}
	if (nrhs==5) {
		if (!mxIsNumeric(N_IN) || mxIsComplex(N_IN) || 
			mxIsSparse(N_IN)	 || !mxIsDouble(N_IN) ) {
			mexErrMsgTxt("RDIFF_INPLACE: N should double");
		}
		np = mxGetPr(N_IN);
		N=*np;
	} else {
		N=1;
	}
	mx=mxGetM(X_IN);		/* rows */
	nx=mxGetN(X_IN);		/* columns */
	my=mxGetM(Y_IN);		
	ny=mxGetN(Y_IN);		
	mr=mxGetM(R_IN);		
	nr=mxGetN(R_IN);		
	ml=mxGetM(L_IN);		
	nl=mxGetN(L_IN);		

	if (nx>1 || ny>1) {
		mexErrMsgTxt("RDIFF_INPLACE: X and Y should be column vectors");
	}
	if (nr !=nl) {
		mexErrMsgTxt("RDIFF_INPLACE: result and lag matrices should have same number of columns");
	}
	if ( ml != 2) {
		mexErrMsgTxt("RDIFF_INPLACE: lags should be a two row matrix");
	}
	if (!nx || !mx || !ny || !my || !nr || !mr || !nl || !ml) {
		mexErrMsgTxt("RDIFF_INPLACE: one of the input matrices is empty");
	}	
		
	xp = mxGetPr(X_IN);
	yp = mxGetPr(Y_IN);
	rp = mxGetPr(R_IN);
	lp = mxGetPr(L_IN);
	
	checkin_matrix((mxArray *) L_IN);
			
	/* find maximum lag for x and y and check that no lag is negative */
	xmax=GET(lp[0]);
	ymax=GET(lp[1]);
	for (j=1;j<nl;j++) {
		xlag = GET(lp[2*j]);
		ylag = GET(lp[2*j+1]);
		if ( xlag<0 || ylag<0 ) {
			mexErrMsgTxt("RDIFF_INPLACE: LAGS should be non-negative");
		}
		if (xlag>xmax ) xmax=xlag;
		if (ylag>ymax ) ymax=ylag;
	}

	if ( (N*mr+ymax)> my) {
		mexPrintf("mr*N: %d, ymax: %d, my: %d\n", mr*N, ymax, my);
		mexErrMsgTxt("RDIFF_INPLACE: Y data too short");
	}
	if ( (N*mr+xmax)> mx) {
		mexPrintf("mr*N: %d, xmax: %d, mx: %d\n", mr*N, xmax, mx);
		mexErrMsgTxt("RDIFF_INPLACE: X data too short");
	}
	
	checkin_matrix((mxArray *) X_IN);
	if (yp!=xp) checkin_matrix((mxArray *) Y_IN);
	checkin_matrix((mxArray *) R_IN);
	
	/* Do the actual computations in a subroutine */
	rdiff_inplace(xp,yp,rp,lp,N,mr,nr);

	checkout_matrix((mxArray *) X_IN);
	if (yp!=xp) checkout_matrix((mxArray *) Y_IN); 
	checkout_matrix((mxArray *) R_IN);
	checkout_matrix((mxArray *) L_IN);
	return;
}

