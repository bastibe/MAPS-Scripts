/* 
* minInterp.c
* determine position of minimum of parabolic fit to local minima
* 
* Alain de Cheveigné, CNRS/Ircam	Dec 2001
* (c) 2001 CNRS
*/

#include <math.h>
#include "mex.h"

/* #define MWRAP */
#include "mwrap_check.h"

/* Input Arguments */
#define	X_IN	prhs[0]


/* Output Arguments */
#define	MIN_OUT	plhs[0]
#define	IND_OUT	plhs[1]


static void mininterp(
            double *xp,
            unsigned int m,
            unsigned int n,
            double *minp,
            double *indp
		   )
{
    double	x1,x2,x3,a,b,min,ind;
    int i,j,k;
  
  	/* each column is interpolated independently */
  	/* each local minimum is replaced by the value of the interpolated minimum, and the 
  	position of the interplolated minimum is recorded in indp */
  	
    for (j=0; j<n; j++) { 
        /* first row: don't touch */
        SET(indp[j*m])=0;         
        SET(minp[j*m])=GET(xp[j*m]);
        /* intermediate: refine */
        for (k=1; k<m-1; k++) {   
            i=k+j*m;
			x1=GET(xp[i-1]); x2=GET(xp[i]); x3=GET(xp[i+1]);  
	        SET(indp[i]) = 0;	/* default */
	        SET(minp[i]) = x2; 
		    if ((x2<=x1) && (x2<=x3)) { 
			    a = (x1 + x3 - 2*x2)/2;
			    b = (x3 - x1)/2;
			    if (a) {                     
	                SET(indp[i]) = - b / (2*a);
	                SET(minp[i]) = x2 +  - (b*b) / (4*a);
	            }
	        }
	   	}
        /* last row: don't touch */
        SET(indp[m-1 + j*m])=0;
        SET(minp[m-1 + j*m])=GET(xp[m-1 + j*m]);
   }
    return;
}

void mexFunction(
                 int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]
		 )
{
  double	*minp, *indp, *xp;
  unsigned int	m,n;
  
  /* Check for proper number of arguments */
  
  if (nrhs != 1) {
    mexErrMsgTxt("minParabolic requires one input argument.");
  } else if (nlhs > 2) {
    mexErrMsgTxt("minParabolic requires one or two output arguments");
  }
  
  
  /* Check type of input */
  
  if (!mxIsNumeric(X_IN) || mxIsComplex(X_IN) || 
      mxIsSparse(X_IN)  || !mxIsDouble(X_IN)) {
    mexErrMsgTxt("minParabolic requires that X be a matrix of doubles");
  }

  m = mxGetM(X_IN); /* rows */
  n = mxGetN(X_IN); /* columns */
  if (m<3) { mexErrMsgTxt("number of rows should be 3 or more"); }
  
  
  /* Create two matrices to return */
  
  MIN_OUT = mxCreateDoubleMatrix(m, n, mxREAL);
  IND_OUT = mxCreateDoubleMatrix(m, n, mxREAL);
  
  
  /* Assign pointers to the various parameters */
  
  minp = mxGetPr(MIN_OUT);
  indp = mxGetPr(IND_OUT);
  xp = mxGetPr(X_IN);  
  
  checkin_matrix((mxArray *) X_IN);
  checkin_matrix(MIN_OUT);
  checkin_matrix(IND_OUT);
  
  /* Do the actual computations in a subroutine */
  

  mininterp(xp,m,n,minp,indp);

  checkout_matrix((mxArray *) X_IN);
  checkout_matrix(MIN_OUT);
  checkout_matrix(IND_OUT);
  return;
}


