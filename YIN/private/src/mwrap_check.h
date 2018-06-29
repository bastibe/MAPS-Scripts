/* 
* stuff to allow checking with mwrap
* 
*/

#include "mwrap.h"
#include "string.h"

#define MACINTOSH	/* compiler (?) converts \n to \r: convert them back */

void checkin_matrix(mxArray *m);
void checkout_matrix(mxArray *m);
void mex_messagefunction(int level, char* message);
void mex_errorfunction(int level);

/* check matrix into mwrap's tree */
void checkin_matrix(mxArray *m) 
{
	char *base, *top;
	
	base = (char *) mxGetPr(m);
	top = base + mxGetN(m) * mxGetM(m) * sizeof(double);
	/* mexPrintf("%d %d\n", base, top); */
	CHECKIN(base, top); 
}

/* check matrix out of mwrap's tree */
void checkout_matrix(mxArray *m) 
{
	char *base;
	
	base = (char *) mxGetPr(m);
	CHECKOUT(base); 
}

/* message function to give to mwrap */
void mex_messagefunction(int level, char* message) {
#ifdef MACINTOSH
	{char *c;
		c = strchr(message, (int) '\r');
	
		while(c) {
			*c = '\n';
			c = strchr(c, (int) '\r');
		}
	}
#endif
	mexPrintf(message); 
}

/* how to exit mex function */
void mex_errorfunction(int level) {
	/* mexPrintf("\n"); */
	mexErrMsgTxt(" "); 
}
