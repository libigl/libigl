
/*
 * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 *
 */
#include <math.h>
#include "slu_sdefs.h"

int sgst04(int n, int nrhs, float *x, int ldx, float *xact,
	      int ldxact, float rcond, float *resid)
{
/*
    Purpose   
    =======   

    SGST04 computes the difference between a computed solution and the   
    true solution to a system of linear equations.   
    RESID =  ( norm(X-XACT) * RCOND ) / ( norm(XACT) * EPS ),   
    where RCOND is the reciprocal of the condition number and EPS is the 
    machine epsilon.   

    Arguments   
    =========   

    N       (input) INT   
            The number of rows of the matrices X and XACT.  N >= 0.   

    NRHS    (input) INT   
            The number of columns of the matrices X and XACT.  NRHS >= 0. 

    X       (input) FLOAT PRECISION array, dimension (LDX,NRHS)   
            The computed solution vectors.  Each vector is stored as a   
            column of the matrix X.   

    LDX     (input) INT   
            The leading dimension of the array X.  LDX >= max(1,N).   

    XACT    (input) FLOAT PRECISION array, dimension( LDX, NRHS )   
            The exact solution vectors.  Each vector is stored as a   
            column of the matrix XACT.   

    LDXACT  (input) INT   
            The leading dimension of the array XACT.  LDXACT >= max(1,N). 

    RCOND   (input) FLOAT PRECISION   
            The reciprocal of the condition number of the coefficient   
            matrix in the system of equations.   

    RESID   (output) FLOAT PRECISION   
            The maximum over the NRHS solution vectors of   
            ( norm(X-XACT) * RCOND ) / ( norm(XACT) * EPS )   

    ===================================================================== 
*/
    /* Table of constant values */
    int c__1 = 1;

    /* System generated locals */
    float d__1, d__2, d__3, d__4;

    /* Local variables */
    int    i, j, n__1;
    int    ix;
    float xnorm;
    float eps;
    float diffnm;

    /* Function prototypes */
    extern int isamax_(int *, float *, int *);

    /* Quick exit if N = 0 or NRHS = 0. */
   if ( n <= 0 || nrhs <= 0 ) {
	*resid = 0.;
	return 0;
    }

    /* Exit with RESID = 1/EPS if RCOND is invalid. */

    eps = slamch_("Epsilon");
    if ( rcond < 0. ) {
	*resid = 1. / eps;
	return 0;
    }

    /* Compute the maximum of norm(X - XACT) / ( norm(XACT) * EPS )   
       over all the vectors X and XACT . */

    *resid = 0.;
    for (j = 0; j < nrhs; ++j) {
	n__1 = n;
	ix = isamax_(&n__1, &xact[j*ldxact], &c__1);
	xnorm = (d__1 = xact[ix-1 + j*ldxact], fabs(d__1));

	diffnm = 0.;
	for (i = 0; i < n; ++i) {
	    /* Computing MAX */
	    d__3 = diffnm;
	    d__4 = (d__1 = x[i+j*ldx]-xact[i+j*ldxact], fabs(d__1));
	    diffnm = SUPERLU_MAX(d__3,d__4);
	}
	if (xnorm <= 0.) {
	    if (diffnm > 0.) {
		*resid = 1. / eps;
	    }
	} else {
	    /* Computing MAX */
	    d__1 = *resid, d__2 = diffnm / xnorm * rcond;
	    *resid = SUPERLU_MAX(d__1,d__2);
	}
    }
    if (*resid * eps < 1.) {
	*resid /= eps;
    }

    return 0;

} /* sgst04_ */
