#include "f2c.h"

/* Subroutine */ int clacgv_(integer *n, complex *x, integer *incx)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    CLACGV conjugates a complex vector of length N.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The length of the vector X.  N >= 0.   

    X       (input/output) COMPLEX array, dimension   
                           (1+(N-1)*abs(INCX))   
            On entry, the vector of length N to be conjugated.   
            On exit, X is overwritten with conjg(X).   

    INCX    (input) INTEGER   
            The spacing between successive elements of X.   

   ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer i__1, i__2;
    complex q__1;
    /* Builtin functions */
    void r_cnjg(complex *, complex *);
    /* Local variables */
    static integer ioff, i;


#define X(I) x[(I)-1]


    if (*incx == 1) {
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    i__2 = i;
	    r_cnjg(&q__1, &X(i));
	    X(i).r = q__1.r, X(i).i = q__1.i;
/* L10: */
	}
    } else {
	ioff = 1;
	if (*incx < 0) {
	    ioff = 1 - (*n - 1) * *incx;
	}
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    i__2 = ioff;
	    r_cnjg(&q__1, &X(ioff));
	    X(ioff).r = q__1.r, X(ioff).i = q__1.i;
	    ioff += *incx;
/* L20: */
	}
    }
    return 0;

/*     End of CLACGV */

} /* clacgv_ */

