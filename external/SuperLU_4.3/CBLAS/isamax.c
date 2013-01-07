
/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

integer isamax_(integer *n, real *sx, integer *incx)
{


    /* System generated locals */
    integer ret_val, i__1;
    real r__1;

    /* Local variables */
    static real smax;
    static integer i, ix;


/*     finds the index of element having max. absolute value.   
       jack dongarra, linpack, 3/11/78.   
       modified 3/93 to return if incx .le. 0.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#define SX(I) sx[(I)-1]


    ret_val = 0;
    if (*n < 1 || *incx <= 0) {
	return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    ix = 1;
    smax = dabs(SX(1));
    ix += *incx;
    i__1 = *n;
    for (i = 2; i <= *n; ++i) {
	if ((r__1 = SX(ix), dabs(r__1)) <= smax) {
	    goto L5;
	}
	ret_val = i;
	smax = (r__1 = SX(ix), dabs(r__1));
L5:
	ix += *incx;
/* L10: */
    }
    return ret_val;

/*        code for increment equal to 1 */

L20:
    smax = dabs(SX(1));
    i__1 = *n;
    for (i = 2; i <= *n; ++i) {
	if ((r__1 = SX(i), dabs(r__1)) <= smax) {
	    goto L30;
	}
	ret_val = i;
	smax = (r__1 = SX(i), dabs(r__1));
L30:
	;
    }
    return ret_val;
} /* isamax_ */

