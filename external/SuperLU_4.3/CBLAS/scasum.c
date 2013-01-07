
/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

real scasum_(integer *n, complex *cx, integer *incx)
{


    /* System generated locals */
    integer i__1, i__2, i__3;
    real ret_val, r__1, r__2;

    /* Builtin functions */
    double r_imag(complex *);

    /* Local variables */
    static integer i, nincx;
    static real stemp;


/*     takes the sum of the absolute values of a complex vector and   
       returns a single precision result.   
       jack dongarra, linpack, 3/11/78.   
       modified 3/93 to return if incx .le. 0.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#define CX(I) cx[(I)-1]


    ret_val = 0.f;
    stemp = 0.f;
    if (*n <= 0 || *incx <= 0) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i = 1; *incx < 0 ? i >= nincx : i <= nincx; i += *incx) {
	i__3 = i;
	stemp = stemp + (r__1 = CX(i).r, dabs(r__1)) + (r__2 = r_imag(&CX(
		i)), dabs(r__2));
/* L10: */
    }
    ret_val = stemp;
    return ret_val;

/*        code for increment equal to 1 */

L20:
    i__2 = *n;
    for (i = 1; i <= *n; ++i) {
	i__1 = i;
	stemp = stemp + (r__1 = CX(i).r, dabs(r__1)) + (r__2 = r_imag(&CX(
		i)), dabs(r__2));
/* L30: */
    }
    ret_val = stemp;
    return ret_val;
} /* scasum_ */

