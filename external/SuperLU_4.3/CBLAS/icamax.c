#include "f2c.h"

integer icamax_(integer *n, complex *cx, integer *incx)
{
    /* System generated locals */
    integer ret_val, i__1, i__2;
    real r__1, r__2;
    /* Builtin functions */
    double r_imag(complex *);
    /* Local variables */
    static real smax;
    static integer i, ix;
/*     finds the index of element having max. absolute value.   
       jack dongarra, linpack, 3/11/78.   
       modified 3/93 to return if incx .le. 0.   
       modified 12/3/93, array(1) declarations changed to array(*)   
    
   Parameter adjustments   
       Function Body */
#define CX(I) cx[(I)-1]
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
    smax = (r__1 = CX(1).r, dabs(r__1)) + (r__2 = r_imag(&CX(1)), dabs(r__2));
    ix += *incx;
    i__1 = *n;
    for (i = 2; i <= *n; ++i) {
	i__2 = ix;
	if ((r__1 = CX(ix).r, dabs(r__1)) + (r__2 = r_imag(&CX(ix)), dabs(
		r__2)) <= smax) {
	    goto L5;
	}
	ret_val = i;
	i__2 = ix;
	smax = (r__1 = CX(ix).r, dabs(r__1)) + (r__2 = r_imag(&CX(ix)), 
		dabs(r__2));
L5:
	ix += *incx;
/* L10: */
    }
    return ret_val;
/*        code for increment equal to 1 */
L20:
    smax = (r__1 = CX(1).r, dabs(r__1)) + (r__2 = r_imag(&CX(1)), dabs(r__2));
    i__1 = *n;
    for (i = 2; i <= *n; ++i) {
	i__2 = i;
	if ((r__1 = CX(i).r, dabs(r__1)) + (r__2 = r_imag(&CX(i)), dabs(
		r__2)) <= smax) {
	    goto L30;
	}
	ret_val = i;
	i__2 = i;
	smax = (r__1 = CX(i).r, dabs(r__1)) + (r__2 = r_imag(&CX(i)), dabs(
		r__2));
L30:
	;
    }
    return ret_val;
} /* icamax_ */


