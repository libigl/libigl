
/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

doublereal dasum_(integer *n, doublereal *dx, integer *incx)
{


    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1, d__2, d__3, d__4, d__5, d__6;

    /* Local variables */
    static integer i, m;
    static doublereal dtemp;
    static integer nincx, mp1;


/*     takes the sum of the absolute values.   
       jack dongarra, linpack, 3/11/78.   
       modified 3/93 to return if incx .le. 0.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#define DX(I) dx[(I)-1]


    ret_val = 0.;
    dtemp = 0.;
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
	dtemp += (d__1 = DX(i), abs(d__1));
/* L10: */
    }
    ret_val = dtemp;
    return ret_val;

/*        code for increment equal to 1   


          clean-up loop */

L20:
    m = *n % 6;
    if (m == 0) {
	goto L40;
    }
    i__2 = m;
    for (i = 1; i <= m; ++i) {
	dtemp += (d__1 = DX(i), abs(d__1));
/* L30: */
    }
    if (*n < 6) {
	goto L60;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i = mp1; i <= *n; i += 6) {
	dtemp = dtemp + (d__1 = DX(i), abs(d__1)) + (d__2 = DX(i + 1), abs(
		d__2)) + (d__3 = DX(i + 2), abs(d__3)) + (d__4 = DX(i + 3), 
		abs(d__4)) + (d__5 = DX(i + 4), abs(d__5)) + (d__6 = DX(i + 5)
		, abs(d__6));
/* L50: */
    }
L60:
    ret_val = dtemp;
    return ret_val;
} /* dasum_ */

