
/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int srot_(integer *n, real *sx, integer *incx, real *sy, 
	integer *incy, real *c, real *s)
{


    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i;
    static real stemp;
    static integer ix, iy;


/*     applies a plane rotation.   
       jack dongarra, linpack, 3/11/78.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#define SY(I) sy[(I)-1]
#define SX(I) sx[(I)-1]


    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*       code for unequal increments or equal increments not equal   
           to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	stemp = *c * SX(ix) + *s * SY(iy);
	SY(iy) = *c * SY(iy) - *s * SX(ix);
	SX(ix) = stemp;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*       code for both increments equal to 1 */

L20:
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	stemp = *c * SX(i) + *s * SY(i);
	SY(i) = *c * SY(i) - *s * SX(i);
	SX(i) = stemp;
/* L30: */
    }
    return 0;
} /* srot_ */

