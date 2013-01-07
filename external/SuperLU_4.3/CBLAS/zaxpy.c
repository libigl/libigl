
/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int zaxpy_(integer *n, doublecomplex *za, doublecomplex *zx, 
	integer *incx, doublecomplex *zy, integer *incy)
{


    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublecomplex z__1, z__2;

    /* Local variables */
    static integer i;
    extern doublereal dcabs1_(doublecomplex *);
    static integer ix, iy;


/*     constant times a vector plus a vector.   
       jack dongarra, 3/11/78.   
       modified 12/3/93, array(1) declarations changed to array(*)   

    
   Parameter adjustments   
       Function Body */
#define ZY(I) zy[(I)-1]
#define ZX(I) zx[(I)-1]


    if (*n <= 0) {
	return 0;
    }
    if (dcabs1_(za) == 0.) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments   
            not equal to 1 */

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
	i__2 = iy;
	i__3 = iy;
	i__4 = ix;
	z__2.r = za->r * ZX(ix).r - za->i * ZX(ix).i, z__2.i = za->r * ZX(
		ix).i + za->i * ZX(ix).r;
	z__1.r = ZY(iy).r + z__2.r, z__1.i = ZY(iy).i + z__2.i;
	ZY(iy).r = z__1.r, ZY(iy).i = z__1.i;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1 */

L20:
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	i__2 = i;
	i__3 = i;
	i__4 = i;
	z__2.r = za->r * ZX(i).r - za->i * ZX(i).i, z__2.i = za->r * ZX(
		i).i + za->i * ZX(i).r;
	z__1.r = ZY(i).r + z__2.r, z__1.i = ZY(i).i + z__2.i;
	ZY(i).r = z__1.r, ZY(i).i = z__1.i;
/* L30: */
    }
    return 0;
} /* zaxpy_ */

