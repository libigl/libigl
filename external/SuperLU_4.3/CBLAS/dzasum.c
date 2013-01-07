
/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

doublereal dzasum_(integer *n, doublecomplex *zx, integer *incx)
{


    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer i;
    static doublereal stemp;
    extern doublereal dcabs1_(doublecomplex *);
    static integer ix;


/*     takes the sum of the absolute values.   
       jack dongarra, 3/11/78.   
       modified 3/93 to return if incx .le. 0.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#define ZX(I) zx[(I)-1]


    ret_val = 0.;
    stemp = 0.;
    if (*n <= 0 || *incx <= 0) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    ix = 1;
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	stemp += dcabs1_(&ZX(ix));
	ix += *incx;
/* L10: */
    }
    ret_val = stemp;
    return ret_val;

/*        code for increment equal to 1 */

L20:
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	stemp += dcabs1_(&ZX(i));
/* L30: */
    }
    ret_val = stemp;
    return ret_val;
} /* dzasum_ */

