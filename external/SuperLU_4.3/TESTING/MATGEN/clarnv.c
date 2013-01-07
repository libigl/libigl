#include "f2c.h"

/* Subroutine */ int clarnv_(integer *idist, integer *iseed, integer *n, 
	complex *x)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    CLARNV returns a vector of n random complex numbers from a uniform or 
  
    normal distribution.   

    Arguments   
    =========   

    IDIST   (input) INTEGER   
            Specifies the distribution of the random numbers:   
            = 1:  real and imaginary parts each uniform (0,1)   
            = 2:  real and imaginary parts each uniform (-1,1)   
            = 3:  real and imaginary parts each normal (0,1)   
            = 4:  uniformly distributed on the disc abs(z) < 1   
            = 5:  uniformly distributed on the circle abs(z) = 1   

    ISEED   (input/output) INTEGER array, dimension (4)   
            On entry, the seed of the random number generator; the array 
  
            elements must be between 0 and 4095, and ISEED(4) must be   
            odd.   
            On exit, the seed is updated.   

    N       (input) INTEGER   
            The number of random numbers to be generated.   

    X       (output) COMPLEX array, dimension (N)   
            The generated random numbers.   

    Further Details   
    ===============   

    This routine calls the auxiliary routine SLARUV to generate random   
    real numbers from a uniform (0,1) distribution, in batches of up to   
    128 using vectorisable code. The Box-Muller method is used to   
    transform numbers from a uniform to a normal distribution.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2;
    complex q__1, q__2, q__3;
    /* Builtin functions */
    double log(doublereal), sqrt(doublereal);
    void c_exp(complex *, complex *);
    /* Local variables */
    static integer i;
    static real u[128];
    static integer il, iv;
    extern /* Subroutine */ int slaruv_(integer *, integer *, real *);


#define X(I) x[(I)-1]
#define ISEED(I) iseed[(I)-1]


    i__1 = *n;
    for (iv = 1; iv <= *n; iv += 64) {
/* Computing MIN */
	i__2 = 64, i__3 = *n - iv + 1;
	il = min(i__2,i__3);

/*        Call SLARUV to generate 2*IL real numbers from a uniform (0,
1)   
          distribution (2*IL <= LV) */

	i__2 = il << 1;
	slaruv_(&ISEED(1), &i__2, u);

	if (*idist == 1) {

/*           Copy generated numbers */

	    i__2 = il;
	    for (i = 1; i <= il; ++i) {
		i__3 = iv + i - 1;
		i__4 = (i << 1) - 2;
		i__5 = (i << 1) - 1;
		q__1.r = u[(i<<1)-2], q__1.i = u[(i<<1)-1];
		X(iv+i-1).r = q__1.r, X(iv+i-1).i = q__1.i;
/* L10: */
	    }
	} else if (*idist == 2) {

/*           Convert generated numbers to uniform (-1,1) distribut
ion */

	    i__2 = il;
	    for (i = 1; i <= il; ++i) {
		i__3 = iv + i - 1;
		d__1 = u[(i << 1) - 2] * 2.f - 1.f;
		d__2 = u[(i << 1) - 1] * 2.f - 1.f;
		q__1.r = d__1, q__1.i = d__2;
		X(iv+i-1).r = q__1.r, X(iv+i-1).i = q__1.i;
/* L20: */
	    }
	} else if (*idist == 3) {

/*           Convert generated numbers to normal (0,1) distributio
n */

	    i__2 = il;
	    for (i = 1; i <= il; ++i) {
		i__3 = iv + i - 1;
		d__1 = sqrt(log(u[(i << 1) - 2]) * -2.f);
		d__2 = u[(i << 1) - 1] * 6.2831853071795864769252867663f;
		q__3.r = 0.f, q__3.i = d__2;
		c_exp(&q__2, &q__3);
		q__1.r = d__1 * q__2.r, q__1.i = d__1 * q__2.i;
		X(iv+i-1).r = q__1.r, X(iv+i-1).i = q__1.i;
/* L30: */
	    }
	} else if (*idist == 4) {

/*           Convert generated numbers to complex numbers uniforml
y   
             distributed on the unit disk */

	    i__2 = il;
	    for (i = 1; i <= il; ++i) {
		i__3 = iv + i - 1;
		d__1 = sqrt(u[(i << 1) - 2]);
		d__2 = u[(i << 1) - 1] * 6.2831853071795864769252867663f;
		q__3.r = 0.f, q__3.i = d__2;
		c_exp(&q__2, &q__3);
		q__1.r = d__1 * q__2.r, q__1.i = d__1 * q__2.i;
		X(iv+i-1).r = q__1.r, X(iv+i-1).i = q__1.i;
/* L40: */
	    }
	} else if (*idist == 5) {

/*           Convert generated numbers to complex numbers uniforml
y   
             distributed on the unit circle */

	    i__2 = il;
	    for (i = 1; i <= il; ++i) {
		i__3 = iv + i - 1;
		d__1 = u[(i << 1) - 1] * 6.2831853071795864769252867663f;
		q__2.r = 0.f, q__2.i = d__1;
		c_exp(&q__1, &q__2);
		X(iv+i-1).r = q__1.r, X(iv+i-1).i = q__1.i;
/* L50: */
	    }
	}
/* L60: */
    }
    return 0;

/*     End of CLARNV */

} /* clarnv_ */

