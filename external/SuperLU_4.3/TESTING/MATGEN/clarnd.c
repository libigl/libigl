/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Complex */ VOID clarnd_(complex * ret_val, integer *idist, integer *iseed)
{
    /* System generated locals */
    doublereal d__1, d__2;
    complex q__1, q__2, q__3;

    /* Builtin functions */
    double log(doublereal), sqrt(doublereal);
    void c_exp(complex *, complex *);

    /* Local variables */
    static real t1, t2;
    extern doublereal slaran_(integer *);


/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    CLARND returns a random complex number from a uniform or normal   
    distribution.   

    Arguments   
    =========   

    IDIST   (input) INTEGER   
            Specifies the distribution of the random numbers:   
            = 1:  real and imaginary parts each uniform (0,1)   
            = 2:  real and imaginary parts each uniform (-1,1)   
            = 3:  real and imaginary parts each normal (0,1)   
            = 4:  uniformly distributed on the disc abs(z) <= 1   
            = 5:  uniformly distributed on the circle abs(z) = 1   

    ISEED   (input/output) INTEGER array, dimension (4)   
            On entry, the seed of the random number generator; the array 
  
            elements must be between 0 and 4095, and ISEED(4) must be   
            odd.   
            On exit, the seed is updated.   

    Further Details   
    ===============   

    This routine calls the auxiliary routine SLARAN to generate a random 
  
    real number from a uniform (0,1) distribution. The Box-Muller method 
  
    is used to transform numbers from a uniform to a normal distribution. 
  

    ===================================================================== 
  


       Generate a pair of real random numbers from a uniform (0,1)   
       distribution   

       Parameter adjustments */
    --iseed;

    /* Function Body */
    t1 = slaran_(&iseed[1]);
    t2 = slaran_(&iseed[1]);

    if (*idist == 1) {

/*        real and imaginary parts each uniform (0,1) */

	q__1.r = t1, q__1.i = t2;
	 ret_val->r = q__1.r,  ret_val->i = q__1.i;
    } else if (*idist == 2) {

/*        real and imaginary parts each uniform (-1,1) */

	d__1 = t1 * 2.f - 1.f;
	d__2 = t2 * 2.f - 1.f;
	q__1.r = d__1, q__1.i = d__2;
	 ret_val->r = q__1.r,  ret_val->i = q__1.i;
    } else if (*idist == 3) {

/*        real and imaginary parts each normal (0,1) */

	d__1 = sqrt(log(t1) * -2.f);
	d__2 = t2 * 6.2831853071795864769252867663f;
	q__3.r = 0.f, q__3.i = d__2;
	c_exp(&q__2, &q__3);
	q__1.r = d__1 * q__2.r, q__1.i = d__1 * q__2.i;
	 ret_val->r = q__1.r,  ret_val->i = q__1.i;
    } else if (*idist == 4) {

/*        uniform distribution on the unit disc abs(z) <= 1 */

	d__1 = sqrt(t1);
	d__2 = t2 * 6.2831853071795864769252867663f;
	q__3.r = 0.f, q__3.i = d__2;
	c_exp(&q__2, &q__3);
	q__1.r = d__1 * q__2.r, q__1.i = d__1 * q__2.i;
	 ret_val->r = q__1.r,  ret_val->i = q__1.i;
    } else if (*idist == 5) {

/*        uniform distribution on the unit circle abs(z) = 1 */

	d__1 = t2 * 6.2831853071795864769252867663f;
	q__2.r = 0.f, q__2.i = d__1;
	c_exp(&q__1, &q__2);
	 ret_val->r = q__1.r,  ret_val->i = q__1.i;
    }
    return ;

/*     End of CLARND */

} /* clarnd_ */

