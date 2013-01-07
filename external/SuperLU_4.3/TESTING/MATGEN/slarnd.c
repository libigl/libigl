/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

doublereal slarnd_(integer *idist, integer *iseed)
{
    /* System generated locals */
    real ret_val;

    /* Builtin functions */
    double log(doublereal), sqrt(doublereal), cos(doublereal);

    /* Local variables */
    static real t1, t2;
    extern doublereal slaran_(integer *);


/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    SLARND returns a random real number from a uniform or normal   
    distribution.   

    Arguments   
    =========   

    IDIST   (input) INTEGER   
            Specifies the distribution of the random numbers:   
            = 1:  uniform (0,1)   
            = 2:  uniform (-1,1)   
            = 3:  normal (0,1)   

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
  


       Generate a real random number from a uniform (0,1) distribution   

       Parameter adjustments */
    --iseed;

    /* Function Body */
    t1 = slaran_(&iseed[1]);

    if (*idist == 1) {

/*        uniform (0,1) */

	ret_val = t1;
    } else if (*idist == 2) {

/*        uniform (-1,1) */

	ret_val = t1 * 2.f - 1.f;
    } else if (*idist == 3) {

/*        normal (0,1) */

	t2 = slaran_(&iseed[1]);
	ret_val = sqrt(log(t1) * -2.f) * cos(t2 * 
		6.2831853071795864769252867663f);
    }
    return ret_val;

/*     End of SLARND */

} /* slarnd_ */

