#include <string.h>
#include "f2c.h"

logical lsamen_(integer *n, char *ca, char *cb)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    LSAMEN  tests if the first N letters of CA are the same as the   
    first N letters of CB, regardless of case.   
    LSAMEN returns .TRUE. if CA and CB are equivalent except for case   
    and .FALSE. otherwise.  LSAMEN also returns .FALSE. if LEN( CA )   
    or LEN( CB ) is less than N.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The number of characters in CA and CB to be compared.   

    CA      (input) CHARACTER*(*)   
    CB      (input) CHARACTER*(*)   
            CA and CB specify two character strings of length at least N. 
  
            Only the first N characters of each string will be accessed. 
  

   ===================================================================== 
*/
    /* System generated locals */
    integer i__1;
    logical ret_val;
    /* Local variables */
    static integer i;
    extern logical lsame_(char *, char *);


    ret_val = FALSE_;
    if (strlen(ca) < *n || strlen(cb) < *n) {
	goto L20;
    }

/*     Do for each character in the two strings. */

    i__1 = *n;
    for (i = 1; i <= *n; ++i) {

/*        Test if the characters are equal using LSAME. */

	if (! lsame_(ca + (i - 1), cb + (i - 1))) {
	    goto L20;
	}

/* L10: */
    }
    ret_val = TRUE_;

L20:
    return ret_val;

/*     End of LSAMEN */

} /* lsamen_ */

