#include "f2c.h"

/* Subroutine */ int claset_(char *uplo, integer *m, integer *n, complex *
	alpha, complex *beta, complex *a, integer *lda)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    CLASET initializes a 2-D array A to BETA on the diagonal and   
    ALPHA on the offdiagonals.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies the part of the matrix A to be set.   
            = 'U':      Upper triangular part is set. The lower triangle 
  
                        is unchanged.   
            = 'L':      Lower triangular part is set. The upper triangle 
  
                        is unchanged.   
            Otherwise:  All of the matrix A is set.   

    M       (input) INTEGER   
            On entry, M specifies the number of rows of A.   

    N       (input) INTEGER   
            On entry, N specifies the number of columns of A.   

    ALPHA   (input) COMPLEX   
            All the offdiagonal array elements are set to ALPHA.   

    BETA    (input) COMPLEX   
            All the diagonal array elements are set to BETA.   

    A       (input/output) COMPLEX array, dimension (LDA,N)   
            On entry, the m by n matrix A.   
            On exit, A(i,j) = ALPHA, 1 <= i <= m, 1 <= j <= n, i.ne.j;   
                     A(i,i) = BETA , 1 <= i <= min(m,n)   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    static integer i, j;
    extern logical lsame_(char *, char *);



#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    if (lsame_(uplo, "U")) {

/*        Set the diagonal to BETA and the strictly upper triangular 
  
          part of the array to ALPHA. */

	i__1 = *n;
	for (j = 2; j <= *n; ++j) {
/* Computing MIN */
	    i__3 = j - 1;
	    i__2 = min(i__3,*m);
	    for (i = 1; i <= min(j-1,*m); ++i) {
		i__3 = i + j * a_dim1;
		A(i,j).r = alpha->r, A(i,j).i = alpha->i;
/* L10: */
	    }
/* L20: */
	}
	i__1 = min(*n,*m);
	for (i = 1; i <= min(*n,*m); ++i) {
	    i__2 = i + i * a_dim1;
	    A(i,i).r = beta->r, A(i,i).i = beta->i;
/* L30: */
	}

    } else if (lsame_(uplo, "L")) {

/*        Set the diagonal to BETA and the strictly lower triangular 
  
          part of the array to ALPHA. */

	i__1 = min(*m,*n);
	for (j = 1; j <= min(*m,*n); ++j) {
	    i__2 = *m;
	    for (i = j + 1; i <= *m; ++i) {
		i__3 = i + j * a_dim1;
		A(i,j).r = alpha->r, A(i,j).i = alpha->i;
/* L40: */
	    }
/* L50: */
	}
	i__1 = min(*n,*m);
	for (i = 1; i <= min(*n,*m); ++i) {
	    i__2 = i + i * a_dim1;
	    A(i,i).r = beta->r, A(i,i).i = beta->i;
/* L60: */
	}

    } else {

/*        Set the array to BETA on the diagonal and ALPHA on the   
          offdiagonal. */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    i__2 = *m;
	    for (i = 1; i <= *m; ++i) {
		i__3 = i + j * a_dim1;
		A(i,j).r = alpha->r, A(i,j).i = alpha->i;
/* L70: */
	    }
/* L80: */
	}
	i__1 = min(*m,*n);
	for (i = 1; i <= min(*m,*n); ++i) {
	    i__2 = i + i * a_dim1;
	    A(i,i).r = beta->r, A(i,i).i = beta->i;
/* L90: */
	}
    }

    return 0;

/*     End of CLASET */

} /* claset_ */

