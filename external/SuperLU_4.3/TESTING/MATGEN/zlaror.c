/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__3 = 3;
static integer c__1 = 1;

/* Subroutine */ int zlaror_(char *side, char *init, integer *m, integer *n, 
	doublecomplex *a, integer *lda, integer *iseed, doublecomplex *x, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer kbeg, jcol;
    static doublereal xabs;
    static integer irow, j;
    extern logical lsame_(char *, char *);
    static doublecomplex csign;
    extern /* Subroutine */ int zgerc_(integer *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    static integer ixfrm;
    extern /* Subroutine */ int zgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *);
    static integer itype, nxfrm;
    static doublereal xnorm;
    extern doublereal dznrm2_(integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static doublereal factor;
    extern /* Subroutine */ int zlacgv_(integer *, doublecomplex *, integer *)
	    ;
    extern /* Double Complex */ VOID zlarnd_(doublecomplex *, integer *, 
	    integer *);
    extern /* Subroutine */ int zlaset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *);
    static doublecomplex xnorms;


/*  -- LAPACK auxiliary test routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

       ZLAROR pre- or post-multiplies an M by N matrix A by a random   
       unitary matrix U, overwriting A. A may optionally be   
       initialized to the identity matrix before multiplying by U.   
       U is generated using the method of G.W. Stewart   
       ( SIAM J. Numer. Anal. 17, 1980, pp. 403-409 ).   
       (BLAS-2 version)   

    Arguments   
    =========   

    SIDE   - CHARACTER*1   
             SIDE specifies whether A is multiplied on the left or right 
  
             by U.   
         SIDE = 'L'   Multiply A on the left (premultiply) by U   
         SIDE = 'R'   Multiply A on the right (postmultiply) by U*   
         SIDE = 'C'   Multiply A on the left by U and the right by U*   
         SIDE = 'T'   Multiply A on the left by U and the right by U'   
             Not modified.   

    INIT   - CHARACTER*1   
             INIT specifies whether or not A should be initialized to   
             the identity matrix.   
                INIT = 'I'   Initialize A to (a section of) the   
                             identity matrix before applying U.   
                INIT = 'N'   No initialization.  Apply U to the   
                             input matrix A.   

             INIT = 'I' may be used to generate square (i.e., unitary)   
             or rectangular orthogonal matrices (orthogonality being   
             in the sense of ZDOTC):   

             For square matrices, M=N, and SIDE many be either 'L' or   
             'R'; the rows will be orthogonal to each other, as will the 
  
             columns.   
             For rectangular matrices where M < N, SIDE = 'R' will   
             produce a dense matrix whose rows will be orthogonal and   
             whose columns will not, while SIDE = 'L' will produce a   
             matrix whose rows will be orthogonal, and whose first M   
             columns will be orthogonal, the remaining columns being   
             zero.   
             For matrices where M > N, just use the previous   
             explaination, interchanging 'L' and 'R' and "rows" and   
             "columns".   

             Not modified.   

    M      - INTEGER   
             Number of rows of A. Not modified.   

    N      - INTEGER   
             Number of columns of A. Not modified.   

    A      - COMPLEX*16 array, dimension ( LDA, N )   
             Input and output array. Overwritten by U A ( if SIDE = 'L' ) 
  
             or by A U ( if SIDE = 'R' )   
             or by U A U* ( if SIDE = 'C')   
             or by U A U' ( if SIDE = 'T') on exit.   

    LDA    - INTEGER   
             Leading dimension of A. Must be at least MAX ( 1, M ).   
             Not modified.   

    ISEED  - INTEGER array, dimension ( 4 )   
             On entry ISEED specifies the seed of the random number   
             generator. The array elements should be between 0 and 4095; 
  
             if not they will be reduced mod 4096.  Also, ISEED(4) must   
             be odd.  The random number generator uses a linear   
             congruential sequence limited to small integers, and so   
             should produce machine independent random numbers. The   
             values of ISEED are changed on exit, and can be used in the 
  
             next call to ZLAROR to continue the same random number   
             sequence.   
             Modified.   

    X      - COMPLEX*16 array, dimension ( 3*MAX( M, N ) )   
             Workspace. Of length:   
                 2*M + N if SIDE = 'L',   
                 2*N + M if SIDE = 'R',   
                 3*N     if SIDE = 'C' or 'T'.   
             Modified.   

    INFO   - INTEGER   
             An error flag.  It is set to:   
              0  if no error.   
              1  if ZLARND returned a bad random number (installation   
                 problem)   
             -1  if SIDE is not L, R, C, or T.   
             -3  if M is negative.   
             -4  if N is negative or if SIDE is C or T and N is not equal 
  
                 to M.   
             -6  if LDA is less than M.   

    ===================================================================== 
  


       Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --iseed;
    --x;

    /* Function Body */
    if (*n == 0 || *m == 0) {
	return 0;
    }

    itype = 0;
    if (lsame_(side, "L")) {
	itype = 1;
    } else if (lsame_(side, "R")) {
	itype = 2;
    } else if (lsame_(side, "C")) {
	itype = 3;
    } else if (lsame_(side, "T")) {
	itype = 4;
    }

/*     Check for argument errors. */

    *info = 0;
    if (itype == 0) {
	*info = -1;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0 || itype == 3 && *n != *m) {
	*info = -4;
    } else if (*lda < *m) {
	*info = -6;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZLAROR", &i__1);
	return 0;
    }

    if (itype == 1) {
	nxfrm = *m;
    } else {
	nxfrm = *n;
    }

/*     Initialize A to the identity matrix if desired */

    if (lsame_(init, "I")) {
	zlaset_("Full", m, n, &c_b1, &c_b2, &a[a_offset], lda);
    }

/*     If no rotation possible, still multiply by   
       a random complex number from the circle |x| = 1   

        2)      Compute Rotation by computing Householder   
                Transformations H(2), H(3), ..., H(n).  Note that the   
                order in which they are computed is irrelevant. */

    i__1 = nxfrm;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	x[i__2].r = 0., x[i__2].i = 0.;
/* L10: */
    }

    i__1 = nxfrm;
    for (ixfrm = 2; ixfrm <= i__1; ++ixfrm) {
	kbeg = nxfrm - ixfrm + 1;

/*        Generate independent normal( 0, 1 ) random numbers */

	i__2 = nxfrm;
	for (j = kbeg; j <= i__2; ++j) {
	    i__3 = j;
	    zlarnd_(&z__1, &c__3, &iseed[1]);
	    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
/* L20: */
	}

/*        Generate a Householder transformation from the random vector
 X */

	xnorm = dznrm2_(&ixfrm, &x[kbeg], &c__1);
	xabs = z_abs(&x[kbeg]);
	if (xabs != 0.) {
	    i__2 = kbeg;
	    z__1.r = x[i__2].r / xabs, z__1.i = x[i__2].i / xabs;
	    csign.r = z__1.r, csign.i = z__1.i;
	} else {
	    csign.r = 1., csign.i = 0.;
	}
	z__1.r = xnorm * csign.r, z__1.i = xnorm * csign.i;
	xnorms.r = z__1.r, xnorms.i = z__1.i;
	i__2 = nxfrm + kbeg;
	z__1.r = -csign.r, z__1.i = -csign.i;
	x[i__2].r = z__1.r, x[i__2].i = z__1.i;
	factor = xnorm * (xnorm + xabs);
	if (abs(factor) < 1e-20) {
	    *info = 1;
	    i__2 = -(*info);
	    xerbla_("ZLAROR", &i__2);
	    return 0;
	} else {
	    factor = 1. / factor;
	}
	i__2 = kbeg;
	i__3 = kbeg;
	z__1.r = x[i__3].r + xnorms.r, z__1.i = x[i__3].i + xnorms.i;
	x[i__2].r = z__1.r, x[i__2].i = z__1.i;

/*        Apply Householder transformation to A */

	if (itype == 1 || itype == 3 || itype == 4) {

/*           Apply H(k) on the left of A */

	    zgemv_("C", &ixfrm, n, &c_b2, &a[kbeg + a_dim1], lda, &x[kbeg], &
		    c__1, &c_b1, &x[(nxfrm << 1) + 1], &c__1);
	    z__2.r = factor, z__2.i = 0.;
	    z__1.r = -z__2.r, z__1.i = -z__2.i;
	    zgerc_(&ixfrm, n, &z__1, &x[kbeg], &c__1, &x[(nxfrm << 1) + 1], &
		    c__1, &a[kbeg + a_dim1], lda);

	}

	if (itype >= 2 && itype <= 4) {

/*           Apply H(k)* (or H(k)') on the right of A */

	    if (itype == 4) {
		zlacgv_(&ixfrm, &x[kbeg], &c__1);
	    }

	    zgemv_("N", m, &ixfrm, &c_b2, &a[kbeg * a_dim1 + 1], lda, &x[kbeg]
		    , &c__1, &c_b1, &x[(nxfrm << 1) + 1], &c__1);
	    z__2.r = factor, z__2.i = 0.;
	    z__1.r = -z__2.r, z__1.i = -z__2.i;
	    zgerc_(m, &ixfrm, &z__1, &x[(nxfrm << 1) + 1], &c__1, &x[kbeg], &
		    c__1, &a[kbeg * a_dim1 + 1], lda);

	}
/* L30: */
    }

    zlarnd_(&z__1, &c__3, &iseed[1]);
    x[1].r = z__1.r, x[1].i = z__1.i;
    xabs = z_abs(&x[1]);
    if (xabs != 0.) {
	z__1.r = x[1].r / xabs, z__1.i = x[1].i / xabs;
	csign.r = z__1.r, csign.i = z__1.i;
    } else {
	csign.r = 1., csign.i = 0.;
    }
    i__1 = nxfrm << 1;
    x[i__1].r = csign.r, x[i__1].i = csign.i;

/*     Scale the matrix A by D. */

    if (itype == 1 || itype == 3 || itype == 4) {
	i__1 = *m;
	for (irow = 1; irow <= i__1; ++irow) {
	    d_cnjg(&z__1, &x[nxfrm + irow]);
	    zscal_(n, &z__1, &a[irow + a_dim1], lda);
/* L40: */
	}
    }

    if (itype == 2 || itype == 3) {
	i__1 = *n;
	for (jcol = 1; jcol <= i__1; ++jcol) {
	    zscal_(m, &x[nxfrm + jcol], &a[jcol * a_dim1 + 1], &c__1);
/* L50: */
	}
    }

    if (itype == 4) {
	i__1 = *n;
	for (jcol = 1; jcol <= i__1; ++jcol) {
	    d_cnjg(&z__1, &x[nxfrm + jcol]);
	    zscal_(m, &z__1, &a[jcol * a_dim1 + 1], &c__1);
/* L60: */
	}
    }
    return 0;

/*     End of ZLAROR */

} /* zlaror_ */

