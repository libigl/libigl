/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static real c_b9 = 0.f;
static real c_b10 = 1.f;
static integer c__3 = 3;
static integer c__1 = 1;

/* Subroutine */ int slaror_(char *side, char *init, integer *m, integer *n, 
	real *a, integer *lda, integer *iseed, real *x, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    real r__1;

    /* Builtin functions */
    double r_sign(real *, real *);

    /* Local variables */
    static integer kbeg, jcol;
    extern /* Subroutine */ int sger_(integer *, integer *, real *, real *, 
	    integer *, real *, integer *, real *, integer *);
    static integer irow;
    extern real snrm2_(integer *, real *, integer *);
    static integer j;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *), 
	    sgemv_(char *, integer *, integer *, real *, real *, integer *, 
	    real *, integer *, real *, real *, integer *);
    static integer ixfrm, itype, nxfrm;
    static real xnorm;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static real factor;
    extern doublereal slarnd_(integer *, integer *);
    extern /* Subroutine */ int slaset_(char *, integer *, integer *, real *, 
	    real *, real *, integer *);
    static real xnorms;


/*  -- LAPACK auxiliary test routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    SLAROR pre- or post-multiplies an M by N matrix A by a random   
    orthogonal matrix U, overwriting A.  A may optionally be initialized 
  
    to the identity matrix before multiplying by U.  U is generated using 
  
    the method of G.W. Stewart (SIAM J. Numer. Anal. 17, 1980, 403-409). 
  

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            Specifies whether A is multiplied on the left or right by U. 
  
            = 'L':         Multiply A on the left (premultiply) by U   
            = 'R':         Multiply A on the right (postmultiply) by U'   
            = 'C' or 'T':  Multiply A on the left by U and the right   
                            by U' (Here, U' means U-transpose.)   

    INIT    (input) CHARACTER*1   
            Specifies whether or not A should be initialized to the   
            identity matrix.   
            = 'I':  Initialize A to (a section of) the identity matrix   
                     before applying U.   
            = 'N':  No initialization.  Apply U to the input matrix A.   

            INIT = 'I' may be used to generate square or rectangular   
            orthogonal matrices:   

            For M = N and SIDE = 'L' or 'R', the rows will be orthogonal 
  
            to each other, as will the columns.   

            If M < N, SIDE = 'R' produces a dense matrix whose rows are   
            orthogonal and whose columns are not, while SIDE = 'L'   
            produces a matrix whose rows are orthogonal, and whose first 
  
            M columns are orthogonal, and whose remaining columns are   
            zero.   

            If M > N, SIDE = 'L' produces a dense matrix whose columns   
            are orthogonal and whose rows are not, while SIDE = 'R'   
            produces a matrix whose columns are orthogonal, and whose   
            first M rows are orthogonal, and whose remaining rows are   
            zero.   

    M       (input) INTEGER   
            The number of rows of A.   

    N       (input) INTEGER   
            The number of columns of A.   

    A       (input/output) REAL array, dimension (LDA, N)   
            On entry, the array A.   
            On exit, overwritten by U A ( if SIDE = 'L' ),   
             or by A U ( if SIDE = 'R' ),   
             or by U A U' ( if SIDE = 'C' or 'T').   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    ISEED   (input/output) INTEGER array, dimension (4)   
            On entry ISEED specifies the seed of the random number   
            generator. The array elements should be between 0 and 4095;   
            if not they will be reduced mod 4096.  Also, ISEED(4) must   
            be odd.  The random number generator uses a linear   
            congruential sequence limited to small integers, and so   
            should produce machine independent random numbers. The   
            values of ISEED are changed on exit, and can be used in the   
            next call to SLAROR to continue the same random number   
            sequence.   

    X       (workspace) REAL array, dimension (3*MAX( M, N ))   
            Workspace of length   
                2*M + N if SIDE = 'L',   
                2*N + M if SIDE = 'R',   
                3*N     if SIDE = 'C' or 'T'.   

    INFO    (output) INTEGER   
            An error flag.  It is set to:   
            = 0:  normal return   
            < 0:  if INFO = -k, the k-th argument had an illegal value   
            = 1:  if the random numbers generated by SLARND are bad.   

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
    } else if (lsame_(side, "C") || lsame_(side, "T")) {
	itype = 3;
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
	xerbla_("SLAROR", &i__1);
	return 0;
    }

    if (itype == 1) {
	nxfrm = *m;
    } else {
	nxfrm = *n;
    }

/*     Initialize A to the identity matrix if desired */

    if (lsame_(init, "I")) {
	slaset_("Full", m, n, &c_b9, &c_b10, &a[a_offset], lda);
    }

/*     If no rotation possible, multiply by random +/-1   

       Compute rotation by computing Householder transformations   
       H(2), H(3), ..., H(nhouse) */

    i__1 = nxfrm;
    for (j = 1; j <= i__1; ++j) {
	x[j] = 0.f;
/* L10: */
    }

    i__1 = nxfrm;
    for (ixfrm = 2; ixfrm <= i__1; ++ixfrm) {
	kbeg = nxfrm - ixfrm + 1;

/*        Generate independent normal( 0, 1 ) random numbers */

	i__2 = nxfrm;
	for (j = kbeg; j <= i__2; ++j) {
	    x[j] = slarnd_(&c__3, &iseed[1]);
/* L20: */
	}

/*        Generate a Householder transformation from the random vector
 X */

	xnorm = snrm2_(&ixfrm, &x[kbeg], &c__1);
	xnorms = r_sign(&xnorm, &x[kbeg]);
	r__1 = -(doublereal)x[kbeg];
	x[kbeg + nxfrm] = r_sign(&c_b10, &r__1);
	factor = xnorms * (xnorms + x[kbeg]);
	if (dabs(factor) < 1e-20f) {
	    *info = 1;
	    xerbla_("SLAROR", info);
	    return 0;
	} else {
	    factor = 1.f / factor;
	}
	x[kbeg] += xnorms;

/*        Apply Householder transformation to A */

	if (itype == 1 || itype == 3) {

/*           Apply H(k) from the left. */

	    sgemv_("T", &ixfrm, n, &c_b10, &a[kbeg + a_dim1], lda, &x[kbeg], &
		    c__1, &c_b9, &x[(nxfrm << 1) + 1], &c__1);
	    r__1 = -(doublereal)factor;
	    sger_(&ixfrm, n, &r__1, &x[kbeg], &c__1, &x[(nxfrm << 1) + 1], &
		    c__1, &a[kbeg + a_dim1], lda);

	}

	if (itype == 2 || itype == 3) {

/*           Apply H(k) from the right. */

	    sgemv_("N", m, &ixfrm, &c_b10, &a[kbeg * a_dim1 + 1], lda, &x[
		    kbeg], &c__1, &c_b9, &x[(nxfrm << 1) + 1], &c__1);
	    r__1 = -(doublereal)factor;
	    sger_(m, &ixfrm, &r__1, &x[(nxfrm << 1) + 1], &c__1, &x[kbeg], &
		    c__1, &a[kbeg * a_dim1 + 1], lda);

	}
/* L30: */
    }

    r__1 = slarnd_(&c__3, &iseed[1]);
    x[nxfrm * 2] = r_sign(&c_b10, &r__1);

/*     Scale the matrix A by D. */

    if (itype == 1 || itype == 3) {
	i__1 = *m;
	for (irow = 1; irow <= i__1; ++irow) {
	    sscal_(n, &x[nxfrm + irow], &a[irow + a_dim1], lda);
/* L40: */
	}
    }

    if (itype == 2 || itype == 3) {
	i__1 = *n;
	for (jcol = 1; jcol <= i__1; ++jcol) {
	    sscal_(m, &x[nxfrm + jcol], &a[jcol * a_dim1 + 1], &c__1);
/* L50: */
	}
    }
    return 0;

/*     End of SLAROR */

} /* slaror_ */

