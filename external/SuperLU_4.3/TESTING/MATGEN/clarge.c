/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static complex c_b1 = {0.f,0.f};
static complex c_b2 = {1.f,0.f};
static integer c__3 = 3;
static integer c__1 = 1;

/* Subroutine */ int clarge_(integer *n, complex *a, integer *lda, integer *
	iseed, complex *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    doublereal d__1;
    complex q__1;

    /* Builtin functions */
    double c_abs(complex *);
    void c_div(complex *, complex *, complex *);

    /* Local variables */
    static integer i;
    extern /* Subroutine */ int cgerc_(integer *, integer *, complex *, 
	    complex *, integer *, complex *, integer *, complex *, integer *),
	     cscal_(integer *, complex *, complex *, integer *), cgemv_(char *
	    , integer *, integer *, complex *, complex *, integer *, complex *
	    , integer *, complex *, complex *, integer *);
    extern real scnrm2_(integer *, complex *, integer *);
    static complex wa, wb;
    static real wn;
    extern /* Subroutine */ int xerbla_(char *, integer *), clarnv_(
	    integer *, integer *, integer *, complex *);
    static complex tau;


/*  -- LAPACK auxiliary test routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    CLARGE pre- and post-multiplies a complex general n by n matrix A   
    with a random unitary matrix: A = U*D*U'.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) COMPLEX array, dimension (LDA,N)   
            On entry, the original n by n matrix A.   
            On exit, A is overwritten by U*A*U' for some random   
            unitary matrix U.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= N.   

    ISEED   (input/output) INTEGER array, dimension (4)   
            On entry, the seed of the random number generator; the array 
  
            elements must be between 0 and 4095, and ISEED(4) must be   
            odd.   
            On exit, the seed is updated.   

    WORK    (workspace) COMPLEX array, dimension (2*N)   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input arguments   

       Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --iseed;
    --work;

    /* Function Body */
    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*lda < max(1,*n)) {
	*info = -3;
    }
    if (*info < 0) {
	i__1 = -(*info);
	xerbla_("CLARGE", &i__1);
	return 0;
    }

/*     pre- and post-multiply A by random unitary matrix */

    for (i = *n; i >= 1; --i) {

/*        generate random reflection */

	i__1 = *n - i + 1;
	clarnv_(&c__3, &iseed[1], &i__1, &work[1]);
	i__1 = *n - i + 1;
	wn = scnrm2_(&i__1, &work[1], &c__1);
	d__1 = wn / c_abs(&work[1]);
	q__1.r = d__1 * work[1].r, q__1.i = d__1 * work[1].i;
	wa.r = q__1.r, wa.i = q__1.i;
	if (wn == 0.f) {
	    tau.r = 0.f, tau.i = 0.f;
	} else {
	    q__1.r = work[1].r + wa.r, q__1.i = work[1].i + wa.i;
	    wb.r = q__1.r, wb.i = q__1.i;
	    i__1 = *n - i;
	    c_div(&q__1, &c_b2, &wb);
	    cscal_(&i__1, &q__1, &work[2], &c__1);
	    work[1].r = 1.f, work[1].i = 0.f;
	    c_div(&q__1, &wb, &wa);
	    d__1 = q__1.r;
	    tau.r = d__1, tau.i = 0.f;
	}

/*        multiply A(i:n,1:n) by random reflection from the left */

	i__1 = *n - i + 1;
	cgemv_("Conjugate transpose", &i__1, n, &c_b2, &a[i + a_dim1], lda, &
		work[1], &c__1, &c_b1, &work[*n + 1], &c__1);
	i__1 = *n - i + 1;
	q__1.r = -(doublereal)tau.r, q__1.i = -(doublereal)tau.i;
	cgerc_(&i__1, n, &q__1, &work[1], &c__1, &work[*n + 1], &c__1, &a[i + 
		a_dim1], lda);

/*        multiply A(1:n,i:n) by random reflection from the right */

	i__1 = *n - i + 1;
	cgemv_("No transpose", n, &i__1, &c_b2, &a[i * a_dim1 + 1], lda, &
		work[1], &c__1, &c_b1, &work[*n + 1], &c__1);
	i__1 = *n - i + 1;
	q__1.r = -(doublereal)tau.r, q__1.i = -(doublereal)tau.i;
	cgerc_(n, &i__1, &q__1, &work[*n + 1], &c__1, &work[1], &c__1, &a[i * 
		a_dim1 + 1], lda);
/* L10: */
    }
    return 0;

/*     End of CLARGE */

} /* clarge_ */

