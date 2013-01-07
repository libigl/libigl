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

/* Subroutine */ int claghe_(integer *n, integer *k, real *d, complex *a, 
	integer *lda, integer *iseed, complex *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;
    complex q__1, q__2, q__3, q__4;

    /* Builtin functions */
    double c_abs(complex *);
    void c_div(complex *, complex *, complex *), r_cnjg(complex *, complex *);

    /* Local variables */
    extern /* Subroutine */ int cher2_(char *, integer *, complex *, complex *
	    , integer *, complex *, integer *, complex *, integer *);
    static integer i, j;
    extern /* Subroutine */ int cgerc_(integer *, integer *, complex *, 
	    complex *, integer *, complex *, integer *, complex *, integer *);
    static complex alpha;
    extern /* Subroutine */ int cscal_(integer *, complex *, complex *, 
	    integer *);
    extern /* Complex */ VOID cdotc_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);
    extern /* Subroutine */ int cgemv_(char *, integer *, integer *, complex *
	    , complex *, integer *, complex *, integer *, complex *, complex *
	    , integer *), chemv_(char *, integer *, complex *, 
	    complex *, integer *, complex *, integer *, complex *, complex *, 
	    integer *), caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *);
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

    CLAGHE generates a complex hermitian matrix A, by pre- and post-   
    multiplying a real diagonal matrix D with a random unitary matrix:   
    A = U*D*U'. The semi-bandwidth may then be reduced to k by additional 
  
    unitary transformations.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    K       (input) INTEGER   
            The number of nonzero subdiagonals within the band of A.   
            0 <= K <= N-1.   

    D       (input) REAL array, dimension (N)   
            The diagonal elements of the diagonal matrix D.   

    A       (output) COMPLEX array, dimension (LDA,N)   
            The generated n by n hermitian matrix A (the full matrix is   
            stored).   

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
    --d;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --iseed;
    --work;

    /* Function Body */
    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*k < 0 || *k > *n - 1) {
	*info = -2;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    }
    if (*info < 0) {
	i__1 = -(*info);
	xerbla_("CLAGHE", &i__1);
	return 0;
    }

/*     initialize lower triangle of A to diagonal matrix */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i = j + 1; i <= i__2; ++i) {
	    i__3 = i + j * a_dim1;
	    a[i__3].r = 0.f, a[i__3].i = 0.f;
/* L10: */
	}
/* L20: */
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	i__2 = i + i * a_dim1;
	i__3 = i;
	a[i__2].r = d[i__3], a[i__2].i = 0.f;
/* L30: */
    }

/*     Generate lower triangle of hermitian matrix */

    for (i = *n - 1; i >= 1; --i) {

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

/*        apply random reflection to A(i:n,i:n) from the left   
          and the right   

          compute  y := tau * A * u */

	i__1 = *n - i + 1;
	chemv_("Lower", &i__1, &tau, &a[i + i * a_dim1], lda, &work[1], &c__1,
		 &c_b1, &work[*n + 1], &c__1);

/*        compute  v := y - 1/2 * tau * ( y, u ) * u */

	q__3.r = -.5f, q__3.i = 0.f;
	q__2.r = q__3.r * tau.r - q__3.i * tau.i, q__2.i = q__3.r * tau.i + 
		q__3.i * tau.r;
	i__1 = *n - i + 1;
	cdotc_(&q__4, &i__1, &work[*n + 1], &c__1, &work[1], &c__1);
	q__1.r = q__2.r * q__4.r - q__2.i * q__4.i, q__1.i = q__2.r * q__4.i 
		+ q__2.i * q__4.r;
	alpha.r = q__1.r, alpha.i = q__1.i;
	i__1 = *n - i + 1;
	caxpy_(&i__1, &alpha, &work[1], &c__1, &work[*n + 1], &c__1);

/*        apply the transformation as a rank-2 update to A(i:n,i:n) */

	i__1 = *n - i + 1;
	q__1.r = -1.f, q__1.i = 0.f;
	cher2_("Lower", &i__1, &q__1, &work[1], &c__1, &work[*n + 1], &c__1, &
		a[i + i * a_dim1], lda);
/* L40: */
    }

/*     Reduce number of subdiagonals to K */

    i__1 = *n - 1 - *k;
    for (i = 1; i <= i__1; ++i) {

/*        generate reflection to annihilate A(k+i+1:n,i) */

	i__2 = *n - *k - i + 1;
	wn = scnrm2_(&i__2, &a[*k + i + i * a_dim1], &c__1);
	d__1 = wn / c_abs(&a[*k + i + i * a_dim1]);
	i__2 = *k + i + i * a_dim1;
	q__1.r = d__1 * a[i__2].r, q__1.i = d__1 * a[i__2].i;
	wa.r = q__1.r, wa.i = q__1.i;
	if (wn == 0.f) {
	    tau.r = 0.f, tau.i = 0.f;
	} else {
	    i__2 = *k + i + i * a_dim1;
	    q__1.r = a[i__2].r + wa.r, q__1.i = a[i__2].i + wa.i;
	    wb.r = q__1.r, wb.i = q__1.i;
	    i__2 = *n - *k - i;
	    c_div(&q__1, &c_b2, &wb);
	    cscal_(&i__2, &q__1, &a[*k + i + 1 + i * a_dim1], &c__1);
	    i__2 = *k + i + i * a_dim1;
	    a[i__2].r = 1.f, a[i__2].i = 0.f;
	    c_div(&q__1, &wb, &wa);
	    d__1 = q__1.r;
	    tau.r = d__1, tau.i = 0.f;
	}

/*        apply reflection to A(k+i:n,i+1:k+i-1) from the left */

	i__2 = *n - *k - i + 1;
	i__3 = *k - 1;
	cgemv_("Conjugate transpose", &i__2, &i__3, &c_b2, &a[*k + i + (i + 1)
		 * a_dim1], lda, &a[*k + i + i * a_dim1], &c__1, &c_b1, &work[
		1], &c__1);
	i__2 = *n - *k - i + 1;
	i__3 = *k - 1;
	q__1.r = -(doublereal)tau.r, q__1.i = -(doublereal)tau.i;
	cgerc_(&i__2, &i__3, &q__1, &a[*k + i + i * a_dim1], &c__1, &work[1], 
		&c__1, &a[*k + i + (i + 1) * a_dim1], lda);

/*        apply reflection to A(k+i:n,k+i:n) from the left and the rig
ht   

          compute  y := tau * A * u */

	i__2 = *n - *k - i + 1;
	chemv_("Lower", &i__2, &tau, &a[*k + i + (*k + i) * a_dim1], lda, &a[*
		k + i + i * a_dim1], &c__1, &c_b1, &work[1], &c__1);

/*        compute  v := y - 1/2 * tau * ( y, u ) * u */

	q__3.r = -.5f, q__3.i = 0.f;
	q__2.r = q__3.r * tau.r - q__3.i * tau.i, q__2.i = q__3.r * tau.i + 
		q__3.i * tau.r;
	i__2 = *n - *k - i + 1;
	cdotc_(&q__4, &i__2, &work[1], &c__1, &a[*k + i + i * a_dim1], &c__1);
	q__1.r = q__2.r * q__4.r - q__2.i * q__4.i, q__1.i = q__2.r * q__4.i 
		+ q__2.i * q__4.r;
	alpha.r = q__1.r, alpha.i = q__1.i;
	i__2 = *n - *k - i + 1;
	caxpy_(&i__2, &alpha, &a[*k + i + i * a_dim1], &c__1, &work[1], &c__1)
		;

/*        apply hermitian rank-2 update to A(k+i:n,k+i:n) */

	i__2 = *n - *k - i + 1;
	q__1.r = -1.f, q__1.i = 0.f;
	cher2_("Lower", &i__2, &q__1, &a[*k + i + i * a_dim1], &c__1, &work[1]
		, &c__1, &a[*k + i + (*k + i) * a_dim1], lda);

	i__2 = *k + i + i * a_dim1;
	q__1.r = -(doublereal)wa.r, q__1.i = -(doublereal)wa.i;
	a[i__2].r = q__1.r, a[i__2].i = q__1.i;
	i__2 = *n;
	for (j = *k + i + 1; j <= i__2; ++j) {
	    i__3 = j + i * a_dim1;
	    a[i__3].r = 0.f, a[i__3].i = 0.f;
/* L50: */
	}
/* L60: */
    }

/*     Store full hermitian matrix */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i = j + 1; i <= i__2; ++i) {
	    i__3 = j + i * a_dim1;
	    r_cnjg(&q__1, &a[i + j * a_dim1]);
	    a[i__3].r = q__1.r, a[i__3].i = q__1.i;
/* L70: */
	}
/* L80: */
    }
    return 0;

/*     End of CLAGHE */

} /* claghe_ */

