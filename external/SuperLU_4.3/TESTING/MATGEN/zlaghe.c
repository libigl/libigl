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

/* Subroutine */ int zlaghe_(integer *n, integer *k, doublereal *d, 
	doublecomplex *a, integer *lda, integer *iseed, doublecomplex *work, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    extern /* Subroutine */ int zher2_(char *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static integer i, j;
    static doublecomplex alpha;
    extern /* Subroutine */ int zgerc_(integer *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int zgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *), 
	    zhemv_(char *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *), zaxpy_(integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *);
    extern doublereal dznrm2_(integer *, doublecomplex *, integer *);
    static doublecomplex wa, wb;
    static doublereal wn;
    extern /* Subroutine */ int xerbla_(char *, integer *), zlarnv_(
	    integer *, integer *, integer *, doublecomplex *);
    static doublecomplex tau;


/*  -- LAPACK auxiliary test routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZLAGHE generates a complex hermitian matrix A, by pre- and post-   
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

    D       (input) DOUBLE PRECISION array, dimension (N)   
            The diagonal elements of the diagonal matrix D.   

    A       (output) COMPLEX*16 array, dimension (LDA,N)   
            The generated n by n hermitian matrix A (the full matrix is   
            stored).   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= N.   

    ISEED   (input/output) INTEGER array, dimension (4)   
            On entry, the seed of the random number generator; the array 
  
            elements must be between 0 and 4095, and ISEED(4) must be   
            odd.   
            On exit, the seed is updated.   

    WORK    (workspace) COMPLEX*16 array, dimension (2*N)   

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
	xerbla_("ZLAGHE", &i__1);
	return 0;
    }

/*     initialize lower triangle of A to diagonal matrix */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i = j + 1; i <= i__2; ++i) {
	    i__3 = i + j * a_dim1;
	    a[i__3].r = 0., a[i__3].i = 0.;
/* L10: */
	}
/* L20: */
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	i__2 = i + i * a_dim1;
	i__3 = i;
	a[i__2].r = d[i__3], a[i__2].i = 0.;
/* L30: */
    }

/*     Generate lower triangle of hermitian matrix */

    for (i = *n - 1; i >= 1; --i) {

/*        generate random reflection */

	i__1 = *n - i + 1;
	zlarnv_(&c__3, &iseed[1], &i__1, &work[1]);
	i__1 = *n - i + 1;
	wn = dznrm2_(&i__1, &work[1], &c__1);
	d__1 = wn / z_abs(&work[1]);
	z__1.r = d__1 * work[1].r, z__1.i = d__1 * work[1].i;
	wa.r = z__1.r, wa.i = z__1.i;
	if (wn == 0.) {
	    tau.r = 0., tau.i = 0.;
	} else {
	    z__1.r = work[1].r + wa.r, z__1.i = work[1].i + wa.i;
	    wb.r = z__1.r, wb.i = z__1.i;
	    i__1 = *n - i;
	    z_div(&z__1, &c_b2, &wb);
	    zscal_(&i__1, &z__1, &work[2], &c__1);
	    work[1].r = 1., work[1].i = 0.;
	    z_div(&z__1, &wb, &wa);
	    d__1 = z__1.r;
	    tau.r = d__1, tau.i = 0.;
	}

/*        apply random reflection to A(i:n,i:n) from the left   
          and the right   

          compute  y := tau * A * u */

	i__1 = *n - i + 1;
	zhemv_("Lower", &i__1, &tau, &a[i + i * a_dim1], lda, &work[1], &c__1,
		 &c_b1, &work[*n + 1], &c__1);

/*        compute  v := y - 1/2 * tau * ( y, u ) * u */

	z__3.r = -.5, z__3.i = 0.;
	z__2.r = z__3.r * tau.r - z__3.i * tau.i, z__2.i = z__3.r * tau.i + 
		z__3.i * tau.r;
	i__1 = *n - i + 1;
	zdotc_(&z__4, &i__1, &work[*n + 1], &c__1, &work[1], &c__1);
	z__1.r = z__2.r * z__4.r - z__2.i * z__4.i, z__1.i = z__2.r * z__4.i 
		+ z__2.i * z__4.r;
	alpha.r = z__1.r, alpha.i = z__1.i;
	i__1 = *n - i + 1;
	zaxpy_(&i__1, &alpha, &work[1], &c__1, &work[*n + 1], &c__1);

/*        apply the transformation as a rank-2 update to A(i:n,i:n) */

	i__1 = *n - i + 1;
	z__1.r = -1., z__1.i = 0.;
	zher2_("Lower", &i__1, &z__1, &work[1], &c__1, &work[*n + 1], &c__1, &
		a[i + i * a_dim1], lda);
/* L40: */
    }

/*     Reduce number of subdiagonals to K */

    i__1 = *n - 1 - *k;
    for (i = 1; i <= i__1; ++i) {

/*        generate reflection to annihilate A(k+i+1:n,i) */

	i__2 = *n - *k - i + 1;
	wn = dznrm2_(&i__2, &a[*k + i + i * a_dim1], &c__1);
	d__1 = wn / z_abs(&a[*k + i + i * a_dim1]);
	i__2 = *k + i + i * a_dim1;
	z__1.r = d__1 * a[i__2].r, z__1.i = d__1 * a[i__2].i;
	wa.r = z__1.r, wa.i = z__1.i;
	if (wn == 0.) {
	    tau.r = 0., tau.i = 0.;
	} else {
	    i__2 = *k + i + i * a_dim1;
	    z__1.r = a[i__2].r + wa.r, z__1.i = a[i__2].i + wa.i;
	    wb.r = z__1.r, wb.i = z__1.i;
	    i__2 = *n - *k - i;
	    z_div(&z__1, &c_b2, &wb);
	    zscal_(&i__2, &z__1, &a[*k + i + 1 + i * a_dim1], &c__1);
	    i__2 = *k + i + i * a_dim1;
	    a[i__2].r = 1., a[i__2].i = 0.;
	    z_div(&z__1, &wb, &wa);
	    d__1 = z__1.r;
	    tau.r = d__1, tau.i = 0.;
	}

/*        apply reflection to A(k+i:n,i+1:k+i-1) from the left */

	i__2 = *n - *k - i + 1;
	i__3 = *k - 1;
	zgemv_("Conjugate transpose", &i__2, &i__3, &c_b2, &a[*k + i + (i + 1)
		 * a_dim1], lda, &a[*k + i + i * a_dim1], &c__1, &c_b1, &work[
		1], &c__1);
	i__2 = *n - *k - i + 1;
	i__3 = *k - 1;
	z__1.r = -tau.r, z__1.i = -tau.i;
	zgerc_(&i__2, &i__3, &z__1, &a[*k + i + i * a_dim1], &c__1, &work[1], 
		&c__1, &a[*k + i + (i + 1) * a_dim1], lda);

/*        apply reflection to A(k+i:n,k+i:n) from the left and the rig
ht   

          compute  y := tau * A * u */

	i__2 = *n - *k - i + 1;
	zhemv_("Lower", &i__2, &tau, &a[*k + i + (*k + i) * a_dim1], lda, &a[*
		k + i + i * a_dim1], &c__1, &c_b1, &work[1], &c__1);

/*        compute  v := y - 1/2 * tau * ( y, u ) * u */

	z__3.r = -.5, z__3.i = 0.;
	z__2.r = z__3.r * tau.r - z__3.i * tau.i, z__2.i = z__3.r * tau.i + 
		z__3.i * tau.r;
	i__2 = *n - *k - i + 1;
	zdotc_(&z__4, &i__2, &work[1], &c__1, &a[*k + i + i * a_dim1], &c__1);
	z__1.r = z__2.r * z__4.r - z__2.i * z__4.i, z__1.i = z__2.r * z__4.i 
		+ z__2.i * z__4.r;
	alpha.r = z__1.r, alpha.i = z__1.i;
	i__2 = *n - *k - i + 1;
	zaxpy_(&i__2, &alpha, &a[*k + i + i * a_dim1], &c__1, &work[1], &c__1)
		;

/*        apply hermitian rank-2 update to A(k+i:n,k+i:n) */

	i__2 = *n - *k - i + 1;
	z__1.r = -1., z__1.i = 0.;
	zher2_("Lower", &i__2, &z__1, &a[*k + i + i * a_dim1], &c__1, &work[1]
		, &c__1, &a[*k + i + (*k + i) * a_dim1], lda);

	i__2 = *k + i + i * a_dim1;
	z__1.r = -wa.r, z__1.i = -wa.i;
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
	i__2 = *n;
	for (j = *k + i + 1; j <= i__2; ++j) {
	    i__3 = j + i * a_dim1;
	    a[i__3].r = 0., a[i__3].i = 0.;
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
	    d_cnjg(&z__1, &a[i + j * a_dim1]);
	    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
/* L70: */
	}
/* L80: */
    }
    return 0;

/*     End of ZLAGHE */

} /* zlaghe_ */

