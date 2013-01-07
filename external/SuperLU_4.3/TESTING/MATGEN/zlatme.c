/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c__0 = 0;
static integer c__5 = 5;

/* Subroutine */ int zlatme_(integer *n, char *dist, integer *iseed, 
	doublecomplex *d, integer *mode, doublereal *cond, doublecomplex *
	dmax__, char *ei, char *rsign, char *upper, char *sim, doublereal *ds,
	 integer *modes, doublereal *conds, integer *kl, integer *ku, 
	doublereal *anorm, doublecomplex *a, integer *lda, doublecomplex *
	work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static logical bads;
    static integer isim;
    static doublereal temp;
    static integer i, j;
    static doublecomplex alpha;
    extern logical lsame_(char *, char *);
    static integer iinfo;
    static doublereal tempa[1];
    static integer icols;
    extern /* Subroutine */ int zgerc_(integer *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static integer idist;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), zgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *);
    static integer irows;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), dlatm1_(integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, integer *, integer 
	    *), zlatm1_(integer *, doublereal *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, integer *);
    static integer ic, jc, ir;
    static doublereal ralpha;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *);
    extern /* Subroutine */ int zdscal_(integer *, doublereal *, 
	    doublecomplex *, integer *), zlarge_(integer *, doublecomplex *, 
	    integer *, integer *, doublecomplex *, integer *), zlarfg_(
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *), zlacgv_(integer *, doublecomplex *, integer *);
    extern /* Double Complex */ void zlarnd_(doublecomplex *, integer *, 
	    integer *);
    static integer irsign;
    extern /* Subroutine */ int zlaset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *);
    static integer iupper;
    extern /* Subroutine */ int zlarnv_(integer *, integer *, integer *, 
	    doublecomplex *);
    static doublecomplex xnorms;
    static integer jcr;
    static doublecomplex tau;


/*  -- LAPACK test routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

       ZLATME generates random non-symmetric square matrices with   
       specified eigenvalues for testing LAPACK programs.   

       ZLATME operates by applying the following sequence of   
       operations:   

       1. Set the diagonal to D, where D may be input or   
            computed according to MODE, COND, DMAX, and RSIGN   
            as described below.   

       2. If UPPER='T', the upper triangle of A is set to random values   
            out of distribution DIST.   

       3. If SIM='T', A is multiplied on the left by a random matrix   
            X, whose singular values are specified by DS, MODES, and   
            CONDS, and on the right by X inverse.   

       4. If KL < N-1, the lower bandwidth is reduced to KL using   
            Householder transformations.  If KU < N-1, the upper   
            bandwidth is reduced to KU.   

       5. If ANORM is not negative, the matrix is scaled to have   
            maximum-element-norm ANORM.   

       (Note: since the matrix cannot be reduced beyond Hessenberg form, 
  
        no packing options are available.)   

    Arguments   
    =========   

    N      - INTEGER   
             The number of columns (or rows) of A. Not modified.   

    DIST   - CHARACTER*1   
             On entry, DIST specifies the type of distribution to be used 
  
             to generate the random eigen-/singular values, and on the   
             upper triangle (see UPPER).   
             'U' => UNIFORM( 0, 1 )  ( 'U' for uniform )   
             'S' => UNIFORM( -1, 1 ) ( 'S' for symmetric )   
             'N' => NORMAL( 0, 1 )   ( 'N' for normal )   
             'D' => uniform on the complex disc |z| < 1.   
             Not modified.   

    ISEED  - INTEGER array, dimension ( 4 )   
             On entry ISEED specifies the seed of the random number   
             generator. They should lie between 0 and 4095 inclusive,   
             and ISEED(4) should be odd. The random number generator   
             uses a linear congruential sequence limited to small   
             integers, and so should produce machine independent   
             random numbers. The values of ISEED are changed on   
             exit, and can be used in the next call to ZLATME   
             to continue the same random number sequence.   
             Changed on exit.   

    D      - COMPLEX*16 array, dimension ( N )   
             This array is used to specify the eigenvalues of A.  If   
             MODE=0, then D is assumed to contain the eigenvalues   
             otherwise they will be computed according to MODE, COND,   
             DMAX, and RSIGN and placed in D.   
             Modified if MODE is nonzero.   

    MODE   - INTEGER   
             On entry this describes how the eigenvalues are to   
             be specified:   
             MODE = 0 means use D as input   
             MODE = 1 sets D(1)=1 and D(2:N)=1.0/COND   
             MODE = 2 sets D(1:N-1)=1 and D(N)=1.0/COND   
             MODE = 3 sets D(I)=COND**(-(I-1)/(N-1))   
             MODE = 4 sets D(i)=1 - (i-1)/(N-1)*(1 - 1/COND)   
             MODE = 5 sets D to random numbers in the range   
                      ( 1/COND , 1 ) such that their logarithms   
                      are uniformly distributed.   
             MODE = 6 set D to random numbers from same distribution   
                      as the rest of the matrix.   
             MODE < 0 has the same meaning as ABS(MODE), except that   
                the order of the elements of D is reversed.   
             Thus if MODE is between 1 and 4, D has entries ranging   
                from 1 to 1/COND, if between -1 and -4, D has entries   
                ranging from 1/COND to 1,   
             Not modified.   

    COND   - DOUBLE PRECISION   
             On entry, this is used as described under MODE above.   
             If used, it must be >= 1. Not modified.   

    DMAX   - COMPLEX*16   
             If MODE is neither -6, 0 nor 6, the contents of D, as   
             computed according to MODE and COND, will be scaled by   
             DMAX / max(abs(D(i))).  Note that DMAX need not be   
             positive or real: if DMAX is negative or complex (or zero), 
  
             D will be scaled by a negative or complex number (or zero). 
  
             If RSIGN='F' then the largest (absolute) eigenvalue will be 
  
             equal to DMAX.   
             Not modified.   

    EI     - CHARACTER*1 (ignored)   
             Not modified.   

    RSIGN  - CHARACTER*1   
             If MODE is not 0, 6, or -6, and RSIGN='T', then the   
             elements of D, as computed according to MODE and COND, will 
  
             be multiplied by a random complex number from the unit   
             circle |z| = 1.  If RSIGN='F', they will not be.  RSIGN may 
  
             only have the values 'T' or 'F'.   
             Not modified.   

    UPPER  - CHARACTER*1   
             If UPPER='T', then the elements of A above the diagonal   
             will be set to random numbers out of DIST.  If UPPER='F',   
             they will not.  UPPER may only have the values 'T' or 'F'.   
             Not modified.   

    SIM    - CHARACTER*1   
             If SIM='T', then A will be operated on by a "similarity   
             transform", i.e., multiplied on the left by a matrix X and   
             on the right by X inverse.  X = U S V, where U and V are   
             random unitary matrices and S is a (diagonal) matrix of   
             singular values specified by DS, MODES, and CONDS.  If   
             SIM='F', then A will not be transformed.   
             Not modified.   

    DS     - DOUBLE PRECISION array, dimension ( N )   
             This array is used to specify the singular values of X,   
             in the same way that D specifies the eigenvalues of A.   
             If MODE=0, the DS contains the singular values, which   
             may not be zero.   
             Modified if MODE is nonzero.   

    MODES  - INTEGER   
    CONDS  - DOUBLE PRECISION   
             Similar to MODE and COND, but for specifying the diagonal   
             of S.  MODES=-6 and +6 are not allowed (since they would   
             result in randomly ill-conditioned eigenvalues.)   

    KL     - INTEGER   
             This specifies the lower bandwidth of the  matrix.  KL=1   
             specifies upper Hessenberg form.  If KL is at least N-1,   
             then A will have full lower bandwidth.   
             Not modified.   

    KU     - INTEGER   
             This specifies the upper bandwidth of the  matrix.  KU=1   
             specifies lower Hessenberg form.  If KU is at least N-1,   
             then A will have full upper bandwidth; if KU and KL   
             are both at least N-1, then A will be dense.  Only one of   
             KU and KL may be less than N-1.   
             Not modified.   

    ANORM  - DOUBLE PRECISION   
             If ANORM is not negative, then A will be scaled by a non-   
             negative real number to make the maximum-element-norm of A   
             to be ANORM.   
             Not modified.   

    A      - COMPLEX*16 array, dimension ( LDA, N )   
             On exit A is the desired test matrix.   
             Modified.   

    LDA    - INTEGER   
             LDA specifies the first dimension of A as declared in the   
             calling program.  LDA must be at least M.   
             Not modified.   

    WORK   - COMPLEX*16 array, dimension ( 3*N )   
             Workspace.   
             Modified.   

    INFO   - INTEGER   
             Error code.  On exit, INFO will be set to one of the   
             following values:   
               0 => normal return   
              -1 => N negative   
              -2 => DIST illegal string   
              -5 => MODE not in range -6 to 6   
              -6 => COND less than 1.0, and MODE neither -6, 0 nor 6   
              -9 => RSIGN is not 'T' or 'F'   
             -10 => UPPER is not 'T' or 'F'   
             -11 => SIM   is not 'T' or 'F'   
             -12 => MODES=0 and DS has a zero singular value.   
             -13 => MODES is not in the range -5 to 5.   
             -14 => MODES is nonzero and CONDS is less than 1.   
             -15 => KL is less than 1.   
             -16 => KU is less than 1, or KL and KU are both less than   
                    N-1.   
             -19 => LDA is less than M.   
              1  => Error return from ZLATM1 (computing D)   
              2  => Cannot scale to DMAX (max. eigenvalue is 0)   
              3  => Error return from DLATM1 (computing DS)   
              4  => Error return from ZLARGE   
              5  => Zero singular value from DLATM1.   

    ===================================================================== 
  


       1)      Decode and Test the input parameters.   
               Initialize flags & seed.   

       Parameter adjustments */
    --iseed;
    --d;
    --ds;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --work;

    /* Function Body */
    *info = 0;

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Decode DIST */

    if (lsame_(dist, "U")) {
	idist = 1;
    } else if (lsame_(dist, "S")) {
	idist = 2;
    } else if (lsame_(dist, "N")) {
	idist = 3;
    } else if (lsame_(dist, "D")) {
	idist = 4;
    } else {
	idist = -1;
    }

/*     Decode RSIGN */

    if (lsame_(rsign, "T")) {
	irsign = 1;
    } else if (lsame_(rsign, "F")) {
	irsign = 0;
    } else {
	irsign = -1;
    }

/*     Decode UPPER */

    if (lsame_(upper, "T")) {
	iupper = 1;
    } else if (lsame_(upper, "F")) {
	iupper = 0;
    } else {
	iupper = -1;
    }

/*     Decode SIM */

    if (lsame_(sim, "T")) {
	isim = 1;
    } else if (lsame_(sim, "F")) {
	isim = 0;
    } else {
	isim = -1;
    }

/*     Check DS, if MODES=0 and ISIM=1 */

    bads = FALSE_;
    if (*modes == 0 && isim == 1) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (ds[j] == 0.) {
		bads = TRUE_;
	    }
/* L10: */
	}
    }

/*     Set INFO if an error */

    if (*n < 0) {
	*info = -1;
    } else if (idist == -1) {
	*info = -2;
    } else if (abs(*mode) > 6) {
	*info = -5;
    } else if (*mode != 0 && abs(*mode) != 6 && *cond < 1.) {
	*info = -6;
    } else if (irsign == -1) {
	*info = -9;
    } else if (iupper == -1) {
	*info = -10;
    } else if (isim == -1) {
	*info = -11;
    } else if (bads) {
	*info = -12;
    } else if (isim == 1 && abs(*modes) > 5) {
	*info = -13;
    } else if (isim == 1 && *modes != 0 && *conds < 1.) {
	*info = -14;
    } else if (*kl < 1) {
	*info = -15;
    } else if (*ku < 1 || *ku < *n - 1 && *kl < *n - 1) {
	*info = -16;
    } else if (*lda < max(1,*n)) {
	*info = -19;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZLATME", &i__1);
	return 0;
    }

/*     Initialize random number generator */

    for (i = 1; i <= 4; ++i) {
	iseed[i] = (i__1 = iseed[i], abs(i__1)) % 4096;
/* L20: */
    }

    if (iseed[4] % 2 != 1) {
	++iseed[4];
    }

/*     2)      Set up diagonal of A   

               Compute D according to COND and MODE */

    zlatm1_(mode, cond, &irsign, &idist, &iseed[1], &d[1], n, &iinfo);
    if (iinfo != 0) {
	*info = 1;
	return 0;
    }
    if (*mode != 0 && abs(*mode) != 6) {

/*        Scale by DMAX */

	temp = z_abs(&d[1]);
	i__1 = *n;
	for (i = 2; i <= i__1; ++i) {
/* Computing MAX */
	    d__1 = temp, d__2 = z_abs(&d[i]);
	    temp = max(d__1,d__2);
/* L30: */
	}

	if (temp > 0.) {
	    z__1.r = dmax__->r / temp, z__1.i = dmax__->i / temp;
	    alpha.r = z__1.r, alpha.i = z__1.i;
	} else {
	    *info = 2;
	    return 0;
	}

	zscal_(n, &alpha, &d[1], &c__1);

    }

    zlaset_("Full", n, n, &c_b1, &c_b1, &a[a_offset], lda);
    i__1 = *lda + 1;
    zcopy_(n, &d[1], &c__1, &a[a_offset], &i__1);

/*     3)      If UPPER='T', set upper triangle of A to random numbers. */

    if (iupper != 0) {
	i__1 = *n;
	for (jc = 2; jc <= i__1; ++jc) {
	    i__2 = jc - 1;
	    zlarnv_(&idist, &iseed[1], &i__2, &a[jc * a_dim1 + 1]);
/* L40: */
	}
    }

/*     4)      If SIM='T', apply similarity transformation.   

                                  -1   
               Transform is  X A X  , where X = U S V, thus   

               it is  U S V A V' (1/S) U' */

    if (isim != 0) {

/*        Compute S (singular values of the eigenvector matrix)   
          according to CONDS and MODES */

	dlatm1_(modes, conds, &c__0, &c__0, &iseed[1], &ds[1], n, &iinfo);
	if (iinfo != 0) {
	    *info = 3;
	    return 0;
	}

/*        Multiply by V and V' */

	zlarge_(n, &a[a_offset], lda, &iseed[1], &work[1], &iinfo);
	if (iinfo != 0) {
	    *info = 4;
	    return 0;
	}

/*        Multiply by S and (1/S) */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    zdscal_(n, &ds[j], &a[j + a_dim1], lda);
	    if (ds[j] != 0.) {
		d__1 = 1. / ds[j];
		zdscal_(n, &d__1, &a[j * a_dim1 + 1], &c__1);
	    } else {
		*info = 5;
		return 0;
	    }
/* L50: */
	}

/*        Multiply by U and U' */

	zlarge_(n, &a[a_offset], lda, &iseed[1], &work[1], &iinfo);
	if (iinfo != 0) {
	    *info = 4;
	    return 0;
	}
    }

/*     5)      Reduce the bandwidth. */

    if (*kl < *n - 1) {

/*        Reduce bandwidth -- kill column */

	i__1 = *n - 1;
	for (jcr = *kl + 1; jcr <= i__1; ++jcr) {
	    ic = jcr - *kl;
	    irows = *n + 1 - jcr;
	    icols = *n + *kl - jcr;

	    zcopy_(&irows, &a[jcr + ic * a_dim1], &c__1, &work[1], &c__1);
	    xnorms.r = work[1].r, xnorms.i = work[1].i;
	    zlarfg_(&irows, &xnorms, &work[2], &c__1, &tau);
	    d_cnjg(&z__1, &tau);
	    tau.r = z__1.r, tau.i = z__1.i;
	    work[1].r = 1., work[1].i = 0.;
	    zlarnd_(&z__1, &c__5, &iseed[1]);
	    alpha.r = z__1.r, alpha.i = z__1.i;

	    zgemv_("C", &irows, &icols, &c_b2, &a[jcr + (ic + 1) * a_dim1], 
		    lda, &work[1], &c__1, &c_b1, &work[irows + 1], &c__1);
	    z__1.r = -tau.r, z__1.i = -tau.i;
	    zgerc_(&irows, &icols, &z__1, &work[1], &c__1, &work[irows + 1], &
		    c__1, &a[jcr + (ic + 1) * a_dim1], lda);

	    zgemv_("N", n, &irows, &c_b2, &a[jcr * a_dim1 + 1], lda, &work[1],
		     &c__1, &c_b1, &work[irows + 1], &c__1);
	    d_cnjg(&z__2, &tau);
	    z__1.r = -z__2.r, z__1.i = -z__2.i;
	    zgerc_(n, &irows, &z__1, &work[irows + 1], &c__1, &work[1], &c__1,
		     &a[jcr * a_dim1 + 1], lda);

	    i__2 = jcr + ic * a_dim1;
	    a[i__2].r = xnorms.r, a[i__2].i = xnorms.i;
	    i__2 = irows - 1;
	    zlaset_("Full", &i__2, &c__1, &c_b1, &c_b1, &a[jcr + 1 + ic * 
		    a_dim1], lda);

	    i__2 = icols + 1;
	    zscal_(&i__2, &alpha, &a[jcr + ic * a_dim1], lda);
	    d_cnjg(&z__1, &alpha);
	    zscal_(n, &z__1, &a[jcr * a_dim1 + 1], &c__1);
/* L60: */
	}
    } else if (*ku < *n - 1) {

/*        Reduce upper bandwidth -- kill a row at a time. */

	i__1 = *n - 1;
	for (jcr = *ku + 1; jcr <= i__1; ++jcr) {
	    ir = jcr - *ku;
	    irows = *n + *ku - jcr;
	    icols = *n + 1 - jcr;

	    zcopy_(&icols, &a[ir + jcr * a_dim1], lda, &work[1], &c__1);
	    xnorms.r = work[1].r, xnorms.i = work[1].i;
	    zlarfg_(&icols, &xnorms, &work[2], &c__1, &tau);
	    d_cnjg(&z__1, &tau);
	    tau.r = z__1.r, tau.i = z__1.i;
	    work[1].r = 1., work[1].i = 0.;
	    i__2 = icols - 1;
	    zlacgv_(&i__2, &work[2], &c__1);
	    zlarnd_(&z__1, &c__5, &iseed[1]);
	    alpha.r = z__1.r, alpha.i = z__1.i;

	    zgemv_("N", &irows, &icols, &c_b2, &a[ir + 1 + jcr * a_dim1], lda,
		     &work[1], &c__1, &c_b1, &work[icols + 1], &c__1);
	    z__1.r = -tau.r, z__1.i = -tau.i;
	    zgerc_(&irows, &icols, &z__1, &work[icols + 1], &c__1, &work[1], &
		    c__1, &a[ir + 1 + jcr * a_dim1], lda);

	    zgemv_("C", &icols, n, &c_b2, &a[jcr + a_dim1], lda, &work[1], &
		    c__1, &c_b1, &work[icols + 1], &c__1);
	    d_cnjg(&z__2, &tau);
	    z__1.r = -z__2.r, z__1.i = -z__2.i;
	    zgerc_(&icols, n, &z__1, &work[1], &c__1, &work[icols + 1], &c__1,
		     &a[jcr + a_dim1], lda);

	    i__2 = ir + jcr * a_dim1;
	    a[i__2].r = xnorms.r, a[i__2].i = xnorms.i;
	    i__2 = icols - 1;
	    zlaset_("Full", &c__1, &i__2, &c_b1, &c_b1, &a[ir + (jcr + 1) * 
		    a_dim1], lda);

	    i__2 = irows + 1;
	    zscal_(&i__2, &alpha, &a[ir + jcr * a_dim1], &c__1);
	    d_cnjg(&z__1, &alpha);
	    zscal_(n, &z__1, &a[jcr + a_dim1], lda);
/* L70: */
	}
    }

/*     Scale the matrix to have norm ANORM */

    if (*anorm >= 0.) {
	temp = zlange_("M", n, n, &a[a_offset], lda, tempa);
	if (temp > 0.) {
	    ralpha = *anorm / temp;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		zdscal_(n, &ralpha, &a[j * a_dim1 + 1], &c__1);
/* L80: */
	    }
	}
    }

    return 0;

/*     End of ZLATME */

} /* zlatme_ */

