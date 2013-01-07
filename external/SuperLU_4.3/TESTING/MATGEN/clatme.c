/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static complex c_b1 = {0.f,0.f};
static complex c_b2 = {1.f,0.f};
static integer c__1 = 1;
static integer c__0 = 0;
static integer c__5 = 5;

/* Subroutine */ int clatme_(integer *n, char *dist, integer *iseed, complex *
	d, integer *mode, real *cond, complex *dmax__, char *ei, char *rsign, 
	char *upper, char *sim, real *ds, integer *modes, real *conds, 
	integer *kl, integer *ku, real *anorm, complex *a, integer *lda, 
	complex *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    real r__1, r__2;
    complex q__1, q__2;

    /* Builtin functions */
    double c_abs(complex *);
    void r_cnjg(complex *, complex *);

    /* Local variables */
    static logical bads;
    static integer isim;
    static real temp;
    static integer i, j;
    extern /* Subroutine */ int cgerc_(integer *, integer *, complex *, 
	    complex *, integer *, complex *, integer *, complex *, integer *);
    static complex alpha;
    extern /* Subroutine */ int cscal_(integer *, complex *, complex *, 
	    integer *);
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int cgemv_(char *, integer *, integer *, complex *
	    , complex *, integer *, complex *, integer *, complex *, complex *
	    , integer *);
    static integer iinfo;
    static real tempa[1];
    static integer icols, idist;
    extern /* Subroutine */ int ccopy_(integer *, complex *, integer *, 
	    complex *, integer *);
    static integer irows;
    extern /* Subroutine */ int clatm1_(integer *, real *, integer *, integer 
	    *, integer *, complex *, integer *, integer *), slatm1_(integer *,
	     real *, integer *, integer *, integer *, real *, integer *, 
	    integer *);
    static integer ic, jc;
    extern doublereal clange_(char *, integer *, integer *, complex *, 
	    integer *, real *);
    static integer ir;
    extern /* Subroutine */ int clarge_(integer *, complex *, integer *, 
	    integer *, complex *, integer *), clarfg_(integer *, complex *, 
	    complex *, integer *, complex *), clacgv_(integer *, complex *, 
	    integer *);
    extern /* Complex */ VOID clarnd_(complex *, integer *, integer *);
    static real ralpha;
    extern /* Subroutine */ int csscal_(integer *, real *, complex *, integer 
	    *), claset_(char *, integer *, integer *, complex *, complex *, 
	    complex *, integer *), xerbla_(char *, integer *),
	     clarnv_(integer *, integer *, integer *, complex *);
    static integer irsign, iupper;
    static complex xnorms;
    static integer jcr;
    static complex tau;


/*  -- LAPACK test routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

       CLATME generates random non-symmetric square matrices with   
       specified eigenvalues for testing LAPACK programs.   

       CLATME operates by applying the following sequence of   
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
             exit, and can be used in the next call to CLATME   
             to continue the same random number sequence.   
             Changed on exit.   

    D      - COMPLEX array, dimension ( N )   
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

    COND   - REAL   
             On entry, this is used as described under MODE above.   
             If used, it must be >= 1. Not modified.   

    DMAX   - COMPLEX   
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

    DS     - REAL array, dimension ( N )   
             This array is used to specify the singular values of X,   
             in the same way that D specifies the eigenvalues of A.   
             If MODE=0, the DS contains the singular values, which   
             may not be zero.   
             Modified if MODE is nonzero.   

    MODES  - INTEGER   
    CONDS  - REAL   
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

    ANORM  - REAL   
             If ANORM is not negative, then A will be scaled by a non-   
             negative real number to make the maximum-element-norm of A   
             to be ANORM.   
             Not modified.   

    A      - COMPLEX array, dimension ( LDA, N )   
             On exit A is the desired test matrix.   
             Modified.   

    LDA    - INTEGER   
             LDA specifies the first dimension of A as declared in the   
             calling program.  LDA must be at least M.   
             Not modified.   

    WORK   - COMPLEX array, dimension ( 3*N )   
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
              1  => Error return from CLATM1 (computing D)   
              2  => Cannot scale to DMAX (max. eigenvalue is 0)   
              3  => Error return from SLATM1 (computing DS)   
              4  => Error return from CLARGE   
              5  => Zero singular value from SLATM1.   

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
	    if (ds[j] == 0.f) {
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
    } else if (*mode != 0 && abs(*mode) != 6 && *cond < 1.f) {
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
    } else if (isim == 1 && *modes != 0 && *conds < 1.f) {
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
	xerbla_("CLATME", &i__1);
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

    clatm1_(mode, cond, &irsign, &idist, &iseed[1], &d[1], n, &iinfo);
    if (iinfo != 0) {
	*info = 1;
	return 0;
    }
    if (*mode != 0 && abs(*mode) != 6) {

/*        Scale by DMAX */

	temp = c_abs(&d[1]);
	i__1 = *n;
	for (i = 2; i <= i__1; ++i) {
/* Computing MAX */
	    r__1 = temp, r__2 = c_abs(&d[i]);
	    temp = dmax(r__1,r__2);
/* L30: */
	}

	if (temp > 0.f) {
	    q__1.r = dmax__->r / temp, q__1.i = dmax__->i / temp;
	    alpha.r = q__1.r, alpha.i = q__1.i;
	} else {
	    *info = 2;
	    return 0;
	}

	cscal_(n, &alpha, &d[1], &c__1);

    }

    claset_("Full", n, n, &c_b1, &c_b1, &a[a_offset], lda);
    i__1 = *lda + 1;
    ccopy_(n, &d[1], &c__1, &a[a_offset], &i__1);

/*     3)      If UPPER='T', set upper triangle of A to random numbers. */

    if (iupper != 0) {
	i__1 = *n;
	for (jc = 2; jc <= i__1; ++jc) {
	    i__2 = jc - 1;
	    clarnv_(&idist, &iseed[1], &i__2, &a[jc * a_dim1 + 1]);
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

	slatm1_(modes, conds, &c__0, &c__0, &iseed[1], &ds[1], n, &iinfo);
	if (iinfo != 0) {
	    *info = 3;
	    return 0;
	}

/*        Multiply by V and V' */

	clarge_(n, &a[a_offset], lda, &iseed[1], &work[1], &iinfo);
	if (iinfo != 0) {
	    *info = 4;
	    return 0;
	}

/*        Multiply by S and (1/S) */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    csscal_(n, &ds[j], &a[j + a_dim1], lda);
	    if (ds[j] != 0.f) {
		r__1 = 1.f / ds[j];
		csscal_(n, &r__1, &a[j * a_dim1 + 1], &c__1);
	    } else {
		*info = 5;
		return 0;
	    }
/* L50: */
	}

/*        Multiply by U and U' */

	clarge_(n, &a[a_offset], lda, &iseed[1], &work[1], &iinfo);
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

	    ccopy_(&irows, &a[jcr + ic * a_dim1], &c__1, &work[1], &c__1);
	    xnorms.r = work[1].r, xnorms.i = work[1].i;
	    clarfg_(&irows, &xnorms, &work[2], &c__1, &tau);
	    r_cnjg(&q__1, &tau);
	    tau.r = q__1.r, tau.i = q__1.i;
	    work[1].r = 1.f, work[1].i = 0.f;
	    clarnd_(&q__1, &c__5, &iseed[1]);
	    alpha.r = q__1.r, alpha.i = q__1.i;

	    cgemv_("C", &irows, &icols, &c_b2, &a[jcr + (ic + 1) * a_dim1], 
		    lda, &work[1], &c__1, &c_b1, &work[irows + 1], &c__1);
	    q__1.r = -(doublereal)tau.r, q__1.i = -(doublereal)tau.i;
	    cgerc_(&irows, &icols, &q__1, &work[1], &c__1, &work[irows + 1], &
		    c__1, &a[jcr + (ic + 1) * a_dim1], lda);

	    cgemv_("N", n, &irows, &c_b2, &a[jcr * a_dim1 + 1], lda, &work[1],
		     &c__1, &c_b1, &work[irows + 1], &c__1);
	    r_cnjg(&q__2, &tau);
	    q__1.r = -(doublereal)q__2.r, q__1.i = -(doublereal)q__2.i;
	    cgerc_(n, &irows, &q__1, &work[irows + 1], &c__1, &work[1], &c__1,
		     &a[jcr * a_dim1 + 1], lda);

	    i__2 = jcr + ic * a_dim1;
	    a[i__2].r = xnorms.r, a[i__2].i = xnorms.i;
	    i__2 = irows - 1;
	    claset_("Full", &i__2, &c__1, &c_b1, &c_b1, &a[jcr + 1 + ic * 
		    a_dim1], lda);

	    i__2 = icols + 1;
	    cscal_(&i__2, &alpha, &a[jcr + ic * a_dim1], lda);
	    r_cnjg(&q__1, &alpha);
	    cscal_(n, &q__1, &a[jcr * a_dim1 + 1], &c__1);
/* L60: */
	}
    } else if (*ku < *n - 1) {

/*        Reduce upper bandwidth -- kill a row at a time. */

	i__1 = *n - 1;
	for (jcr = *ku + 1; jcr <= i__1; ++jcr) {
	    ir = jcr - *ku;
	    irows = *n + *ku - jcr;
	    icols = *n + 1 - jcr;

	    ccopy_(&icols, &a[ir + jcr * a_dim1], lda, &work[1], &c__1);
	    xnorms.r = work[1].r, xnorms.i = work[1].i;
	    clarfg_(&icols, &xnorms, &work[2], &c__1, &tau);
	    r_cnjg(&q__1, &tau);
	    tau.r = q__1.r, tau.i = q__1.i;
	    work[1].r = 1.f, work[1].i = 0.f;
	    i__2 = icols - 1;
	    clacgv_(&i__2, &work[2], &c__1);
	    clarnd_(&q__1, &c__5, &iseed[1]);
	    alpha.r = q__1.r, alpha.i = q__1.i;

	    cgemv_("N", &irows, &icols, &c_b2, &a[ir + 1 + jcr * a_dim1], lda,
		     &work[1], &c__1, &c_b1, &work[icols + 1], &c__1);
	    q__1.r = -(doublereal)tau.r, q__1.i = -(doublereal)tau.i;
	    cgerc_(&irows, &icols, &q__1, &work[icols + 1], &c__1, &work[1], &
		    c__1, &a[ir + 1 + jcr * a_dim1], lda);

	    cgemv_("C", &icols, n, &c_b2, &a[jcr + a_dim1], lda, &work[1], &
		    c__1, &c_b1, &work[icols + 1], &c__1);
	    r_cnjg(&q__2, &tau);
	    q__1.r = -(doublereal)q__2.r, q__1.i = -(doublereal)q__2.i;
	    cgerc_(&icols, n, &q__1, &work[1], &c__1, &work[icols + 1], &c__1,
		     &a[jcr + a_dim1], lda);

	    i__2 = ir + jcr * a_dim1;
	    a[i__2].r = xnorms.r, a[i__2].i = xnorms.i;
	    i__2 = icols - 1;
	    claset_("Full", &c__1, &i__2, &c_b1, &c_b1, &a[ir + (jcr + 1) * 
		    a_dim1], lda);

	    i__2 = irows + 1;
	    cscal_(&i__2, &alpha, &a[ir + jcr * a_dim1], &c__1);
	    r_cnjg(&q__1, &alpha);
	    cscal_(n, &q__1, &a[jcr + a_dim1], lda);
/* L70: */
	}
    }

/*     Scale the matrix to have norm ANORM */

    if (*anorm >= 0.f) {
	temp = clange_("M", n, n, &a[a_offset], lda, tempa);
	if (temp > 0.f) {
	    ralpha = *anorm / temp;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		csscal_(n, &ralpha, &a[j * a_dim1 + 1], &c__1);
/* L80: */
	    }
	}
    }

    return 0;

/*     End of CLATME */

} /* clatme_ */

