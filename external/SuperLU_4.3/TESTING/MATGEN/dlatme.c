/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b23 = 0.;
static integer c__0 = 0;
static doublereal c_b39 = 1.;

/* Subroutine */ int dlatme_(integer *n, char *dist, integer *iseed, 
	doublereal *d, integer *mode, doublereal *cond, doublereal *dmax__, 
	char *ei, char *rsign, char *upper, char *sim, doublereal *ds, 
	integer *modes, doublereal *conds, integer *kl, integer *ku, 
	doublereal *anorm, doublereal *a, integer *lda, doublereal *work, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static logical bads;
    extern /* Subroutine */ int dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer isim;
    static doublereal temp;
    static logical badei;
    static integer i, j;
    static doublereal alpha;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *);
    static integer iinfo;
    static doublereal tempa[1];
    static integer icols;
    static logical useei;
    static integer idist;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer irows;
    extern /* Subroutine */ int dlatm1_(integer *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *);
    static integer ic, jc;
    extern doublereal dlange_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *);
    static integer ir, jr;
    extern /* Subroutine */ int dlarge_(integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *), dlarfg_(integer *, 
	    doublereal *, doublereal *, integer *, doublereal *);
    extern doublereal dlaran_(integer *);
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *), 
	    xerbla_(char *, integer *), dlarnv_(integer *, integer *, 
	    integer *, doublereal *);
    static integer irsign, iupper;
    static doublereal xnorms;
    static integer jcr;
    static doublereal tau;


/*  -- LAPACK test routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

       DLATME generates random non-symmetric square matrices with   
       specified eigenvalues for testing LAPACK programs.   

       DLATME operates by applying the following sequence of   
       operations:   

       1. Set the diagonal to D, where D may be input or   
            computed according to MODE, COND, DMAX, and RSIGN   
            as described below.   

       2. If complex conjugate pairs are desired (MODE=0 and EI(1)='R',   
            or MODE=5), certain pairs of adjacent elements of D are   
            interpreted as the real and complex parts of a complex   
            conjugate pair; A thus becomes block diagonal, with 1x1   
            and 2x2 blocks.   

       3. If UPPER='T', the upper triangle of A is set to random values   
            out of distribution DIST.   

       4. If SIM='T', A is multiplied on the left by a random matrix   
            X, whose singular values are specified by DS, MODES, and   
            CONDS, and on the right by X inverse.   

       5. If KL < N-1, the lower bandwidth is reduced to KL using   
            Householder transformations.  If KU < N-1, the upper   
            bandwidth is reduced to KU.   

       6. If ANORM is not negative, the matrix is scaled to have   
            maximum-element-norm ANORM.   

       (Note: since the matrix cannot be reduced beyond Hessenberg form, 
  
        no packing options are available.)   

    Arguments   
    =========   

    N      - INTEGER   
             The number of columns (or rows) of A. Not modified.   

    DIST   - CHARACTER*1   
             On entry, DIST specifies the type of distribution to be used 
  
             to generate the random eigen-/singular values, and for the   
             upper triangle (see UPPER).   
             'U' => UNIFORM( 0, 1 )  ( 'U' for uniform )   
             'S' => UNIFORM( -1, 1 ) ( 'S' for symmetric )   
             'N' => NORMAL( 0, 1 )   ( 'N' for normal )   
             Not modified.   

    ISEED  - INTEGER array, dimension ( 4 )   
             On entry ISEED specifies the seed of the random number   
             generator. They should lie between 0 and 4095 inclusive,   
             and ISEED(4) should be odd. The random number generator   
             uses a linear congruential sequence limited to small   
             integers, and so should produce machine independent   
             random numbers. The values of ISEED are changed on   
             exit, and can be used in the next call to DLATME   
             to continue the same random number sequence.   
             Changed on exit.   

    D      - DOUBLE PRECISION array, dimension ( N )   
             This array is used to specify the eigenvalues of A.  If   
             MODE=0, then D is assumed to contain the eigenvalues (but   
             see the description of EI), otherwise they will be   
             computed according to MODE, COND, DMAX, and RSIGN and   
             placed in D.   
             Modified if MODE is nonzero.   

    MODE   - INTEGER   
             On entry this describes how the eigenvalues are to   
             be specified:   
             MODE = 0 means use D (with EI) as input   
             MODE = 1 sets D(1)=1 and D(2:N)=1.0/COND   
             MODE = 2 sets D(1:N-1)=1 and D(N)=1.0/COND   
             MODE = 3 sets D(I)=COND**(-(I-1)/(N-1))   
             MODE = 4 sets D(i)=1 - (i-1)/(N-1)*(1 - 1/COND)   
             MODE = 5 sets D to random numbers in the range   
                      ( 1/COND , 1 ) such that their logarithms   
                      are uniformly distributed.  Each odd-even pair   
                      of elements will be either used as two real   
                      eigenvalues or as the real and imaginary part   
                      of a complex conjugate pair of eigenvalues;   
                      the choice of which is done is random, with   
                      50-50 probability, for each pair.   
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

    DMAX   - DOUBLE PRECISION   
             If MODE is neither -6, 0 nor 6, the contents of D, as   
             computed according to MODE and COND, will be scaled by   
             DMAX / max(abs(D(i))).  Note that DMAX need not be   
             positive: if DMAX is negative (or zero), D will be   
             scaled by a negative number (or zero).   
             Not modified.   

    EI     - CHARACTER*1 array, dimension ( N )   
             If MODE is 0, and EI(1) is not ' ' (space character),   
             this array specifies which elements of D (on input) are   
             real eigenvalues and which are the real and imaginary parts 
  
             of a complex conjugate pair of eigenvalues.  The elements   
             of EI may then only have the values 'R' and 'I'.  If   
             EI(j)='R' and EI(j+1)='I', then the j-th eigenvalue is   
             CMPLX( D(j) , D(j+1) ), and the (j+1)-th is the complex   
             conjugate thereof.  If EI(j)=EI(j+1)='R', then the j-th   
             eigenvalue is D(j) (i.e., real).  EI(1) may not be 'I',   
             nor may two adjacent elements of EI both have the value 'I'. 
  
             If MODE is not 0, then EI is ignored.  If MODE is 0 and   
             EI(1)=' ', then the eigenvalues will all be real.   
             Not modified.   

    RSIGN  - CHARACTER*1   
             If MODE is not 0, 6, or -6, and RSIGN='T', then the   
             elements of D, as computed according to MODE and COND, will 
  
             be multiplied by a random sign (+1 or -1).  If RSIGN='F',   
             they will not be.  RSIGN may only have the values 'T' or   
             'F'.   
             Not modified.   

    UPPER  - CHARACTER*1   
             If UPPER='T', then the elements of A above the diagonal   
             (and above the 2x2 diagonal blocks, if A has complex   
             eigenvalues) will be set to random numbers out of DIST.   
             If UPPER='F', they will not.  UPPER may only have the   
             values 'T' or 'F'.   
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
             Same as MODE and COND, but for specifying the diagonal   
             of S.  MODES=-6 and +6 are not allowed (since they would   
             result in randomly ill-conditioned eigenvalues.)   

    KL     - INTEGER   
             This specifies the lower bandwidth of the  matrix.  KL=1   
             specifies upper Hessenberg form.  If KL is at least N-1,   
             then A will have full lower bandwidth.  KL must be at   
             least 1.   
             Not modified.   

    KU     - INTEGER   
             This specifies the upper bandwidth of the  matrix.  KU=1   
             specifies lower Hessenberg form.  If KU is at least N-1,   
             then A will have full upper bandwidth; if KU and KL   
             are both at least N-1, then A will be dense.  Only one of   
             KU and KL may be less than N-1.  KU must be at least 1.   
             Not modified.   

    ANORM  - DOUBLE PRECISION   
             If ANORM is not negative, then A will be scaled by a non-   
             negative real number to make the maximum-element-norm of A   
             to be ANORM.   
             Not modified.   

    A      - DOUBLE PRECISION array, dimension ( LDA, N )   
             On exit A is the desired test matrix.   
             Modified.   

    LDA    - INTEGER   
             LDA specifies the first dimension of A as declared in the   
             calling program.  LDA must be at least N.   
             Not modified.   

    WORK   - DOUBLE PRECISION array, dimension ( 3*N )   
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
              -8 => EI(1) is not ' ' or 'R', EI(j) is not 'R' or 'I', or 
  
                    two adjacent elements of EI are 'I'.   
              -9 => RSIGN is not 'T' or 'F'   
             -10 => UPPER is not 'T' or 'F'   
             -11 => SIM   is not 'T' or 'F'   
             -12 => MODES=0 and DS has a zero singular value.   
             -13 => MODES is not in the range -5 to 5.   
             -14 => MODES is nonzero and CONDS is less than 1.   
             -15 => KL is less than 1.   
             -16 => KU is less than 1, or KL and KU are both less than   
                    N-1.   
             -19 => LDA is less than N.   
              1  => Error return from DLATM1 (computing D)   
              2  => Cannot scale to DMAX (max. eigenvalue is 0)   
              3  => Error return from DLATM1 (computing DS)   
              4  => Error return from DLARGE   
              5  => Zero singular value from DLATM1.   

    ===================================================================== 
  


       1)      Decode and Test the input parameters.   
               Initialize flags & seed.   

       Parameter adjustments */
    --iseed;
    --d;
    --ei;
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
    } else {
	idist = -1;
    }

/*     Check EI */

    useei = TRUE_;
    badei = FALSE_;
    if (lsame_(ei + 1, " ") || *mode != 0) {
	useei = FALSE_;
    } else {
	if (lsame_(ei + 1, "R")) {
	    i__1 = *n;
	    for (j = 2; j <= i__1; ++j) {
		if (lsame_(ei + j, "I")) {
		    if (lsame_(ei + (j - 1), "I")) {
			badei = TRUE_;
		    }
		} else {
		    if (! lsame_(ei + j, "R")) {
			badei = TRUE_;
		    }
		}
/* L10: */
	    }
	} else {
	    badei = TRUE_;
	}
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
/* L20: */
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
    } else if (badei) {
	*info = -8;
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
	xerbla_("DLATME", &i__1);
	return 0;
    }

/*     Initialize random number generator */

    for (i = 1; i <= 4; ++i) {
	iseed[i] = (i__1 = iseed[i], abs(i__1)) % 4096;
/* L30: */
    }

    if (iseed[4] % 2 != 1) {
	++iseed[4];
    }

/*     2)      Set up diagonal of A   

               Compute D according to COND and MODE */

    dlatm1_(mode, cond, &irsign, &idist, &iseed[1], &d[1], n, &iinfo);
    if (iinfo != 0) {
	*info = 1;
	return 0;
    }
    if (*mode != 0 && abs(*mode) != 6) {

/*        Scale by DMAX */

	temp = abs(d[1]);
	i__1 = *n;
	for (i = 2; i <= i__1; ++i) {
/* Computing MAX */
	    d__2 = temp, d__3 = (d__1 = d[i], abs(d__1));
	    temp = max(d__2,d__3);
/* L40: */
	}

	if (temp > 0.) {
	    alpha = *dmax__ / temp;
	} else if (*dmax__ != 0.) {
	    *info = 2;
	    return 0;
	} else {
	    alpha = 0.;
	}

	dscal_(n, &alpha, &d[1], &c__1);

    }

    dlaset_("Full", n, n, &c_b23, &c_b23, &a[a_offset], lda);
    i__1 = *lda + 1;
    dcopy_(n, &d[1], &c__1, &a[a_offset], &i__1);

/*     Set up complex conjugate pairs */

    if (*mode == 0) {
	if (useei) {
	    i__1 = *n;
	    for (j = 2; j <= i__1; ++j) {
		if (lsame_(ei + j, "I")) {
		    a[j - 1 + j * a_dim1] = a[j + j * a_dim1];
		    a[j + (j - 1) * a_dim1] = -a[j + j * a_dim1];
		    a[j + j * a_dim1] = a[j - 1 + (j - 1) * a_dim1];
		}
/* L50: */
	    }
	}

    } else if (abs(*mode) == 5) {

	i__1 = *n;
	for (j = 2; j <= i__1; j += 2) {
	    if (dlaran_(&iseed[1]) > .5) {
		a[j - 1 + j * a_dim1] = a[j + j * a_dim1];
		a[j + (j - 1) * a_dim1] = -a[j + j * a_dim1];
		a[j + j * a_dim1] = a[j - 1 + (j - 1) * a_dim1];
	    }
/* L60: */
	}
    }

/*     3)      If UPPER='T', set upper triangle of A to random numbers.   
               (but don't modify the corners of 2x2 blocks.) */

    if (iupper != 0) {
	i__1 = *n;
	for (jc = 2; jc <= i__1; ++jc) {
	    if (a[jc - 1 + jc * a_dim1] != 0.) {
		jr = jc - 2;
	    } else {
		jr = jc - 1;
	    }
	    dlarnv_(&idist, &iseed[1], &jr, &a[jc * a_dim1 + 1]);
/* L70: */
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

	dlarge_(n, &a[a_offset], lda, &iseed[1], &work[1], &iinfo);
	if (iinfo != 0) {
	    *info = 4;
	    return 0;
	}

/*        Multiply by S and (1/S) */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    dscal_(n, &ds[j], &a[j + a_dim1], lda);
	    if (ds[j] != 0.) {
		d__1 = 1. / ds[j];
		dscal_(n, &d__1, &a[j * a_dim1 + 1], &c__1);
	    } else {
		*info = 5;
		return 0;
	    }
/* L80: */
	}

/*        Multiply by U and U' */

	dlarge_(n, &a[a_offset], lda, &iseed[1], &work[1], &iinfo);
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

	    dcopy_(&irows, &a[jcr + ic * a_dim1], &c__1, &work[1], &c__1);
	    xnorms = work[1];
	    dlarfg_(&irows, &xnorms, &work[2], &c__1, &tau);
	    work[1] = 1.;

	    dgemv_("T", &irows, &icols, &c_b39, &a[jcr + (ic + 1) * a_dim1], 
		    lda, &work[1], &c__1, &c_b23, &work[irows + 1], &c__1)
		    ;
	    d__1 = -tau;
	    dger_(&irows, &icols, &d__1, &work[1], &c__1, &work[irows + 1], &
		    c__1, &a[jcr + (ic + 1) * a_dim1], lda);

	    dgemv_("N", n, &irows, &c_b39, &a[jcr * a_dim1 + 1], lda, &work[1]
		    , &c__1, &c_b23, &work[irows + 1], &c__1);
	    d__1 = -tau;
	    dger_(n, &irows, &d__1, &work[irows + 1], &c__1, &work[1], &c__1, 
		    &a[jcr * a_dim1 + 1], lda);

	    a[jcr + ic * a_dim1] = xnorms;
	    i__2 = irows - 1;
	    dlaset_("Full", &i__2, &c__1, &c_b23, &c_b23, &a[jcr + 1 + ic * 
		    a_dim1], lda);
/* L90: */
	}
    } else if (*ku < *n - 1) {

/*        Reduce upper bandwidth -- kill a row at a time. */

	i__1 = *n - 1;
	for (jcr = *ku + 1; jcr <= i__1; ++jcr) {
	    ir = jcr - *ku;
	    irows = *n + *ku - jcr;
	    icols = *n + 1 - jcr;

	    dcopy_(&icols, &a[ir + jcr * a_dim1], lda, &work[1], &c__1);
	    xnorms = work[1];
	    dlarfg_(&icols, &xnorms, &work[2], &c__1, &tau);
	    work[1] = 1.;

	    dgemv_("N", &irows, &icols, &c_b39, &a[ir + 1 + jcr * a_dim1], 
		    lda, &work[1], &c__1, &c_b23, &work[icols + 1], &c__1)
		    ;
	    d__1 = -tau;
	    dger_(&irows, &icols, &d__1, &work[icols + 1], &c__1, &work[1], &
		    c__1, &a[ir + 1 + jcr * a_dim1], lda);

	    dgemv_("C", &icols, n, &c_b39, &a[jcr + a_dim1], lda, &work[1], &
		    c__1, &c_b23, &work[icols + 1], &c__1);
	    d__1 = -tau;
	    dger_(&icols, n, &d__1, &work[1], &c__1, &work[icols + 1], &c__1, 
		    &a[jcr + a_dim1], lda);

	    a[ir + jcr * a_dim1] = xnorms;
	    i__2 = icols - 1;
	    dlaset_("Full", &c__1, &i__2, &c_b23, &c_b23, &a[ir + (jcr + 1) * 
		    a_dim1], lda);
/* L100: */
	}
    }

/*     Scale the matrix to have norm ANORM */

    if (*anorm >= 0.) {
	temp = dlange_("M", n, n, &a[a_offset], lda, tempa);
	if (temp > 0.) {
	    alpha = *anorm / temp;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		dscal_(n, &alpha, &a[j * a_dim1 + 1], &c__1);
/* L110: */
	    }
	}
    }

    return 0;

/*     End of DLATME */

} /* dlatme_ */

