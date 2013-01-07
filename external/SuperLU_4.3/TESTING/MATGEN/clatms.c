/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static complex c_b1 = {0.f,0.f};
static integer c__1 = 1;
static integer c__5 = 5;
static logical c_true = TRUE_;
static logical c_false = FALSE_;

/* Subroutine */ int clatms_(integer *m, integer *n, char *dist, integer *
	iseed, char *sym, real *d, integer *mode, real *cond, real *dmax__, 
	integer *kl, integer *ku, char *pack, complex *a, integer *lda, 
	complex *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    real r__1, r__2, r__3;
    doublereal d__1;
    complex q__1, q__2, q__3;
    logical L__1;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);
    void r_cnjg(complex *, complex *);

    /* Local variables */
    static integer ilda, icol;
    static real temp;
    static logical csym;
    static integer irow, isym;
    static complex c;
    static integer i, j, k;
    static complex s;
    static real alpha, angle;
    static integer ipack;
    static real realc;
    static integer ioffg;
    extern logical lsame_(char *, char *);
    static integer iinfo;
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *);
    static complex ctemp;
    static integer idist, mnmin, iskew;
    static complex extra, dummy;
    extern /* Subroutine */ int slatm1_(integer *, real *, integer *, integer 
	    *, integer *, real *, integer *, integer *);
    static integer ic, jc, nc;
    extern /* Subroutine */ int clagge_(integer *, integer *, integer *, 
	    integer *, real *, complex *, integer *, integer *, complex *, 
	    integer *), claghe_(integer *, integer *, real *, complex *, 
	    integer *, integer *, complex *, integer *);
    static integer il;
    static complex ct;
    static integer iendch, ir, jr, ipackg, mr;
    extern /* Complex */ VOID clarnd_(complex *, integer *, integer *);
    static integer minlda;
    static complex st;
    extern /* Subroutine */ int claset_(char *, integer *, integer *, complex 
	    *, complex *, complex *, integer *), clartg_(complex *, 
	    complex *, real *, complex *, complex *), xerbla_(char *, integer 
	    *), clagsy_(integer *, integer *, real *, complex *, 
	    integer *, integer *, complex *, integer *);
    extern doublereal slarnd_(integer *, integer *);
    extern /* Subroutine */ int clarot_(logical *, logical *, logical *, 
	    integer *, complex *, complex *, complex *, integer *, complex *, 
	    complex *);
    static logical iltemp, givens;
    static integer ioffst, irsign;
    static logical ilextr, topdwn;
    static integer ir1, ir2, isympk, jch, llb, jkl, jku, uub;


/*  -- LAPACK test routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

       CLATMS generates random matrices with specified singular values   
       (or hermitian with specified eigenvalues)   
       for testing LAPACK programs.   

       CLATMS operates by applying the following sequence of   
       operations:   

         Set the diagonal to D, where D may be input or   
            computed according to MODE, COND, DMAX, and SYM   
            as described below.   

         Generate a matrix with the appropriate band structure, by one   
            of two methods:   

         Method A:   
             Generate a dense M x N matrix by multiplying D on the left   
                 and the right by random unitary matrices, then:   

             Reduce the bandwidth according to KL and KU, using   
                 Householder transformations.   

         Method B:   
             Convert the bandwidth-0 (i.e., diagonal) matrix to a   
                 bandwidth-1 matrix using Givens rotations, "chasing"   
                 out-of-band elements back, much as in QR; then convert   
                 the bandwidth-1 to a bandwidth-2 matrix, etc.  Note   
                 that for reasonably small bandwidths (relative to M and 
  
                 N) this requires less storage, as a dense matrix is not 
  
                 generated.  Also, for hermitian or symmetric matrices,   
                 only one triangle is generated.   

         Method A is chosen if the bandwidth is a large fraction of the   
             order of the matrix, and LDA is at least M (so a dense   
             matrix can be stored.)  Method B is chosen if the bandwidth 
  
             is small (< 1/2 N for hermitian or symmetric, < .3 N+M for   
             non-symmetric), or LDA is less than M and not less than the 
  
             bandwidth.   

         Pack the matrix if desired. Options specified by PACK are:   
            no packing   
            zero out upper half (if hermitian)   
            zero out lower half (if hermitian)   
            store the upper half columnwise (if hermitian or upper   
                  triangular)   
            store the lower half columnwise (if hermitian or lower   
                  triangular)   
            store the lower triangle in banded format (if hermitian or   
                  lower triangular)   
            store the upper triangle in banded format (if hermitian or   
                  upper triangular)   
            store the entire matrix in banded format   
         If Method B is chosen, and band format is specified, then the   
            matrix will be generated in the band format, so no repacking 
  
            will be necessary.   

    Arguments   
    =========   

    M      - INTEGER   
             The number of rows of A. Not modified.   

    N      - INTEGER   
             The number of columns of A. N must equal M if the matrix   
             is symmetric or hermitian (i.e., if SYM is not 'N')   
             Not modified.   

    DIST   - CHARACTER*1   
             On entry, DIST specifies the type of distribution to be used 
  
             to generate the random eigen-/singular values.   
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
             exit, and can be used in the next call to CLATMS   
             to continue the same random number sequence.   
             Changed on exit.   

    SYM    - CHARACTER*1   
             If SYM='H', the generated matrix is hermitian, with   
               eigenvalues specified by D, COND, MODE, and DMAX; they   
               may be positive, negative, or zero.   
             If SYM='P', the generated matrix is hermitian, with   
               eigenvalues (= singular values) specified by D, COND,   
               MODE, and DMAX; they will not be negative.   
             If SYM='N', the generated matrix is nonsymmetric, with   
               singular values specified by D, COND, MODE, and DMAX;   
               they will not be negative.   
             If SYM='S', the generated matrix is (complex) symmetric,   
               with singular values specified by D, COND, MODE, and   
               DMAX; they will not be negative.   
             Not modified.   

    D      - REAL array, dimension ( MIN( M, N ) )   
             This array is used to specify the singular values or   
             eigenvalues of A (see SYM, above.)  If MODE=0, then D is   
             assumed to contain the singular/eigenvalues, otherwise   
             they will be computed according to MODE, COND, and DMAX,   
             and placed in D.   
             Modified if MODE is nonzero.   

    MODE   - INTEGER   
             On entry this describes how the singular/eigenvalues are to 
  
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
             Thus if MODE is positive, D has entries ranging from   
                1 to 1/COND, if negative, from 1/COND to 1,   
             If SYM='H', and MODE is neither 0, 6, nor -6, then   
                the elements of D will also be multiplied by a random   
                sign (i.e., +1 or -1.)   
             Not modified.   

    COND   - REAL   
             On entry, this is used as described under MODE above.   
             If used, it must be >= 1. Not modified.   

    DMAX   - REAL   
             If MODE is neither -6, 0 nor 6, the contents of D, as   
             computed according to MODE and COND, will be scaled by   
             DMAX / max(abs(D(i))); thus, the maximum absolute eigen- or 
  
             singular value (which is to say the norm) will be abs(DMAX). 
  
             Note that DMAX need not be positive: if DMAX is negative   
             (or zero), D will be scaled by a negative number (or zero). 
  
             Not modified.   

    KL     - INTEGER   
             This specifies the lower bandwidth of the  matrix. For   
             example, KL=0 implies upper triangular, KL=1 implies upper   
             Hessenberg, and KL being at least M-1 means that the matrix 
  
             has full lower bandwidth.  KL must equal KU if the matrix   
             is symmetric or hermitian.   
             Not modified.   

    KU     - INTEGER   
             This specifies the upper bandwidth of the  matrix. For   
             example, KU=0 implies lower triangular, KU=1 implies lower   
             Hessenberg, and KU being at least N-1 means that the matrix 
  
             has full upper bandwidth.  KL must equal KU if the matrix   
             is symmetric or hermitian.   
             Not modified.   

    PACK   - CHARACTER*1   
             This specifies packing of matrix as follows:   
             'N' => no packing   
             'U' => zero out all subdiagonal entries (if symmetric   
                    or hermitian)   
             'L' => zero out all superdiagonal entries (if symmetric   
                    or hermitian)   
             'C' => store the upper triangle columnwise (only if the   
                    matrix is symmetric, hermitian, or upper triangular) 
  
             'R' => store the lower triangle columnwise (only if the   
                    matrix is symmetric, hermitian, or lower triangular) 
  
             'B' => store the lower triangle in band storage scheme   
                    (only if the matrix is symmetric, hermitian, or   
                    lower triangular)   
             'Q' => store the upper triangle in band storage scheme   
                    (only if the matrix is symmetric, hermitian, or   
                    upper triangular)   
             'Z' => store the entire matrix in band storage scheme   
                        (pivoting can be provided for by using this   
                        option to store A in the trailing rows of   
                        the allocated storage)   

             Using these options, the various LAPACK packed and banded   
             storage schemes can be obtained:   
             GB                    - use 'Z'   
             PB, SB, HB, or TB     - use 'B' or 'Q'   
             PP, SP, HB, or TP     - use 'C' or 'R'   

             If two calls to CLATMS differ only in the PACK parameter,   
             they will generate mathematically equivalent matrices.   
             Not modified.   

    A      - COMPLEX array, dimension ( LDA, N )   
             On exit A is the desired test matrix.  A is first generated 
  
             in full (unpacked) form, and then packed, if so specified   
             by PACK.  Thus, the first M elements of the first N   
             columns will always be modified.  If PACK specifies a   
             packed or banded storage scheme, all LDA elements of the   
             first N columns will be modified; the elements of the   
             array which do not correspond to elements of the generated   
             matrix are set to zero.   
             Modified.   

    LDA    - INTEGER   
             LDA specifies the first dimension of A as declared in the   
             calling program.  If PACK='N', 'U', 'L', 'C', or 'R', then   
             LDA must be at least M.  If PACK='B' or 'Q', then LDA must   
             be at least MIN( KL, M-1) (which is equal to MIN(KU,N-1)).   
             If PACK='Z', LDA must be large enough to hold the packed   
             array: MIN( KU, N-1) + MIN( KL, M-1) + 1.   
             Not modified.   

    WORK   - COMPLEX array, dimension ( 3*MAX( N, M ) )   
             Workspace.   
             Modified.   

    INFO   - INTEGER   
             Error code.  On exit, INFO will be set to one of the   
             following values:   
               0 => normal return   
              -1 => M negative or unequal to N and SYM='S', 'H', or 'P'   
              -2 => N negative   
              -3 => DIST illegal string   
              -5 => SYM illegal string   
              -7 => MODE not in range -6 to 6   
              -8 => COND less than 1.0, and MODE neither -6, 0 nor 6   
             -10 => KL negative   
             -11 => KU negative, or SYM is not 'N' and KU is not equal to 
  
                    KL   
             -12 => PACK illegal string, or PACK='U' or 'L', and SYM='N'; 
  
                    or PACK='C' or 'Q' and SYM='N' and KL is not zero;   
                    or PACK='R' or 'B' and SYM='N' and KU is not zero;   
                    or PACK='U', 'L', 'C', 'R', 'B', or 'Q', and M is not 
  
                    N.   
             -14 => LDA is less than M, or PACK='Z' and LDA is less than 
  
                    MIN(KU,N-1) + MIN(KL,M-1) + 1.   
              1  => Error return from SLATM1   
              2  => Cannot scale to DMAX (max. sing. value is 0)   
              3  => Error return from CLAGGE, CLAGHE or CLAGSY   

    ===================================================================== 
  


       1)      Decode and Test the input parameters.   
               Initialize flags & seed.   

       Parameter adjustments */
    --iseed;
    --d;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --work;

    /* Function Body */
    *info = 0;

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
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

/*     Decode SYM */

    if (lsame_(sym, "N")) {
	isym = 1;
	irsign = 0;
	csym = FALSE_;
    } else if (lsame_(sym, "P")) {
	isym = 2;
	irsign = 0;
	csym = FALSE_;
    } else if (lsame_(sym, "S")) {
	isym = 2;
	irsign = 0;
	csym = TRUE_;
    } else if (lsame_(sym, "H")) {
	isym = 2;
	irsign = 1;
	csym = FALSE_;
    } else {
	isym = -1;
    }

/*     Decode PACK */

    isympk = 0;
    if (lsame_(pack, "N")) {
	ipack = 0;
    } else if (lsame_(pack, "U")) {
	ipack = 1;
	isympk = 1;
    } else if (lsame_(pack, "L")) {
	ipack = 2;
	isympk = 1;
    } else if (lsame_(pack, "C")) {
	ipack = 3;
	isympk = 2;
    } else if (lsame_(pack, "R")) {
	ipack = 4;
	isympk = 3;
    } else if (lsame_(pack, "B")) {
	ipack = 5;
	isympk = 3;
    } else if (lsame_(pack, "Q")) {
	ipack = 6;
	isympk = 2;
    } else if (lsame_(pack, "Z")) {
	ipack = 7;
    } else {
	ipack = -1;
    }

/*     Set certain internal parameters */

    mnmin = min(*m,*n);
/* Computing MIN */
    i__1 = *kl, i__2 = *m - 1;
    llb = min(i__1,i__2);
/* Computing MIN */
    i__1 = *ku, i__2 = *n - 1;
    uub = min(i__1,i__2);
/* Computing MIN */
    i__1 = *m, i__2 = *n + llb;
    mr = min(i__1,i__2);
/* Computing MIN */
    i__1 = *n, i__2 = *m + uub;
    nc = min(i__1,i__2);

    if (ipack == 5 || ipack == 6) {
	minlda = uub + 1;
    } else if (ipack == 7) {
	minlda = llb + uub + 1;
    } else {
	minlda = *m;
    }

/*     Use Givens rotation method if bandwidth small enough,   
       or if LDA is too small to store the matrix unpacked. */

    givens = FALSE_;
    if (isym == 1) {
/* Computing MAX */
	i__1 = 1, i__2 = mr + nc;
	if ((real) (llb + uub) < (real) max(i__1,i__2) * .3f) {
	    givens = TRUE_;
	}
    } else {
	if (llb << 1 < *m) {
	    givens = TRUE_;
	}
    }
    if (*lda < *m && *lda >= minlda) {
	givens = TRUE_;
    }

/*     Set INFO if an error */

    if (*m < 0) {
	*info = -1;
    } else if (*m != *n && isym != 1) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (idist == -1) {
	*info = -3;
    } else if (isym == -1) {
	*info = -5;
    } else if (abs(*mode) > 6) {
	*info = -7;
    } else if (*mode != 0 && abs(*mode) != 6 && *cond < 1.f) {
	*info = -8;
    } else if (*kl < 0) {
	*info = -10;
    } else if (*ku < 0 || isym != 1 && *kl != *ku) {
	*info = -11;
    } else if (ipack == -1 || isympk == 1 && isym == 1 || isympk == 2 && isym 
	    == 1 && *kl > 0 || isympk == 3 && isym == 1 && *ku > 0 || isympk 
	    != 0 && *m != *n) {
	*info = -12;
    } else if (*lda < max(1,minlda)) {
	*info = -14;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CLATMS", &i__1);
	return 0;
    }

/*     Initialize random number generator */

    for (i = 1; i <= 4; ++i) {
	iseed[i] = (i__1 = iseed[i], abs(i__1)) % 4096;
/* L10: */
    }

    if (iseed[4] % 2 != 1) {
	++iseed[4];
    }

/*     2)      Set up D  if indicated.   

               Compute D according to COND and MODE */

    slatm1_(mode, cond, &irsign, &idist, &iseed[1], &d[1], &mnmin, &iinfo);
    if (iinfo != 0) {
	*info = 1;
	return 0;
    }

/*     Choose Top-Down if D is (apparently) increasing,   
       Bottom-Up if D is (apparently) decreasing. */

    if (dabs(d[1]) <= (r__1 = d[mnmin], dabs(r__1))) {
	topdwn = TRUE_;
    } else {
	topdwn = FALSE_;
    }

    if (*mode != 0 && abs(*mode) != 6) {

/*        Scale by DMAX */

	temp = dabs(d[1]);
	i__1 = mnmin;
	for (i = 2; i <= i__1; ++i) {
/* Computing MAX */
	    r__2 = temp, r__3 = (r__1 = d[i], dabs(r__1));
	    temp = dmax(r__2,r__3);
/* L20: */
	}

	if (temp > 0.f) {
	    alpha = *dmax__ / temp;
	} else {
	    *info = 2;
	    return 0;
	}

	sscal_(&mnmin, &alpha, &d[1], &c__1);

    }

    claset_("Full", lda, n, &c_b1, &c_b1, &a[a_offset], lda);

/*     3)      Generate Banded Matrix using Givens rotations.   
               Also the special case of UUB=LLB=0   

                 Compute Addressing constants to cover all   
                 storage formats.  Whether GE, HE, SY, GB, HB, or SB,   
                 upper or lower triangle or both,   
                 the (i,j)-th element is in   
                 A( i - ISKEW*j + IOFFST, j ) */

    if (ipack > 4) {
	ilda = *lda - 1;
	iskew = 1;
	if (ipack > 5) {
	    ioffst = uub + 1;
	} else {
	    ioffst = 1;
	}
    } else {
	ilda = *lda;
	iskew = 0;
	ioffst = 0;
    }

/*     IPACKG is the format that the matrix is generated in. If this is   
       different from IPACK, then the matrix must be repacked at the   
       end.  It also signals how to compute the norm, for scaling. */

    ipackg = 0;

/*     Diagonal Matrix -- We are done, unless it   
       is to be stored HP/SP/PP/TP (PACK='R' or 'C') */

    if (llb == 0 && uub == 0) {
	i__1 = mnmin;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = (1 - iskew) * j + ioffst + j * a_dim1;
	    i__3 = j;
	    q__1.r = d[i__3], q__1.i = 0.f;
	    a[i__2].r = q__1.r, a[i__2].i = q__1.i;
/* L30: */
	}

	if (ipack <= 2 || ipack >= 5) {
	    ipackg = ipack;
	}

    } else if (givens) {

/*        Check whether to use Givens rotations,   
          Householder transformations, or nothing. */

	if (isym == 1) {

/*           Non-symmetric -- A = U D V */

	    if (ipack > 4) {
		ipackg = ipack;
	    } else {
		ipackg = 0;
	    }

	    i__1 = mnmin;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = (1 - iskew) * j + ioffst + j * a_dim1;
		i__3 = j;
		q__1.r = d[i__3], q__1.i = 0.f;
		a[i__2].r = q__1.r, a[i__2].i = q__1.i;
/* L40: */
	    }

	    if (topdwn) {
		jkl = 0;
		i__1 = uub;
		for (jku = 1; jku <= i__1; ++jku) {

/*                 Transform from bandwidth JKL, JKU-1 to 
JKL, JKU   

                   Last row actually rotated is M   
                   Last column actually rotated is MIN( M+
JKU, N )   

   Computing MIN */
		    i__3 = *m + jku;
		    i__2 = min(i__3,*n) + jkl - 1;
		    for (jr = 1; jr <= i__2; ++jr) {
			extra.r = 0.f, extra.i = 0.f;
			angle = slarnd_(&c__1, &iseed[1]) * 
				6.2831853071795864769252867663f;
			d__1 = cos(angle);
			clarnd_(&q__2, &c__5, &iseed[1]);
			q__1.r = d__1 * q__2.r, q__1.i = d__1 * q__2.i;
			c.r = q__1.r, c.i = q__1.i;
			d__1 = sin(angle);
			clarnd_(&q__2, &c__5, &iseed[1]);
			q__1.r = d__1 * q__2.r, q__1.i = d__1 * q__2.i;
			s.r = q__1.r, s.i = q__1.i;
/* Computing MAX */
			i__3 = 1, i__4 = jr - jkl;
			icol = max(i__3,i__4);
			if (jr < *m) {
/* Computing MIN */
			    i__3 = *n, i__4 = jr + jku;
			    il = min(i__3,i__4) + 1 - icol;
			    L__1 = jr > jkl;
			    clarot_(&c_true, &L__1, &c_false, &il, &c, &s, &a[
				    jr - iskew * icol + ioffst + icol * 
				    a_dim1], &ilda, &extra, &dummy);
			}

/*                    Chase "EXTRA" back up */

			ir = jr;
			ic = icol;
			i__3 = -jkl - jku;
			for (jch = jr - jkl; i__3 < 0 ? jch >= 1 : jch <= 1; 
				jch += i__3) {
			    if (ir < *m) {
				clartg_(&a[ir + 1 - iskew * (ic + 1) + ioffst 
					+ (ic + 1) * a_dim1], &extra, &realc, 
					&s, &dummy);
				clarnd_(&q__1, &c__5, &iseed[1]);
				dummy.r = q__1.r, dummy.i = q__1.i;
				q__2.r = realc * dummy.r, q__2.i = realc * 
					dummy.i;
				r_cnjg(&q__1, &q__2);
				c.r = q__1.r, c.i = q__1.i;
				q__3.r = -(doublereal)s.r, q__3.i = -(
					doublereal)s.i;
				q__2.r = q__3.r * dummy.r - q__3.i * dummy.i, 
					q__2.i = q__3.r * dummy.i + q__3.i * 
					dummy.r;
				r_cnjg(&q__1, &q__2);
				s.r = q__1.r, s.i = q__1.i;
			    }
/* Computing MAX */
			    i__4 = 1, i__5 = jch - jku;
			    irow = max(i__4,i__5);
			    il = ir + 2 - irow;
			    ctemp.r = 0.f, ctemp.i = 0.f;
			    iltemp = jch > jku;
			    clarot_(&c_false, &iltemp, &c_true, &il, &c, &s, &
				    a[irow - iskew * ic + ioffst + ic * 
				    a_dim1], &ilda, &ctemp, &extra);
			    if (iltemp) {
				clartg_(&a[irow + 1 - iskew * (ic + 1) + 
					ioffst + (ic + 1) * a_dim1], &ctemp, &
					realc, &s, &dummy);
				clarnd_(&q__1, &c__5, &iseed[1]);
				dummy.r = q__1.r, dummy.i = q__1.i;
				q__2.r = realc * dummy.r, q__2.i = realc * 
					dummy.i;
				r_cnjg(&q__1, &q__2);
				c.r = q__1.r, c.i = q__1.i;
				q__3.r = -(doublereal)s.r, q__3.i = -(
					doublereal)s.i;
				q__2.r = q__3.r * dummy.r - q__3.i * dummy.i, 
					q__2.i = q__3.r * dummy.i + q__3.i * 
					dummy.r;
				r_cnjg(&q__1, &q__2);
				s.r = q__1.r, s.i = q__1.i;

/* Computing MAX */
				i__4 = 1, i__5 = jch - jku - jkl;
				icol = max(i__4,i__5);
				il = ic + 2 - icol;
				extra.r = 0.f, extra.i = 0.f;
				L__1 = jch > jku + jkl;
				clarot_(&c_true, &L__1, &c_true, &il, &c, &s, 
					&a[irow - iskew * icol + ioffst + 
					icol * a_dim1], &ilda, &extra, &ctemp)
					;
				ic = icol;
				ir = irow;
			    }
/* L50: */
			}
/* L60: */
		    }
/* L70: */
		}

		jku = uub;
		i__1 = llb;
		for (jkl = 1; jkl <= i__1; ++jkl) {

/*                 Transform from bandwidth JKL-1, JKU to 
JKL, JKU   

   Computing MIN */
		    i__3 = *n + jkl;
		    i__2 = min(i__3,*m) + jku - 1;
		    for (jc = 1; jc <= i__2; ++jc) {
			extra.r = 0.f, extra.i = 0.f;
			angle = slarnd_(&c__1, &iseed[1]) * 
				6.2831853071795864769252867663f;
			d__1 = cos(angle);
			clarnd_(&q__2, &c__5, &iseed[1]);
			q__1.r = d__1 * q__2.r, q__1.i = d__1 * q__2.i;
			c.r = q__1.r, c.i = q__1.i;
			d__1 = sin(angle);
			clarnd_(&q__2, &c__5, &iseed[1]);
			q__1.r = d__1 * q__2.r, q__1.i = d__1 * q__2.i;
			s.r = q__1.r, s.i = q__1.i;
/* Computing MAX */
			i__3 = 1, i__4 = jc - jku;
			irow = max(i__3,i__4);
			if (jc < *n) {
/* Computing MIN */
			    i__3 = *m, i__4 = jc + jkl;
			    il = min(i__3,i__4) + 1 - irow;
			    L__1 = jc > jku;
			    clarot_(&c_false, &L__1, &c_false, &il, &c, &s, &
				    a[irow - iskew * jc + ioffst + jc * 
				    a_dim1], &ilda, &extra, &dummy);
			}

/*                    Chase "EXTRA" back up */

			ic = jc;
			ir = irow;
			i__3 = -jkl - jku;
			for (jch = jc - jku; i__3 < 0 ? jch >= 1 : jch <= 1; 
				jch += i__3) {
			    if (ic < *n) {
				clartg_(&a[ir + 1 - iskew * (ic + 1) + ioffst 
					+ (ic + 1) * a_dim1], &extra, &realc, 
					&s, &dummy);
				clarnd_(&q__1, &c__5, &iseed[1]);
				dummy.r = q__1.r, dummy.i = q__1.i;
				q__2.r = realc * dummy.r, q__2.i = realc * 
					dummy.i;
				r_cnjg(&q__1, &q__2);
				c.r = q__1.r, c.i = q__1.i;
				q__3.r = -(doublereal)s.r, q__3.i = -(
					doublereal)s.i;
				q__2.r = q__3.r * dummy.r - q__3.i * dummy.i, 
					q__2.i = q__3.r * dummy.i + q__3.i * 
					dummy.r;
				r_cnjg(&q__1, &q__2);
				s.r = q__1.r, s.i = q__1.i;
			    }
/* Computing MAX */
			    i__4 = 1, i__5 = jch - jkl;
			    icol = max(i__4,i__5);
			    il = ic + 2 - icol;
			    ctemp.r = 0.f, ctemp.i = 0.f;
			    iltemp = jch > jkl;
			    clarot_(&c_true, &iltemp, &c_true, &il, &c, &s, &
				    a[ir - iskew * icol + ioffst + icol * 
				    a_dim1], &ilda, &ctemp, &extra);
			    if (iltemp) {
				clartg_(&a[ir + 1 - iskew * (icol + 1) + 
					ioffst + (icol + 1) * a_dim1], &ctemp,
					 &realc, &s, &dummy);
				clarnd_(&q__1, &c__5, &iseed[1]);
				dummy.r = q__1.r, dummy.i = q__1.i;
				q__2.r = realc * dummy.r, q__2.i = realc * 
					dummy.i;
				r_cnjg(&q__1, &q__2);
				c.r = q__1.r, c.i = q__1.i;
				q__3.r = -(doublereal)s.r, q__3.i = -(
					doublereal)s.i;
				q__2.r = q__3.r * dummy.r - q__3.i * dummy.i, 
					q__2.i = q__3.r * dummy.i + q__3.i * 
					dummy.r;
				r_cnjg(&q__1, &q__2);
				s.r = q__1.r, s.i = q__1.i;
/* Computing MAX */
				i__4 = 1, i__5 = jch - jkl - jku;
				irow = max(i__4,i__5);
				il = ir + 2 - irow;
				extra.r = 0.f, extra.i = 0.f;
				L__1 = jch > jkl + jku;
				clarot_(&c_false, &L__1, &c_true, &il, &c, &s,
					 &a[irow - iskew * icol + ioffst + 
					icol * a_dim1], &ilda, &extra, &ctemp)
					;
				ic = icol;
				ir = irow;
			    }
/* L80: */
			}
/* L90: */
		    }
/* L100: */
		}

	    } else {

/*              Bottom-Up -- Start at the bottom right. */

		jkl = 0;
		i__1 = uub;
		for (jku = 1; jku <= i__1; ++jku) {

/*                 Transform from bandwidth JKL, JKU-1 to 
JKL, JKU   

                   First row actually rotated is M   
                   First column actually rotated is MIN( M
+JKU, N )   

   Computing MIN */
		    i__2 = *m, i__3 = *n + jkl;
		    iendch = min(i__2,i__3) - 1;
/* Computing MIN */
		    i__2 = *m + jku;
		    i__3 = 1 - jkl;
		    for (jc = min(i__2,*n) - 1; jc >= i__3; --jc) {
			extra.r = 0.f, extra.i = 0.f;
			angle = slarnd_(&c__1, &iseed[1]) * 
				6.2831853071795864769252867663f;
			d__1 = cos(angle);
			clarnd_(&q__2, &c__5, &iseed[1]);
			q__1.r = d__1 * q__2.r, q__1.i = d__1 * q__2.i;
			c.r = q__1.r, c.i = q__1.i;
			d__1 = sin(angle);
			clarnd_(&q__2, &c__5, &iseed[1]);
			q__1.r = d__1 * q__2.r, q__1.i = d__1 * q__2.i;
			s.r = q__1.r, s.i = q__1.i;
/* Computing MAX */
			i__2 = 1, i__4 = jc - jku + 1;
			irow = max(i__2,i__4);
			if (jc > 0) {
/* Computing MIN */
			    i__2 = *m, i__4 = jc + jkl + 1;
			    il = min(i__2,i__4) + 1 - irow;
			    L__1 = jc + jkl < *m;
			    clarot_(&c_false, &c_false, &L__1, &il, &c, &s, &
				    a[irow - iskew * jc + ioffst + jc * 
				    a_dim1], &ilda, &dummy, &extra);
			}

/*                    Chase "EXTRA" back down */

			ic = jc;
			i__2 = iendch;
			i__4 = jkl + jku;
			for (jch = jc + jkl; i__4 < 0 ? jch >= i__2 : jch <= 
				i__2; jch += i__4) {
			    ilextr = ic > 0;
			    if (ilextr) {
				clartg_(&a[jch - iskew * ic + ioffst + ic * 
					a_dim1], &extra, &realc, &s, &dummy);
				clarnd_(&q__1, &c__5, &iseed[1]);
				dummy.r = q__1.r, dummy.i = q__1.i;
				q__1.r = realc * dummy.r, q__1.i = realc * 
					dummy.i;
				c.r = q__1.r, c.i = q__1.i;
				q__1.r = s.r * dummy.r - s.i * dummy.i, 
					q__1.i = s.r * dummy.i + s.i * 
					dummy.r;
				s.r = q__1.r, s.i = q__1.i;
			    }
			    ic = max(1,ic);
/* Computing MIN */
			    i__5 = *n - 1, i__6 = jch + jku;
			    icol = min(i__5,i__6);
			    iltemp = jch + jku < *n;
			    ctemp.r = 0.f, ctemp.i = 0.f;
			    i__5 = icol + 2 - ic;
			    clarot_(&c_true, &ilextr, &iltemp, &i__5, &c, &s, 
				    &a[jch - iskew * ic + ioffst + ic * 
				    a_dim1], &ilda, &extra, &ctemp);
			    if (iltemp) {
				clartg_(&a[jch - iskew * icol + ioffst + icol 
					* a_dim1], &ctemp, &realc, &s, &dummy)
					;
				clarnd_(&q__1, &c__5, &iseed[1]);
				dummy.r = q__1.r, dummy.i = q__1.i;
				q__1.r = realc * dummy.r, q__1.i = realc * 
					dummy.i;
				c.r = q__1.r, c.i = q__1.i;
				q__1.r = s.r * dummy.r - s.i * dummy.i, 
					q__1.i = s.r * dummy.i + s.i * 
					dummy.r;
				s.r = q__1.r, s.i = q__1.i;
/* Computing MIN */
				i__5 = iendch, i__6 = jch + jkl + jku;
				il = min(i__5,i__6) + 2 - jch;
				extra.r = 0.f, extra.i = 0.f;
				L__1 = jch + jkl + jku <= iendch;
				clarot_(&c_false, &c_true, &L__1, &il, &c, &s,
					 &a[jch - iskew * icol + ioffst + 
					icol * a_dim1], &ilda, &ctemp, &extra)
					;
				ic = icol;
			    }
/* L110: */
			}
/* L120: */
		    }
/* L130: */
		}

		jku = uub;
		i__1 = llb;
		for (jkl = 1; jkl <= i__1; ++jkl) {

/*                 Transform from bandwidth JKL-1, JKU to 
JKL, JKU   

                   First row actually rotated is MIN( N+JK
L, M )   
                   First column actually rotated is N   

   Computing MIN */
		    i__3 = *n, i__4 = *m + jku;
		    iendch = min(i__3,i__4) - 1;
/* Computing MIN */
		    i__3 = *n + jkl;
		    i__4 = 1 - jku;
		    for (jr = min(i__3,*m) - 1; jr >= i__4; --jr) {
			extra.r = 0.f, extra.i = 0.f;
			angle = slarnd_(&c__1, &iseed[1]) * 
				6.2831853071795864769252867663f;
			d__1 = cos(angle);
			clarnd_(&q__2, &c__5, &iseed[1]);
			q__1.r = d__1 * q__2.r, q__1.i = d__1 * q__2.i;
			c.r = q__1.r, c.i = q__1.i;
			d__1 = sin(angle);
			clarnd_(&q__2, &c__5, &iseed[1]);
			q__1.r = d__1 * q__2.r, q__1.i = d__1 * q__2.i;
			s.r = q__1.r, s.i = q__1.i;
/* Computing MAX */
			i__3 = 1, i__2 = jr - jkl + 1;
			icol = max(i__3,i__2);
			if (jr > 0) {
/* Computing MIN */
			    i__3 = *n, i__2 = jr + jku + 1;
			    il = min(i__3,i__2) + 1 - icol;
			    L__1 = jr + jku < *n;
			    clarot_(&c_true, &c_false, &L__1, &il, &c, &s, &a[
				    jr - iskew * icol + ioffst + icol * 
				    a_dim1], &ilda, &dummy, &extra);
			}

/*                    Chase "EXTRA" back down */

			ir = jr;
			i__3 = iendch;
			i__2 = jkl + jku;
			for (jch = jr + jku; i__2 < 0 ? jch >= i__3 : jch <= 
				i__3; jch += i__2) {
			    ilextr = ir > 0;
			    if (ilextr) {
				clartg_(&a[ir - iskew * jch + ioffst + jch * 
					a_dim1], &extra, &realc, &s, &dummy);
				clarnd_(&q__1, &c__5, &iseed[1]);
				dummy.r = q__1.r, dummy.i = q__1.i;
				q__1.r = realc * dummy.r, q__1.i = realc * 
					dummy.i;
				c.r = q__1.r, c.i = q__1.i;
				q__1.r = s.r * dummy.r - s.i * dummy.i, 
					q__1.i = s.r * dummy.i + s.i * 
					dummy.r;
				s.r = q__1.r, s.i = q__1.i;
			    }
			    ir = max(1,ir);
/* Computing MIN */
			    i__5 = *m - 1, i__6 = jch + jkl;
			    irow = min(i__5,i__6);
			    iltemp = jch + jkl < *m;
			    ctemp.r = 0.f, ctemp.i = 0.f;
			    i__5 = irow + 2 - ir;
			    clarot_(&c_false, &ilextr, &iltemp, &i__5, &c, &s,
				     &a[ir - iskew * jch + ioffst + jch * 
				    a_dim1], &ilda, &extra, &ctemp);
			    if (iltemp) {
				clartg_(&a[irow - iskew * jch + ioffst + jch *
					 a_dim1], &ctemp, &realc, &s, &dummy);
				clarnd_(&q__1, &c__5, &iseed[1]);
				dummy.r = q__1.r, dummy.i = q__1.i;
				q__1.r = realc * dummy.r, q__1.i = realc * 
					dummy.i;
				c.r = q__1.r, c.i = q__1.i;
				q__1.r = s.r * dummy.r - s.i * dummy.i, 
					q__1.i = s.r * dummy.i + s.i * 
					dummy.r;
				s.r = q__1.r, s.i = q__1.i;
/* Computing MIN */
				i__5 = iendch, i__6 = jch + jkl + jku;
				il = min(i__5,i__6) + 2 - jch;
				extra.r = 0.f, extra.i = 0.f;
				L__1 = jch + jkl + jku <= iendch;
				clarot_(&c_true, &c_true, &L__1, &il, &c, &s, 
					&a[irow - iskew * jch + ioffst + jch *
					 a_dim1], &ilda, &ctemp, &extra);
				ir = irow;
			    }
/* L140: */
			}
/* L150: */
		    }
/* L160: */
		}

	    }

	} else {

/*           Symmetric -- A = U D U'   
             Hermitian -- A = U D U* */

	    ipackg = ipack;
	    ioffg = ioffst;

	    if (topdwn) {

/*              Top-Down -- Generate Upper triangle only */

		if (ipack >= 5) {
		    ipackg = 6;
		    ioffg = uub + 1;
		} else {
		    ipackg = 1;
		}

		i__1 = mnmin;
		for (j = 1; j <= i__1; ++j) {
		    i__4 = (1 - iskew) * j + ioffg + j * a_dim1;
		    i__2 = j;
		    q__1.r = d[i__2], q__1.i = 0.f;
		    a[i__4].r = q__1.r, a[i__4].i = q__1.i;
/* L170: */
		}

		i__1 = uub;
		for (k = 1; k <= i__1; ++k) {
		    i__4 = *n - 1;
		    for (jc = 1; jc <= i__4; ++jc) {
/* Computing MAX */
			i__2 = 1, i__3 = jc - k;
			irow = max(i__2,i__3);
/* Computing MIN */
			i__2 = jc + 1, i__3 = k + 2;
			il = min(i__2,i__3);
			extra.r = 0.f, extra.i = 0.f;
			i__2 = jc - iskew * (jc + 1) + ioffg + (jc + 1) * 
				a_dim1;
			ctemp.r = a[i__2].r, ctemp.i = a[i__2].i;
			angle = slarnd_(&c__1, &iseed[1]) * 
				6.2831853071795864769252867663f;
			d__1 = cos(angle);
			clarnd_(&q__2, &c__5, &iseed[1]);
			q__1.r = d__1 * q__2.r, q__1.i = d__1 * q__2.i;
			c.r = q__1.r, c.i = q__1.i;
			d__1 = sin(angle);
			clarnd_(&q__2, &c__5, &iseed[1]);
			q__1.r = d__1 * q__2.r, q__1.i = d__1 * q__2.i;
			s.r = q__1.r, s.i = q__1.i;
			if (csym) {
			    ct.r = c.r, ct.i = c.i;
			    st.r = s.r, st.i = s.i;
			} else {
			    r_cnjg(&q__1, &ctemp);
			    ctemp.r = q__1.r, ctemp.i = q__1.i;
			    r_cnjg(&q__1, &c);
			    ct.r = q__1.r, ct.i = q__1.i;
			    r_cnjg(&q__1, &s);
			    st.r = q__1.r, st.i = q__1.i;
			}
			L__1 = jc > k;
			clarot_(&c_false, &L__1, &c_true, &il, &c, &s, &a[
				irow - iskew * jc + ioffg + jc * a_dim1], &
				ilda, &extra, &ctemp);
/* Computing MIN */
			i__3 = k, i__5 = *n - jc;
			i__2 = min(i__3,i__5) + 1;
			clarot_(&c_true, &c_true, &c_false, &i__2, &ct, &st, &
				a[(1 - iskew) * jc + ioffg + jc * a_dim1], &
				ilda, &ctemp, &dummy);

/*                    Chase EXTRA back up the matrix 
*/

			icol = jc;
			i__2 = -k;
			for (jch = jc - k; i__2 < 0 ? jch >= 1 : jch <= 1; 
				jch += i__2) {
			    clartg_(&a[jch + 1 - iskew * (icol + 1) + ioffg + 
				    (icol + 1) * a_dim1], &extra, &realc, &s, 
				    &dummy);
			    clarnd_(&q__1, &c__5, &iseed[1]);
			    dummy.r = q__1.r, dummy.i = q__1.i;
			    q__2.r = realc * dummy.r, q__2.i = realc * 
				    dummy.i;
			    r_cnjg(&q__1, &q__2);
			    c.r = q__1.r, c.i = q__1.i;
			    q__3.r = -(doublereal)s.r, q__3.i = -(doublereal)
				    s.i;
			    q__2.r = q__3.r * dummy.r - q__3.i * dummy.i, 
				    q__2.i = q__3.r * dummy.i + q__3.i * 
				    dummy.r;
			    r_cnjg(&q__1, &q__2);
			    s.r = q__1.r, s.i = q__1.i;
			    i__3 = jch - iskew * (jch + 1) + ioffg + (jch + 1)
				     * a_dim1;
			    ctemp.r = a[i__3].r, ctemp.i = a[i__3].i;
			    if (csym) {
				ct.r = c.r, ct.i = c.i;
				st.r = s.r, st.i = s.i;
			    } else {
				r_cnjg(&q__1, &ctemp);
				ctemp.r = q__1.r, ctemp.i = q__1.i;
				r_cnjg(&q__1, &c);
				ct.r = q__1.r, ct.i = q__1.i;
				r_cnjg(&q__1, &s);
				st.r = q__1.r, st.i = q__1.i;
			    }
			    i__3 = k + 2;
			    clarot_(&c_true, &c_true, &c_true, &i__3, &c, &s, 
				    &a[(1 - iskew) * jch + ioffg + jch * 
				    a_dim1], &ilda, &ctemp, &extra);
/* Computing MAX */
			    i__3 = 1, i__5 = jch - k;
			    irow = max(i__3,i__5);
/* Computing MIN */
			    i__3 = jch + 1, i__5 = k + 2;
			    il = min(i__3,i__5);
			    extra.r = 0.f, extra.i = 0.f;
			    L__1 = jch > k;
			    clarot_(&c_false, &L__1, &c_true, &il, &ct, &st, &
				    a[irow - iskew * jch + ioffg + jch * 
				    a_dim1], &ilda, &extra, &ctemp);
			    icol = jch;
/* L180: */
			}
/* L190: */
		    }
/* L200: */
		}

/*              If we need lower triangle, copy from upper. No
te that   
                the order of copying is chosen to work for 'q'
 -> 'b' */

		if (ipack != ipackg && ipack != 3) {
		    i__1 = *n;
		    for (jc = 1; jc <= i__1; ++jc) {
			irow = ioffst - iskew * jc;
			if (csym) {
/* Computing MIN */
			    i__2 = *n, i__3 = jc + uub;
			    i__4 = min(i__2,i__3);
			    for (jr = jc; jr <= i__4; ++jr) {
				i__2 = jr + irow + jc * a_dim1;
				i__3 = jc - iskew * jr + ioffg + jr * a_dim1;
				a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
/* L210: */
			    }
			} else {
/* Computing MIN */
			    i__2 = *n, i__3 = jc + uub;
			    i__4 = min(i__2,i__3);
			    for (jr = jc; jr <= i__4; ++jr) {
				i__2 = jr + irow + jc * a_dim1;
				r_cnjg(&q__1, &a[jc - iskew * jr + ioffg + jr 
					* a_dim1]);
				a[i__2].r = q__1.r, a[i__2].i = q__1.i;
/* L220: */
			    }
			}
/* L230: */
		    }
		    if (ipack == 5) {
			i__1 = *n;
			for (jc = *n - uub + 1; jc <= i__1; ++jc) {
			    i__4 = uub + 1;
			    for (jr = *n + 2 - jc; jr <= i__4; ++jr) {
				i__2 = jr + jc * a_dim1;
				a[i__2].r = 0.f, a[i__2].i = 0.f;
/* L240: */
			    }
/* L250: */
			}
		    }
		    if (ipackg == 6) {
			ipackg = ipack;
		    } else {
			ipackg = 0;
		    }
		}
	    } else {

/*              Bottom-Up -- Generate Lower triangle only */

		if (ipack >= 5) {
		    ipackg = 5;
		    if (ipack == 6) {
			ioffg = 1;
		    }
		} else {
		    ipackg = 2;
		}

		i__1 = mnmin;
		for (j = 1; j <= i__1; ++j) {
		    i__4 = (1 - iskew) * j + ioffg + j * a_dim1;
		    i__2 = j;
		    q__1.r = d[i__2], q__1.i = 0.f;
		    a[i__4].r = q__1.r, a[i__4].i = q__1.i;
/* L260: */
		}

		i__1 = uub;
		for (k = 1; k <= i__1; ++k) {
		    for (jc = *n - 1; jc >= 1; --jc) {
/* Computing MIN */
			i__4 = *n + 1 - jc, i__2 = k + 2;
			il = min(i__4,i__2);
			extra.r = 0.f, extra.i = 0.f;
			i__4 = (1 - iskew) * jc + 1 + ioffg + jc * a_dim1;
			ctemp.r = a[i__4].r, ctemp.i = a[i__4].i;
			angle = slarnd_(&c__1, &iseed[1]) * 
				6.2831853071795864769252867663f;
			d__1 = cos(angle);
			clarnd_(&q__2, &c__5, &iseed[1]);
			q__1.r = d__1 * q__2.r, q__1.i = d__1 * q__2.i;
			c.r = q__1.r, c.i = q__1.i;
			d__1 = sin(angle);
			clarnd_(&q__2, &c__5, &iseed[1]);
			q__1.r = d__1 * q__2.r, q__1.i = d__1 * q__2.i;
			s.r = q__1.r, s.i = q__1.i;
			if (csym) {
			    ct.r = c.r, ct.i = c.i;
			    st.r = s.r, st.i = s.i;
			} else {
			    r_cnjg(&q__1, &ctemp);
			    ctemp.r = q__1.r, ctemp.i = q__1.i;
			    r_cnjg(&q__1, &c);
			    ct.r = q__1.r, ct.i = q__1.i;
			    r_cnjg(&q__1, &s);
			    st.r = q__1.r, st.i = q__1.i;
			}
			L__1 = *n - jc > k;
			clarot_(&c_false, &c_true, &L__1, &il, &c, &s, &a[(1 
				- iskew) * jc + ioffg + jc * a_dim1], &ilda, &
				ctemp, &extra);
/* Computing MAX */
			i__4 = 1, i__2 = jc - k + 1;
			icol = max(i__4,i__2);
			i__4 = jc + 2 - icol;
			clarot_(&c_true, &c_false, &c_true, &i__4, &ct, &st, &
				a[jc - iskew * icol + ioffg + icol * a_dim1], 
				&ilda, &dummy, &ctemp);

/*                    Chase EXTRA back down the matrix
 */

			icol = jc;
			i__4 = *n - 1;
			i__2 = k;
			for (jch = jc + k; i__2 < 0 ? jch >= i__4 : jch <= 
				i__4; jch += i__2) {
			    clartg_(&a[jch - iskew * icol + ioffg + icol * 
				    a_dim1], &extra, &realc, &s, &dummy);
			    clarnd_(&q__1, &c__5, &iseed[1]);
			    dummy.r = q__1.r, dummy.i = q__1.i;
			    q__1.r = realc * dummy.r, q__1.i = realc * 
				    dummy.i;
			    c.r = q__1.r, c.i = q__1.i;
			    q__1.r = s.r * dummy.r - s.i * dummy.i, q__1.i = 
				    s.r * dummy.i + s.i * dummy.r;
			    s.r = q__1.r, s.i = q__1.i;
			    i__3 = (1 - iskew) * jch + 1 + ioffg + jch * 
				    a_dim1;
			    ctemp.r = a[i__3].r, ctemp.i = a[i__3].i;
			    if (csym) {
				ct.r = c.r, ct.i = c.i;
				st.r = s.r, st.i = s.i;
			    } else {
				r_cnjg(&q__1, &ctemp);
				ctemp.r = q__1.r, ctemp.i = q__1.i;
				r_cnjg(&q__1, &c);
				ct.r = q__1.r, ct.i = q__1.i;
				r_cnjg(&q__1, &s);
				st.r = q__1.r, st.i = q__1.i;
			    }
			    i__3 = k + 2;
			    clarot_(&c_true, &c_true, &c_true, &i__3, &c, &s, 
				    &a[jch - iskew * icol + ioffg + icol * 
				    a_dim1], &ilda, &extra, &ctemp);
/* Computing MIN */
			    i__3 = *n + 1 - jch, i__5 = k + 2;
			    il = min(i__3,i__5);
			    extra.r = 0.f, extra.i = 0.f;
			    L__1 = *n - jch > k;
			    clarot_(&c_false, &c_true, &L__1, &il, &ct, &st, &
				    a[(1 - iskew) * jch + ioffg + jch * 
				    a_dim1], &ilda, &ctemp, &extra);
			    icol = jch;
/* L270: */
			}
/* L280: */
		    }
/* L290: */
		}

/*              If we need upper triangle, copy from lower. No
te that   
                the order of copying is chosen to work for 'b'
 -> 'q' */

		if (ipack != ipackg && ipack != 4) {
		    for (jc = *n; jc >= 1; --jc) {
			irow = ioffst - iskew * jc;
			if (csym) {
/* Computing MAX */
			    i__2 = 1, i__4 = jc - uub;
			    i__1 = max(i__2,i__4);
			    for (jr = jc; jr >= i__1; --jr) {
				i__2 = jr + irow + jc * a_dim1;
				i__4 = jc - iskew * jr + ioffg + jr * a_dim1;
				a[i__2].r = a[i__4].r, a[i__2].i = a[i__4].i;
/* L300: */
			    }
			} else {
/* Computing MAX */
			    i__2 = 1, i__4 = jc - uub;
			    i__1 = max(i__2,i__4);
			    for (jr = jc; jr >= i__1; --jr) {
				i__2 = jr + irow + jc * a_dim1;
				r_cnjg(&q__1, &a[jc - iskew * jr + ioffg + jr 
					* a_dim1]);
				a[i__2].r = q__1.r, a[i__2].i = q__1.i;
/* L310: */
			    }
			}
/* L320: */
		    }
		    if (ipack == 6) {
			i__1 = uub;
			for (jc = 1; jc <= i__1; ++jc) {
			    i__2 = uub + 1 - jc;
			    for (jr = 1; jr <= i__2; ++jr) {
				i__4 = jr + jc * a_dim1;
				a[i__4].r = 0.f, a[i__4].i = 0.f;
/* L330: */
			    }
/* L340: */
			}
		    }
		    if (ipackg == 5) {
			ipackg = ipack;
		    } else {
			ipackg = 0;
		    }
		}
	    }

/*           Ensure that the diagonal is real if Hermitian */

	    if (! csym) {
		i__1 = *n;
		for (jc = 1; jc <= i__1; ++jc) {
		    irow = ioffst + (1 - iskew) * jc;
		    i__2 = irow + jc * a_dim1;
		    i__4 = irow + jc * a_dim1;
		    d__1 = a[i__4].r;
		    q__1.r = d__1, q__1.i = 0.f;
		    a[i__2].r = q__1.r, a[i__2].i = q__1.i;
/* L350: */
		}
	    }

	}

    } else {

/*        4)      Generate Banded Matrix by first   
                  Rotating by random Unitary matrices,   
                  then reducing the bandwidth using Householder   
                  transformations.   

                  Note: we should get here only if LDA .ge. N */

	if (isym == 1) {

/*           Non-symmetric -- A = U D V */

	    clagge_(&mr, &nc, &llb, &uub, &d[1], &a[a_offset], lda, &iseed[1],
		     &work[1], &iinfo);
	} else {

/*           Symmetric -- A = U D U' or   
             Hermitian -- A = U D U* */

	    if (csym) {
		clagsy_(m, &llb, &d[1], &a[a_offset], lda, &iseed[1], &work[1]
			, &iinfo);
	    } else {
		claghe_(m, &llb, &d[1], &a[a_offset], lda, &iseed[1], &work[1]
			, &iinfo);
	    }
	}

	if (iinfo != 0) {
	    *info = 3;
	    return 0;
	}
    }

/*     5)      Pack the matrix */

    if (ipack != ipackg) {
	if (ipack == 1) {

/*           'U' -- Upper triangular, not packed */

	    i__1 = *m;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i = j + 1; i <= i__2; ++i) {
		    i__4 = i + j * a_dim1;
		    a[i__4].r = 0.f, a[i__4].i = 0.f;
/* L360: */
		}
/* L370: */
	    }

	} else if (ipack == 2) {

/*           'L' -- Lower triangular, not packed */

	    i__1 = *m;
	    for (j = 2; j <= i__1; ++j) {
		i__2 = j - 1;
		for (i = 1; i <= i__2; ++i) {
		    i__4 = i + j * a_dim1;
		    a[i__4].r = 0.f, a[i__4].i = 0.f;
/* L380: */
		}
/* L390: */
	    }

	} else if (ipack == 3) {

/*           'C' -- Upper triangle packed Columnwise. */

	    icol = 1;
	    irow = 0;
	    i__1 = *m;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = j;
		for (i = 1; i <= i__2; ++i) {
		    ++irow;
		    if (irow > *lda) {
			irow = 1;
			++icol;
		    }
		    i__4 = irow + icol * a_dim1;
		    i__3 = i + j * a_dim1;
		    a[i__4].r = a[i__3].r, a[i__4].i = a[i__3].i;
/* L400: */
		}
/* L410: */
	    }

	} else if (ipack == 4) {

/*           'R' -- Lower triangle packed Columnwise. */

	    icol = 1;
	    irow = 0;
	    i__1 = *m;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i = j; i <= i__2; ++i) {
		    ++irow;
		    if (irow > *lda) {
			irow = 1;
			++icol;
		    }
		    i__4 = irow + icol * a_dim1;
		    i__3 = i + j * a_dim1;
		    a[i__4].r = a[i__3].r, a[i__4].i = a[i__3].i;
/* L420: */
		}
/* L430: */
	    }

	} else if (ipack >= 5) {

/*           'B' -- The lower triangle is packed as a band matrix.
   
             'Q' -- The upper triangle is packed as a band matrix.
   
             'Z' -- The whole matrix is packed as a band matrix. 
*/

	    if (ipack == 5) {
		uub = 0;
	    }
	    if (ipack == 6) {
		llb = 0;
	    }

	    i__1 = uub;
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
		i__2 = j + llb;
		for (i = min(i__2,*m); i >= 1; --i) {
		    i__2 = i - j + uub + 1 + j * a_dim1;
		    i__4 = i + j * a_dim1;
		    a[i__2].r = a[i__4].r, a[i__2].i = a[i__4].i;
/* L440: */
		}
/* L450: */
	    }

	    i__1 = *n;
	    for (j = uub + 2; j <= i__1; ++j) {
/* Computing MIN */
		i__4 = j + llb;
		i__2 = min(i__4,*m);
		for (i = j - uub; i <= i__2; ++i) {
		    i__4 = i - j + uub + 1 + j * a_dim1;
		    i__3 = i + j * a_dim1;
		    a[i__4].r = a[i__3].r, a[i__4].i = a[i__3].i;
/* L460: */
		}
/* L470: */
	    }
	}

/*        If packed, zero out extraneous elements.   

          Symmetric/Triangular Packed --   
          zero out everything after A(IROW,ICOL) */

	if (ipack == 3 || ipack == 4) {
	    i__1 = *m;
	    for (jc = icol; jc <= i__1; ++jc) {
		i__2 = *lda;
		for (jr = irow + 1; jr <= i__2; ++jr) {
		    i__4 = jr + jc * a_dim1;
		    a[i__4].r = 0.f, a[i__4].i = 0.f;
/* L480: */
		}
		irow = 0;
/* L490: */
	    }

	} else if (ipack >= 5) {

/*           Packed Band --   
                1st row is now in A( UUB+2-j, j), zero above it   
                m-th row is now in A( M+UUB-j,j), zero below it   
                last non-zero diagonal is now in A( UUB+LLB+1,j ),
   
                   zero below it, too. */

	    ir1 = uub + llb + 2;
	    ir2 = uub + *m + 2;
	    i__1 = *n;
	    for (jc = 1; jc <= i__1; ++jc) {
		i__2 = uub + 1 - jc;
		for (jr = 1; jr <= i__2; ++jr) {
		    i__4 = jr + jc * a_dim1;
		    a[i__4].r = 0.f, a[i__4].i = 0.f;
/* L500: */
		}
/* Computing MAX   
   Computing MIN */
		i__3 = ir1, i__5 = ir2 - jc;
		i__2 = 1, i__4 = min(i__3,i__5);
		i__6 = *lda;
		for (jr = max(i__2,i__4); jr <= i__6; ++jr) {
		    i__2 = jr + jc * a_dim1;
		    a[i__2].r = 0.f, a[i__2].i = 0.f;
/* L510: */
		}
/* L520: */
	    }
	}
    }

    return 0;

/*     End of CLATMS */

} /* clatms_ */

