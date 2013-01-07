/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__4 = 4;
static integer c__8 = 8;
static integer c__1 = 1;

/* Subroutine */ int dlarot_(logical *lrows, logical *lleft, logical *lright, 
	integer *nl, doublereal *c, doublereal *s, doublereal *a, integer *
	lda, doublereal *xleft, doublereal *xright)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer iinc;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer inext, ix, iy, nt;
    static doublereal xt[2], yt[2];
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static integer iyt;


/*  -- LAPACK auxiliary test routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

       DLAROT applies a (Givens) rotation to two adjacent rows or   
       columns, where one element of the first and/or last column/row   
       may be a separate variable.  This is specifically indended   
       for use on matrices stored in some format other than GE, so   
       that elements of the matrix may be used or modified for which   
       no array element is provided.   

       One example is a symmetric matrix in SB format (bandwidth=4), for 
  
       which UPLO='L':  Two adjacent rows will have the format:   

       row j:     *  *  *  *  *  .  .  .  .   
       row j+1:      *  *  *  *  *  .  .  .  .   

       '*' indicates elements for which storage is provided,   
       '.' indicates elements for which no storage is provided, but   
       are not necessarily zero; their values are determined by   
       symmetry.  ' ' indicates elements which are necessarily zero,   
        and have no storage provided.   

       Those columns which have two '*'s can be handled by DROT.   
       Those columns which have no '*'s can be ignored, since as long   
       as the Givens rotations are carefully applied to preserve   
       symmetry, their values are determined.   
       Those columns which have one '*' have to be handled separately,   
       by using separate variables "p" and "q":   

       row j:     *  *  *  *  *  p  .  .  .   
       row j+1:   q  *  *  *  *  *  .  .  .  .   

       The element p would have to be set correctly, then that column   
       is rotated, setting p to its new value.  The next call to   
       DLAROT would rotate columns j and j+1, using p, and restore   
       symmetry.  The element q would start out being zero, and be   
       made non-zero by the rotation.  Later, rotations would presumably 
  
       be chosen to zero q out.   

       Typical Calling Sequences: rotating the i-th and (i+1)-st rows.   
       ------- ------- ---------   

         General dense matrix:   

                 CALL DLAROT(.TRUE.,.FALSE.,.FALSE., N, C,S,   
                         A(i,1),LDA, DUMMY, DUMMY)   

         General banded matrix in GB format:   

                 j = MAX(1, i-KL )   
                 NL = MIN( N, i+KU+1 ) + 1-j   
                 CALL DLAROT( .TRUE., i-KL.GE.1, i+KU.LT.N, NL, C,S,   
                         A(KU+i+1-j,j),LDA-1, XLEFT, XRIGHT )   

                 [ note that i+1-j is just MIN(i,KL+1) ]   

         Symmetric banded matrix in SY format, bandwidth K,   
         lower triangle only:   

                 j = MAX(1, i-K )   
                 NL = MIN( K+1, i ) + 1   
                 CALL DLAROT( .TRUE., i-K.GE.1, .TRUE., NL, C,S,   
                         A(i,j), LDA, XLEFT, XRIGHT )   

         Same, but upper triangle only:   

                 NL = MIN( K+1, N-i ) + 1   
                 CALL DLAROT( .TRUE., .TRUE., i+K.LT.N, NL, C,S,   
                         A(i,i), LDA, XLEFT, XRIGHT )   

         Symmetric banded matrix in SB format, bandwidth K,   
         lower triangle only:   

                 [ same as for SY, except:]   
                     . . . .   
                         A(i+1-j,j), LDA-1, XLEFT, XRIGHT )   

                 [ note that i+1-j is just MIN(i,K+1) ]   

         Same, but upper triangle only:   
                      . . .   
                         A(K+1,i), LDA-1, XLEFT, XRIGHT )   

         Rotating columns is just the transpose of rotating rows, except 
  
         for GB and SB: (rotating columns i and i+1)   

         GB:   
                 j = MAX(1, i-KU )   
                 NL = MIN( N, i+KL+1 ) + 1-j   
                 CALL DLAROT( .TRUE., i-KU.GE.1, i+KL.LT.N, NL, C,S,   
                         A(KU+j+1-i,i),LDA-1, XTOP, XBOTTM )   

                 [note that KU+j+1-i is just MAX(1,KU+2-i)]   

         SB: (upper triangle)   

                      . . . . . .   
                         A(K+j+1-i,i),LDA-1, XTOP, XBOTTM )   

         SB: (lower triangle)   

                      . . . . . .   
                         A(1,i),LDA-1, XTOP, XBOTTM )   

    Arguments   
    =========   

    LROWS  - LOGICAL   
             If .TRUE., then DLAROT will rotate two rows.  If .FALSE.,   
             then it will rotate two columns.   
             Not modified.   

    LLEFT  - LOGICAL   
             If .TRUE., then XLEFT will be used instead of the   
             corresponding element of A for the first element in the   
             second row (if LROWS=.FALSE.) or column (if LROWS=.TRUE.)   
             If .FALSE., then the corresponding element of A will be   
             used.   
             Not modified.   

    LRIGHT - LOGICAL   
             If .TRUE., then XRIGHT will be used instead of the   
             corresponding element of A for the last element in the   
             first row (if LROWS=.FALSE.) or column (if LROWS=.TRUE.) If 
  
             .FALSE., then the corresponding element of A will be used.   
             Not modified.   

    NL     - INTEGER   
             The length of the rows (if LROWS=.TRUE.) or columns (if   
             LROWS=.FALSE.) to be rotated.  If XLEFT and/or XRIGHT are   
             used, the columns/rows they are in should be included in   
             NL, e.g., if LLEFT = LRIGHT = .TRUE., then NL must be at   
             least 2.  The number of rows/columns to be rotated   
             exclusive of those involving XLEFT and/or XRIGHT may   
             not be negative, i.e., NL minus how many of LLEFT and   
             LRIGHT are .TRUE. must be at least zero; if not, XERBLA   
             will be called.   
             Not modified.   

    C, S   - DOUBLE PRECISION   
             Specify the Givens rotation to be applied.  If LROWS is   
             true, then the matrix ( c  s )   
                                   (-s  c )  is applied from the left;   
             if false, then the transpose thereof is applied from the   
             right.  For a Givens rotation, C**2 + S**2 should be 1,   
             but this is not checked.   
             Not modified.   

    A      - DOUBLE PRECISION array.   
             The array containing the rows/columns to be rotated.  The   
             first element of A should be the upper left element to   
             be rotated.   
             Read and modified.   

    LDA    - INTEGER   
             The "effective" leading dimension of A.  If A contains   
             a matrix stored in GE or SY format, then this is just   
             the leading dimension of A as dimensioned in the calling   
             routine.  If A contains a matrix stored in band (GB or SB)   
             format, then this should be *one less* than the leading   
             dimension used in the calling routine.  Thus, if   
             A were dimensioned A(LDA,*) in DLAROT, then A(1,j) would   
             be the j-th element in the first of the two rows   
             to be rotated, and A(2,j) would be the j-th in the second,   
             regardless of how the array may be stored in the calling   
             routine.  [A cannot, however, actually be dimensioned thus, 
  
             since for band format, the row number may exceed LDA, which 
  
             is not legal FORTRAN.]   
             If LROWS=.TRUE., then LDA must be at least 1, otherwise   
             it must be at least NL minus the number of .TRUE. values   
             in XLEFT and XRIGHT.   
             Not modified.   

    XLEFT  - DOUBLE PRECISION   
             If LLEFT is .TRUE., then XLEFT will be used and modified   
             instead of A(2,1) (if LROWS=.TRUE.) or A(1,2)   
             (if LROWS=.FALSE.).   
             Read and modified.   

    XRIGHT - DOUBLE PRECISION   
             If LRIGHT is .TRUE., then XRIGHT will be used and modified   
             instead of A(1,NL) (if LROWS=.TRUE.) or A(NL,1)   
             (if LROWS=.FALSE.).   
             Read and modified.   

    ===================================================================== 
  


       Set up indices, arrays for ends   

       Parameter adjustments */
    --a;

    /* Function Body */
    if (*lrows) {
	iinc = *lda;
	inext = 1;
    } else {
	iinc = 1;
	inext = *lda;
    }

    if (*lleft) {
	nt = 1;
	ix = iinc + 1;
	iy = *lda + 2;
	xt[0] = a[1];
	yt[0] = *xleft;
    } else {
	nt = 0;
	ix = 1;
	iy = inext + 1;
    }

    if (*lright) {
	iyt = inext + 1 + (*nl - 1) * iinc;
	++nt;
	xt[nt - 1] = *xright;
	yt[nt - 1] = a[iyt];
    }

/*     Check for errors */

    if (*nl < nt) {
	xerbla_("DLAROT", &c__4);
	return 0;
    }
    if (*lda <= 0 || ! (*lrows) && *lda < *nl - nt) {
	xerbla_("DLAROT", &c__8);
	return 0;
    }

/*     Rotate */

    i__1 = *nl - nt;
    drot_(&i__1, &a[ix], &iinc, &a[iy], &iinc, c, s);
    drot_(&nt, xt, &c__1, yt, &c__1, c, s);

/*     Stuff values back into XLEFT, XRIGHT, etc. */

    if (*lleft) {
	a[1] = xt[0];
	*xleft = yt[0];
    }

    if (*lright) {
	*xright = xt[nt - 1];
	a[iyt] = yt[nt - 1];
    }

    return 0;

/*     End of DLAROT */

} /* dlarot_ */

