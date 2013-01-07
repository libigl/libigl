/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__4 = 4;
static integer c__8 = 8;

/* Subroutine */ int clarot_(logical *lrows, logical *lleft, logical *lright, 
	integer *nl, complex *c, complex *s, complex *a, integer *lda, 
	complex *xleft, complex *xright)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    complex q__1, q__2, q__3, q__4, q__5, q__6;

    /* Builtin functions */
    void r_cnjg(complex *, complex *);

    /* Local variables */
    static integer iinc, j, inext;
    static complex tempx;
    static integer ix, iy, nt;
    static complex xt[2], yt[2];
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static integer iyt;


/*  -- LAPACK auxiliary test routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

       CLAROT applies a (Givens) rotation to two adjacent rows or   
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

       Those columns which have two '*'s can be handled by SROT.   
       Those columns which have no '*'s can be ignored, since as long   
       as the Givens rotations are carefully applied to preserve   
       symmetry, their values are determined.   
       Those columns which have one '*' have to be handled separately,   
       by using separate variables "p" and "q":   

       row j:     *  *  *  *  *  p  .  .  .   
       row j+1:   q  *  *  *  *  *  .  .  .  .   

       The element p would have to be set correctly, then that column   
       is rotated, setting p to its new value.  The next call to   
       CLAROT would rotate columns j and j+1, using p, and restore   
       symmetry.  The element q would start out being zero, and be   
       made non-zero by the rotation.  Later, rotations would presumably 
  
       be chosen to zero q out.   

       Typical Calling Sequences: rotating the i-th and (i+1)-st rows.   
       ------- ------- ---------   

         General dense matrix:   

                 CALL CLAROT(.TRUE.,.FALSE.,.FALSE., N, C,S,   
                         A(i,1),LDA, DUMMY, DUMMY)   

         General banded matrix in GB format:   

                 j = MAX(1, i-KL )   
                 NL = MIN( N, i+KU+1 ) + 1-j   
                 CALL CLAROT( .TRUE., i-KL.GE.1, i+KU.LT.N, NL, C,S,   
                         A(KU+i+1-j,j),LDA-1, XLEFT, XRIGHT )   

                 [ note that i+1-j is just MIN(i,KL+1) ]   

         Symmetric banded matrix in SY format, bandwidth K,   
         lower triangle only:   

                 j = MAX(1, i-K )   
                 NL = MIN( K+1, i ) + 1   
                 CALL CLAROT( .TRUE., i-K.GE.1, .TRUE., NL, C,S,   
                         A(i,j), LDA, XLEFT, XRIGHT )   

         Same, but upper triangle only:   

                 NL = MIN( K+1, N-i ) + 1   
                 CALL CLAROT( .TRUE., .TRUE., i+K.LT.N, NL, C,S,   
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
                 CALL CLAROT( .TRUE., i-KU.GE.1, i+KL.LT.N, NL, C,S,   
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
             If .TRUE., then CLAROT will rotate two rows.  If .FALSE.,   
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

    C, S   - COMPLEX   
             Specify the Givens rotation to be applied.  If LROWS is   
             true, then the matrix ( c  s )   
                                   ( _  _ )   
                                   (-s  c )  is applied from the left;   
             if false, then the transpose (not conjugated) thereof is   
             applied from the right.  Note that in contrast to the   
             output of CROTG or to most versions of CROT, both C and S   
             are complex.  For a Givens rotation, |C|**2 + |S|**2 should 
  
             be 1, but this is not checked.   
             Not modified.   

    A      - COMPLEX array.   
             The array containing the rows/columns to be rotated.  The   
             first element of A should be the upper left element to   
             be rotated.   
             Read and modified.   

    LDA    - INTEGER   
             The "effective" leading dimension of A.  If A contains   
             a matrix stored in GE, HE, or SY format, then this is just   
             the leading dimension of A as dimensioned in the calling   
             routine.  If A contains a matrix stored in band (GB, HB, or 
  
             SB) format, then this should be *one less* than the leading 
  
             dimension used in the calling routine.  Thus, if A were   
             dimensioned A(LDA,*) in CLAROT, then A(1,j) would be the   
             j-th element in the first of the two rows to be rotated,   
             and A(2,j) would be the j-th in the second, regardless of   
             how the array may be stored in the calling routine.  [A   
             cannot, however, actually be dimensioned thus, since for   
             band format, the row number may exceed LDA, which is not   
             legal FORTRAN.]   
             If LROWS=.TRUE., then LDA must be at least 1, otherwise   
             it must be at least NL minus the number of .TRUE. values   
             in XLEFT and XRIGHT.   
             Not modified.   

    XLEFT  - COMPLEX   
             If LLEFT is .TRUE., then XLEFT will be used and modified   
             instead of A(2,1) (if LROWS=.TRUE.) or A(1,2)   
             (if LROWS=.FALSE.).   
             Read and modified.   

    XRIGHT - COMPLEX   
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
	xt[0].r = a[1].r, xt[0].i = a[1].i;
	yt[0].r = xleft->r, yt[0].i = xleft->i;
    } else {
	nt = 0;
	ix = 1;
	iy = inext + 1;
    }

    if (*lright) {
	iyt = inext + 1 + (*nl - 1) * iinc;
	++nt;
	i__1 = nt - 1;
	xt[i__1].r = xright->r, xt[i__1].i = xright->i;
	i__1 = nt - 1;
	i__2 = iyt;
	yt[i__1].r = a[i__2].r, yt[i__1].i = a[i__2].i;
    }

/*     Check for errors */

    if (*nl < nt) {
	xerbla_("CLAROT", &c__4);
	return 0;
    }
    if (*lda <= 0 || ! (*lrows) && *lda < *nl - nt) {
	xerbla_("CLAROT", &c__8);
	return 0;
    }

/*     Rotate   

       CROT( NL-NT, A(IX),IINC, A(IY),IINC, C, S ) with complex C, S */

    i__1 = *nl - nt - 1;
    for (j = 0; j <= i__1; ++j) {
	i__2 = ix + j * iinc;
	q__2.r = c->r * a[i__2].r - c->i * a[i__2].i, q__2.i = c->r * a[i__2]
		.i + c->i * a[i__2].r;
	i__3 = iy + j * iinc;
	q__3.r = s->r * a[i__3].r - s->i * a[i__3].i, q__3.i = s->r * a[i__3]
		.i + s->i * a[i__3].r;
	q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
	tempx.r = q__1.r, tempx.i = q__1.i;
	i__2 = iy + j * iinc;
	r_cnjg(&q__4, s);
	q__3.r = -(doublereal)q__4.r, q__3.i = -(doublereal)q__4.i;
	i__3 = ix + j * iinc;
	q__2.r = q__3.r * a[i__3].r - q__3.i * a[i__3].i, q__2.i = q__3.r * a[
		i__3].i + q__3.i * a[i__3].r;
	r_cnjg(&q__6, c);
	i__4 = iy + j * iinc;
	q__5.r = q__6.r * a[i__4].r - q__6.i * a[i__4].i, q__5.i = q__6.r * a[
		i__4].i + q__6.i * a[i__4].r;
	q__1.r = q__2.r + q__5.r, q__1.i = q__2.i + q__5.i;
	a[i__2].r = q__1.r, a[i__2].i = q__1.i;
	i__2 = ix + j * iinc;
	a[i__2].r = tempx.r, a[i__2].i = tempx.i;
/* L10: */
    }

/*     CROT( NT, XT,1, YT,1, C, S ) with complex C, S */

    i__1 = nt;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j - 1;
	q__2.r = c->r * xt[i__2].r - c->i * xt[i__2].i, q__2.i = c->r * xt[
		i__2].i + c->i * xt[i__2].r;
	i__3 = j - 1;
	q__3.r = s->r * yt[i__3].r - s->i * yt[i__3].i, q__3.i = s->r * yt[
		i__3].i + s->i * yt[i__3].r;
	q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
	tempx.r = q__1.r, tempx.i = q__1.i;
	i__2 = j - 1;
	r_cnjg(&q__4, s);
	q__3.r = -(doublereal)q__4.r, q__3.i = -(doublereal)q__4.i;
	i__3 = j - 1;
	q__2.r = q__3.r * xt[i__3].r - q__3.i * xt[i__3].i, q__2.i = q__3.r * 
		xt[i__3].i + q__3.i * xt[i__3].r;
	r_cnjg(&q__6, c);
	i__4 = j - 1;
	q__5.r = q__6.r * yt[i__4].r - q__6.i * yt[i__4].i, q__5.i = q__6.r * 
		yt[i__4].i + q__6.i * yt[i__4].r;
	q__1.r = q__2.r + q__5.r, q__1.i = q__2.i + q__5.i;
	yt[i__2].r = q__1.r, yt[i__2].i = q__1.i;
	i__2 = j - 1;
	xt[i__2].r = tempx.r, xt[i__2].i = tempx.i;
/* L20: */
    }

/*     Stuff values back into XLEFT, XRIGHT, etc. */

    if (*lleft) {
	a[1].r = xt[0].r, a[1].i = xt[0].i;
	xleft->r = yt[0].r, xleft->i = yt[0].i;
    }

    if (*lright) {
	i__1 = nt - 1;
	xright->r = xt[i__1].r, xright->i = xt[i__1].i;
	i__1 = iyt;
	i__2 = nt - 1;
	a[i__1].r = yt[i__2].r, a[i__1].i = yt[i__2].i;
    }

    return 0;

/*     End of CLAROT */

} /* clarot_ */

