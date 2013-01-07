/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Complex */ VOID clatm2_(complex * ret_val, integer *m, integer *n, integer 
	*i, integer *j, integer *kl, integer *ku, integer *idist, integer *
	iseed, complex *d, integer *igrade, complex *dl, complex *dr, integer 
	*ipvtng, integer *iwork, real *sparse)
{
    /* System generated locals */
    integer i__1, i__2;
    complex q__1, q__2, q__3;

    /* Builtin functions */
    void c_div(complex *, complex *, complex *), r_cnjg(complex *, complex *);

    /* Local variables */
    static integer isub, jsub;
    static complex ctemp;
    extern /* Complex */ VOID clarnd_(complex *, integer *, integer *);
    extern doublereal slaran_(integer *);


/*  -- LAPACK auxiliary test routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   





    Purpose   
    =======   

       CLATM2 returns the (I,J) entry of a random matrix of dimension   
       (M, N) described by the other paramters. It is called by the   
       CLATMR routine in order to build random test matrices. No error   
       checking on parameters is done, because this routine is called in 
  
       a tight loop by CLATMR which has already checked the parameters.   

       Use of CLATM2 differs from CLATM3 in the order in which the random 
  
       number generator is called to fill in random matrix entries.   
       With CLATM2, the generator is called to fill in the pivoted matrix 
  
       columnwise. With CLATM3, the generator is called to fill in the   
       matrix columnwise, after which it is pivoted. Thus, CLATM3 can   
       be used to construct random matrices which differ only in their   
       order of rows and/or columns. CLATM2 is used to construct band   
       matrices while avoiding calling the random number generator for   
       entries outside the band (and therefore generating random numbers 
  

       The matrix whose (I,J) entry is returned is constructed as   
       follows (this routine only computes one entry):   

         If I is outside (1..M) or J is outside (1..N), return zero   
            (this is convenient for generating matrices in band format). 
  

         Generate a matrix A with random entries of distribution IDIST.   

         Set the diagonal to D.   

         Grade the matrix, if desired, from the left (by DL) and/or   
            from the right (by DR or DL) as specified by IGRADE.   

         Permute, if desired, the rows and/or columns as specified by   
            IPVTNG and IWORK.   

         Band the matrix to have lower bandwidth KL and upper   
            bandwidth KU.   

         Set random entries to zero as specified by SPARSE.   

    Arguments   
    =========   

    M      - INTEGER   
             Number of rows of matrix. Not modified.   

    N      - INTEGER   
             Number of columns of matrix. Not modified.   

    I      - INTEGER   
             Row of entry to be returned. Not modified.   

    J      - INTEGER   
             Column of entry to be returned. Not modified.   

    KL     - INTEGER   
             Lower bandwidth. Not modified.   

    KU     - INTEGER   
             Upper bandwidth. Not modified.   

    IDIST  - INTEGER   
             On entry, IDIST specifies the type of distribution to be   
             used to generate a random matrix .   
             1 => real and imaginary parts each UNIFORM( 0, 1 )   
             2 => real and imaginary parts each UNIFORM( -1, 1 )   
             3 => real and imaginary parts each NORMAL( 0, 1 )   
             4 => complex number uniform in DISK( 0 , 1 )   
             Not modified.   

    ISEED  - INTEGER            array of dimension ( 4 )   
             Seed for random number generator.   
             Changed on exit.   

    D      - COMPLEX            array of dimension ( MIN( I , J ) )   
             Diagonal entries of matrix. Not modified.   

    IGRADE - INTEGER   
             Specifies grading of matrix as follows:   
             0  => no grading   
             1  => matrix premultiplied by diag( DL )   
             2  => matrix postmultiplied by diag( DR )   
             3  => matrix premultiplied by diag( DL ) and   
                           postmultiplied by diag( DR )   
             4  => matrix premultiplied by diag( DL ) and   
                           postmultiplied by inv( diag( DL ) )   
             5  => matrix premultiplied by diag( DL ) and   
                           postmultiplied by diag( CONJG(DL) )   
             6  => matrix premultiplied by diag( DL ) and   
                           postmultiplied by diag( DL )   
             Not modified.   

    DL     - COMPLEX            array ( I or J, as appropriate )   
             Left scale factors for grading matrix.  Not modified.   

    DR     - COMPLEX            array ( I or J, as appropriate )   
             Right scale factors for grading matrix.  Not modified.   

    IPVTNG - INTEGER   
             On entry specifies pivoting permutations as follows:   
             0 => none.   
             1 => row pivoting.   
             2 => column pivoting.   
             3 => full pivoting, i.e., on both sides.   
             Not modified.   

    IWORK  - INTEGER            array ( I or J, as appropriate )   
             This array specifies the permutation used. The   
             row (or column) in position K was originally in   
             position IWORK( K ).   
             This differs from IWORK for CLATM3. Not modified.   

    SPARSE - REAL               between 0. and 1.   
             On entry specifies the sparsity of the matrix   
             if sparse matix is to be generated.   
             SPARSE should lie between 0 and 1.   
             A uniform ( 0, 1 ) random number x is generated and   
             compared to SPARSE; if x is larger the matrix entry   
             is unchanged and if x is smaller the entry is set   
             to zero. Thus on the average a fraction SPARSE of the   
             entries will be set to zero.   
             Not modified.   

    ===================================================================== 
  









   -----------------------------------------------------------------------
   



       Check for I and J in range   

       Parameter adjustments */
    --iwork;
    --dr;
    --dl;
    --d;
    --iseed;

    /* Function Body */
    if (*i < 1 || *i > *m || *j < 1 || *j > *n) {
	 ret_val->r = 0.f,  ret_val->i = 0.f;
	return ;
    }

/*     Check for banding */

    if (*j > *i + *ku || *j < *i - *kl) {
	 ret_val->r = 0.f,  ret_val->i = 0.f;
	return ;
    }

/*     Check for sparsity */

    if (*sparse > 0.f) {
	if (slaran_(&iseed[1]) < *sparse) {
	     ret_val->r = 0.f,  ret_val->i = 0.f;
	    return ;
	}
    }

/*     Compute subscripts depending on IPVTNG */

    if (*ipvtng == 0) {
	isub = *i;
	jsub = *j;
    } else if (*ipvtng == 1) {
	isub = iwork[*i];
	jsub = *j;
    } else if (*ipvtng == 2) {
	isub = *i;
	jsub = iwork[*j];
    } else if (*ipvtng == 3) {
	isub = iwork[*i];
	jsub = iwork[*j];
    }

/*     Compute entry and grade it according to IGRADE */

    if (isub == jsub) {
	i__1 = isub;
	ctemp.r = d[i__1].r, ctemp.i = d[i__1].i;
    } else {
	clarnd_(&q__1, idist, &iseed[1]);
	ctemp.r = q__1.r, ctemp.i = q__1.i;
    }
    if (*igrade == 1) {
	i__1 = isub;
	q__1.r = ctemp.r * dl[i__1].r - ctemp.i * dl[i__1].i, q__1.i = 
		ctemp.r * dl[i__1].i + ctemp.i * dl[i__1].r;
	ctemp.r = q__1.r, ctemp.i = q__1.i;
    } else if (*igrade == 2) {
	i__1 = jsub;
	q__1.r = ctemp.r * dr[i__1].r - ctemp.i * dr[i__1].i, q__1.i = 
		ctemp.r * dr[i__1].i + ctemp.i * dr[i__1].r;
	ctemp.r = q__1.r, ctemp.i = q__1.i;
    } else if (*igrade == 3) {
	i__1 = isub;
	q__2.r = ctemp.r * dl[i__1].r - ctemp.i * dl[i__1].i, q__2.i = 
		ctemp.r * dl[i__1].i + ctemp.i * dl[i__1].r;
	i__2 = jsub;
	q__1.r = q__2.r * dr[i__2].r - q__2.i * dr[i__2].i, q__1.i = q__2.r * 
		dr[i__2].i + q__2.i * dr[i__2].r;
	ctemp.r = q__1.r, ctemp.i = q__1.i;
    } else if (*igrade == 4 && isub != jsub) {
	i__1 = isub;
	q__2.r = ctemp.r * dl[i__1].r - ctemp.i * dl[i__1].i, q__2.i = 
		ctemp.r * dl[i__1].i + ctemp.i * dl[i__1].r;
	c_div(&q__1, &q__2, &dl[jsub]);
	ctemp.r = q__1.r, ctemp.i = q__1.i;
    } else if (*igrade == 5) {
	i__1 = isub;
	q__2.r = ctemp.r * dl[i__1].r - ctemp.i * dl[i__1].i, q__2.i = 
		ctemp.r * dl[i__1].i + ctemp.i * dl[i__1].r;
	r_cnjg(&q__3, &dl[jsub]);
	q__1.r = q__2.r * q__3.r - q__2.i * q__3.i, q__1.i = q__2.r * q__3.i 
		+ q__2.i * q__3.r;
	ctemp.r = q__1.r, ctemp.i = q__1.i;
    } else if (*igrade == 6) {
	i__1 = isub;
	q__2.r = ctemp.r * dl[i__1].r - ctemp.i * dl[i__1].i, q__2.i = 
		ctemp.r * dl[i__1].i + ctemp.i * dl[i__1].r;
	i__2 = jsub;
	q__1.r = q__2.r * dl[i__2].r - q__2.i * dl[i__2].i, q__1.i = q__2.r * 
		dl[i__2].i + q__2.i * dl[i__2].r;
	ctemp.r = q__1.r, ctemp.i = q__1.i;
    }
     ret_val->r = ctemp.r,  ret_val->i = ctemp.i;
    return ;

/*     End of CLATM2 */

} /* clatm2_ */

