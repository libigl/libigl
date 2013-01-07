/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Double Complex */ VOID zlatm3_(doublecomplex * ret_val, integer *m, 
	integer *n, integer *i, integer *j, integer *isub, integer *jsub, 
	integer *kl, integer *ku, integer *idist, integer *iseed, 
	doublecomplex *d, integer *igrade, doublecomplex *dl, doublecomplex *
	dr, integer *ipvtng, integer *iwork, doublereal *sparse)
{
    /* System generated locals */
    integer i__1, i__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static doublecomplex ctemp;
    extern doublereal dlaran_(integer *);
    extern /* Double Complex */ VOID zlarnd_(doublecomplex *, integer *, 
	    integer *);


/*  -- LAPACK auxiliary test routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   





    Purpose   
    =======   

       ZLATM3 returns the (ISUB,JSUB) entry of a random matrix of   
       dimension (M, N) described by the other paramters. (ISUB,JSUB)   
       is the final position of the (I,J) entry after pivoting   
       according to IPVTNG and IWORK. ZLATM3 is called by the   
       ZLATMR routine in order to build random test matrices. No error   
       checking on parameters is done, because this routine is called in 
  
       a tight loop by ZLATMR which has already checked the parameters.   

       Use of ZLATM3 differs from CLATM2 in the order in which the random 
  
       number generator is called to fill in random matrix entries.   
       With ZLATM2, the generator is called to fill in the pivoted matrix 
  
       columnwise. With ZLATM3, the generator is called to fill in the   
       matrix columnwise, after which it is pivoted. Thus, ZLATM3 can   
       be used to construct random matrices which differ only in their   
       order of rows and/or columns. ZLATM2 is used to construct band   
       matrices while avoiding calling the random number generator for   
       entries outside the band (and therefore generating random numbers 
  
       in different orders for different pivot orders).   

       The matrix whose (ISUB,JSUB) entry is returned is constructed as   
       follows (this routine only computes one entry):   

         If ISUB is outside (1..M) or JSUB is outside (1..N), return zero 
  
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
             Row of unpivoted entry to be returned. Not modified.   

    J      - INTEGER   
             Column of unpivoted entry to be returned. Not modified.   

    ISUB   - INTEGER   
             Row of pivoted entry to be returned. Changed on exit.   

    JSUB   - INTEGER   
             Column of pivoted entry to be returned. Changed on exit.   

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

    D      - COMPLEX*16            array of dimension ( MIN( I , J ) )   
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

    DL     - COMPLEX*16            array ( I or J, as appropriate )   
             Left scale factors for grading matrix.  Not modified.   

    DR     - COMPLEX*16            array ( I or J, as appropriate )   
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
             row (or column) originally in position K is in   
             position IWORK( K ) after pivoting.   
             This differs from IWORK for ZLATM2. Not modified.   

    SPARSE - DOUBLE PRECISION               between 0. and 1.   
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
	*isub = *i;
	*jsub = *j;
	 ret_val->r = 0.,  ret_val->i = 0.;
	return ;
    }

/*     Compute subscripts depending on IPVTNG */

    if (*ipvtng == 0) {
	*isub = *i;
	*jsub = *j;
    } else if (*ipvtng == 1) {
	*isub = iwork[*i];
	*jsub = *j;
    } else if (*ipvtng == 2) {
	*isub = *i;
	*jsub = iwork[*j];
    } else if (*ipvtng == 3) {
	*isub = iwork[*i];
	*jsub = iwork[*j];
    }

/*     Check for banding */

    if (*jsub > *isub + *ku || *jsub < *isub - *kl) {
	 ret_val->r = 0.,  ret_val->i = 0.;
	return ;
    }

/*     Check for sparsity */

    if (*sparse > 0.) {
	if (dlaran_(&iseed[1]) < *sparse) {
	     ret_val->r = 0.,  ret_val->i = 0.;
	    return ;
	}
    }

/*     Compute entry and grade it according to IGRADE */

    if (*i == *j) {
	i__1 = *i;
	ctemp.r = d[i__1].r, ctemp.i = d[i__1].i;
    } else {
	zlarnd_(&z__1, idist, &iseed[1]);
	ctemp.r = z__1.r, ctemp.i = z__1.i;
    }
    if (*igrade == 1) {
	i__1 = *i;
	z__1.r = ctemp.r * dl[i__1].r - ctemp.i * dl[i__1].i, z__1.i = 
		ctemp.r * dl[i__1].i + ctemp.i * dl[i__1].r;
	ctemp.r = z__1.r, ctemp.i = z__1.i;
    } else if (*igrade == 2) {
	i__1 = *j;
	z__1.r = ctemp.r * dr[i__1].r - ctemp.i * dr[i__1].i, z__1.i = 
		ctemp.r * dr[i__1].i + ctemp.i * dr[i__1].r;
	ctemp.r = z__1.r, ctemp.i = z__1.i;
    } else if (*igrade == 3) {
	i__1 = *i;
	z__2.r = ctemp.r * dl[i__1].r - ctemp.i * dl[i__1].i, z__2.i = 
		ctemp.r * dl[i__1].i + ctemp.i * dl[i__1].r;
	i__2 = *j;
	z__1.r = z__2.r * dr[i__2].r - z__2.i * dr[i__2].i, z__1.i = z__2.r * 
		dr[i__2].i + z__2.i * dr[i__2].r;
	ctemp.r = z__1.r, ctemp.i = z__1.i;
    } else if (*igrade == 4 && *i != *j) {
	i__1 = *i;
	z__2.r = ctemp.r * dl[i__1].r - ctemp.i * dl[i__1].i, z__2.i = 
		ctemp.r * dl[i__1].i + ctemp.i * dl[i__1].r;
	z_div(&z__1, &z__2, &dl[*j]);
	ctemp.r = z__1.r, ctemp.i = z__1.i;
    } else if (*igrade == 5) {
	i__1 = *i;
	z__2.r = ctemp.r * dl[i__1].r - ctemp.i * dl[i__1].i, z__2.i = 
		ctemp.r * dl[i__1].i + ctemp.i * dl[i__1].r;
	d_cnjg(&z__3, &dl[*j]);
	z__1.r = z__2.r * z__3.r - z__2.i * z__3.i, z__1.i = z__2.r * z__3.i 
		+ z__2.i * z__3.r;
	ctemp.r = z__1.r, ctemp.i = z__1.i;
    } else if (*igrade == 6) {
	i__1 = *i;
	z__2.r = ctemp.r * dl[i__1].r - ctemp.i * dl[i__1].i, z__2.i = 
		ctemp.r * dl[i__1].i + ctemp.i * dl[i__1].r;
	i__2 = *j;
	z__1.r = z__2.r * dl[i__2].r - z__2.i * dl[i__2].i, z__1.i = z__2.r * 
		dl[i__2].i + z__2.i * dl[i__2].r;
	ctemp.r = z__1.r, ctemp.i = z__1.i;
    }
     ret_val->r = ctemp.r,  ret_val->i = ctemp.i;
    return ;

/*     End of ZLATM3 */

} /* zlatm3_ */

