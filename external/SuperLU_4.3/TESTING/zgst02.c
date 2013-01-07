
/*
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 */
#include "slu_zdefs.h"

int zgst02(trans_t trans, int m, int n, int nrhs, SuperMatrix *A,
	      doublecomplex *x, int ldx, doublecomplex *b, int ldb, double *resid)
{
/*  
    Purpose   
    =======   

    ZGST02 computes the residual for a solution of a system of linear   
    equations  A*x = b  or  A'*x = b:   
       RESID = norm(B - A*X) / ( norm(A) * norm(X) * EPS ),   
    where EPS is the machine epsilon.   

    Arguments   
    =========   

    TRANS   (input) trans_t
            Specifies the form of the system of equations:   
            = NOTRANS:  A *x = b   
            = TRANS  :  A'*x = b, where A' is the transpose of A   
            = CONJ   :  A'*x = b, where A' is the transpose of A   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of columns of B, the matrix of right hand sides.   
            NRHS >= 0.
	    
    A       (input) SuperMatrix*, dimension (LDA,N)   
            The original M x N sparse matrix A.   

    X       (input) DOUBLE COMPLEX PRECISION array, dimension (LDX,NRHS)   
            The computed solution vectors for the system of linear   
            equations.   

    LDX     (input) INTEGER   
            The leading dimension of the array X.  If TRANS = NOTRANS,   
            LDX >= max(1,N); if TRANS = TRANS or CONJ, LDX >= max(1,M).   

    B       (input/output) DOUBLE COMPLEX PRECISION array, dimension (LDB,NRHS)   
            On entry, the right hand side vectors for the system of   
            linear equations.   
            On exit, B is overwritten with the difference B - A*X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  IF TRANS = NOTRANS,
            LDB >= max(1,M); if TRANS = TRANS or CONJ, LDB >= max(1,N).
	    
    RESID   (output) DOUBLE PRECISION   
            The maximum over the number of right hand sides of   
            norm(B - A*X) / ( norm(A) * norm(X) * EPS ).   

    =====================================================================
*/

    /* Table of constant values */
    doublecomplex alpha = {-1., 0.0};
    doublecomplex beta  = {1., 0.0};
    int    c__1  = 1;
    
    /* System generated locals */
    double d__1, d__2;

    /* Local variables */
    int j;
    int n1, n2;
    double anorm, bnorm;
    double xnorm;
    double eps;
    char transc[1];

    /* Function prototypes */
    extern int lsame_(char *, char *);
    extern double zlangs(char *, SuperMatrix *);
    extern double dzasum_(int *, doublecomplex *, int *);
    
    /* Function Body */
    if ( m <= 0 || n <= 0 || nrhs == 0) {
	*resid = 0.;
	return 0;
    }

    if ( (trans == TRANS) || (trans == CONJ) ) {
	n1 = n;
	n2 = m;
        *transc = 'T';
    } else {
	n1 = m;
	n2 = n;
	*transc = 'N';
    }

    /* Exit with RESID = 1/EPS if ANORM = 0. */

    eps = dlamch_("Epsilon");
    anorm = zlangs("1", A);
    if (anorm <= 0.) {
	*resid = 1. / eps;
	return 0;
    }

    /* Compute  B - A*X  (or  B - A'*X ) and store in B. */

    sp_zgemm(transc, "N", n1, nrhs, n2, alpha, A, x, ldx, beta, b, ldb);

    /* Compute the maximum over the number of right hand sides of   
       norm(B - A*X) / ( norm(A) * norm(X) * EPS ) . */

    *resid = 0.;
    for (j = 0; j < nrhs; ++j) {
	bnorm = dzasum_(&n1, &b[j*ldb], &c__1);
	xnorm = dzasum_(&n2, &x[j*ldx], &c__1);
	if (xnorm <= 0.) {
	    *resid = 1. / eps;
	} else {
	    /* Computing MAX */
	    d__1 = *resid, d__2 = bnorm / anorm / xnorm / eps;
	    *resid = SUPERLU_MAX(d__1, d__2);
	}
    }

    return 0;

} /* zgst02 */

