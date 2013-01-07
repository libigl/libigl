
/*
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 */
#include <math.h>
#include "slu_zdefs.h"

int zgst07(trans_t trans, int n, int nrhs, SuperMatrix *A, doublecomplex *b, 
	      int ldb, doublecomplex *x, int ldx, doublecomplex *xact, 
              int ldxact, double *ferr, double *berr, double *reslts)
{
/*
    Purpose   
    =======   

    ZGST07 tests the error bounds from iterative refinement for the   
    computed solution to a system of equations op(A)*X = B, where A is a 
    general n by n matrix and op(A) = A or A**T, depending on TRANS.
    
    RESLTS(1) = test of the error bound   
              = norm(X - XACT) / ( norm(X) * FERR )   
    A large value is returned if this ratio is not less than one.   

    RESLTS(2) = residual from the iterative refinement routine   
              = the maximum of BERR / ( (n+1)*EPS + (*) ), where   
                (*) = (n+1)*UNFL / (min_i (abs(op(A))*abs(X) +abs(b))_i ) 

    Arguments   
    =========   

    TRANS   (input) trans_t
            Specifies the form of the system of equations.   
            = NOTRANS:  A *x = b   
            = TRANS  :  A'*x = b, where A' is the transpose of A   
            = CONJ   :  A'*x = b, where A' is the transpose of A   

    N       (input) INT
            The number of rows of the matrices X and XACT.  N >= 0.   

    NRHS    (input) INT   
            The number of columns of the matrices X and XACT.  NRHS >= 0. 
  

    A       (input) SuperMatrix *, dimension (A->nrow, A->ncol)
            The original n by n matrix A.   

    B       (input) DOUBLE COMPLEX PRECISION array, dimension (LDB,NRHS)   
            The right hand side vectors for the system of linear   
            equations.   

    LDB     (input) INT   
            The leading dimension of the array B.  LDB >= max(1,N).   

    X       (input) DOUBLE COMPLEX PRECISION array, dimension (LDX,NRHS)   
            The computed solution vectors.  Each vector is stored as a   
            column of the matrix X.   

    LDX     (input) INT   
            The leading dimension of the array X.  LDX >= max(1,N).   

    XACT    (input) DOUBLE COMPLEX PRECISION array, dimension (LDX,NRHS)   
            The exact solution vectors.  Each vector is stored as a   
            column of the matrix XACT.   

    LDXACT  (input) INT   
            The leading dimension of the array XACT.  LDXACT >= max(1,N). 
  

    FERR    (input) DOUBLE COMPLEX PRECISION array, dimension (NRHS)   
            The estimated forward error bounds for each solution vector   
            X.  If XTRUE is the true solution, FERR bounds the magnitude 
            of the largest entry in (X - XTRUE) divided by the magnitude 
            of the largest entry in X.   

    BERR    (input) DOUBLE COMPLEX PRECISION array, dimension (NRHS)   
            The componentwise relative backward error of each solution   
            vector (i.e., the smallest relative change in any entry of A 
  
            or B that makes X an exact solution).   

    RESLTS  (output) DOUBLE PRECISION array, dimension (2)   
            The maximum over the NRHS solution vectors of the ratios:   
            RESLTS(1) = norm(X - XACT) / ( norm(X) * FERR )   
            RESLTS(2) = BERR / ( (n+1)*EPS + (*) )   

    ===================================================================== 
*/
    
    /* Table of constant values */
    int c__1 = 1;

    /* System generated locals */
    double d__1, d__2;
    double d__3, d__4;

    /* Local variables */
    double diff, axbi;
    int    imax, irow, n__1;
    int    i, j, k;
    double unfl, ovfl;
    double xnorm;
    double errbnd;
    int    notran;
    double eps, tmp;
    double *rwork;
    doublecomplex *Aval;
    NCformat *Astore;

    /* Function prototypes */
    extern int    lsame_(char *, char *);
    extern int    izamax_(int *, doublecomplex *, int *);


    /* Quick exit if N = 0 or NRHS = 0. */
    if ( n <= 0 || nrhs <= 0 ) {
	reslts[0] = 0.;
	reslts[1] = 0.;
	return 0;
    }

    eps = dlamch_("Epsilon");
    unfl = dlamch_("Safe minimum");
    ovfl   = 1. / unfl;
    notran = (trans == NOTRANS);

    rwork  = (double *) SUPERLU_MALLOC(n*sizeof(double));
    if ( !rwork ) ABORT("SUPERLU_MALLOC fails for rwork");
    Astore = A->Store;
    Aval   = (doublecomplex *) Astore->nzval;
    
    /* Test 1:  Compute the maximum of   
       norm(X - XACT) / ( norm(X) * FERR )   
       over all the vectors X and XACT using the infinity-norm. */

    errbnd = 0.;
    for (j = 0; j < nrhs; ++j) {
	n__1 = n;
	imax = izamax_(&n__1, &x[j*ldx], &c__1);
	d__1 = (d__2 = x[imax-1 + j*ldx].r, fabs(d__2)) + 
               (d__3 = x[imax-1 + j*ldx].i, fabs(d__3));
	xnorm = SUPERLU_MAX(d__1,unfl);
	diff = 0.;
	for (i = 0; i < n; ++i) {
	    d__1 = (d__2 = x[i+j*ldx].r - xact[i+j*ldxact].r, fabs(d__2)) +
                   (d__3 = x[i+j*ldx].i - xact[i+j*ldxact].i, fabs(d__3));
	    diff = SUPERLU_MAX(diff, d__1);
	}

	if (xnorm > 1.) {
	    goto L20;
	} else if (diff <= ovfl * xnorm) {
	    goto L20;
	} else {
	    errbnd = 1. / eps;
	    goto L30;
	}

L20:
#if 0	
	if (diff / xnorm <= ferr[j]) {
	    d__1 = diff / xnorm / ferr[j];
	    errbnd = SUPERLU_MAX(errbnd,d__1);
	} else {
	    errbnd = 1. / eps;
	}
#endif
	d__1 = diff / xnorm / ferr[j];
	errbnd = SUPERLU_MAX(errbnd,d__1);
	/*printf("Ferr: %f\n", errbnd);*/
L30:
	;
    }
    reslts[0] = errbnd;

    /* Test 2: Compute the maximum of BERR / ( (n+1)*EPS + (*) ), where 
       (*) = (n+1)*UNFL / (min_i (abs(op(A))*abs(X) + abs(b))_i ) */

    for (k = 0; k < nrhs; ++k) {
	for (i = 0; i < n; ++i) 
            rwork[i] = (d__1 = b[i + k*ldb].r, fabs(d__1)) +
                       (d__2 = b[i + k*ldb].i, fabs(d__2));
	if ( notran ) {
	    for (j = 0; j < n; ++j) {
		tmp = (d__1 = x[j + k*ldx].r, fabs(d__1)) +
                      (d__2 = x[j + k*ldx].i, fabs(d__2));
		for (i = Astore->colptr[j]; i < Astore->colptr[j+1]; ++i) {
		    d__1 = (d__2 = Aval[i].r, fabs(d__2)) +
                           (d__3 = Aval[i].i, fabs(d__3));
		    rwork[Astore->rowind[i]] += d__1 * tmp;
                }
	    }
	} else {
	    for (j = 0; j < n; ++j) {
		tmp = 0.;
		for (i = Astore->colptr[j]; i < Astore->colptr[j+1]; ++i) {
		    irow = Astore->rowind[i];
		    d__1 = (d__2 = x[irow + k*ldx].r, fabs(d__2)) +
                           (d__3 = x[irow + k*ldx].i, fabs(d__3));
                    d__2 = (d__3 = Aval[i].r, fabs(d__3)) +
                           (d__4 = Aval[i].i, fabs(d__4));
		    tmp += d__2 * d__1;
		}
		rwork[j] += tmp;
	    }
	}

	axbi = rwork[0];
	for (i = 1; i < n; ++i) axbi = SUPERLU_MIN(axbi, rwork[i]);
	
	/* Computing MAX */
	d__1 = axbi, d__2 = (n + 1) * unfl;
	tmp = berr[k] / ((n + 1) * eps + (n + 1) * unfl / SUPERLU_MAX(d__1,d__2));
	
	if (k == 0) {
	    reslts[1] = tmp;
	} else {
	    reslts[1] = SUPERLU_MAX(reslts[1],tmp);
	}
    }

    SUPERLU_FREE(rwork);
    return 0;

} /* zgst07 */
