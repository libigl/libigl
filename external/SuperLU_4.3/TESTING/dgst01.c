
/*
 * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 *
 */
#include <math.h>
#include "slu_ddefs.h"

int dgst01(int m, int n, SuperMatrix *A, SuperMatrix *L, 
		SuperMatrix *U, int *perm_c, int *perm_r, double *resid)
{
/* 
    Purpose   
    =======   

    DGST01 reconstructs a matrix A from its L*U factorization and   
    computes the residual   
       norm(L*U - A) / ( N * norm(A) * EPS ),   
    where EPS is the machine epsilon.   

    Arguments   
    ==========   

    M       (input) INT   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INT   
            The number of columns of the matrix A.  N >= 0.   

    A       (input) SuperMatrix *, dimension (A->nrow, A->ncol)
            The original M x N matrix A.   

    L       (input) SuperMatrix *, dimension (L->nrow, L->ncol)
            The factor matrix L.

    U       (input) SuperMatrix *, dimension (U->nrow, U->ncol)
            The factor matrix U.

    perm_c (input) INT array, dimension (N)
            The column permutation from DGSTRF.   

    perm_r  (input) INT array, dimension (M)
            The pivot indices from DGSTRF.   

    RESID   (output) DOUBLE*
            norm(L*U - A) / ( N * norm(A) * EPS )   

    ===================================================================== 
*/  

    /* Local variables */
    double zero = 0.0;
    int i, j, k, arow, lptr,isub,  urow, superno, fsupc, u_part;
    double utemp, comp_temp;
    double anorm, tnorm, cnorm;
    double eps;
    double *work;
    SCformat *Lstore;
    NCformat *Astore, *Ustore;
    double *Aval, *Lval, *Uval;
    int *colbeg, *colend;

    /* Function prototypes */
    extern double dlangs(char *, SuperMatrix *);

    /* Quick exit if M = 0 or N = 0. */

    if (m <= 0 || n <= 0) {
	*resid = 0.f;
	return 0;
    }

    work = (double *)doubleCalloc(m);

    Astore = A->Store;
    Aval = Astore->nzval;
    Lstore = L->Store;
    Lval = Lstore->nzval;
    Ustore = U->Store;
    Uval = Ustore->nzval;

    /* Determine EPS and the norm of A. */
    eps = dlamch_("Epsilon");
    anorm = dlangs("1", A);
    cnorm = 0.;

    /* Compute the product L*U, one column at a time */
    for (k = 0; k < n; ++k) {

	/* The U part outside the rectangular supernode */
        for (i = U_NZ_START(k); i < U_NZ_START(k+1); ++i) {
	    urow = U_SUB(i);
	    utemp = Uval[i];
            superno = Lstore->col_to_sup[urow];
	    fsupc = L_FST_SUPC(superno);
	    u_part = urow - fsupc + 1;
	    lptr = L_SUB_START(fsupc) + u_part;
            work[L_SUB(lptr-1)] -= utemp;   /* L_ii = 1 */
	    for (j = L_NZ_START(urow) + u_part; j < L_NZ_START(urow+1); ++j) {
                isub = L_SUB(lptr);
	        work[isub] -= Lval[j] * utemp;
	        ++lptr;
	    }
	}

	/* The U part inside the rectangular supernode */
	superno = Lstore->col_to_sup[k];
	fsupc = L_FST_SUPC(superno);
	urow = L_NZ_START(k);
	for (i = fsupc; i <= k; ++i) {
	    utemp = Lval[urow++];
	    u_part = i - fsupc + 1;
	    lptr = L_SUB_START(fsupc) + u_part;
            work[L_SUB(lptr-1)] -= utemp;   /* L_ii = 1 */
	    for (j = L_NZ_START(i)+u_part; j < L_NZ_START(i+1); ++j) {
                isub = L_SUB(lptr);
	        work[isub] -= Lval[j] * utemp;
	        ++lptr;
	    }
	}

	/* Now compute A[k] - (L*U)[k] (Both matrices may be permuted.) */

	colbeg = intMalloc(n);
	colend = intMalloc(n);
	for (i = 0; i < n; i++) {
	    colbeg[perm_c[i]] = Astore->colptr[i]; 
	    colend[perm_c[i]] = Astore->colptr[i+1];
	}
	
	for (i = colbeg[k]; i < colend[k]; ++i) {
	    arow = Astore->rowind[i];
	    work[perm_r[arow]] += Aval[i];
        }

	/* Now compute the 1-norm of the column vector work */
        tnorm = 0.;
	for (i = 0; i < m; ++i) {
	    tnorm += fabs(work[i]);
	    work[i] = zero;
	}
	cnorm = SUPERLU_MAX(tnorm, cnorm);
    }

    *resid = cnorm;

    if (anorm <= 0.f) {
	if (*resid != 0.f) {
	    *resid = 1.f / eps;
	}
    } else {
	*resid = *resid / (float) n / anorm / eps;
    }

    SUPERLU_FREE(work);
    SUPERLU_FREE(colbeg);
    SUPERLU_FREE(colend);
    return 0;

/*     End of DGST01 */

} /* dgst01_ */

