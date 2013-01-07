
/*
 * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 *
 */

#include "slu_sdefs.h"

/*
 * Convert a full matrix into a sparse matrix format. 
 */
int
sp_sconvert(int m, int n, float *A, int lda, int kl, int ku,
	   float *a, int *asub, int *xa, int *nnz)
{
    int     lasta = 0;
    int     i, j, ilow, ihigh;
    int     *row;
    float  *val;

    for (j = 0; j < n; ++j) {
	xa[j] = lasta;
	val = &a[xa[j]];
	row = &asub[xa[j]];

	ilow = SUPERLU_MAX(0, j - ku);
	ihigh = SUPERLU_MIN(n-1, j + kl);
	for (i = ilow; i <= ihigh; ++i) {
	    val[i-ilow] = A[i + j*lda];
	    row[i-ilow] = i;
	}
	lasta += ihigh - ilow + 1;
    }

    xa[n] = *nnz = lasta;
    return 0;
}


