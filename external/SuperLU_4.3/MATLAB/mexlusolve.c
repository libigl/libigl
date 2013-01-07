/*
 * -- SuperLU routine (version 4.2) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 * Modified: September 25, 2011,  compatible with 64-bit integer in R2006b
 */
#include <stdio.h>
#include "mex.h"
#include "slu_ddefs.h"

#define  MatlabMatrix mxArray

#if 0   /* V4 */
#define  MatlabMatrix Matrix
#endif


/* Aliases for input and output arguments */
#define A_in		prhs[0]
#define b_in    	prhs[1]
#define Pc_in		prhs[2]
#define x_out		plhs[0]

#define verbose (SPUMONI>0)
#define babble  (SPUMONI>1)
#define burble  (SPUMONI>2)

void mexFunction(
    int          nlhs,           /* number of expected outputs */
    MatlabMatrix *plhs[],        /* matrix pointer array returning outputs */
    int          nrhs,           /* number of inputs */
    const MatlabMatrix *prhs[]   /* matrix pointer array for inputs. */
#if 0 /* V4 */
    MatlabMatrix *prhs[]         /* matrix pointer array for inputs */
#endif
    )
{
    int SPUMONI;             /* ... as should the sparse monitor flag */
    double FlopsInSuperLU;   /* ... as should the flop counter. */
    extern flops_t LUFactFlops(SuperLUStat_t*);
    extern flops_t LUSolveFlops(SuperLUStat_t*);
    
    /* Arguments to dgssv(). */
    SuperMatrix A;
    SuperMatrix B;
    SuperMatrix L, U;
    int	   	m, n, nnz;
    int         numrhs;
    double 	*vb, *x;
    double      *val;
    int       	*rowind;
    int		*colptr;
    int    	*perm_r, *perm_c;
    mwSize      *perm_c_64;
    mwSize      *rowind_64;
    mwSize      *colptr_64;
    int		info;
    MatlabMatrix *X, *Y;            /* args to calls back to Matlab */
    int         i, mexerr;
    superlu_options_t options;
    SuperLUStat_t stat;

    /* Check number of arguments passed from Matlab. */
    if (nrhs != 3) {
	mexErrMsgTxt("LUSOLVE requires 3 input arguments.");
    } else if (nlhs != 1) {
      	mexErrMsgTxt("LUSOLVE requires 1 output argument.");
    }   

    /* Read the Sparse Monitor Flag */
    X = mxCreateString("spumoni");
    mexerr = mexCallMATLAB(1, &Y, 1, &X, "sparsfun");
    SPUMONI = mxGetScalar(Y);
    mxDestroyArray(Y);
    mxDestroyArray(X);
#if 0 /* V4 */
    mxFreeMatrix(Y);
    mxFreeMatrix(X);
#endif

    m = mxGetM(A_in);
    n = mxGetN(A_in);
    numrhs = mxGetN(b_in);
    if ( babble ) printf("m=%d, n=%d, numrhs=%d\n", m, n, numrhs);
    vb = mxGetPr(b_in);
    x_out = mxCreateDoubleMatrix(m, numrhs, mxREAL);

    x = mxGetPr(x_out);
    val = mxGetPr(A_in);
    rowind_64 = mxGetIr(A_in);
    colptr_64 = mxGetJc(A_in);
    perm_c_64 = mxGetIr(Pc_in); 
    nnz = colptr_64[n];
    perm_r = (int *) mxCalloc(m, sizeof(int));
    perm_c = (int *) mxMalloc(n * sizeof(int));
    rowind = (int *) mxMalloc(nnz * sizeof(int));
    colptr = (int *) mxMalloc((n+1) * sizeof(int));
    for (i = 0; i < n; ++i) {
	perm_c[i] = perm_c_64[i];
	colptr[i] = colptr_64[i];
	/*printf("perm_c[%d] %d\n", i, perm_c[i]);*/
    }
    colptr[n] = colptr_64[n];
    for (i = 0; i < nnz; ++i) rowind[i] = rowind_64[i];

    dCreate_CompCol_Matrix(&A, m, n, nnz, val, rowind, colptr,
			   SLU_NC, SLU_D, SLU_GE);
    dCopy_Dense_Matrix(m, numrhs, vb, m, x, m);
    dCreate_Dense_Matrix(&B, m, numrhs, x, m, SLU_DN, SLU_D, SLU_GE);

    FlopsInSuperLU = 0;

    set_default_options(&options);
    options.ColPerm = MY_PERMC;
    StatInit(&stat);

    /* Call simple driver */
    if ( verbose )
      mexPrintf("Call LUSOLVE, use SUPERLU to factor first ...\n");
    dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

#if 0 /* FLOPS is not available in the new Matlab. */
    /* Tell Matlab how many flops we did. */
    FlopsInSuperLU += LUFactFlops(&stat) + LUSolveFlops(&stat);
    if ( verbose ) mexPrintf("LUSOLVE flops: %.f\n", FlopsInSuperLU);
    mexerr = mexCallMATLAB(1, &X, 0, NULL, "flops");
    *(mxGetPr(X)) += FlopsInSuperLU;
    mexerr = mexCallMATLAB(1, &Y, 1, &X, "flops");
    mxDestroyArray(Y);
    mxDestroyArray(X);
#endif

    /* Construct Matlab solution matrix. */
    if ( !info ) {
        Destroy_SuperNode_Matrix(&L);
        Destroy_CompCol_Matrix(&U);
        if ( babble ) printf("Destroy L & U from SuperLU...\n");
    } else {
        printf("dgssv info = %d\n", info);
	mexErrMsgTxt("Error returned from C dgssv().");
    }

    mxFree(perm_r);
    mxFree(perm_c);
    mxFree(rowind);
    mxFree(colptr);
    StatFree(&stat);

    return;
 
}
