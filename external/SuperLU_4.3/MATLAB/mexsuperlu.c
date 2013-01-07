/*
 * -- SuperLU routine (version 4.2) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * June 30, 2009
 *
 * Modified: September 25, 2011,  compatible with 64-bit integer in R2006b
 */
#include <stdio.h>
#include "mex.h"
#include "matrix.h"

#include "slu_ddefs.h"

#define  MatlabMatrix mxArray


/* Aliases for input and output arguments */
#define A_in		prhs[0]
#define Pc_in		prhs[1]
#define L_out    	plhs[0]
#define U_out          	plhs[1]
#define Pr_out     	plhs[2]
#define Pc_out   	plhs[3]

void LUextract(SuperMatrix *, SuperMatrix *, double *, mwIndex *, mwIndex *, 
	       double *, mwIndex *, mwIndex *, int *, int*);

#define verbose (SPUMONI>0)
#define babble  (SPUMONI>1)
#define burble  (SPUMONI>2)

void mexFunction(
    int          nlhs,           /* number of expected outputs */
    MatlabMatrix *plhs[],        /* matrix pointer array returning outputs */
    int          nrhs,           /* number of inputs */
    const MatlabMatrix *prhs[]   /* matrix pointer array for inputs */
    )
{
    int SPUMONI;             /* ... as should the sparse monitor flag */
    double FlopsInSuperLU;   /* ... as should the flop counter */
    extern flops_t LUFactFlops(SuperLUStat_t *);
    
    /* Arguments to C dgstrf(). */
    SuperMatrix A;
    SuperMatrix Ac;        /* Matrix postmultiplied by Pc */
    SuperMatrix L, U;
    int	   	m, n, nnz;
    double      *val;
    int       	*rowind;
    int		*colptr;
    int    	*etree, *perm_r, *perm_c;
    mwSize      *perm_c_64;
    mwSize      *rowind_64;
    mwSize      *colptr_64;
    int         panel_size, relax;
    double      thresh = 1.0;       /* diagonal pivoting threshold */
    int		info;
    MatlabMatrix *X, *Y;            /* args to calls back to Matlab */
    int         i, mexerr;
    double      *dp;
    double      *Lval, *Uval;
    int         *Lrow, *Urow;
    int         *Lcol, *Ucol;
    mwIndex     *Lrow_64, *Lcol_64, *Urow_64, *Ucol_64;
    int         nnzL, nnzU, snnzL, snnzU;
    superlu_options_t options;
    SuperLUStat_t stat;

    /* Check number of arguments passed from Matlab. */
    if (nrhs != 2) {
	mexErrMsgTxt("SUPERLU requires 2 input arguments.");
    } else if (nlhs != 4) {
      	mexErrMsgTxt("SUPERLU requires 4 output arguments.");
    }   

    /* Read the Sparse Monitor Flag */
    X = mxCreateString("spumoni");
    mexerr = mexCallMATLAB(1, &Y, 1, &X, "sparsfun");
    SPUMONI = mxGetScalar(Y);
    mxDestroyArray(Y);
    mxDestroyArray(X);

    m = mxGetM(A_in);
    n = mxGetN(A_in);
    val = mxGetPr(A_in);
    rowind_64 = mxGetIr(A_in);
    colptr_64 = mxGetJc(A_in);
    perm_c_64 = mxGetIr(Pc_in); 
    nnz = colptr_64[n];
    if ( verbose ) mexPrintf("m = %d, n = %d, nnz  = %d\n", m, n, nnz);

    etree = (int *) mxCalloc(n, sizeof(int));
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
    panel_size = sp_ienv(1);
    relax      = sp_ienv(2);
    thresh     = 1.0;
    FlopsInSuperLU      = 0;

    set_default_options(&options);
    StatInit(&stat);

    if ( verbose ) mexPrintf("Apply column perm to A and compute etree...\n");
    sp_preorder(&options, &A, perm_c, etree, &Ac);

    if ( verbose ) {
	mexPrintf("LU factorization...\n");
	mexPrintf("\tpanel_size %d, relax %d, diag_pivot_thresh %.2g\n",
		  panel_size, relax, thresh);
    }

    dgstrf(&options, &Ac, relax, panel_size, etree,
	   NULL, 0, perm_c, perm_r, &L, &U, &stat, &info);

    if ( verbose ) mexPrintf("INFO from dgstrf %d\n", info);

#if 0 /* FLOPS is not available in the new Matlab. */
    /* Tell Matlab how many flops we did. */
    FlopsInSuperLU += LUFactFlops(&stat);
    if (verbose) mexPrintf("SUPERLU flops: %.f\n", FlopsInSuperLU);
    mexerr = mexCallMATLAB(1, &X, 0, NULL, "flops");
    *(mxGetPr(X)) += FlopsInSuperLU;
    mexerr = mexCallMATLAB(1, &Y, 1, &X, "flops");
    mxDestroyArray(Y);
    mxDestroyArray(X);
#endif
	
    /* Construct output arguments for Matlab. */
    if ( info >= 0 && info <= n ) {
	Pr_out = mxCreateDoubleMatrix(m, 1, mxREAL);   /* output row perm */
	dp = mxGetPr(Pr_out);
	for (i = 0; i < m; *dp++ = (double) perm_r[i++]+1);

	Pc_out = mxCreateDoubleMatrix(n, 1, mxREAL);   /* output col perm */
	dp = mxGetPr(Pc_out);
	for (i = 0; i < m; ++i) dp[i] = (double) perm_c[i]+1;
	
	/* Now for L and U */
	nnzL = ((SCformat*)L.Store)->nnz; /* count diagonals */
   	nnzU = ((NCformat*)U.Store)->nnz;
	L_out = mxCreateSparse(m, n, nnzL, mxREAL);
	Lval = mxGetPr(L_out);
	Lrow_64 = mxGetIr(L_out);
	Lcol_64 = mxGetJc(L_out);
	U_out = mxCreateSparse(m, n, nnzU, mxREAL);
	Uval = mxGetPr(U_out);
	Urow_64 = mxGetIr(U_out);
	Ucol_64 = mxGetJc(U_out);

	LUextract(&L, &U, Lval, Lrow_64, Lcol_64, Uval, Urow_64, Ucol_64,
		  &snnzL, &snnzU);

	if ( babble ) 
	  for (i = 0; i <= n; ++i) printf("Lcol_64[%d] %d\n", i, Lcol_64[i]);

	printf("nnzL = %d, nnzU = %d\n", nnzL, nnzU);
	if ( babble ) {
	  for (i=0; i < nnzL; ++i) 
	    mexPrintf("Lrow_64[%d] %d\n", i, Lrow_64[i]);
	  for (i = 0; i < snnzU; ++i) 
	    mexPrintf("Urow_64[%d] = %d\n", i, Urow_64[i]);
	}

        Destroy_CompCol_Permuted(&Ac);
	Destroy_SuperNode_Matrix(&L);
	Destroy_CompCol_Matrix(&U);

	if (verbose) mexPrintf("factor nonzeros: %d unsqueezed, %d squeezed.\n",
			      nnzL + nnzU, snnzL + snnzU);
    } else {
	mexErrMsgTxt("Error returned from C dgstrf().");
    }

    mxFree(etree);
    mxFree(perm_r);
    mxFree(perm_c);
    mxFree(rowind);
    mxFree(colptr);

    StatFree(&stat);

    return;
}

void
LUextract(SuperMatrix *L, SuperMatrix *U, double *Lval, mwIndex *Lrow,
	  mwIndex *Lcol, double *Uval, mwIndex *Urow, mwIndex *Ucol,
	  int *snnzL, int *snnzU)
{
    int         i, j, k;
    int         upper;
    int         fsupc, istart, nsupr;
    int         lastl = 0, lastu = 0;
    SCformat    *Lstore;
    NCformat    *Ustore;
    double      *SNptr;

    Lstore = L->Store;
    Ustore = U->Store;
    Lcol[0] = 0;
    Ucol[0] = 0;
    
    /* for each supernode */
    for (k = 0; k <= Lstore->nsuper; ++k) {
	
	fsupc = L_FST_SUPC(k);
	istart = L_SUB_START(fsupc);
	nsupr = L_SUB_START(fsupc+1) - istart;
	upper = 1;
	
	/* for each column in the supernode */
	for (j = fsupc; j < L_FST_SUPC(k+1); ++j) {
	    SNptr = &((double*)Lstore->nzval)[L_NZ_START(j)];

	    /* Extract U */
	    for (i = U_NZ_START(j); i < U_NZ_START(j+1); ++i) {
		Uval[lastu] = ((double*)Ustore->nzval)[i];
 		/* Matlab doesn't like explicit zero. */
		if (Uval[lastu] != 0.0) Urow[lastu++] = (mwIndex) U_SUB(i);
	    }
	    for (i = 0; i < upper; ++i) { /* upper triangle in the supernode */
		Uval[lastu] = SNptr[i];
 		/* Matlab doesn't like explicit zero. */
		if (Uval[lastu] != 0.0) Urow[lastu++] = (mwIndex)L_SUB(istart+i);
	    }
	    Ucol[j+1] = lastu;

	    /* Extract L */
	    Lval[lastl] = 1.0; /* unit diagonal */
	    Lrow[lastl++] = L_SUB(istart + upper - 1);
	    for (i = upper; i < nsupr; ++i) {
		Lval[lastl] = SNptr[i];
 		/* Matlab doesn't like explicit zero. */
		if (Lval[lastl] != 0.0) Lrow[lastl++] = (mwIndex)L_SUB(istart+i);
	    }
	    Lcol[j+1] = lastl;

	    ++upper;
	    
	} /* for j ... */
	
    } /* for k ... */

    *snnzL = lastl;
    *snnzU = lastu;
}
