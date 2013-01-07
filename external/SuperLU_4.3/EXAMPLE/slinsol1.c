
/*
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 */
#include "slu_sdefs.h"

int main(int argc, char *argv[])
{
    SuperMatrix A;
    NCformat *Astore;
    float   *a;
    int      *asub, *xa;
    int      *perm_c; /* column permutation vector */
    int      *perm_r; /* row permutations from partial pivoting */
    SuperMatrix L;      /* factor L */
    SCformat *Lstore;
    SuperMatrix U;      /* factor U */
    NCformat *Ustore;
    SuperMatrix B;
    int      nrhs, ldx, info, m, n, nnz;
    float   *xact, *rhs;
    mem_usage_t   mem_usage;
    superlu_options_t options;
    SuperLUStat_t stat;
    
#if ( DEBUGlevel>=1 )
    CHECK_MALLOC("Enter main()");
#endif

    /* Set the default input options:
	options.Fact = DOFACT;
        options.Equil = YES;
    	options.ColPerm = COLAMD;
	options.DiagPivotThresh = 1.0;
    	options.Trans = NOTRANS;
    	options.IterRefine = NOREFINE;
    	options.SymmetricMode = NO;
    	options.PivotGrowth = NO;
    	options.ConditionNumber = NO;
    	options.PrintStat = YES;
     */
    set_default_options(&options);

    /* Now we modify the default options to use the symmetric mode. */
    options.SymmetricMode = YES;
    options.ColPerm = MMD_AT_PLUS_A;
    options.DiagPivotThresh = 0.001;

    /* Read the matrix in Harwell-Boeing format. */
    sreadhb(&m, &n, &nnz, &a, &asub, &xa);

    sCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_S, SLU_GE);
    Astore = A.Store;
    printf("Dimension %dx%d; # nonzeros %d\n", A.nrow, A.ncol, Astore->nnz);
    
    nrhs   = 1;
    if ( !(rhs = floatMalloc(m * nrhs)) ) ABORT("Malloc fails for rhs[].");
    sCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_S, SLU_GE);
    xact = floatMalloc(n * nrhs);
    ldx = n;
    sGenXtrue(n, nrhs, xact, ldx);
    sFillRHS(options.Trans, nrhs, xact, ldx, &A, &B);

    if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
    if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");

    /* Initialize the statistics variables. */
    StatInit(&stat);
    
    sgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);
    
    if ( info == 0 ) {

	/* This is how you could access the solution matrix. */
        float *sol = (float*) ((DNformat*) B.Store)->nzval; 

	 /* Compute the infinity norm of the error. */
	sinf_norm_error(nrhs, &B, xact);

	Lstore = (SCformat *) L.Store;
	Ustore = (NCformat *) U.Store;
    	printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
    	printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
    	printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - n);
    	printf("FILL ratio = %.1f\n", (float)(Lstore->nnz + Ustore->nnz - n)/nnz);
	
	sQuerySpace(&L, &U, &mem_usage);
	printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
	       mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
	
    } else {
	printf("sgssv() error returns INFO= %d\n", info);
	if ( info <= n ) { /* factorization completes */
	    sQuerySpace(&L, &U, &mem_usage);
	    printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
		   mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
	}
    }

    if ( options.PrintStat ) StatPrint(&stat);
    StatFree(&stat);

    SUPERLU_FREE (rhs);
    SUPERLU_FREE (xact);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC("Exit main()");
#endif
}

