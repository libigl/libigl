
/*
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 */

#include "slu_sdefs.h"

#define HANDLE_SIZE  8
/* kind of integer to hold a pointer.  Use int.
   This might need to be changed on 64-bit systems. */
typedef int fptr;  /* 32-bit by default */

typedef struct {
    SuperMatrix *L;
    SuperMatrix *U;
    int *perm_c;
    int *perm_r;
} factors_t;

void
c_fortran_sgssv_(int *iopt, int *n, int *nnz, int *nrhs, 
                 float *values, int *rowind, int *colptr,
                 float *b, int *ldb,
		 fptr *f_factors, /* a handle containing the address
				     pointing to the factored matrices */
		 int *info)

{
/* 
 * This routine can be called from Fortran.
 *
 * iopt (input) int
 *      Specifies the operation:
 *      = 1, performs LU decomposition for the first time
 *      = 2, performs triangular solve
 *      = 3, free all the storage in the end
 *
 * f_factors (input/output) fptr* 
 *      If iopt == 1, it is an output and contains the pointer pointing to
 *                    the structure of the factored matrices.
 *      Otherwise, it it an input.
 *
 */
 
    SuperMatrix A, AC, B;
    SuperMatrix *L, *U;
    int *perm_r; /* row permutations from partial pivoting */
    int *perm_c; /* column permutation vector */
    int *etree;  /* column elimination tree */
    SCformat *Lstore;
    NCformat *Ustore;
    int      i, panel_size, permc_spec, relax;
    trans_t  trans;
    mem_usage_t   mem_usage;
    superlu_options_t options;
    SuperLUStat_t stat;
    factors_t *LUfactors;

    trans = NOTRANS;

    if ( *iopt == 1 ) { /* LU decomposition */

        /* Set the default input options. */
        set_default_options(&options);

	/* Initialize the statistics variables. */
	StatInit(&stat);

	/* Adjust to 0-based indexing */
	for (i = 0; i < *nnz; ++i) --rowind[i];
	for (i = 0; i <= *n; ++i) --colptr[i];

	sCreate_CompCol_Matrix(&A, *n, *n, *nnz, values, rowind, colptr,
			       SLU_NC, SLU_S, SLU_GE);
	L = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
	U = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
	if ( !(perm_r = intMalloc(*n)) ) ABORT("Malloc fails for perm_r[].");
	if ( !(perm_c = intMalloc(*n)) ) ABORT("Malloc fails for perm_c[].");
	if ( !(etree = intMalloc(*n)) ) ABORT("Malloc fails for etree[].");

	/*
	 * Get column permutation vector perm_c[], according to permc_spec:
	 *   permc_spec = 0: natural ordering 
	 *   permc_spec = 1: minimum degree on structure of A'*A
	 *   permc_spec = 2: minimum degree on structure of A'+A
	 *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
	 */    	
	permc_spec = options.ColPerm;
	get_perm_c(permc_spec, &A, perm_c);
	
	sp_preorder(&options, &A, perm_c, etree, &AC);

	panel_size = sp_ienv(1);
	relax = sp_ienv(2);

	sgstrf(&options, &AC, relax, panel_size, etree,
                NULL, 0, perm_c, perm_r, L, U, &stat, info);

	if ( *info == 0 ) {
	    Lstore = (SCformat *) L->Store;
	    Ustore = (NCformat *) U->Store;
	    printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
	    printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
	    printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz);
	    sQuerySpace(L, U, &mem_usage);
	    printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
		   mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
	} else {
	    printf("sgstrf() error returns INFO= %d\n", *info);
	    if ( *info <= *n ) { /* factorization completes */
		sQuerySpace(L, U, &mem_usage);
		printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
		       mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
	    }
	}
	
	/* Restore to 1-based indexing */
	for (i = 0; i < *nnz; ++i) ++rowind[i];
	for (i = 0; i <= *n; ++i) ++colptr[i];

	/* Save the LU factors in the factors handle */
	LUfactors = (factors_t*) SUPERLU_MALLOC(sizeof(factors_t));
	LUfactors->L = L;
	LUfactors->U = U;
	LUfactors->perm_c = perm_c;
	LUfactors->perm_r = perm_r;
	*f_factors = (fptr) LUfactors;

	/* Free un-wanted storage */
	SUPERLU_FREE(etree);
	Destroy_SuperMatrix_Store(&A);
	Destroy_CompCol_Permuted(&AC);
	StatFree(&stat);

    } else if ( *iopt == 2 ) { /* Triangular solve */
	/* Initialize the statistics variables. */
	StatInit(&stat);

	/* Extract the LU factors in the factors handle */
	LUfactors = (factors_t*) *f_factors;
	L = LUfactors->L;
	U = LUfactors->U;
	perm_c = LUfactors->perm_c;
	perm_r = LUfactors->perm_r;

	sCreate_Dense_Matrix(&B, *n, *nrhs, b, *ldb, SLU_DN, SLU_S, SLU_GE);

        /* Solve the system A*X=B, overwriting B with X. */
        sgstrs (trans, L, U, perm_c, perm_r, &B, &stat, info);

	Destroy_SuperMatrix_Store(&B);
	StatFree(&stat);

    } else if ( *iopt == 3 ) { /* Free storage */
	/* Free the LU factors in the factors handle */
	LUfactors = (factors_t*) *f_factors;
	SUPERLU_FREE (LUfactors->perm_r);
	SUPERLU_FREE (LUfactors->perm_c);
	Destroy_SuperNode_Matrix(LUfactors->L);
	Destroy_CompCol_Matrix(LUfactors->U);
        SUPERLU_FREE (LUfactors->L);
        SUPERLU_FREE (LUfactors->U);
	SUPERLU_FREE (LUfactors);
    } else {
	fprintf(stderr,"Invalid iopt=%d passed to c_fortran_sgssv()\n",*iopt);
	exit(-1);
    }
}


