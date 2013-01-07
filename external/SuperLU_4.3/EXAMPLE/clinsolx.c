
/*
 * -- SuperLU routine (version 3.1) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 1, 2008
 *
 */
#include "slu_cdefs.h"

int main(int argc, char *argv[])
{
    char           equed[1];
    yes_no_t       equil;
    trans_t        trans;
    SuperMatrix    A, L, U;
    SuperMatrix    B, X;
    NCformat       *Astore;
    NCformat       *Ustore;
    SCformat       *Lstore;
    complex         *a;
    int            *asub, *xa;
    int            *perm_r; /* row permutations from partial pivoting */
    int            *perm_c; /* column permutation vector */
    int            *etree;
    void           *work;
    int            info, lwork, nrhs, ldx;
    int            i, m, n, nnz;
    complex         *rhsb, *rhsx, *xact;
    float         *R, *C;
    float         *ferr, *berr;
    float         u, rpg, rcond;
    mem_usage_t    mem_usage;
    superlu_options_t options;
    SuperLUStat_t stat;
    extern void  parse_command_line();

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC("Enter main()");
#endif

    /* Defaults */
    lwork = 0;
    nrhs  = 1;
    equil = YES;	
    u     = 1.0;
    trans = NOTRANS;
    
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

    /* Can use command line input to modify the defaults. */
    parse_command_line(argc, argv, &lwork, &u, &equil, &trans);
    options.Equil = equil;
    options.DiagPivotThresh = u;
    options.Trans = trans;

    /* Add more functionalities that the defaults. */
    options.PivotGrowth = YES;    /* Compute reciprocal pivot growth */
    options.ConditionNumber = YES;/* Compute reciprocal condition number */
    options.IterRefine = SLU_SINGLE;  /* Perform single-precision refinement */
    
    if ( lwork > 0 ) {
	work = SUPERLU_MALLOC(lwork);
	if ( !work ) {
	    ABORT("CLINSOLX: cannot allocate work[]");
	}
    }

    /* Read matrix A from a file in Harwell-Boeing format.*/
    creadhb(&m, &n, &nnz, &a, &asub, &xa);
    
    cCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_C, SLU_GE);
    Astore = A.Store;
    printf("Dimension %dx%d; # nonzeros %d\n", A.nrow, A.ncol, Astore->nnz);
    
    if ( !(rhsb = complexMalloc(m * nrhs)) ) ABORT("Malloc fails for rhsb[].");
    if ( !(rhsx = complexMalloc(m * nrhs)) ) ABORT("Malloc fails for rhsx[].");
    cCreate_Dense_Matrix(&B, m, nrhs, rhsb, m, SLU_DN, SLU_C, SLU_GE);
    cCreate_Dense_Matrix(&X, m, nrhs, rhsx, m, SLU_DN, SLU_C, SLU_GE);
    xact = complexMalloc(n * nrhs);
    ldx = n;
    cGenXtrue(n, nrhs, xact, ldx);
    cFillRHS(trans, nrhs, xact, ldx, &A, &B);
    
    if ( !(etree = intMalloc(n)) ) ABORT("Malloc fails for etree[].");
    if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");
    if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
    if ( !(R = (float *) SUPERLU_MALLOC(A.nrow * sizeof(float))) ) 
        ABORT("SUPERLU_MALLOC fails for R[].");
    if ( !(C = (float *) SUPERLU_MALLOC(A.ncol * sizeof(float))) )
        ABORT("SUPERLU_MALLOC fails for C[].");
    if ( !(ferr = (float *) SUPERLU_MALLOC(nrhs * sizeof(float))) )
        ABORT("SUPERLU_MALLOC fails for ferr[].");
    if ( !(berr = (float *) SUPERLU_MALLOC(nrhs * sizeof(float))) ) 
        ABORT("SUPERLU_MALLOC fails for berr[].");

    
    /* Initialize the statistics variables. */
    StatInit(&stat);
    
    /* Solve the system and compute the condition number
       and error bounds using dgssvx.      */
    
    cgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C,
           &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
           &mem_usage, &stat, &info);

    printf("cgssvx(): info %d\n", info);

    if ( info == 0 || info == n+1 ) {

        /* This is how you could access the solution matrix. */
        complex *sol = (complex*) ((DNformat*) X.Store)->nzval; 

	if ( options.PivotGrowth == YES )
            printf("Recip. pivot growth = %e\n", rpg);
	if ( options.ConditionNumber == YES )
	    printf("Recip. condition number = %e\n", rcond);
	if ( options.IterRefine != NOREFINE ) {
            printf("Iterative Refinement:\n");
	    printf("%8s%8s%16s%16s\n", "rhs", "Steps", "FERR", "BERR");
	    for (i = 0; i < nrhs; ++i)
	      printf("%8d%8d%16e%16e\n", i+1, stat.RefineSteps, ferr[i], berr[i]);
	}
        Lstore = (SCformat *) L.Store;
        Ustore = (NCformat *) U.Store;
	printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
    	printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
    	printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - n);
    	printf("FILL ratio = %.1f\n", (float)(Lstore->nnz + Ustore->nnz - n)/nnz);

	printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
	       mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
	     
	fflush(stdout);

    } else if ( info > 0 && lwork == -1 ) {
        printf("** Estimated memory: %d bytes\n", info - n);
    }

    if ( options.PrintStat ) StatPrint(&stat);
    StatFree(&stat);

    SUPERLU_FREE (rhsb);
    SUPERLU_FREE (rhsx);
    SUPERLU_FREE (xact);
    SUPERLU_FREE (etree);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    SUPERLU_FREE (R);
    SUPERLU_FREE (C);
    SUPERLU_FREE (ferr);
    SUPERLU_FREE (berr);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperMatrix_Store(&X);
    if ( lwork == 0 ) {
        Destroy_SuperNode_Matrix(&L);
        Destroy_CompCol_Matrix(&U);
    } else if ( lwork > 0 ) {
        SUPERLU_FREE(work);
    }

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC("Exit main()");
#endif
}


/*  
 * Parse command line inputs.
 */
void
parse_command_line(int argc, char *argv[], int *lwork,
                   float *u, yes_no_t *equil, trans_t *trans )
{
    int c;
    extern char *optarg;

    while ( (c = getopt(argc, argv, "hl:w:r:u:f:t:p:e:")) != EOF ) {
	switch (c) {
	  case 'h':
	    printf("Options:\n");
	    printf("\t-l <int> - length of work[*] array\n");
	    printf("\t-u <int> - pivoting threshold\n");
	    printf("\t-e <0 or 1> - equilibrate or not\n");
	    printf("\t-t <0 or 1> - solve transposed system or not\n");
	    exit(1);
	    break;
	  case 'l': *lwork = atoi(optarg);
	            break;
	  case 'u': *u = atof(optarg); 
	            break;
	  case 'e': *equil = atoi(optarg); 
	            break;
	  case 't': *trans = atoi(optarg);
	            break;
  	}
    }
}
