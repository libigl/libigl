
/*
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 */
/*
 * File name:		cdrive.c
 * Purpose:             MAIN test program
 */
#include <string.h>
#include "slu_cdefs.h"

#define NTESTS    5      /* Number of test types */
#define NTYPES    11     /* Number of matrix types */
#define NTRAN     2    
#define THRESH    20.0
#define FMT1      "%10s:n=%d, test(%d)=%12.5g\n"
#define	FMT2      "%10s:fact=%4d, trans=%4d, equed=%c, n=%d, imat=%d, test(%d)=%12.5g\n"
#define FMT3      "%10s:info=%d, izero=%d, n=%d, nrhs=%d, imat=%d, nfail=%d\n"

static void
parse_command_line(int argc, char *argv[], char *matrix_type,
		   int *n, int *w, int *relax, int *nrhs, int *maxsuper,
		   int *rowblk, int *colblk, int *lwork, double *u);

main(int argc, char *argv[])
{
/* 
 * Purpose
 * =======
 *
 * CDRIVE is the main test program for the COMPLEX linear 
 * equation driver routines CGSSV and CGSSVX.
 * 
 * The program is invoked by a shell script file -- ctest.csh.
 * The output from the tests are written into a file -- ctest.out.
 *
 * =====================================================================
 */
    complex         *a, *a_save;
    int            *asub, *asub_save;
    int            *xa, *xa_save;
    SuperMatrix  A, B, X, L, U;
    SuperMatrix  ASAV, AC;
    mem_usage_t    mem_usage;
    int            *perm_r; /* row permutation from partial pivoting */
    int            *perm_c, *pc_save; /* column permutation */
    int            *etree;
    complex  zero = {0.0, 0.0};
    float         *R, *C;
    float         *ferr, *berr;
    float         *rwork;
    complex	   *wwork;
    void           *work;
    int            info, lwork, nrhs, panel_size, relax;
    int            m, n, nnz;
    complex         *xact;
    complex         *rhsb, *solx, *bsav;
    int            ldb, ldx;
    float         rpg, rcond;
    int            i, j, k1;
    float         rowcnd, colcnd, amax;
    int            maxsuper, rowblk, colblk;
    int            prefact, nofact, equil, iequed;
    int            nt, nrun, nfail, nerrs, imat, fimat, nimat;
    int            nfact, ifact, itran;
    int            kl, ku, mode, lda;
    int            zerot, izero, ioff;
    double         u;
    float         anorm, cndnum;
    complex         *Afull;
    float         result[NTESTS];
    superlu_options_t options;
    fact_t         fact;
    trans_t        trans;
    SuperLUStat_t  stat;
    static char    matrix_type[8];
    static char    equed[1], path[4], sym[1], dist[1];

    /* Fixed set of parameters */
    int            iseed[]  = {1988, 1989, 1990, 1991};
    static char    equeds[]  = {'N', 'R', 'C', 'B'};
    static fact_t  facts[] = {FACTORED, DOFACT, SamePattern,
			      SamePattern_SameRowPerm};
    static trans_t transs[]  = {NOTRANS, TRANS, CONJ};

    /* Some function prototypes */ 
    extern int cgst01(int, int, SuperMatrix *, SuperMatrix *, 
		      SuperMatrix *, int *, int *, float *);
    extern int cgst02(trans_t, int, int, int, SuperMatrix *, complex *,
                      int, complex *, int, float *resid);
    extern int cgst04(int, int, complex *, int, 
                      complex *, int, float rcond, float *resid);
    extern int cgst07(trans_t, int, int, SuperMatrix *, complex *, int,
                         complex *, int, complex *, int, 
                         float *, float *, float *);
    extern int clatb4_(char *, int *, int *, int *, char *, int *, int *, 
	               float *, int *, float *, char *);
    extern int clatms_(int *, int *, char *, int *, char *, float *d,
                       int *, float *, float *, int *, int *,
                       char *, complex *, int *, complex *, int *);
    extern int sp_cconvert(int, int, complex *, int, int, int,
	                   complex *a, int *, int *, int *);


    /* Executable statements */

    strcpy(path, "CGE");
    nrun  = 0;
    nfail = 0;
    nerrs = 0;

    /* Defaults */
    lwork      = 0;
    n          = 1;
    nrhs       = 1;
    panel_size = sp_ienv(1);
    relax      = sp_ienv(2);
    u          = 1.0;
    strcpy(matrix_type, "LA");
    parse_command_line(argc, argv, matrix_type, &n,
		       &panel_size, &relax, &nrhs, &maxsuper,
		       &rowblk, &colblk, &lwork, &u);
    if ( lwork > 0 ) {
	work = SUPERLU_MALLOC(lwork);
	if ( !work ) {
	    fprintf(stderr, "expert: cannot allocate %d bytes\n", lwork);
	    exit (-1);
	}
    }

    /* Set the default input options. */
    set_default_options(&options);
    options.DiagPivotThresh = u;
    options.PrintStat = NO;
    options.PivotGrowth = YES;
    options.ConditionNumber = YES;
    options.IterRefine = SLU_SINGLE;
    
    if ( strcmp(matrix_type, "LA") == 0 ) {
	/* Test LAPACK matrix suite. */
	m = n;
	lda = SUPERLU_MAX(n, 1);
	nnz = n * n;        /* upper bound */
	fimat = 1;
	nimat = NTYPES;
	Afull = complexCalloc(lda * n);
	callocateA(n, nnz, &a, &asub, &xa);
    } else {
	/* Read a sparse matrix */
	fimat = nimat = 0;
	creadhb(&m, &n, &nnz, &a, &asub, &xa);
    }

    callocateA(n, nnz, &a_save, &asub_save, &xa_save);
    rhsb = complexMalloc(m * nrhs);
    bsav = complexMalloc(m * nrhs);
    solx = complexMalloc(n * nrhs);
    ldb  = m;
    ldx  = n;
    cCreate_Dense_Matrix(&B, m, nrhs, rhsb, ldb, SLU_DN, SLU_C, SLU_GE);
    cCreate_Dense_Matrix(&X, n, nrhs, solx, ldx, SLU_DN, SLU_C, SLU_GE);
    xact = complexMalloc(n * nrhs);
    etree   = intMalloc(n);
    perm_r  = intMalloc(n);
    perm_c  = intMalloc(n);
    pc_save = intMalloc(n);
    R       = (float *) SUPERLU_MALLOC(m*sizeof(float));
    C       = (float *) SUPERLU_MALLOC(n*sizeof(float));
    ferr    = (float *) SUPERLU_MALLOC(nrhs*sizeof(float));
    berr    = (float *) SUPERLU_MALLOC(nrhs*sizeof(float));
    j = SUPERLU_MAX(m,n) * SUPERLU_MAX(4,nrhs);    
    rwork   = (float *) SUPERLU_MALLOC(j*sizeof(float));
    for (i = 0; i < j; ++i) rwork[i] = 0.;
    if ( !R ) ABORT("SUPERLU_MALLOC fails for R");
    if ( !C ) ABORT("SUPERLU_MALLOC fails for C");
    if ( !ferr ) ABORT("SUPERLU_MALLOC fails for ferr");
    if ( !berr ) ABORT("SUPERLU_MALLOC fails for berr");
    if ( !rwork ) ABORT("SUPERLU_MALLOC fails for rwork");
    wwork   = complexCalloc( SUPERLU_MAX(m,n) * SUPERLU_MAX(4,nrhs) );

    for (i = 0; i < n; ++i) perm_c[i] = pc_save[i] = i;
    options.ColPerm = MY_PERMC;

    for (imat = fimat; imat <= nimat; ++imat) { /* All matrix types */
	
	if ( imat ) {

	    /* Skip types 5, 6, or 7 if the matrix size is too small. */
	    zerot = (imat >= 5 && imat <= 7);
	    if ( zerot && n < imat-4 )
		continue;
	    
	    /* Set up parameters with CLATB4 and generate a test matrix
	       with CLATMS.  */
	    clatb4_(path, &imat, &n, &n, sym, &kl, &ku, &anorm, &mode,
		    &cndnum, dist);

	    clatms_(&n, &n, dist, iseed, sym, &rwork[0], &mode, &cndnum,
		    &anorm, &kl, &ku, "No packing", Afull, &lda,
		    &wwork[0], &info);

	    if ( info ) {
		printf(FMT3, "CLATMS", info, izero, n, nrhs, imat, nfail);
		continue;
	    }

	    /* For types 5-7, zero one or more columns of the matrix
	       to test that INFO is returned correctly.   */
	    if ( zerot ) {
		if ( imat == 5 ) izero = 1;
		else if ( imat == 6 ) izero = n;
		else izero = n / 2 + 1;
		ioff = (izero - 1) * lda;
		if ( imat < 7 ) {
		    for (i = 0; i < n; ++i) Afull[ioff + i] = zero;
		} else {
		    for (j = 0; j < n - izero + 1; ++j)
			for (i = 0; i < n; ++i)
			    Afull[ioff + i + j*lda] = zero;
		}
	    } else {
		izero = 0;
	    }

	    /* Convert to sparse representation. */
	    sp_cconvert(n, n, Afull, lda, kl, ku, a, asub, xa, &nnz);

	} else {
	    izero = 0;
	    zerot = 0;
	}
	
	cCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_C, SLU_GE);

	/* Save a copy of matrix A in ASAV */
	cCreate_CompCol_Matrix(&ASAV, m, n, nnz, a_save, asub_save, xa_save,
			      SLU_NC, SLU_C, SLU_GE);
	cCopy_CompCol_Matrix(&A, &ASAV);
	
	/* Form exact solution. */
	cGenXtrue(n, nrhs, xact, ldx);
	
	StatInit(&stat);

	for (iequed = 0; iequed < 4; ++iequed) {
	    *equed = equeds[iequed];
	    if (iequed == 0) nfact = 4;
	    else nfact = 1; /* Only test factored, pre-equilibrated matrix */

	    for (ifact = 0; ifact < nfact; ++ifact) {
		fact = facts[ifact];
		options.Fact = fact;

		for (equil = 0; equil < 2; ++equil) {
		    options.Equil = equil;
		    prefact   = ( options.Fact == FACTORED ||
				  options.Fact == SamePattern_SameRowPerm );
                                /* Need a first factor */
		    nofact    = (options.Fact != FACTORED);  /* Not factored */

		    /* Restore the matrix A. */
		    cCopy_CompCol_Matrix(&ASAV, &A);
			
		    if ( zerot ) {
                        if ( prefact ) continue;
		    } else if ( options.Fact == FACTORED ) {
                        if ( equil || iequed ) {
			    /* Compute row and column scale factors to
			       equilibrate matrix A.    */
			    cgsequ(&A, R, C, &rowcnd, &colcnd, &amax, &info);

			    /* Force equilibration. */
			    if ( !info && n > 0 ) {
				if ( lsame_(equed, "R") ) {
				    rowcnd = 0.;
				    colcnd = 1.;
				} else if ( lsame_(equed, "C") ) {
				    rowcnd = 1.;
				    colcnd = 0.;
				} else if ( lsame_(equed, "B") ) {
				    rowcnd = 0.;
				    colcnd = 0.;
				}
			    }
			
			    /* Equilibrate the matrix. */
			    claqgs(&A, R, C, rowcnd, colcnd, amax, equed);
			}
		    }
		    
		    if ( prefact ) { /* Need a factor for the first time */
			
		        /* Save Fact option. */
		        fact = options.Fact;
			options.Fact = DOFACT;

			/* Preorder the matrix, obtain the column etree. */
			sp_preorder(&options, &A, perm_c, etree, &AC);

			/* Factor the matrix AC. */
			cgstrf(&options, &AC, relax, panel_size,
                               etree, work, lwork, perm_c, perm_r, &L, &U,
                               &stat, &info);

			if ( info ) { 
                            printf("** First factor: info %d, equed %c\n",
				   info, *equed);
                            if ( lwork == -1 ) {
                                printf("** Estimated memory: %d bytes\n",
                                        info - n);
                                exit(0);
                            }
                        }
	
                        Destroy_CompCol_Permuted(&AC);
			
		        /* Restore Fact option. */
			options.Fact = fact;
		    } /* if .. first time factor */
		    
		    for (itran = 0; itran < NTRAN; ++itran) {
			trans = transs[itran];
                        options.Trans = trans;

			/* Restore the matrix A. */
			cCopy_CompCol_Matrix(&ASAV, &A);
			
 			/* Set the right hand side. */
			cFillRHS(trans, nrhs, xact, ldx, &A, &B);
			cCopy_Dense_Matrix(m, nrhs, rhsb, ldb, bsav, ldb);

			/*----------------
			 * Test cgssv
			 *----------------*/
			if ( options.Fact == DOFACT && itran == 0) {
                            /* Not yet factored, and untransposed */
	
			    cCopy_Dense_Matrix(m, nrhs, rhsb, ldb, solx, ldx);
			    cgssv(&options, &A, perm_c, perm_r, &L, &U, &X,
                                  &stat, &info);
			    
			    if ( info && info != izero ) {
                                printf(FMT3, "cgssv",
				       info, izero, n, nrhs, imat, nfail);
			    } else {
                                /* Reconstruct matrix from factors and
	                           compute residual. */
                                cgst01(m, n, &A, &L, &U, perm_c, perm_r,
                                         &result[0]);
				nt = 1;
				if ( izero == 0 ) {
				    /* Compute residual of the computed
				       solution. */
				    cCopy_Dense_Matrix(m, nrhs, rhsb, ldb,
						       wwork, ldb);
				    cgst02(trans, m, n, nrhs, &A, solx,
                                              ldx, wwork,ldb, &result[1]);
				    nt = 2;
				}
				
				/* Print information about the tests that
				   did not pass the threshold.      */
				for (i = 0; i < nt; ++i) {
				    if ( result[i] >= THRESH ) {
					printf(FMT1, "cgssv", n, i,
					       result[i]);
					++nfail;
				    }
				}
				nrun += nt;
			    } /* else .. info == 0 */

			    /* Restore perm_c. */
			    for (i = 0; i < n; ++i) perm_c[i] = pc_save[i];

		            if (lwork == 0) {
			        Destroy_SuperNode_Matrix(&L);
			        Destroy_CompCol_Matrix(&U);
			    }
			} /* if .. end of testing cgssv */
    
			/*----------------
			 * Test cgssvx
			 *----------------*/
    
			/* Equilibrate the matrix if fact = FACTORED and
			   equed = 'R', 'C', or 'B'.   */
			if ( options.Fact == FACTORED &&
			     (equil || iequed) && n > 0 ) {
			    claqgs(&A, R, C, rowcnd, colcnd, amax, equed);
			}
			
			/* Solve the system and compute the condition number
			   and error bounds using cgssvx.      */
			cgssvx(&options, &A, perm_c, perm_r, etree,
                               equed, R, C, &L, &U, work, lwork, &B, &X, &rpg,
                               &rcond, ferr, berr, &mem_usage, &stat, &info);

			if ( info && info != izero ) {
			    printf(FMT3, "cgssvx",
				   info, izero, n, nrhs, imat, nfail);
                            if ( lwork == -1 ) {
                                printf("** Estimated memory: %.0f bytes\n",
                                        mem_usage.total_needed);
                                exit(0);
                            }
			} else {
			    if ( !prefact ) {
			    	/* Reconstruct matrix from factors and
	 			   compute residual. */
                                cgst01(m, n, &A, &L, &U, perm_c, perm_r,
                                         &result[0]);
				k1 = 0;
			    } else {
			   	k1 = 1;
			    }

			    if ( !info ) {
				/* Compute residual of the computed solution.*/
				cCopy_Dense_Matrix(m, nrhs, bsav, ldb,
						  wwork, ldb);
				cgst02(trans, m, n, nrhs, &ASAV, solx, ldx,
					  wwork, ldb, &result[1]);

				/* Check solution from generated exact
				   solution. */
				cgst04(n, nrhs, solx, ldx, xact, ldx, rcond,
					  &result[2]);

				/* Check the error bounds from iterative
				   refinement. */
				cgst07(trans, n, nrhs, &ASAV, bsav, ldb,
					  solx, ldx, xact, ldx, ferr, berr,
					  &result[3]);

				/* Print information about the tests that did
				   not pass the threshold.    */
				for (i = k1; i < NTESTS; ++i) {
				    if ( result[i] >= THRESH ) {
					printf(FMT2, "cgssvx",
					       options.Fact, trans, *equed,
					       n, imat, i, result[i]);
					++nfail;
				    }
				}
				nrun += NTESTS;
			    } /* if .. info == 0 */
			} /* else .. end of testing cgssvx */

		    } /* for itran ... */

		    if ( lwork == 0 ) {
			Destroy_SuperNode_Matrix(&L);
			Destroy_CompCol_Matrix(&U);
		    }

		} /* for equil ... */
	    } /* for ifact ... */
	} /* for iequed ... */
#if 0    
    if ( !info ) {
	PrintPerf(&L, &U, &mem_usage, rpg, rcond, ferr, berr, equed);
    }
#endif    

    } /* for imat ... */

    /* Print a summary of the results. */
    PrintSumm("CGE", nfail, nrun, nerrs);

    SUPERLU_FREE (rhsb);
    SUPERLU_FREE (bsav);
    SUPERLU_FREE (solx);    
    SUPERLU_FREE (xact);
    SUPERLU_FREE (etree);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    SUPERLU_FREE (pc_save);
    SUPERLU_FREE (R);
    SUPERLU_FREE (C);
    SUPERLU_FREE (ferr);
    SUPERLU_FREE (berr);
    SUPERLU_FREE (rwork);
    SUPERLU_FREE (wwork);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperMatrix_Store(&X);
    Destroy_CompCol_Matrix(&A);
    Destroy_CompCol_Matrix(&ASAV);
    if ( lwork > 0 ) {
	SUPERLU_FREE (work);
	Destroy_SuperMatrix_Store(&L);
	Destroy_SuperMatrix_Store(&U);
    }
    StatFree(&stat);

    return 0;
}

/*  
 * Parse command line options to get relaxed snode size, panel size, etc.
 */
static void
parse_command_line(int argc, char *argv[], char *matrix_type,
		   int *n, int *w, int *relax, int *nrhs, int *maxsuper,
		   int *rowblk, int *colblk, int *lwork, double *u)
{
    int c;
    extern char *optarg;

    while ( (c = getopt(argc, argv, "ht:n:w:r:s:m:b:c:l:")) != EOF ) {
	switch (c) {
	  case 'h':
	    printf("Options:\n");
	    printf("\t-w <int> - panel size\n");
	    printf("\t-r <int> - granularity of relaxed supernodes\n");
	    exit(1);
	    break;
	  case 't': strcpy(matrix_type, optarg);
	            break;
	  case 'n': *n = atoi(optarg);
	            break;
	  case 'w': *w = atoi(optarg);
	            break;
	  case 'r': *relax = atoi(optarg); 
	            break;
	  case 's': *nrhs = atoi(optarg); 
	            break;
	  case 'm': *maxsuper = atoi(optarg); 
	            break;
	  case 'b': *rowblk = atoi(optarg); 
	            break;
	  case 'c': *colblk = atoi(optarg); 
	            break;
	  case 'l': *lwork = atoi(optarg); 
	            break;
	  case 'u': *u = atof(optarg); 
	            break;
  	}
    }
}
