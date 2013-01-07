#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "slu_ddefs.h"

#define  epsmac  1.0e-16
extern double ddot_(int *, double [], int *, double [], int *);
extern double dnrm2_(int *, double [], int *);
extern void daxpy_(int *, double *, double [], int *, double [], int *);
extern double dcopy_(int *, double [], int *, double [], int *);

int fgmr(int n,
     void (*matvec) (double, double[], double, double[]),
     void (*psolve) (int, double[], double[]),
     double *rhs, double *sol, double tol, int im, int *itmax, FILE * fits)
{
/*----------------------------------------------------------------------
|                 *** Preconditioned FGMRES ***
+-----------------------------------------------------------------------
| This is a simple version of the ARMS preconditioned FGMRES algorithm.
+-----------------------------------------------------------------------
| Y. S. Dec. 2000. -- Apr. 2008
+-----------------------------------------------------------------------
| on entry:
|----------
|
| rhs     = real vector of length n containing the right hand side.
| sol     = real vector of length n containing an initial guess to the
|           solution on input.
| tol     = tolerance for stopping iteration
| im      = Krylov subspace dimension
| (itmax) = max number of iterations allowed.
| fits    = NULL: no output
|        != NULL: file handle to output " resid vs time and its"
|
| on return:
|----------
| fgmr      int =  0 --> successful return.
|           int =  1 --> convergence not achieved in itmax iterations.
| sol     = contains an approximate solution (upon successful return).
| itmax   = has changed. It now contains the number of steps required
|           to converge --
+-----------------------------------------------------------------------
| internal work arrays:
|----------
| vv      = work array of length [im+1][n] (used to store the Arnoldi
|           basis)
| hh      = work array of length [im][im+1] (Householder matrix)
| z       = work array of length [im][n] to store preconditioned vectors
+-----------------------------------------------------------------------
| subroutines called :
| matvec - matrix-vector multiplication operation
| psolve - (right) preconditionning operation
|	   psolve can be a NULL pointer (GMRES without preconditioner)
+---------------------------------------------------------------------*/

    int maxits = *itmax;
    int i, i1, ii, j, k, k1, its, retval, i_1 = 1;
    double **hh, *c, *s, *rs, t, t0;
    double beta, eps1 = 0.0, gam, **vv, **z;

    its = 0;
    vv = (double **)SUPERLU_MALLOC((im + 1) * sizeof(double *));
    for (i = 0; i <= im; i++)
	vv[i] = doubleMalloc(n);
    z = (double **)SUPERLU_MALLOC(im * sizeof(double *));
    hh = (double **)SUPERLU_MALLOC(im * sizeof(double *));
    for (i = 0; i < im; i++)
    {
	hh[i] = doubleMalloc(i + 2);
	z[i] = doubleMalloc(n);
    }
    c = doubleMalloc(im);
    s = doubleMalloc(im);
    rs = doubleMalloc(im + 1);

    /*---- outer loop starts here ----*/
    do
    {
	/*---- compute initial residual vector ----*/
	matvec(1.0, sol, 0.0, vv[0]);
	for (j = 0; j < n; j++)
	    vv[0][j] = rhs[j] - vv[0][j];	/* vv[0]= initial residual */
	beta = dnrm2_(&n, vv[0], &i_1);

	/*---- print info if fits != null ----*/
	if (fits != NULL && its == 0)
	    fprintf(fits, "%8d   %10.2e\n", its, beta);
	/*if ( beta < tol * dnrm2_(&n, rhs, &i_1) )*/
	if ( !(beta >= tol * dnrm2_(&n, rhs, &i_1)) )
	    break;
	t = 1.0 / beta;

	/*---- normalize: vv[0] = vv[0] / beta ----*/
	for (j = 0; j < n; j++)
	    vv[0][j] = vv[0][j] * t;
	if (its == 0)
	    eps1 = tol * beta;

	/*---- initialize 1-st term of rhs of hessenberg system ----*/
	rs[0] = beta;
	for (i = 0; i < im; i++)
	{
	    its++;
	    i1 = i + 1;

	    /*------------------------------------------------------------
	    |  (Right) Preconditioning Operation   z_{j} = M^{-1} v_{j}
	    +-----------------------------------------------------------*/
	    if (psolve)
		psolve(n, z[i], vv[i]);
	    else
		dcopy_(&n, vv[i], &i_1, z[i], &i_1);

	    /*---- matvec operation w = A z_{j} = A M^{-1} v_{j} ----*/
	    matvec(1.0, z[i], 0.0, vv[i1]);

	    /*------------------------------------------------------------
	    |     modified gram - schmidt...
	    |     h_{i,j} = (w,v_{i})
	    |     w  = w - h_{i,j} v_{i}
	    +------------------------------------------------------------*/
	    t0 = dnrm2_(&n, vv[i1], &i_1);
	    for (j = 0; j <= i; j++)
	    {
		double negt;
		t = ddot_(&n, vv[j], &i_1, vv[i1], &i_1);
		hh[i][j] = t;
		negt = -t;
		daxpy_(&n, &negt, vv[j], &i_1, vv[i1], &i_1);
	    }

	    /*---- h_{j+1,j} = ||w||_{2} ----*/
	    t = dnrm2_(&n, vv[i1], &i_1);
	    while (t < 0.5 * t0)
	    {
		t0 = t;
		for (j = 0; j <= i; j++)
		{
		    double negt;
		    t = ddot_(&n, vv[j], &i_1, vv[i1], &i_1);
		    hh[i][j] += t;
		    negt = -t;
		    daxpy_(&n, &negt, vv[j], &i_1, vv[i1], &i_1);
		}
		t = dnrm2_(&n, vv[i1], &i_1);
	    }
	    hh[i][i1] = t;
	    if (t != 0.0)
	    {
		/*---- v_{j+1} = w / h_{j+1,j} ----*/
		t = 1.0 / t;
		for (k = 0; k < n; k++)
		    vv[i1][k] = vv[i1][k] * t;
	    }
	    /*---------------------------------------------------
	    |     done with modified gram schimdt and arnoldi step
	    |     now  update factorization of hh
	    +--------------------------------------------------*/

	    /*--------------------------------------------------------
	    |   perform previous transformations  on i-th column of h
	    +-------------------------------------------------------*/
	    for (k = 1; k <= i; k++)
	    {
		k1 = k - 1;
		t = hh[i][k1];
		hh[i][k1] = c[k1] * t + s[k1] * hh[i][k];
		hh[i][k] = -s[k1] * t + c[k1] * hh[i][k];
	    }
	    gam = sqrt(pow(hh[i][i], 2) + pow(hh[i][i1], 2));

	    /*---------------------------------------------------
	    |     if gamma is zero then any small value will do
	    |     affect only residual estimate
	    +--------------------------------------------------*/
	    /* if (gam == 0.0) gam = epsmac; */

	    /*---- get next plane rotation ---*/
	    if (gam > 0.0)
	    {
		c[i] = hh[i][i] / gam;
		s[i] = hh[i][i1] / gam;
	    }
	    else
	    {
		c[i] = 1.0;
		s[i] = 0.0;
	    }
	    rs[i1] = -s[i] * rs[i];
	    rs[i] = c[i] * rs[i];

	    /*----------------------------------------------------
	    |   determine residual norm and test for convergence
	    +---------------------------------------------------*/
	    hh[i][i] = c[i] * hh[i][i] + s[i] * hh[i][i1];
	    beta = fabs(rs[i1]);
	    if (fits != NULL)
		fprintf(fits, "%8d   %10.2e\n", its, beta);
	    if (beta <= eps1 || its >= maxits)
		break;
	}

	if (i == im) i--;
	/*---- now compute solution. 1st, solve upper triangular system ----*/
	rs[i] = rs[i] / hh[i][i];
	for (ii = 1; ii <= i; ii++)
	{
	    k = i - ii;
	    k1 = k + 1;
	    t = rs[k];
	    for (j = k1; j <= i; j++)
		t = t - hh[j][k] * rs[j];
	    rs[k] = t / hh[k][k];
	}

	/*---- linear combination of v[i]'s to get sol. ----*/
	for (j = 0; j <= i; j++)
	{
	    t = rs[j];
	    for (k = 0; k < n; k++)
		sol[k] += t * z[j][k];
	}

	/* calculate the residual and output */
	matvec(1.0, sol, 0.0, vv[0]);
	for (j = 0; j < n; j++)
	    vv[0][j] = rhs[j] - vv[0][j];	/* vv[0]= initial residual */

	/*---- print info if fits != null ----*/
	beta = dnrm2_(&n, vv[0], &i_1);

	/*---- restart outer loop if needed ----*/
	/*if (beta >= eps1 / tol)*/
	if ( !(beta < eps1 / tol) )
	{
	    its = maxits + 10;
	    break;
	}
	if (beta <= eps1)
	    break;
    } while(its < maxits);

    retval = (its >= maxits);
    for (i = 0; i <= im; i++)
	SUPERLU_FREE(vv[i]);
    SUPERLU_FREE(vv);
    for (i = 0; i < im; i++)
    {
	SUPERLU_FREE(hh[i]);
	SUPERLU_FREE(z[i]);
    }
    SUPERLU_FREE(hh);
    SUPERLU_FREE(z);
    SUPERLU_FREE(c);
    SUPERLU_FREE(s);
    SUPERLU_FREE(rs);

    *itmax = its;

    return retval;
} /*----end of fgmr ----*/
