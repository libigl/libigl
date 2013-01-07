
/*! @file cfgmr.c
 * \brief flexible GMRES from ITSOL developed by Yousef Saad.
 */

/* ITSOL COPYRIGHT

Copyright (C) 2006, the University of Minnesota 

ITSOL is free software; you can redistribute it and/or modify it under
the terms of  the GNU General Public License as  published by the Free
Software Foundation [version 2 of  the License, or any later version]
For details, see 

http://www.gnu.org/copyleft/gpl.html

A copy of the GNU licencing agreement is attached to the ITSOL package
in the file GNU.  For additional information contact the Free Software
Foundation Inc., 65 Mass Ave, Cambridge, MA 02139, USA.

DISCLAIMER
----------

This program  is distributed in the  hope that it will  be useful, but
WITHOUT   ANY  WARRANTY;   without  even   the  implied   warranty  of
MERCHANTABILITY  or FITNESS  FOR A  PARTICULAR PURPOSE.   See  the GNU
General Public License for more details. 

For information on ITSOL contact saad@cs.umn.edu
*/

#include "slu_cdefs.h"

#define  epsmac  1.0e-16

extern void cdotc_(complex *, int *, complex [], int *, complex [], int *);
extern float scnrm2_(int *, complex [], int *);


int cfgmr(int n,
     void (*cmatvec) (complex, complex[], complex, complex[]),
     void (*cpsolve) (int, complex[], complex[]),
     complex *rhs, complex *sol, double tol, int im, int *itmax, FILE * fits)
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
    int i, i1, ii, j, k, k1, its, retval, i_1 = 1, i_2 = 2;
    float beta, eps1 = 0.0, t, t0, gam;
    complex **hh, *c, *s, *rs;
    complex **vv, **z, tt;
    complex zero = {0.0, 0.0};
    complex one = {1.0, 0.0};
    complex tt1, tt2;

    its = 0;
    vv = (complex **)SUPERLU_MALLOC((im + 1) * sizeof(complex *));
    for (i = 0; i <= im; i++) vv[i] = complexMalloc(n);
    z = (complex **)SUPERLU_MALLOC(im * sizeof(complex *));
    hh = (complex **)SUPERLU_MALLOC(im * sizeof(complex *));
    for (i = 0; i < im; i++)
    {
	hh[i] = complexMalloc(i + 2);
	z[i] = complexMalloc(n);
    }
    c = complexMalloc(im);
    s = complexMalloc(im);
    rs = complexMalloc(im + 1);

    /*---- outer loop starts here ----*/
    do
    {
	/*---- compute initial residual vector ----*/
	cmatvec(one, sol, zero, vv[0]);
	for (j = 0; j < n; j++)
	    c_sub(&vv[0][j], &rhs[j], &vv[0][j]);	/* vv[0]= initial residual */
	beta = scnrm2_(&n, vv[0], &i_1);

	/*---- print info if fits != null ----*/
	if (fits != NULL && its == 0)
	    fprintf(fits, "%8d   %10.2e\n", its, beta);
	/*if ( beta <= tol * dnrm2_(&n, rhs, &i_1) )*/
	if ( !(beta > tol * scnrm2_(&n, rhs, &i_1)) )
	    break;
	t = 1.0 / beta;

	/*---- normalize: vv[0] = vv[0] / beta ----*/
	for (j = 0; j < n; j++)
	    cs_mult(&vv[0][j], &vv[0][j], t);
	if (its == 0)
	    eps1 = tol * beta;

	/*---- initialize 1-st term of rhs of hessenberg system ----*/
	rs[0].r = beta;
	rs[0].i = 0.0;
	for (i = 0; i < im; i++)
	{
	    its++;
	    i1 = i + 1;

	    /*------------------------------------------------------------
	    |  (Right) Preconditioning Operation   z_{j} = M^{-1} v_{j}
	    +-----------------------------------------------------------*/
	    if (cpsolve)
		cpsolve(n, z[i], vv[i]);
	    else
		ccopy_(&n, vv[i], &i_1, z[i], &i_1);

	    /*---- matvec operation w = A z_{j} = A M^{-1} v_{j} ----*/
	    cmatvec(one, z[i], zero, vv[i1]);

	    /*------------------------------------------------------------
	    |     modified gram - schmidt...
	    |     h_{i,j} = (w,v_{i})
	    |     w  = w - h_{i,j} v_{i}
	    +------------------------------------------------------------*/
	    t0 = scnrm2_(&n, vv[i1], &i_1);
	    for (j = 0; j <= i; j++)
	    {
		complex negt;
#if 0
		cdotc_(&tt, &n, vv[j], &i_1, vv[i1], &i_1);
#else
		tt = zero;
		for (k = 0; k < n; ++k) {
		    cc_conj(&tt1, &vv[j][k]);
		    cc_mult(&tt2, &tt1, &vv[i1][k]);
		    c_add(&tt, &tt, &tt2);
		}
#endif
		hh[i][j] = tt;
		negt.r = -tt.r;
		negt.i = -tt.i;
		caxpy_(&n, &negt, vv[j], &i_1, vv[i1], &i_1);
	    }

	    /*---- h_{j+1,j} = ||w||_{2} ----*/
	    t = scnrm2_(&n, vv[i1], &i_1);
	    while (t < 0.5 * t0)
	    {
		t0 = t;
		for (j = 0; j <= i; j++)
		{
		    complex negt;
#if 0
		    cdotc_(&tt, &n, vv[j], &i_1, vv[i1], &i_1);
#else
   	            tt = zero;
		    for (k = 0; k < n; ++k) {
		        cc_conj(&tt1, &vv[j][k]);
		        cc_mult(&tt2, &tt1, &vv[i1][k]);
		        c_add(&tt, &tt, &tt2);
		    }
#endif
		    c_add(&hh[i][j], &hh[i][j], &tt);
		    negt.r = -tt.r;
		    negt.i = -tt.i;
		    caxpy_(&n, &negt, vv[j], &i_1, vv[i1], &i_1);
		}
		t = scnrm2_(&n, vv[i1], &i_1);
	    }

	    hh[i][i1].r = t;
	    hh[i][i1].i = 0.0;

	    if (t != 0.0)
	    {
		/*---- v_{j+1} = w / h_{j+1,j} ----*/
		t = 1.0 / t;
		for (k = 0; k < n; k++)
	            cs_mult(&vv[i1][k], &vv[i1][k], t);
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
		tt = hh[i][k1];
                cc_mult(&tt1, &c[k1], &tt);
                cc_mult(&tt2, &s[k1], &hh[i][k]);
                c_add(&hh[i][k1], &tt1, &tt2);

                cc_mult(&tt1, &s[k1], &tt);
                cc_mult(&tt2, &c[k1], &hh[i][k]);
                c_sub(&hh[i][k], &tt2, &tt1);
	    }

	    gam = scnrm2_(&i_2, &hh[i][i], &i_1);

	    /*---------------------------------------------------
	    |     if gamma is zero then any small value will do
	    |     affect only residual estimate
	    +--------------------------------------------------*/
	    /* if (gam == 0.0) gam = epsmac; */

	    /*---- get next plane rotation ---*/
	    if (gam == 0.0)
	    {
		c[i] = one;
		s[i] = zero;
	    }
            else
	    {
                gam = 1.0 / gam;
		cs_mult(&c[i], &hh[i][i], gam);
		cs_mult(&s[i], &hh[i][i1], gam);
	    }

	    cc_mult(&rs[i1], &s[i], &rs[i]);
            rs[i1].r = -rs[i1].r;  rs[i1].i = -rs[i1].i;
	    cc_mult(&rs[i], &c[i], &rs[i]);

	    /*----------------------------------------------------
	    |   determine residual norm and test for convergence
	    +---------------------------------------------------*/
            cc_mult(&tt1, &c[i], &hh[i][i]);
            cc_mult(&tt2, &s[i], &hh[i][i1]);
            c_add(&hh[i][i], &tt1, &tt2);
            beta = c_abs1(&rs[i1]);
	    if (fits != NULL)
		fprintf(fits, "%8d   %10.2e\n", its, beta);
	    if (beta <= eps1 || its >= maxits)
		break;
	}

	if (i == im) i--;

	/*---- now compute solution. 1st, solve upper triangular system ----*/
	c_div(&rs[i], &rs[i], &hh[i][i]);

	for (ii = 1; ii <= i; ii++)
	{
	    k = i - ii;
	    k1 = k + 1;
	    tt = rs[k];
	    for (j = k1; j <= i; j++) {
                cc_mult(&tt1, &hh[j][k], &rs[j]);
		c_sub(&tt, &tt, &tt1);
            }
            c_div(&rs[k], &tt, &hh[k][k]);
	}

	/*---- linear combination of v[i]'s to get sol. ----*/
	for (j = 0; j <= i; j++)
	{
	    tt = rs[j];
	    for (k = 0; k < n; k++) {
                cc_mult(&tt1, &tt, &z[j][k]);
		c_add(&sol[k], &sol[k], &tt1);
            }
	}

	/* calculate the residual and output */
	cmatvec(one, sol, zero, vv[0]);
	for (j = 0; j < n; j++)
	    c_sub(&vv[0][j], &rhs[j], &vv[0][j]);/* vv[0]= initial residual */

	/*---- print info if fits != null ----*/
	beta = scnrm2_(&n, vv[0], &i_1);

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
