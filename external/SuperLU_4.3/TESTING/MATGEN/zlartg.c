#include "f2c.h"

/* Subroutine */ int zlartg_(doublecomplex *f, doublecomplex *g, doublereal *
	cs, doublecomplex *sn, doublecomplex *r)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZLARTG generates a plane rotation so that   

       [  CS  SN  ]     [ F ]     [ R ]   
       [  __      ]  .  [   ]  =  [   ]   where CS**2 + |SN|**2 = 1.   
       [ -SN  CS  ]     [ G ]     [ 0 ]   

    This is a faster version of the BLAS1 routine ZROTG, except for   
    the following differences:   
       F and G are unchanged on return.   
       If G=0, then CS=1 and SN=0.   
       If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any   
          floating point operations.   

    Arguments   
    =========   

    F       (input) COMPLEX*16   
            The first component of vector to be rotated.   

    G       (input) COMPLEX*16   
            The second component of vector to be rotated.   

    CS      (output) DOUBLE PRECISION   
            The cosine of the rotation.   

    SN      (output) COMPLEX*16   
            The sine of the rotation.   

    R       (output) COMPLEX*16   
            The nonzero component of the rotated vector.   

    ===================================================================== 
  


       [ 25 or 38 ops for main paths ] */
    /* System generated locals */
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3;
    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);
    double z_abs(doublecomplex *), d_imag(doublecomplex *), sqrt(doublereal);
    /* Local variables */
    static doublereal d, f1, f2, g1, g2, fa, ga, di;
    static doublecomplex fs, gs, ss;


    if (g->r == 0. && g->i == 0.) {
	*cs = 1.;
	sn->r = 0., sn->i = 0.;
	r->r = f->r, r->i = f->i;
    } else if (f->r == 0. && f->i == 0.) {
	*cs = 0.;

	d_cnjg(&z__2, g);
	d__1 = z_abs(g);
	z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
	sn->r = z__1.r, sn->i = z__1.i;
	d__1 = z_abs(g);
	r->r = d__1, r->i = 0.;

/*         SN = ONE   
           R = G */

    } else {
	f1 = (d__1 = f->r, abs(d__1)) + (d__2 = d_imag(f), abs(d__2));
	g1 = (d__1 = g->r, abs(d__1)) + (d__2 = d_imag(g), abs(d__2));
	if (f1 >= g1) {
	    z__1.r = g->r / f1, z__1.i = g->i / f1;
	    gs.r = z__1.r, gs.i = z__1.i;
/* Computing 2nd power */
	    d__1 = gs.r;
/* Computing 2nd power */
	    d__2 = d_imag(&gs);
	    g2 = d__1 * d__1 + d__2 * d__2;
	    z__1.r = f->r / f1, z__1.i = f->i / f1;
	    fs.r = z__1.r, fs.i = z__1.i;
/* Computing 2nd power */
	    d__1 = fs.r;
/* Computing 2nd power */
	    d__2 = d_imag(&fs);
	    f2 = d__1 * d__1 + d__2 * d__2;
	    d = sqrt(g2 / f2 + 1.);
	    *cs = 1. / d;
	    d_cnjg(&z__3, &gs);
	    z__2.r = z__3.r * fs.r - z__3.i * fs.i, z__2.i = z__3.r * fs.i + 
		    z__3.i * fs.r;
	    d__1 = *cs / f2;
	    z__1.r = d__1 * z__2.r, z__1.i = d__1 * z__2.i;
	    sn->r = z__1.r, sn->i = z__1.i;
	    z__1.r = d * f->r, z__1.i = d * f->i;
	    r->r = z__1.r, r->i = z__1.i;
	} else {
	    z__1.r = f->r / g1, z__1.i = f->i / g1;
	    fs.r = z__1.r, fs.i = z__1.i;
/* Computing 2nd power */
	    d__1 = fs.r;
/* Computing 2nd power */
	    d__2 = d_imag(&fs);
	    f2 = d__1 * d__1 + d__2 * d__2;
	    fa = sqrt(f2);
	    z__1.r = g->r / g1, z__1.i = g->i / g1;
	    gs.r = z__1.r, gs.i = z__1.i;
/* Computing 2nd power */
	    d__1 = gs.r;
/* Computing 2nd power */
	    d__2 = d_imag(&gs);
	    g2 = d__1 * d__1 + d__2 * d__2;
	    ga = sqrt(g2);
	    d = sqrt(f2 / g2 + 1.);
	    di = 1. / d;
	    *cs = fa / ga * di;
	    d_cnjg(&z__3, &gs);
	    z__2.r = z__3.r * fs.r - z__3.i * fs.i, z__2.i = z__3.r * fs.i + 
		    z__3.i * fs.r;
	    d__1 = fa * ga;
	    z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
	    ss.r = z__1.r, ss.i = z__1.i;
	    z__1.r = di * ss.r, z__1.i = di * ss.i;
	    sn->r = z__1.r, sn->i = z__1.i;
	    z__2.r = g->r * ss.r - g->i * ss.i, z__2.i = g->r * ss.i + g->i * 
		    ss.r;
	    z__1.r = d * z__2.r, z__1.i = d * z__2.i;
	    r->r = z__1.r, r->i = z__1.i;
	}
    }
    return 0;

/*     End of ZLARTG */

} /* zlartg_ */

