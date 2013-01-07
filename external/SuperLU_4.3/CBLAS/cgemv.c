
/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int cgemv_(char *trans, integer *m, integer *n, complex *
	alpha, complex *a, integer *lda, complex *x, integer *incx, complex *
	beta, complex *y, integer *incy)
{


    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    complex q__1, q__2, q__3;

    /* Builtin functions */
    void r_cnjg(complex *, complex *);

    /* Local variables */
    static integer info;
    static complex temp;
    static integer lenx, leny, i, j;
    extern logical lsame_(char *, char *);
    static integer ix, iy, jx, jy, kx, ky;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static logical noconj;


/*  Purpose   
    =======   

    CGEMV  performs one of the matrix-vector operations   

       y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   or   

       y := alpha*conjg( A' )*x + beta*y,   

    where alpha and beta are scalars, x and y are vectors and A is an   
    m by n matrix.   

    Parameters   
    ==========   

    TRANS  - CHARACTER*1.   
             On entry, TRANS specifies the operation to be performed as   
             follows:   

                TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.   

                TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.   

                TRANS = 'C' or 'c'   y := alpha*conjg( A' )*x + beta*y.   

             Unchanged on exit.   

    M      - INTEGER.   
             On entry, M specifies the number of rows of the matrix A.   
             M must be at least zero.   
             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the number of columns of the matrix A. 
  
             N must be at least zero.   
             Unchanged on exit.   

    ALPHA  - COMPLEX         .   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    A      - COMPLEX          array of DIMENSION ( LDA, n ).   
             Before entry, the leading m by n part of the array A must   
             contain the matrix of coefficients.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             max( 1, m ).   
             Unchanged on exit.   

    X      - COMPLEX          array of DIMENSION at least   
             ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'   
             and at least   
             ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.   
             Before entry, the incremented array X must contain the   
             vector x.   
             Unchanged on exit.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   

    BETA   - COMPLEX         .   
             On entry, BETA specifies the scalar beta. When BETA is   
             supplied as zero then Y need not be set on input.   
             Unchanged on exit.   

    Y      - COMPLEX          array of DIMENSION at least   
             ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'   
             and at least   
             ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.   
             Before entry with BETA non-zero, the incremented array Y   
             must contain the vector y. On exit, Y is overwritten by the 
  
             updated vector y.   

    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   


    Level 2 Blas routine.   

    -- Written on 22-October-1986.   
       Jack Dongarra, Argonne National Lab.   
       Jeremy Du Croz, Nag Central Office.   
       Sven Hammarling, Nag Central Office.   
       Richard Hanson, Sandia National Labs.   



       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]
#define Y(I) y[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! lsame_(trans, "N") && ! lsame_(trans, "T") && ! 
	    lsame_(trans, "C")) {
	info = 1;
    } else if (*m < 0) {
	info = 2;
    } else if (*n < 0) {
	info = 3;
    } else if (*lda < max(1,*m)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    } else if (*incy == 0) {
	info = 11;
    }
    if (info != 0) {
	xerbla_("CGEMV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || alpha->r == 0.f && alpha->i == 0.f && (beta->r 
	    == 1.f && beta->i == 0.f)) {
	return 0;
    }

    noconj = lsame_(trans, "T");

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set 
  
       up the start points in  X  and  Y. */

    if (lsame_(trans, "N")) {
	lenx = *n;
	leny = *m;
    } else {
	lenx = *m;
	leny = *n;
    }
    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (lenx - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (leny - 1) * *incy;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through A.   

       First form  y := beta*y. */

    if (beta->r != 1.f || beta->i != 0.f) {
	if (*incy == 1) {
	    if (beta->r == 0.f && beta->i == 0.f) {
		i__1 = leny;
		for (i = 1; i <= leny; ++i) {
		    i__2 = i;
		    Y(i).r = 0.f, Y(i).i = 0.f;
/* L10: */
		}
	    } else {
		i__1 = leny;
		for (i = 1; i <= leny; ++i) {
		    i__2 = i;
		    i__3 = i;
		    q__1.r = beta->r * Y(i).r - beta->i * Y(i).i, 
			    q__1.i = beta->r * Y(i).i + beta->i * Y(i)
			    .r;
		    Y(i).r = q__1.r, Y(i).i = q__1.i;
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (beta->r == 0.f && beta->i == 0.f) {
		i__1 = leny;
		for (i = 1; i <= leny; ++i) {
		    i__2 = iy;
		    Y(iy).r = 0.f, Y(iy).i = 0.f;
		    iy += *incy;
/* L30: */
		}
	    } else {
		i__1 = leny;
		for (i = 1; i <= leny; ++i) {
		    i__2 = iy;
		    i__3 = iy;
		    q__1.r = beta->r * Y(iy).r - beta->i * Y(iy).i, 
			    q__1.i = beta->r * Y(iy).i + beta->i * Y(iy)
			    .r;
		    Y(iy).r = q__1.r, Y(iy).i = q__1.i;
		    iy += *incy;
/* L40: */
		}
	    }
	}
    }
    if (alpha->r == 0.f && alpha->i == 0.f) {
	return 0;
    }
    if (lsame_(trans, "N")) {

/*        Form  y := alpha*A*x + y. */

	jx = kx;
	if (*incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = jx;
		if (X(jx).r != 0.f || X(jx).i != 0.f) {
		    i__2 = jx;
		    q__1.r = alpha->r * X(jx).r - alpha->i * X(jx).i, 
			    q__1.i = alpha->r * X(jx).i + alpha->i * X(jx)
			    .r;
		    temp.r = q__1.r, temp.i = q__1.i;
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			i__3 = i;
			i__4 = i;
			i__5 = i + j * a_dim1;
			q__2.r = temp.r * A(i,j).r - temp.i * A(i,j).i, 
				q__2.i = temp.r * A(i,j).i + temp.i * A(i,j)
				.r;
			q__1.r = Y(i).r + q__2.r, q__1.i = Y(i).i + 
				q__2.i;
			Y(i).r = q__1.r, Y(i).i = q__1.i;
/* L50: */
		    }
		}
		jx += *incx;
/* L60: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = jx;
		if (X(jx).r != 0.f || X(jx).i != 0.f) {
		    i__2 = jx;
		    q__1.r = alpha->r * X(jx).r - alpha->i * X(jx).i, 
			    q__1.i = alpha->r * X(jx).i + alpha->i * X(jx)
			    .r;
		    temp.r = q__1.r, temp.i = q__1.i;
		    iy = ky;
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			i__3 = iy;
			i__4 = iy;
			i__5 = i + j * a_dim1;
			q__2.r = temp.r * A(i,j).r - temp.i * A(i,j).i, 
				q__2.i = temp.r * A(i,j).i + temp.i * A(i,j)
				.r;
			q__1.r = Y(iy).r + q__2.r, q__1.i = Y(iy).i + 
				q__2.i;
			Y(iy).r = q__1.r, Y(iy).i = q__1.i;
			iy += *incy;
/* L70: */
		    }
		}
		jx += *incx;
/* L80: */
	    }
	}
    } else {

/*        Form  y := alpha*A'*x + y  or  y := alpha*conjg( A' )*x + y.
 */

	jy = ky;
	if (*incx == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		temp.r = 0.f, temp.i = 0.f;
		if (noconj) {
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			i__3 = i + j * a_dim1;
			i__4 = i;
			q__2.r = A(i,j).r * X(i).r - A(i,j).i * X(i)
				.i, q__2.i = A(i,j).r * X(i).i + A(i,j)
				.i * X(i).r;
			q__1.r = temp.r + q__2.r, q__1.i = temp.i + q__2.i;
			temp.r = q__1.r, temp.i = q__1.i;
/* L90: */
		    }
		} else {
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			r_cnjg(&q__3, &A(i,j));
			i__3 = i;
			q__2.r = q__3.r * X(i).r - q__3.i * X(i).i, 
				q__2.i = q__3.r * X(i).i + q__3.i * X(i)
				.r;
			q__1.r = temp.r + q__2.r, q__1.i = temp.i + q__2.i;
			temp.r = q__1.r, temp.i = q__1.i;
/* L100: */
		    }
		}
		i__2 = jy;
		i__3 = jy;
		q__2.r = alpha->r * temp.r - alpha->i * temp.i, q__2.i = 
			alpha->r * temp.i + alpha->i * temp.r;
		q__1.r = Y(jy).r + q__2.r, q__1.i = Y(jy).i + q__2.i;
		Y(jy).r = q__1.r, Y(jy).i = q__1.i;
		jy += *incy;
/* L110: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		temp.r = 0.f, temp.i = 0.f;
		ix = kx;
		if (noconj) {
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			i__3 = i + j * a_dim1;
			i__4 = ix;
			q__2.r = A(i,j).r * X(ix).r - A(i,j).i * X(ix)
				.i, q__2.i = A(i,j).r * X(ix).i + A(i,j)
				.i * X(ix).r;
			q__1.r = temp.r + q__2.r, q__1.i = temp.i + q__2.i;
			temp.r = q__1.r, temp.i = q__1.i;
			ix += *incx;
/* L120: */
		    }
		} else {
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			r_cnjg(&q__3, &A(i,j));
			i__3 = ix;
			q__2.r = q__3.r * X(ix).r - q__3.i * X(ix).i, 
				q__2.i = q__3.r * X(ix).i + q__3.i * X(ix)
				.r;
			q__1.r = temp.r + q__2.r, q__1.i = temp.i + q__2.i;
			temp.r = q__1.r, temp.i = q__1.i;
			ix += *incx;
/* L130: */
		    }
		}
		i__2 = jy;
		i__3 = jy;
		q__2.r = alpha->r * temp.r - alpha->i * temp.i, q__2.i = 
			alpha->r * temp.i + alpha->i * temp.r;
		q__1.r = Y(jy).r + q__2.r, q__1.i = Y(jy).i + q__2.i;
		Y(jy).r = q__1.r, Y(jy).i = q__1.i;
		jy += *incy;
/* L140: */
	    }
	}
    }

    return 0;

/*     End of CGEMV . */

} /* cgemv_ */

