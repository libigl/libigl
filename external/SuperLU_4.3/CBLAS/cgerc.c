
/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int cgerc_(integer *m, integer *n, complex *alpha, complex *
	x, integer *incx, complex *y, integer *incy, complex *a, integer *lda)
{


    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    complex q__1, q__2;

    /* Builtin functions */
    void r_cnjg(complex *, complex *);

    /* Local variables */
    static integer info;
    static complex temp;
    static integer i, j, ix, jy, kx;
    extern /* Subroutine */ int xerbla_(char *, integer *);


/*  Purpose   
    =======   

    CGERC  performs the rank 1 operation   

       A := alpha*x*conjg( y' ) + A,   

    where alpha is a scalar, x is an m element vector, y is an n element 
  
    vector and A is an m by n matrix.   

    Parameters   
    ==========   

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

    X      - COMPLEX          array of dimension at least   
             ( 1 + ( m - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the m   
             element vector x.   
             Unchanged on exit.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   

    Y      - COMPLEX          array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCY ) ).   
             Before entry, the incremented array Y must contain the n   
             element vector y.   
             Unchanged on exit.   

    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   

    A      - COMPLEX          array of DIMENSION ( LDA, n ).   
             Before entry, the leading m by n part of the array A must   
             contain the matrix of coefficients. On exit, A is   
             overwritten by the updated matrix.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             max( 1, m ).   
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
    if (*m < 0) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    } else if (*incy == 0) {
	info = 7;
    } else if (*lda < max(1,*m)) {
	info = 9;
    }
    if (info != 0) {
	xerbla_("CGERC ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || alpha->r == 0.f && alpha->i == 0.f) {
	return 0;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through A. */

    if (*incy > 0) {
	jy = 1;
    } else {
	jy = 1 - (*n - 1) * *incy;
    }
    if (*incx == 1) {
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    i__2 = jy;
	    if (Y(jy).r != 0.f || Y(jy).i != 0.f) {
		r_cnjg(&q__2, &Y(jy));
		q__1.r = alpha->r * q__2.r - alpha->i * q__2.i, q__1.i = 
			alpha->r * q__2.i + alpha->i * q__2.r;
		temp.r = q__1.r, temp.i = q__1.i;
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    i__3 = i + j * a_dim1;
		    i__4 = i + j * a_dim1;
		    i__5 = i;
		    q__2.r = X(i).r * temp.r - X(i).i * temp.i, q__2.i =
			     X(i).r * temp.i + X(i).i * temp.r;
		    q__1.r = A(i,j).r + q__2.r, q__1.i = A(i,j).i + q__2.i;
		    A(i,j).r = q__1.r, A(i,j).i = q__1.i;
/* L10: */
		}
	    }
	    jy += *incy;
/* L20: */
	}
    } else {
	if (*incx > 0) {
	    kx = 1;
	} else {
	    kx = 1 - (*m - 1) * *incx;
	}
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    i__2 = jy;
	    if (Y(jy).r != 0.f || Y(jy).i != 0.f) {
		r_cnjg(&q__2, &Y(jy));
		q__1.r = alpha->r * q__2.r - alpha->i * q__2.i, q__1.i = 
			alpha->r * q__2.i + alpha->i * q__2.r;
		temp.r = q__1.r, temp.i = q__1.i;
		ix = kx;
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    i__3 = i + j * a_dim1;
		    i__4 = i + j * a_dim1;
		    i__5 = ix;
		    q__2.r = X(ix).r * temp.r - X(ix).i * temp.i, q__2.i =
			     X(ix).r * temp.i + X(ix).i * temp.r;
		    q__1.r = A(i,j).r + q__2.r, q__1.i = A(i,j).i + q__2.i;
		    A(i,j).r = q__1.r, A(i,j).i = q__1.i;
		    ix += *incx;
/* L30: */
		}
	    }
	    jy += *incy;
/* L40: */
	}
    }

    return 0;

/*     End of CGERC . */

} /* cgerc_ */

