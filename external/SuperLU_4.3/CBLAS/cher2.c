
/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int cher2_(char *uplo, integer *n, complex *alpha, complex *
	x, integer *incx, complex *y, integer *incy, complex *a, integer *lda)
{


    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1;
    complex q__1, q__2, q__3, q__4;

    /* Builtin functions */
    void r_cnjg(complex *, complex *);

    /* Local variables */
    static integer info;
    static complex temp1, temp2;
    static integer i, j;
    extern logical lsame_(char *, char *);
    static integer ix, iy, jx, jy, kx, ky;
    extern /* Subroutine */ int xerbla_(char *, integer *);


/*  Purpose   
    =======   

    CHER2  performs the hermitian rank 2 operation   

       A := alpha*x*conjg( y' ) + conjg( alpha )*y*conjg( x' ) + A,   

    where alpha is a scalar, x and y are n element vectors and A is an n 
  
    by n hermitian matrix.   

    Parameters   
    ==========   

    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the upper or lower   
             triangular part of the array A is to be referenced as   
             follows:   

                UPLO = 'U' or 'u'   Only the upper triangular part of A   
                                    is to be referenced.   

                UPLO = 'L' or 'l'   Only the lower triangular part of A   
                                    is to be referenced.   

             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   

    ALPHA  - COMPLEX         .   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    X      - COMPLEX          array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the n   
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
             Before entry with  UPLO = 'U' or 'u', the leading n by n   
             upper triangular part of the array A must contain the upper 
  
             triangular part of the hermitian matrix and the strictly   
             lower triangular part of A is not referenced. On exit, the   
             upper triangular part of the array A is overwritten by the   
             upper triangular part of the updated matrix.   
             Before entry with UPLO = 'L' or 'l', the leading n by n   
             lower triangular part of the array A must contain the lower 
  
             triangular part of the hermitian matrix and the strictly   
             upper triangular part of A is not referenced. On exit, the   
             lower triangular part of the array A is overwritten by the   
             lower triangular part of the updated matrix.   
             Note that the imaginary parts of the diagonal elements need 
  
             not be set, they are assumed to be zero, and on exit they   
             are set to zero.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             max( 1, n ).   
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
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    } else if (*incy == 0) {
	info = 7;
    } else if (*lda < max(1,*n)) {
	info = 9;
    }
    if (info != 0) {
	xerbla_("CHER2 ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || alpha->r == 0.f && alpha->i == 0.f) {
	return 0;
    }

/*     Set up the start points in X and Y if the increments are not both 
  
       unity. */

    if (*incx != 1 || *incy != 1) {
	if (*incx > 0) {
	    kx = 1;
	} else {
	    kx = 1 - (*n - 1) * *incx;
	}
	if (*incy > 0) {
	    ky = 1;
	} else {
	    ky = 1 - (*n - 1) * *incy;
	}
	jx = kx;
	jy = ky;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through the triangular part   
       of A. */

    if (lsame_(uplo, "U")) {

/*        Form  A  when A is stored in the upper triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = j;
		i__3 = j;
		if (X(j).r != 0.f || X(j).i != 0.f || (Y(j).r != 0.f 
			|| Y(j).i != 0.f)) {
		    r_cnjg(&q__2, &Y(j));
		    q__1.r = alpha->r * q__2.r - alpha->i * q__2.i, q__1.i = 
			    alpha->r * q__2.i + alpha->i * q__2.r;
		    temp1.r = q__1.r, temp1.i = q__1.i;
		    i__2 = j;
		    q__2.r = alpha->r * X(j).r - alpha->i * X(j).i, 
			    q__2.i = alpha->r * X(j).i + alpha->i * X(j)
			    .r;
		    r_cnjg(&q__1, &q__2);
		    temp2.r = q__1.r, temp2.i = q__1.i;
		    i__2 = j - 1;
		    for (i = 1; i <= j-1; ++i) {
			i__3 = i + j * a_dim1;
			i__4 = i + j * a_dim1;
			i__5 = i;
			q__3.r = X(i).r * temp1.r - X(i).i * temp1.i, 
				q__3.i = X(i).r * temp1.i + X(i).i * 
				temp1.r;
			q__2.r = A(i,j).r + q__3.r, q__2.i = A(i,j).i + 
				q__3.i;
			i__6 = i;
			q__4.r = Y(i).r * temp2.r - Y(i).i * temp2.i, 
				q__4.i = Y(i).r * temp2.i + Y(i).i * 
				temp2.r;
			q__1.r = q__2.r + q__4.r, q__1.i = q__2.i + q__4.i;
			A(i,j).r = q__1.r, A(i,j).i = q__1.i;
/* L10: */
		    }
		    i__2 = j + j * a_dim1;
		    i__3 = j + j * a_dim1;
		    i__4 = j;
		    q__2.r = X(j).r * temp1.r - X(j).i * temp1.i, 
			    q__2.i = X(j).r * temp1.i + X(j).i * 
			    temp1.r;
		    i__5 = j;
		    q__3.r = Y(j).r * temp2.r - Y(j).i * temp2.i, 
			    q__3.i = Y(j).r * temp2.i + Y(j).i * 
			    temp2.r;
		    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
		    d__1 = A(j,j).r + q__1.r;
		    A(j,j).r = d__1, A(j,j).i = 0.f;
		} else {
		    i__2 = j + j * a_dim1;
		    i__3 = j + j * a_dim1;
		    d__1 = A(j,j).r;
		    A(j,j).r = d__1, A(j,j).i = 0.f;
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = jx;
		i__3 = jy;
		if (X(jx).r != 0.f || X(jx).i != 0.f || (Y(jy).r != 0.f 
			|| Y(jy).i != 0.f)) {
		    r_cnjg(&q__2, &Y(jy));
		    q__1.r = alpha->r * q__2.r - alpha->i * q__2.i, q__1.i = 
			    alpha->r * q__2.i + alpha->i * q__2.r;
		    temp1.r = q__1.r, temp1.i = q__1.i;
		    i__2 = jx;
		    q__2.r = alpha->r * X(jx).r - alpha->i * X(jx).i, 
			    q__2.i = alpha->r * X(jx).i + alpha->i * X(jx)
			    .r;
		    r_cnjg(&q__1, &q__2);
		    temp2.r = q__1.r, temp2.i = q__1.i;
		    ix = kx;
		    iy = ky;
		    i__2 = j - 1;
		    for (i = 1; i <= j-1; ++i) {
			i__3 = i + j * a_dim1;
			i__4 = i + j * a_dim1;
			i__5 = ix;
			q__3.r = X(ix).r * temp1.r - X(ix).i * temp1.i, 
				q__3.i = X(ix).r * temp1.i + X(ix).i * 
				temp1.r;
			q__2.r = A(i,j).r + q__3.r, q__2.i = A(i,j).i + 
				q__3.i;
			i__6 = iy;
			q__4.r = Y(iy).r * temp2.r - Y(iy).i * temp2.i, 
				q__4.i = Y(iy).r * temp2.i + Y(iy).i * 
				temp2.r;
			q__1.r = q__2.r + q__4.r, q__1.i = q__2.i + q__4.i;
			A(i,j).r = q__1.r, A(i,j).i = q__1.i;
			ix += *incx;
			iy += *incy;
/* L30: */
		    }
		    i__2 = j + j * a_dim1;
		    i__3 = j + j * a_dim1;
		    i__4 = jx;
		    q__2.r = X(jx).r * temp1.r - X(jx).i * temp1.i, 
			    q__2.i = X(jx).r * temp1.i + X(jx).i * 
			    temp1.r;
		    i__5 = jy;
		    q__3.r = Y(jy).r * temp2.r - Y(jy).i * temp2.i, 
			    q__3.i = Y(jy).r * temp2.i + Y(jy).i * 
			    temp2.r;
		    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
		    d__1 = A(j,j).r + q__1.r;
		    A(j,j).r = d__1, A(j,j).i = 0.f;
		} else {
		    i__2 = j + j * a_dim1;
		    i__3 = j + j * a_dim1;
		    d__1 = A(j,j).r;
		    A(j,j).r = d__1, A(j,j).i = 0.f;
		}
		jx += *incx;
		jy += *incy;
/* L40: */
	    }
	}
    } else {

/*        Form  A  when A is stored in the lower triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = j;
		i__3 = j;
		if (X(j).r != 0.f || X(j).i != 0.f || (Y(j).r != 0.f 
			|| Y(j).i != 0.f)) {
		    r_cnjg(&q__2, &Y(j));
		    q__1.r = alpha->r * q__2.r - alpha->i * q__2.i, q__1.i = 
			    alpha->r * q__2.i + alpha->i * q__2.r;
		    temp1.r = q__1.r, temp1.i = q__1.i;
		    i__2 = j;
		    q__2.r = alpha->r * X(j).r - alpha->i * X(j).i, 
			    q__2.i = alpha->r * X(j).i + alpha->i * X(j)
			    .r;
		    r_cnjg(&q__1, &q__2);
		    temp2.r = q__1.r, temp2.i = q__1.i;
		    i__2 = j + j * a_dim1;
		    i__3 = j + j * a_dim1;
		    i__4 = j;
		    q__2.r = X(j).r * temp1.r - X(j).i * temp1.i, 
			    q__2.i = X(j).r * temp1.i + X(j).i * 
			    temp1.r;
		    i__5 = j;
		    q__3.r = Y(j).r * temp2.r - Y(j).i * temp2.i, 
			    q__3.i = Y(j).r * temp2.i + Y(j).i * 
			    temp2.r;
		    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
		    d__1 = A(j,j).r + q__1.r;
		    A(j,j).r = d__1, A(j,j).i = 0.f;
		    i__2 = *n;
		    for (i = j + 1; i <= *n; ++i) {
			i__3 = i + j * a_dim1;
			i__4 = i + j * a_dim1;
			i__5 = i;
			q__3.r = X(i).r * temp1.r - X(i).i * temp1.i, 
				q__3.i = X(i).r * temp1.i + X(i).i * 
				temp1.r;
			q__2.r = A(i,j).r + q__3.r, q__2.i = A(i,j).i + 
				q__3.i;
			i__6 = i;
			q__4.r = Y(i).r * temp2.r - Y(i).i * temp2.i, 
				q__4.i = Y(i).r * temp2.i + Y(i).i * 
				temp2.r;
			q__1.r = q__2.r + q__4.r, q__1.i = q__2.i + q__4.i;
			A(i,j).r = q__1.r, A(i,j).i = q__1.i;
/* L50: */
		    }
		} else {
		    i__2 = j + j * a_dim1;
		    i__3 = j + j * a_dim1;
		    d__1 = A(j,j).r;
		    A(j,j).r = d__1, A(j,j).i = 0.f;
		}
/* L60: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = jx;
		i__3 = jy;
		if (X(jx).r != 0.f || X(jx).i != 0.f || (Y(jy).r != 0.f 
			|| Y(jy).i != 0.f)) {
		    r_cnjg(&q__2, &Y(jy));
		    q__1.r = alpha->r * q__2.r - alpha->i * q__2.i, q__1.i = 
			    alpha->r * q__2.i + alpha->i * q__2.r;
		    temp1.r = q__1.r, temp1.i = q__1.i;
		    i__2 = jx;
		    q__2.r = alpha->r * X(jx).r - alpha->i * X(jx).i, 
			    q__2.i = alpha->r * X(jx).i + alpha->i * X(jx)
			    .r;
		    r_cnjg(&q__1, &q__2);
		    temp2.r = q__1.r, temp2.i = q__1.i;
		    i__2 = j + j * a_dim1;
		    i__3 = j + j * a_dim1;
		    i__4 = jx;
		    q__2.r = X(jx).r * temp1.r - X(jx).i * temp1.i, 
			    q__2.i = X(jx).r * temp1.i + X(jx).i * 
			    temp1.r;
		    i__5 = jy;
		    q__3.r = Y(jy).r * temp2.r - Y(jy).i * temp2.i, 
			    q__3.i = Y(jy).r * temp2.i + Y(jy).i * 
			    temp2.r;
		    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
		    d__1 = A(j,j).r + q__1.r;
		    A(j,j).r = d__1, A(j,j).i = 0.f;
		    ix = jx;
		    iy = jy;
		    i__2 = *n;
		    for (i = j + 1; i <= *n; ++i) {
			ix += *incx;
			iy += *incy;
			i__3 = i + j * a_dim1;
			i__4 = i + j * a_dim1;
			i__5 = ix;
			q__3.r = X(ix).r * temp1.r - X(ix).i * temp1.i, 
				q__3.i = X(ix).r * temp1.i + X(ix).i * 
				temp1.r;
			q__2.r = A(i,j).r + q__3.r, q__2.i = A(i,j).i + 
				q__3.i;
			i__6 = iy;
			q__4.r = Y(iy).r * temp2.r - Y(iy).i * temp2.i, 
				q__4.i = Y(iy).r * temp2.i + Y(iy).i * 
				temp2.r;
			q__1.r = q__2.r + q__4.r, q__1.i = q__2.i + q__4.i;
			A(i,j).r = q__1.r, A(i,j).i = q__1.i;
/* L70: */
		    }
		} else {
		    i__2 = j + j * a_dim1;
		    i__3 = j + j * a_dim1;
		    d__1 = A(j,j).r;
		    A(j,j).r = d__1, A(j,j).i = 0.f;
		}
		jx += *incx;
		jy += *incy;
/* L80: */
	    }
	}
    }

    return 0;

/*     End of CHER2 . */

} /* cher2_ */

