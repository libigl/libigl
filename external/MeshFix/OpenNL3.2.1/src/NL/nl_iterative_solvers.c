/*
 *  Copyright (c) 2004-2010, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     levy@loria.fr
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 */

#include <NL/nl_iterative_solvers.h>
#include <NL/nl_blas.h>
#include <NL/nl_matrix.h>
#include <NL/nl_context.h>

/************************************************************************/
/* Solvers */

/*
 * The implementation of the solvers is inspired by 
 * the lsolver library, by Christian Badura, available from:
 * http://www.mathematik.uni-freiburg.de
 * /IAM/Research/projectskr/lin_solver/
 *
 * About the Conjugate Gradient, details can be found in:
 *  Ashby, Manteuffel, Saylor
 *     A taxononmy for conjugate gradient methods
 *     SIAM J Numer Anal 27, 1542-1568 (1990)
 */

NLuint nlSolve_CG() {
    NLdouble* b        = nlCurrentContext->b ;
    NLdouble* x        = nlCurrentContext->x ;
    NLdouble  eps      = nlCurrentContext->threshold ;
    NLuint    max_iter = nlCurrentContext->max_iterations ;
    NLint     N        = nlCurrentContext->n ;

    NLdouble *g = NL_NEW_ARRAY(NLdouble, N) ;
    NLdouble *r = NL_NEW_ARRAY(NLdouble, N) ; 
    NLdouble *p = NL_NEW_ARRAY(NLdouble, N) ;
    NLuint its=0;
    NLint i;
    NLdouble t, tau, sig, rho, gam;
    NLdouble b_square=ddot(N,b,1,b,1);
    NLdouble err=eps*eps*b_square;
    NLdouble accu =0.0;
    NLdouble * Ax=NL_NEW_ARRAY(NLdouble,nlCurrentContext->n);
    NLdouble curr_err;
    
    nlCurrentContext->matrix_vector_prod(x,g);
    daxpy(N,-1.,b,1,g,1);
    dscal(N,-1.,g,1);
    dcopy(N,g,1,r,1);
    curr_err = ddot(N,g,1,g,1);
    while ( curr_err >err && its < max_iter) {
	    if(!(its % 100)) {
	        printf ( "%d : %.10e -- %.10e\n", its, curr_err, err ) ;
	    }
        nlCurrentContext->matrix_vector_prod(r,p);
        rho=ddot(N,p,1,p,1);
        sig=ddot(N,r,1,p,1);
        tau=ddot(N,g,1,r,1);
        t=tau/sig;
        daxpy(N,t,r,1,x,1);
        daxpy(N,-t,p,1,g,1);
        gam=(t*t*rho-tau)/tau;
        dscal(N,gam,r,1);
        daxpy(N,1.,g,1,r,1);
        ++its;
        curr_err = ddot(N,g,1,g,1); 
    }
    nlCurrentContext->matrix_vector_prod(x,Ax);
    for(i = 0 ; i < N ; ++i)
        accu+=(Ax[i]-b[i])*(Ax[i]-b[i]);
    printf("in OpenNL : ||Ax-b||/||b|| = %e\n",sqrt(accu)/sqrt(b_square));
    NL_DELETE_ARRAY(Ax);
    NL_DELETE_ARRAY(g) ;
    NL_DELETE_ARRAY(r) ;
    NL_DELETE_ARRAY(p) ;
    return its;
} 


NLuint nlSolve_CG_precond()  {
    NLdouble* b        = nlCurrentContext->b ;
    NLdouble* x        = nlCurrentContext->x ;
    NLdouble  eps      = nlCurrentContext->threshold ;
    NLuint    max_iter = nlCurrentContext->max_iterations ;
    NLint     N        = nlCurrentContext->n ;

    NLdouble* r = NL_NEW_ARRAY(NLdouble, N) ;
    NLdouble* d = NL_NEW_ARRAY(NLdouble, N) ;
    NLdouble* h = NL_NEW_ARRAY(NLdouble, N) ;
    NLdouble *Ad = h;
    NLuint its=0;
    NLdouble rh, alpha, beta;
    NLdouble b_square = ddot(N,b,1,b,1);
    NLdouble err=eps*eps*b_square;
    NLint i;
    NLdouble * Ax=NL_NEW_ARRAY(NLdouble,nlCurrentContext->n);
    NLdouble accu =0.0;
    NLdouble curr_err;
    

    nlCurrentContext->matrix_vector_prod(x,r);
    daxpy(N,-1.,b,1,r,1);
    nlCurrentContext->precond_vector_prod(r,d);
    dcopy(N,d,1,h,1);
    rh=ddot(N,r,1,h,1);
    curr_err = ddot(N,r,1,r,1);
    while ( curr_err >err && its < max_iter) {
	if(!(its % 100)) {
	   printf ( "%d : %.10e -- %.10e\n", its, curr_err, err ) ;
	}
        nlCurrentContext->matrix_vector_prod(d,Ad);
        alpha=rh/ddot(N,d,1,Ad,1);
        daxpy(N,-alpha,d,1,x,1);
        daxpy(N,-alpha,Ad,1,r,1);
        nlCurrentContext->precond_vector_prod(r,h);
        beta=1./rh; rh=ddot(N,r,1,h,1); beta*=rh;
        dscal(N,beta,d,1);
        daxpy(N,1.,h,1,d,1);
        ++its;
        // calcul de l'erreur courante
        curr_err = ddot(N,r,1,r,1);

    }
    nlCurrentContext->matrix_vector_prod(x,Ax);
    for(i = 0 ; i < N ; ++i)
        accu+=(Ax[i]-b[i])*(Ax[i]-b[i]);
    printf("in OpenNL : ||Ax-b||/||b|| = %e\n",sqrt(accu)/sqrt(b_square));
    NL_DELETE_ARRAY(Ax);
    NL_DELETE_ARRAY(r) ;
    NL_DELETE_ARRAY(d) ;
    NL_DELETE_ARRAY(h) ;
    
    return its;
}

NLuint nlSolve_BICGSTAB() {
    NLdouble* b        = nlCurrentContext->b ;
    NLdouble* x        = nlCurrentContext->x ;
    NLdouble  eps      = nlCurrentContext->threshold ;
    NLuint    max_iter = nlCurrentContext->max_iterations ;
    NLint     N        = nlCurrentContext->n ;
    NLint     i ;

    NLdouble *rT  = NL_NEW_ARRAY(NLdouble, N) ; 
    NLdouble *d   = NL_NEW_ARRAY(NLdouble, N) ; 
    NLdouble *h   = NL_NEW_ARRAY(NLdouble, N) ; 
    NLdouble *u   = NL_NEW_ARRAY(NLdouble, N) ; 
    NLdouble *Ad  = NL_NEW_ARRAY(NLdouble, N) ; 
    NLdouble *t   = NL_NEW_ARRAY(NLdouble, N) ; 
    NLdouble *s   = h;
    NLdouble rTh, rTAd, rTr, alpha, beta, omega, st, tt ;
    NLuint its=0;
    NLdouble b_square = ddot(N,b,1,b,1);
    NLdouble err=eps*eps*b_square;
    NLdouble *r = NL_NEW_ARRAY(NLdouble, N) ;
    NLdouble * Ax=NL_NEW_ARRAY(NLdouble,nlCurrentContext->n);
    NLdouble accu =0.0;

    nlCurrentContext->matrix_vector_prod(x,r);
    daxpy(N,-1.,b,1,r,1);
    dcopy(N,r,1,d,1);
    dcopy(N,d,1,h,1);
    dcopy(N,h,1,rT,1);
    nl_assert( ddot(N,rT,1,rT,1)>1e-40 );
    rTh=ddot(N,rT,1,h,1);
    rTr=ddot(N,r,1,r,1);
    while ( rTr>err && its < max_iter) {
	    if(!(its % 100)) {
	        printf ( "%d : %.10e -- %.10e\n", its, rTr, err ) ;
	    }
        nlCurrentContext->matrix_vector_prod(d,Ad);
        rTAd=ddot(N,rT,1,Ad,1);
        nl_assert( fabs(rTAd)>1e-40 );
        alpha=rTh/rTAd;
        daxpy(N,-alpha,Ad,1,r,1);
        dcopy(N,h,1,s,1);
        daxpy(N,-alpha,Ad,1,s,1);
        nlCurrentContext->matrix_vector_prod(s,t);
        daxpy(N,1.,t,1,u,1);
        dscal(N,alpha,u,1);
        st=ddot(N,s,1,t,1);
        tt=ddot(N,t,1,t,1);
        if ( fabs(st)<1e-40 || fabs(tt)<1e-40 ) {
            omega = 0.;
        } else {
            omega = st/tt;
        }
        daxpy(N,-omega,t,1,r,1);
        daxpy(N,-alpha,d,1,x,1);
        daxpy(N,-omega,s,1,x,1);
        dcopy(N,s,1,h,1);
        daxpy(N,-omega,t,1,h,1);
        beta=(alpha/omega)/rTh; rTh=ddot(N,rT,1,h,1); beta*=rTh;
        dscal(N,beta,d,1);
        daxpy(N,1.,h,1,d,1);
        daxpy(N,-beta*omega,Ad,1,d,1);
        rTr=ddot(N,r,1,r,1);
        ++its;
    }

    nlCurrentContext->matrix_vector_prod(x,Ax);
    for(i = 0 ; i < N ; ++i){
        accu+=(Ax[i]-b[i])*(Ax[i]-b[i]);
    }
    printf("in OpenNL : ||Ax-b||/||b|| = %e\n",sqrt(accu)/sqrt(b_square));
    NL_DELETE_ARRAY(Ax);
    NL_DELETE_ARRAY(r) ;
    NL_DELETE_ARRAY(rT) ;
    NL_DELETE_ARRAY(d) ;
    NL_DELETE_ARRAY(h) ;
    NL_DELETE_ARRAY(u) ;
    NL_DELETE_ARRAY(Ad) ;
    NL_DELETE_ARRAY(t) ;
    
    return its;
}


NLuint nlSolve_BICGSTAB_precond() {

    NLdouble* b        = nlCurrentContext->b ;
    NLdouble* x        = nlCurrentContext->x ;
    NLdouble  eps      = nlCurrentContext->threshold ;
    NLuint    max_iter = nlCurrentContext->max_iterations ;
    NLint     N        = nlCurrentContext->n ;
    NLint     i;

    NLdouble *rT  = NL_NEW_ARRAY(NLdouble, N) ;
    NLdouble *d   = NL_NEW_ARRAY(NLdouble, N) ;
    NLdouble *h   = NL_NEW_ARRAY(NLdouble, N) ;
    NLdouble *u   = NL_NEW_ARRAY(NLdouble, N) ;
    NLdouble *Sd  = NL_NEW_ARRAY(NLdouble, N) ;
    NLdouble *t   = NL_NEW_ARRAY(NLdouble, N) ;
    NLdouble *aux = NL_NEW_ARRAY(NLdouble, N) ;
    NLdouble *s   = h;
    NLdouble rTh, rTSd, rTr, alpha, beta, omega, st, tt;
    NLuint its=0;
    NLdouble b_square = ddot(N,b,1,b,1);
    NLdouble err  = eps*eps*b_square;
    NLdouble *r   = NL_NEW_ARRAY(NLdouble, N);
    NLdouble * Ax = NL_NEW_ARRAY(NLdouble,nlCurrentContext->n);
    NLdouble accu =0.0;

    nlCurrentContext->matrix_vector_prod(x,r);
    daxpy(N,-1.,b,1,r,1);
    nlCurrentContext->precond_vector_prod(r,d);
    dcopy(N,d,1,h,1);
    dcopy(N,h,1,rT,1);
    nl_assert( ddot(N,rT,1,rT,1)>1e-40 );
    rTh=ddot(N,rT,1,h,1);
    rTr=ddot(N,r,1,r,1);
    while ( rTr>err && its < max_iter) {
	    if(!(its % 100)) {
	       printf ( "%d : %.10e -- %.10e\n", its, rTr, err ) ;
	    }
        nlCurrentContext->matrix_vector_prod(d,aux);
        nlCurrentContext->precond_vector_prod(aux,Sd);
        rTSd=ddot(N,rT,1,Sd,1);
        nl_assert( fabs(rTSd)>1e-40 );
        alpha=rTh/rTSd;
        daxpy(N,-alpha,aux,1,r,1);
        dcopy(N,h,1,s,1);
        daxpy(N,-alpha,Sd,1,s,1);
        nlCurrentContext->matrix_vector_prod(s,aux);
        nlCurrentContext->precond_vector_prod(aux,t);
        daxpy(N,1.,t,1,u,1);
        dscal(N,alpha,u,1);
        st=ddot(N,s,1,t,1);
        tt=ddot(N,t,1,t,1);
        if ( fabs(st)<1e-40 || fabs(tt)<1e-40 ) {
            omega = 0.;
        } else {
            omega = st/tt;
        }
        daxpy(N,-omega,aux,1,r,1);
        daxpy(N,-alpha,d,1,x,1);
        daxpy(N,-omega,s,1,x,1);
        dcopy(N,s,1,h,1);
        daxpy(N,-omega,t,1,h,1);
        beta=(alpha/omega)/rTh; rTh=ddot(N,rT,1,h,1); beta*=rTh;
        dscal(N,beta,d,1);
        daxpy(N,1.,h,1,d,1);
        daxpy(N,-beta*omega,Sd,1,d,1);
        rTr=ddot(N,r,1,r,1);
        ++its;
    }

    nlCurrentContext->matrix_vector_prod(x,Ax);
    for(i = 0 ; i < N ; ++i){
        accu+=(Ax[i]-b[i])*(Ax[i]-b[i]);
    }
    printf("in OpenNL : ||Ax-b||/||b|| = %e\n",sqrt(accu)/sqrt(b_square));
    NL_DELETE_ARRAY(Ax);
    NL_DELETE_ARRAY(r);
    NL_DELETE_ARRAY(rT);
    NL_DELETE_ARRAY(d);
    NL_DELETE_ARRAY(h);
    NL_DELETE_ARRAY(u);
    NL_DELETE_ARRAY(Sd);
    NL_DELETE_ARRAY(t);
    NL_DELETE_ARRAY(aux);

    return its;
}

NLuint nlSolve_GMRES() {

    NLdouble* b        = nlCurrentContext->b ;
    NLdouble* x        = nlCurrentContext->x ;
    NLdouble  eps      = nlCurrentContext->threshold ;
    NLint    max_iter  = nlCurrentContext->max_iterations ;
    NLint    n         = nlCurrentContext->n ;
    NLint    m         = nlCurrentContext->inner_iterations ;

    typedef NLdouble *NLdoubleP;
    NLdouble *V   = NL_NEW_ARRAY(NLdouble, n*(m+1)   ) ;
    NLdouble *U   = NL_NEW_ARRAY(NLdouble, m*(m+1)/2 ) ;
    NLdouble *r   = NL_NEW_ARRAY(NLdouble, n         ) ;
    NLdouble *y   = NL_NEW_ARRAY(NLdouble, m+1       ) ;
    NLdouble *c   = NL_NEW_ARRAY(NLdouble, m         ) ;
    NLdouble *s   = NL_NEW_ARRAY(NLdouble, m         ) ;
    NLdouble **v  = NL_NEW_ARRAY(NLdoubleP, m+1      ) ;
    NLdouble * Ax = NL_NEW_ARRAY(NLdouble,nlCurrentContext->n);
    NLdouble accu =0.0;
    NLint i, j, io, uij, u0j ; 
    NLint its = -1 ;
    NLdouble beta, h, rd, dd, nrm2b ;

    for ( i=0; i<=m; ++i ){
        v[i]=V+i*n ;
    }
    
    nrm2b=dnrm2(n,b,1);
    io=0;
    do  { /* outer loop */
        ++io;
        nlCurrentContext->matrix_vector_prod(x,r);
        daxpy(n,-1.,b,1,r,1);
        beta=dnrm2(n,r,1);
        dcopy(n,r,1,v[0],1);
        dscal(n,1./beta,v[0],1);
        
        y[0]=beta;
        j=0;
        uij=0;
        do { /* inner loop: j=0,...,m-1 */
            u0j=uij;
            nlCurrentContext->matrix_vector_prod(v[j],v[j+1]);
            dgemv(
                Transpose,n,j+1,1.,V,n,v[j+1],1,0.,U+u0j,1
            );
            dgemv(
                NoTranspose,n,j+1,-1.,V,n,U+u0j,1,1.,v[j+1],1
            );
            h=dnrm2(n,v[j+1],1);
            dscal(n,1./h,v[j+1],1);
            for (i=0; i<j; ++i ) { /* rotiere neue Spalte */
                double tmp = c[i]*U[uij]-s[i]*U[uij+1];
                U[uij+1]   = s[i]*U[uij]+c[i]*U[uij+1];
                U[uij]     = tmp;
                ++uij;
            }
            { /* berechne neue Rotation */
                rd     = U[uij];
                dd     = sqrt(rd*rd+h*h);
                c[j]   = rd/dd;
                s[j]   = -h/dd;
                U[uij] = dd;
                ++uij;
            }
            { /* rotiere rechte Seite y (vorher: y[j+1]=0) */
                y[j+1] = s[j]*y[j];
                y[j]   = c[j]*y[j];
            }
            ++j;
        } while ( 
            j<m && fabs(y[j])>=eps*nrm2b 
        ) ;
        { /* minimiere bzgl Y */
            dtpsv(
                UpperTriangle,
                NoTranspose,
                NotUnitTriangular,
                j,U,y,1
            );
            /* correct X */
            dgemv(NoTranspose,n,j,-1.,V,n,y,1,1.,x,1);
        }
    } while ( fabs(y[j])>=eps*nrm2b && (m*(io-1)+j) < max_iter);
    
    /* Count the inner iterations */
    its = m*(io-1)+j;

    nlCurrentContext->matrix_vector_prod(x,Ax);
    for(i = 0 ; i < n ; ++i)
        accu+=(Ax[i]-b[i])*(Ax[i]-b[i]);
    printf("in OpenNL : ||Ax-b||/||b|| = %e\n",sqrt(accu)/nrm2b);
    NL_DELETE_ARRAY(Ax);
    NL_DELETE_ARRAY(V) ;
    NL_DELETE_ARRAY(U) ;
    NL_DELETE_ARRAY(r) ;
    NL_DELETE_ARRAY(y) ;
    NL_DELETE_ARRAY(c) ;
    NL_DELETE_ARRAY(s) ;
    NL_DELETE_ARRAY(v) ;
    
    return its;
}





