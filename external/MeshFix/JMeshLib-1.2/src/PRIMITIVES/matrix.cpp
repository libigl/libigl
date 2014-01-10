/****************************************************************************
* JMeshLib                                                                  *
*                                                                           *
* Consiglio Nazionale delle Ricerche                                        *
* Istituto di Matematica Applicata e Tecnologie Informatiche                *
* Sezione di Genova                                                         *
* IMATI-GE / CNR                                                            *
*                                                                           *
* Authors: Marco Attene                                                     *
*                                                                           *
* Copyright(C) 2006: IMATI-GE / CNR                                         *
*                                                                           *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

#include <math.h>
#include "matrix.h"

#define FABS(a) (((a)<0)?(-(a)):((a)))


//////////////////////////////////////////////////////////////////////////
//
// Generic 3x3 matrix
//
//////////////////////////////////////////////////////////////////////////

// Plain constructor

Matrix3x3::Matrix3x3(const double& a11, const double& a12, const double& a13,
                     const double& a21, const double& a22, const double& a23, 
                     const double& a31, const double& a32, const double& a33)
{
 M[0] = a11; M[1] = a12; M[2] = a13;
 M[3] = a21; M[4] = a22; M[5] = a23;
 M[6] = a31; M[7] = a32; M[8] = a33;
}


// Matrix T(v1,v2,v3)*(w1,w2,w3).

Matrix3x3::Matrix3x3(const double& v1, const double& v2, const double& v3,
                     const double& w1, const double& w2, const double& w3)
{
 M[0] = v1*w1; M[1] = v1*w2; M[2] = v1*w3;
 M[3] = v2*w1; M[4] = v2*w2; M[5] = v2*w3;
 M[6] = v3*w1; M[7] = v3*w2; M[8] = v3*w3;
}


// Symmetric matrix T(x,y,z)*(x,y,z)

Matrix3x3::Matrix3x3(const double& x, const double& y, const double& z)
{
 M[0] = x*x;  M[1] = x*y;  M[2] = x*z;
 M[3] = M[1]; M[4] = y*y;  M[5] = y*z;
 M[6] = M[2]; M[7] = M[5]; M[8] = z*z;
}


// Self-sum

void Matrix3x3::operator+=(const Matrix3x3& s)
{
 M[0] += s.M[0]; M[1] += s.M[1]; M[2] += s.M[2];
 M[3] += s.M[3]; M[4] += s.M[4]; M[5] += s.M[5];
 M[6] += s.M[6]; M[7] += s.M[7]; M[8] += s.M[8];
}


// Self-subtraction

void Matrix3x3::operator-=(const Matrix3x3& s)
{
 M[0] -= s.M[0]; M[1] -= s.M[1]; M[2] -= s.M[2];
 M[3] -= s.M[3]; M[4] -= s.M[4]; M[5] -= s.M[5];
 M[6] -= s.M[6]; M[7] -= s.M[7]; M[8] -= s.M[8];
}


// Self-multiplication

void Matrix3x3::operator*=(const double& d)
{
 M[0] *= d; M[1] *= d; M[2] *= d;
 M[3] *= d; M[4] *= d; M[5] *= d;
 M[6] *= d; M[7] *= d; M[8] *= d;
}


// Sum

Matrix3x3 Matrix3x3::operator+(const Matrix3x3& s) const
{
 return Matrix3x3(M[0]+s.M[0], M[1]+s.M[1], M[2]+s.M[2], 
                  M[3]+s.M[3], M[4]+s.M[4], M[5]+s.M[5], 
                  M[6]+s.M[6], M[7]+s.M[7], M[8]+s.M[8]);
}


// Scalar Multiplication

Matrix3x3 Matrix3x3::operator*(const double& d) const
{
 return Matrix3x3(M[0]*d, M[1]*d, M[2]*d, 
                  M[3]*d, M[4]*d, M[5]*d, 
                  M[6]*d, M[7]*d, M[8]*d);
}


// Matrix multiplication

Matrix3x3 Matrix3x3::operator*(const Matrix3x3& q) const
{
 return Matrix3x3(
		  M[0]*q.M[0]+M[1]*q.M[3]+M[2]*q.M[6], M[0]*q.M[1]+M[1]*q.M[4]+M[2]*q.M[7], M[0]*q.M[2]+M[1]*q.M[5]+M[2]*q.M[8], 
		  M[3]*q.M[0]+M[4]*q.M[3]+M[5]*q.M[6], M[3]*q.M[1]+M[4]*q.M[4]+M[5]*q.M[7], M[3]*q.M[2]+M[4]*q.M[5]+M[5]*q.M[8], 
		  M[6]*q.M[0]+M[7]*q.M[3]+M[8]*q.M[6], M[6]*q.M[1]+M[7]*q.M[4]+M[8]*q.M[7], M[6]*q.M[2]+M[7]*q.M[5]+M[8]*q.M[8]
		 );
}


// Matrix transpose

Matrix3x3 Matrix3x3::operator~() const
{
 return Matrix3x3(M[0],M[3],M[6],M[1],M[4],M[7],M[2],M[5],M[8]);
}


// Computes (x,y,z)*M*(x,y,z)

double Matrix3x3::lrMultiply(const double& x, const double& y, const double& z) const
{
 return (x*(x*M[0] + y*M[3] + z*M[6]) + y*(x*M[1] + y*M[4] + z*M[7]) + z*(x*M[2] + y*M[5] + z*M[8]));
}


// Computes v*M*w
double Matrix3x3::lrMultiply(const double& v1, const double& v2, const double& v3,
                             const double& w1, const double& w2, const double& w3) const
{
 return (w1*(v1*M[0] + v2*M[3] + v3*M[6]) + w2*(v1*M[1] + v2*M[4] + v3*M[7]) + w3*(v1*M[2] + v2*M[5] + v3*M[8]));
}

Matrix3x3 Matrix3x3::transpose() const
{
 return Matrix3x3(M[0], M[3], M[6], 
                  M[1], M[4], M[7], 
                  M[2], M[5], M[8]);
}

//////////////////////////////////////////////////////////////////////////
//
// Symmetric 3x3 matrix
//
//////////////////////////////////////////////////////////////////////////

// Plain constructor

SymMatrix3x3::SymMatrix3x3(const double&a11, const double&a12, const double&a22,
			   const double&a13, const double&a23, const double&a33)
{
 M[0] = a11; M[1] = a12; M[3] = a13;
             M[2] = a22; M[4] = a23;
                         M[5] = a33;
}


// Symmetric matrix T(x,y,z)*(x,y,z)

SymMatrix3x3::SymMatrix3x3(const double& x, const double& y, const double& z)
{
 M[0] = x*x;  M[1] = x*y;  M[3] = x*z;
              M[2] = y*y;  M[4] = y*z;
                           M[5] = z*z;
}


// Constructor from generic matrix

SymMatrix3x3::SymMatrix3x3(const Matrix3x3& q)
{
 M[0] = q.M[0]; M[1] = q.M[1]; M[3] = q.M[2];
                M[2] = q.M[4]; M[4] = q.M[5];
                               M[5] = q.M[8];
}


// Self-sum

void SymMatrix3x3::operator+=(const SymMatrix3x3& s)
{
 M[0] += s.M[0]; M[1] += s.M[1]; M[2] += s.M[2];
 M[3] += s.M[3]; M[4] += s.M[4]; M[5] += s.M[5];
}


// Self-subtraction

void SymMatrix3x3::operator-=(const SymMatrix3x3& s)
{
 M[0] -= s.M[0]; M[1] -= s.M[1]; M[2] -= s.M[2];
 M[3] -= s.M[3]; M[4] -= s.M[4]; M[5] -= s.M[5];
}


// Self-multiplication

void SymMatrix3x3::operator*=(const double& d)
{
 M[0] *= d; M[1] *= d; M[2] *= d;
 M[3] *= d; M[4] *= d; M[5] *= d;
}


// Sum

SymMatrix3x3 SymMatrix3x3::operator+(const SymMatrix3x3& s) const
{
 return SymMatrix3x3(M[0]+s.M[0], M[1]+s.M[1], M[2]+s.M[2], 
                     M[3]+s.M[3], M[4]+s.M[4], M[5]+s.M[5]);
}


// Multiplication

SymMatrix3x3 SymMatrix3x3::operator*(const double& d) const
{
 return SymMatrix3x3(M[0]*d, M[1]*d, M[2]*d, M[3]*d, M[4]*d, M[5]*d);
}


// Computes (x,y,z)*M*(x,y,z)

double SymMatrix3x3::lrMultiply(const double& x, const double& y, const double& z) const
{
 double a,b,c;
 a = x*M[0] + y*M[1] + z*M[3];
 b = x*M[1] + y*M[2] + z*M[4];
 c = x*M[3] + y*M[4] + z*M[5];
 return (x*a + y*b + z*c);
}



// Computes v*M*w
double SymMatrix3x3::lrMultiply(const double& v1, const double& v2, const double& v3,
                             const double& w1, const double& w2, const double& w3) const
{
 return (w1*(v1*M[0] + v2*M[1] + v3*M[3]) + w2*(v1*M[1] + v2*M[2] + v3*M[4]) + w3*(v1*M[3] + v2*M[4] + v3*M[5]));
}

// Invert the matrix. If singular return FALSE

bool SymMatrix3x3::invert()
{
 double det, pos, neg, t, out[6];

 pos = neg = 0.0;
 t =  M[0]*M[2]*M[5]; ((t >= 0.0)?(pos):(neg)) += t;
 t =  M[1]*M[4]*M[3]; ((t >= 0.0)?(pos):(neg)) += t;
 t =  M[3]*M[1]*M[4]; ((t >= 0.0)?(pos):(neg)) += t;
 t = -M[3]*M[2]*M[3]; ((t >= 0.0)?(pos):(neg)) += t;
 t = -M[1]*M[1]*M[5]; ((t >= 0.0)?(pos):(neg)) += t;
 t = -M[0]*M[4]*M[4]; ((t >= 0.0)?(pos):(neg)) += t;
 det = pos+neg;

 t = det/(pos-neg);
 if (FABS(t) >= 1.0e-15)
 {
  out[0] =  (M[2] * M[5] - M[4] * M[4]) /det;
  out[1] = -(M[1] * M[5] - M[4] * M[3]) /det;
  out[2] =  (M[0] * M[5] - M[3] * M[3]) /det;
  out[3] =  (M[1] * M[4] - M[2] * M[3]) /det;
  out[4] = -(M[0] * M[4] - M[1] * M[3]) /det;
  out[5] =  (M[0] * M[2] - M[1] * M[1]) /det;
  M[0] = out[0]; M[1] = out[1]; M[2] = out[2];
  M[3] = out[3]; M[4] = out[4]; M[5] = out[5];
  return 1;
 }

 return 0;
}


// Compute eigenvalues and eigenvectors of the matrix (JACOBI method).
// The calling function is responsible of verifying that the matrix
// is diagonalizable.
//
// This routine was inspired from a software to estimate curvature
// tensors developed at INRIA. Visit the following link for details:
// http://www-sop.inria.fr/geometrica/team/Pierre.Alliez/demos/curvature/
//
// This version has been slightly optimized.
//

void SymMatrix3x3::diagonalize(double *eigen_val, double *eigen_vec) const
{
 static const double EPS = 0.00001;
 static const double cos_pi_4 = 0.70710678;
 static const int MAX_ITER = 100;
 double a_P[6], v_P[9];
 double *a = a_P, *v = v_P;
 double a_norm,a_normEPS,thr,thr_nn;
 int nb_iter = 0;
 int i,j,k,ij,jj,ik,l,m,lm,mq,lq,ll,mm,imv,im,iq,ilv,il;
 int index_P[3];
 int *index = index_P;
 double a_ij,a_lm,a_ll,a_mm,a_im,a_il,a_lm_2,v_ilv,v_imv,x;
 double sinx,sinx_2,cosx,cosx_2,sincos,delta;

 a[0] = M[0]; a[1] = M[1]; a[2] = M[2];
 a[3] = M[3]; a[4] = M[4]; a[5] = M[5];
 a--;
    
 // Step 2 : Init diagonalization matrix as the unit matrix
    
 for (ij=0, i=0; i<3; i++) for (j=0; j<3; j++) v[ij++] = (i==j)?(1.0):(0.0); 
 v--;
    
 // Step 3 : compute the weight of the non diagonal terms

 a_norm = 0.0;
 for (i=1, ij=1; i<=3; i++) for (j=1; j<=i; j++, ij++) if( i!=j ) {a_ij = a[ij]; a_norm += a_ij*a_ij;}
    
 if( a_norm != 0.0 )
 {   
   a_normEPS = a_norm*EPS;
   thr   = a_norm  ;
  
   // Step 4 : rotations
   while (thr > a_normEPS && nb_iter < MAX_ITER)
   { 
   nb_iter++;
   thr_nn = thr / 6;
    
   for (l=1; l<3; l++)
   for (m=l+1; m<=3; m++)
   {     
     // compute sinx and cosx 
        
     lq = (l*l-l)/2; mq = (m*m-m)/2;    
     lm = l+mq; a_lm = a[lm];
     a_lm_2 = a_lm*a_lm;
        
     if( a_lm_2 < thr_nn ) continue;
        
     ll = l+lq; mm = m+mq;
     a_ll = a[ll]; a_mm = a[mm];  
     delta = a_ll - a_mm;
        
     if (delta==0.0) {sinx = -cos_pi_4; cosx = cos_pi_4;}
     else {x = -atan( (a_lm+a_lm) / delta )/2.0; sinx=sin(x); cosx=cos(x);}

     sinx_2  = sinx*sinx;
     cosx_2  = cosx*cosx;
     sincos  = sinx*cosx;
        
     // rotate L and M columns 
    
     ilv = 3*(l-1); imv = 3*(m-1);
        
     for( i=1; i<=3;i++ )
     {
     if( (i!=l) && (i!=m) )
     {
     iq = (i*i-i)/2;
     im = (i<m)?(i+mq):(m+iq);   
     a_im = a[im];
     il = (i<l)?(i+lq):(l+iq);   
     a_il = a[il];
     a[il] =  a_il*cosx - a_im*sinx;
     a[im] =  a_il*sinx + a_im*cosx;
     }
        
     ilv++; imv++;  
     v_ilv = v[ilv]; v_imv = v[imv];
     v[ilv] = cosx*v_ilv - sinx*v_imv;
     v[imv] = sinx*v_ilv + cosx*v_imv;
     } 
        
     x = a_lm*sincos; x+=x;
        
     a[ll] =  a_ll*cosx_2 + a_mm*sinx_2 - x;
     a[mm] =  a_ll*sinx_2 + a_mm*cosx_2 + x;
     a[lm] =  0.0;
     thr = FABS( thr - a_lm_2 );
   }
   }   
 }
    
 // Step 5: index conversion and copy eigen values 
    
 // back from Fortran to C++
 a++;
    
 for( i=0; i<3; i++ ) { k = i + (i*(i+1))/2; eigen_val[i] = a[k];}
      
 // Step 6: sort the eigen values and eigen vectors 
    
 for( i=0; i<3; i++ ) index[i] = i;
    
 for( i=0; i<(3-1); i++ ) {
   x = eigen_val[i];
   k = i;
    
   for( j=i+1; j<3; j++ ) if( x < eigen_val[j] ) { k = j; x = eigen_val[j];}
    
   eigen_val[k] = eigen_val[i];
   eigen_val[i] = x;
    
   jj   = index[k];
   index[k] = index[i];
   index[i] = jj;
 }

 // Step 7: save the eigen vectors 
  
 v++; // back from Fortran to to C++
    
 ij = 0;
 for( k=0; k<3; k++ ) for(ik = index[k]*3, i=0; i<3; i++ ) eigen_vec[ij++] = v[ik++];
}


// Computes the eigenvalues of the matrix

void SymMatrix3x3::getEigenvalues(double *L1, double *L2, double *L3) const
{
 double a11 = M[0], a12 = M[1], a22 = M[2], a13 = M[3], a23 = M[4], a33 = M[5];
 double c0 = (a11*a22*a33)+(2.0*a12*a13*a23)-(a11*a23*a23)-(a22*a13*a13)-(a33*a12*a12); // Mat. Determinant
 double c1 = a11*a22-a12*a12+a11*a33-a13*a13+a22*a33-a23*a23;
 double c2 = a11+a22+a33;
 double a = (3*c1-c2*c2)/3.0;
 double b = (9*c1*c2-2*c2*c2*c2-27*c0)/27.0;
 double Q = ((b*b)/4.0)+((a*a*a)/27.0);

 if (Q>1.0e-12) {*L1 = *L2 = *L3 = a11; return;} // Evecs = (1,0,0), (0,1,0), (0,0,1)

 double l1, l2, l3; // Eigenvalues to be computed

 if (Q>=0) {double p=(b>0)?(pow(b/2, 1.0/3.0)):(0); l1=l2=((c2/3.0)+p); l3=((c2/3.0)-2.0*p);}
 else
 {
  double t = atan2(sqrt(-Q), -b/2.0)/3.0, r = pow(((b*b)/4.0)-Q, 1.0/6.0);
  double cos_t = cos(t), sin_t = sin(t);
  const double sq3 = sqrt(3.0);
  l1 = l2 = l3 = (c2/3.0);
  l1 += (2*r*cos_t);
  l2 -= r*(cos_t+sq3*sin_t);
  l3 -= r*(cos_t-sq3*sin_t);
 }

 if (l1<=l2 && l1<=l3) {*L1=l1; *L2=(l2<l3)?(l2):(l3); *L3=(l2<l3)?(l3):(l2);}
 else if (l2<=l1 && l2<=l3) {*L1=l2; *L2=(l1<l3)?(l1):(l3); *L3=(l1<l3)?(l3):(l1);}
 else {*L1=l3; *L2=(l1<l2)?(l1):(l2); *L3=(l1<l2)?(l2):(l1);}
}


// Computes the eigenvector corresp. to the minimum eigenvalue

void SymMatrix3x3::getMinEigenvector(double *x, double *y, double *z) const
{
 double a11 = M[0], a12 = M[1], a22 = M[2], a13 = M[3], a23 = M[4], a33 = M[5];
 double l, l1, l2, l3, c0, c1, c2;
 getEigenvalues(&l, &l2, &l3);

 if (l==l3 && l==l2) {*x = 1; *y = *z = 0; return;}

 a11-=l; a22-=l; a33-=l;
 double u11 = a22*a33-a23*a23, u12 = a13*a23-a12*a33, u13 = a12*a23-a13*a22;
 double u22 = a11*a33-a13*a13, u23 = a12*a13-a23*a11, u33 = a11*a22-a12*a12;
 l1 = u11*u11+u12*u12+u13*u13;
 l2 = u12*u12+u22*u22+u23*u23;
 l3 = u13*u13+u23*u23+u33*u33;

 if (l1>=l2 && l1>=l3) {c0=u11; c1=u12; c2=u13; l=l1;}
 else if (l2>=l1 && l2>=l3) {c0=u12; c1=u22; c2=u23; l=l2;}
 else {c0=u13; c1=u23; c2=u33; l=l3;}

 l = sqrt(l); *x = c0/l; *y = c1/l; *z = c2/l;
}


// Computes the eigenvector corresp. to the maximum eigenvalue

void SymMatrix3x3::getMaxEigenvector(double *x, double *y, double *z) const
{
 SymMatrix3x3(-M[0], -M[1], -M[2], -M[3], -M[4], -M[5]).getMinEigenvector(x, y, z);
}

// Prints the matrix values

void SymMatrix3x3::print(FILE *fp) const
{
 fprintf(fp,"%e %e %e\n",M[0],M[1],M[3]);
 fprintf(fp,"%e %e %e\n",M[1],M[2],M[4]);
 fprintf(fp,"%e %e %e\n",M[3],M[4],M[5]);
}


//////////////////////////////////////////////////////////////////////////
//
// Symmetric 4x4 matrix
//
//////////////////////////////////////////////////////////////////////////

// Extend a 3x3 matrix

SymMatrix4x4::SymMatrix4x4(const SymMatrix3x3& q)
{
 a2 = q.M[0]; ab = q.M[1]; ac = q.M[3]; ad = 0;
 b2 = q.M[2]; bc = q.M[4]; bd = 0;
 c2 = q.M[5]; cd = 0; d2 = 1;
}


// Build a quadric

SymMatrix4x4::SymMatrix4x4(const double& a, const double& b, const double& c, const double& d)
{
 a2 = a*a; ab = a*b; ac = a*c; ad = a*d;
 b2 = b*b; bc = b*c; bd = b*d;
 c2 = c*c; cd = c*d; d2 = d*d;
}

// True iff equal

bool SymMatrix4x4::operator==(const SymMatrix4x4& q)
{
 return (
 a2 == q.a2 && ab == q.ab && ac == q.ac && ad == q.ad &&\
 b2 == q.b2 && bc == q.bc && bd == q.bd &&\
 c2 == q.c2 && cd == q.cd && d2 == q.d2);
}

// True iff different

bool SymMatrix4x4::operator!=(const SymMatrix4x4& q)
{
 return (
 a2 != q.a2 || ab != q.ab || ac != q.ac || ad != q.ad ||\
 b2 != q.b2 || bc != q.bc || bd != q.bd ||\
 c2 != q.c2 || cd != q.cd || d2 != q.d2);
}

// Sum the matrix 'q' to the current matrix

void SymMatrix4x4::operator+=(const SymMatrix4x4& q)
{
 a2 += q.a2; ab += q.ab; ac += q.ac; ad += q.ad;
 b2 += q.b2; bc += q.bc; bd += q.bd;
 c2 += q.c2; cd += q.cd; d2 += q.d2;
}

// Returns the sum of the matrix with another matrix 'q'

SymMatrix4x4 SymMatrix4x4::operator+(const SymMatrix4x4& q) const
{
 SymMatrix4x4 n;

 n.a2=a2+q.a2; n.ab=ab+q.ab; n.ac=ac+q.ac; n.ad=ad+q.ad;
 n.b2=b2+q.b2; n.bc=bc+q.bc; n.bd=bd+q.bd;
 n.c2=c2+q.c2; n.cd=cd+q.cd; n.d2=d2+q.d2;

 return n;
}

// Returns the product of the matrix with a scalar

SymMatrix4x4 SymMatrix4x4::operator*(const double& d) const
{
 SymMatrix4x4 n;

 n.a2=a2*d; n.ab=ab*d; n.ac=ac*d; n.ad=ad*d;
 n.b2=b2*d; n.bc=bc*d; n.bd=bd*d;
 n.c2=c2*d; n.cd=cd*d; n.d2=d2*d;

 return n;
}

void SymMatrix4x4::add(const double& a, const double& b, const double& c, const double& d)
{
 a2 += a*a; ab += a*b; ac += a*c; ad += a*d;
 b2 += b*b; bc += b*c; bd += b*d;
 c2 += c*c; cd += c*d; d2 += d*d;
}


// Computes (a,b,c,w)*M*(a,b,c,w)

double SymMatrix4x4::lrMultiply(const double& x, const double& y, const double& z, const double& w) const
{
 double a,b,c,d;
 a = x*a2 + y*ab + z*ac + w*ad;
 b = x*ab + y*b2 + z*bc + w*bd;
 c = x*ac + y*bc + z*c2 + w*cd;
 d = x*ad + y*bd + z*cd + w*d2;
 return (x*a + y*b + z*c + w*d);
}

// Computes (a,b,c) s.t. (a,b,c,1)*M*(a,b,c,1) is minimized
// Returns FALSE if the minimizer is not unique

bool SymMatrix4x4::getMinimizer(double *a, double *b, double *c) const
{
 double det, pos, neg, t;

 pos = neg = 0.0;
 t =  a2*b2*c2; ((t >= 0.0)?(pos):(neg)) += t;
 t =  ab*bc*ac; ((t >= 0.0)?(pos):(neg)) += t;
 t =  ac*ab*bc; ((t >= 0.0)?(pos):(neg)) += t;
 t = -ac*b2*ac; ((t >= 0.0)?(pos):(neg)) += t;
 t = -ab*ab*c2; ((t >= 0.0)?(pos):(neg)) += t;
 t = -a2*bc*bc; ((t >= 0.0)?(pos):(neg)) += t;
 det = pos+neg;

 t = det/(pos-neg);
 if (FABS(t) >= 1.0e-15)
 {
  *a = -(ad*( (b2*c2 - bc*bc)) + bd*(-(ab*c2 - bc*ac)) + cd*( (ab*bc - b2*ac)))/det;
  *b = -(ad*(-(ab*c2 - ac*bc)) + bd*( (a2*c2 - ac*ac)) + cd*(-(a2*bc - ab*ac)))/det;
  *c = -(ad*( (ab*bc - ac*b2)) + bd*(-(a2*bc - ac*ab)) + cd*( (a2*b2 - ab*ab)))/det;
  return 1;
 }

 return 0;
}


// Invert a symmetric 4x4 matrix using L*D*L^T decomposition.
// The calling function is responsible of verifying that the matrix
// is positive definite.
//

bool SymMatrix4x4::invert()
{
 if (a2 <= 0) return false;

 double d00 = 1.0/a2; 
 double L10 = ab;
 double l10 = ab*d00;
 double L20 = ac;
 double l20 = ac*d00;
 double L30 = ad;
 double l30 = ad*d00;
 double d11 = b2-(L10*l10);

 if (d11 <= 0) return false; else d11 = 1.0/d11;

 double L21 = (bc-(L10*l20));
 double l21 = L21*d11;
 double L31 = (bd-(L10*l30));
 double l31 = L31*d11;
 double d22 = c2-(L20*l20)-(L21*l21);

 if (d22 <= 0) return false; else d22 = 1.0/d22;

 double L32 = (cd-(L20*l30)-(L21*l31));
 double l32 = L32*d22;
 double d33 = d2-(L30*l30)-(L31*l31)-(L32*l32);

 if (d33 <= 0) return false; else d33 = 1.0/d33;

 L20 = l10*l21-l20;
 L31 = l21*l32-l31;
 L30 = l32*l20-L31*l10-l30;

 a2 = d00+(-l10)*((-l10)*d11)+L20*(L20*d22)+L30*(L30*d33);
 ab = ((-l10)*d11)+L20*((-l21)*d22)+L30*(L31*d33);
 b2 = d11+(-l21)*((-l21)*d22)+L31*(L31*d33);
 ac = (L20*d22)+L30*((-l32)*d33);
 bc = ((-l21)*d22)+L31*((-l32)*d33);
 c2 = d22+(-l32)*((-l32)*d33);
 ad = (L30*d33);
 bd = (L31*d33);
 cd = ((-l32)*d33);
 d2 = d33;

 return true;
}


//////////////////////////////////////////////////////////////////////////
//
// Generic 4x4 matrix
//
//////////////////////////////////////////////////////////////////////////

Matrix4x4::Matrix4x4() {}

Matrix4x4::Matrix4x4(const double& d)
{
 matrix[0][0] = matrix[1][1] = matrix[2][2] = matrix[3][3] = d;
 matrix[1][0] = matrix[0][1] = matrix[1][2] = matrix[1][3] = 0;
 matrix[2][0] = matrix[2][1] = matrix[0][2] = matrix[2][3] = 0;
 matrix[3][0] = matrix[3][1] = matrix[3][2] = matrix[0][3] = 0;
}

Matrix4x4::Matrix4x4(const double& a11, const double& a12, const double& a13, const double& a14,
		     const double& a21, const double& a22, const double& a23, const double& a24, 
		     const double& a31, const double& a32, const double& a33, const double& a34, 
		     const double& a41, const double& a42, const double& a43, const double& a44)
{
 matrix[0][0] = a11; matrix[0][1] = a12; matrix[0][2] = a13; matrix[0][3] = a14;
 matrix[1][0] = a21; matrix[1][1] = a22; matrix[1][2] = a23; matrix[1][3] = a24;
 matrix[2][0] = a31; matrix[2][1] = a32; matrix[2][2] = a33; matrix[2][3] = a34;
 matrix[3][0] = a41; matrix[3][1] = a42; matrix[3][2] = a43; matrix[3][3] = a44;
}

void Matrix4x4::setRotation(const double& rx, const double& ry, const double& rz, const double& rw)
{
  matrix[0][0] = rw*rw + rx*rx - ry*ry - rz*rz;
  matrix[0][1] = 2*rx*ry + 2*rw*rz;
  matrix[0][2] = 2*rx*rz - 2*rw*ry;
  matrix[0][3] = 0.0;

  matrix[1][0] = 2*rx*ry-2*rw*rz;
  matrix[1][1] = rw*rw - rx*rx + ry*ry - rz*rz;
  matrix[1][2] = 2*ry*rz + 2*rw*rx;
  matrix[1][3] = 0.0;

  matrix[2][0] = 2*rx*rz + 2*rw*ry;
  matrix[2][1] = 2*ry*rz - 2*rw*rx;
  matrix[2][2] = rw*rw - rx*rx - ry*ry + rz*rz;
  matrix[2][3] = 0.0;

  matrix[3][0] = 0.0;
  matrix[3][1] = 0.0;
  matrix[3][2] = 0.0;
  matrix[3][3] = rw*rw + rx*rx + ry*ry + rz*rz;
}

void Matrix4x4::setTranslation(const double& x, const double& y, const double& z)
{
 matrix[0][0] = matrix[1][1] = matrix[2][2] = matrix[3][3] = 1;
 matrix[1][0] = matrix[0][1] = matrix[1][2] = 0;
 matrix[2][0] = matrix[2][1] = matrix[0][2] = 0;
 matrix[3][0] = matrix[3][1] = matrix[3][2] = 0;
 matrix[0][3] = x;
 matrix[1][3] = y;
 matrix[2][3] = z;
}

Matrix4x4 Matrix4x4::operator*(const Matrix4x4& q) const
{
 int i, j;
 Matrix4x4 m;

 for (i=0; i<4; i++) for (j=0; j<4; j++) 
  m.matrix[i][j] = matrix[i][0]*q.matrix[0][j] + matrix[i][1]*q.matrix[1][j] + matrix[i][2]*q.matrix[2][j] + matrix[i][3]*q.matrix[3][j];

 return m;
}

void Matrix4x4::transform(double *x, double *y, double *z)
{
 double w, a = *x, b = *y, c = *z;

 *x = matrix[0][0]*a + matrix[0][1]*b + matrix[0][2]*c + matrix[0][3];
 *y = matrix[1][0]*a + matrix[1][1]*b + matrix[1][2]*c + matrix[1][3];
 *z = matrix[2][0]*a + matrix[2][1]*b + matrix[2][2]*c + matrix[2][3];
 w  = matrix[3][0]*a + matrix[3][1]*b + matrix[3][2]*c + matrix[3][3];

 (*x) /= w; (*y) /= w; (*z) /= w;
}
