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

#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <float.h>
#include "list.h"

//////////////////////////////////////////////////////////////////////////
//
// Generic 3x3 matrix
//
//////////////////////////////////////////////////////////////////////////

//! Generic 3x3 matrix.

//! Elements are stored in a row-dominant order, thus
//! for example, M[4] is the first element of the second row.

class Matrix3x3
{
 public:
 double M[9];	//!< Actual values of the matrix

 Matrix3x3() {M[0]=M[1]=M[2]=M[3]=M[4]=M[5]=M[6]=M[7]=M[8]=0.0;} //!< Contructs a null matrix

 //! Constructs a fully initialized matrix.
 Matrix3x3(const double& a11, const double& a12, const double& a13,
           const double& a21, const double& a22, const double& a23, 
           const double& a31, const double& a32, const double& a33);

 //! Constructs a 3x3 matrix as the product of Transpose(v1,v2,v3) and (w1,w2,w3).
 Matrix3x3(const double& v1, const double& v2, const double& v3,
           const double& w1, const double& w2, const double& w3); 

 //! Constructs a 3x3 matrix as the product of Transpose(a,b,c) and (a,b,c).
 Matrix3x3(const double& a, const double& b, const double& c);

 //! Returns TRUE if the matrix is symmetric
 bool isSymmetric() const {return (M[2]==M[4] && M[3]==M[7] && M[6]==M[8]);}


 //! Initializes all elements to 'd'
 void operator=(const double& d) {M[0]=M[1]=M[2]=M[3]=M[4]=M[5]=M[6]=M[7]=M[8]=d;}

 void operator+=(const Matrix3x3&); //!< Sum another matrix
 void operator-=(const Matrix3x3&); //!< Subtract another matrix
 void operator*=(const double&);    //!< Multiply by a scalar
 void operator/=(const double& d) {operator *=(1.0/d);}       //!< Divide by a scalar
 Matrix3x3 operator+(const Matrix3x3&) const; //!< Returns the sum of this and another matrix
 Matrix3x3 operator*(const double&) const;    //!< Returns the product of this matrix with a scalar
 Matrix3x3 operator*(const Matrix3x3&) const; //!< Returns the product of this and another matrix (rows by columns)
 Matrix3x3 operator~() const;		 //!< Returns the transpose of this matrix

 //! Return the matrix transpose
 Matrix3x3 transpose() const;

 //! Returns Transpose(a,b,c)*M*(a,b,c)

 //! Returns the (scalar) result of multiplying the matrix on
 //! the left and on the right by the vector (a,b,c).
 double lrMultiply(const double& a, const double& b, const double& c) const;

 //! Returns the (scalar) result of v*M*w
 double lrMultiply(const double& v1, const double& v2, const double& v3,
                   const double& w1, const double& w2, const double& w3) const;
};


//////////////////////////////////////////////////////////////////////////
//
// Symmetric 3x3 matrix
//
//////////////////////////////////////////////////////////////////////////

//! Symmetric 3x3 matrix

//! Compact storage:	\n
//! M[0] M[1] M[3]	\n
//! M[1] M[2] M[4]	\n
//! M[3] M[4] M[5]	\n

class SymMatrix3x3
{
 public:
 double M[6];	//!< Actual values of the matrix

 SymMatrix3x3() {M[0]=M[1]=M[2]=M[3]=M[4]=M[5]=0.0;} //!< Constructs a null matrix.

 //! Constructs a fully initialized matrix.
 SymMatrix3x3(const double&a11, const double&a12, const double&a22,
              const double&a13, const double&a23, const double&a33);

 //! Constructs a symmetric matrix as the product of Transpose(a,b,c) and (a,b,c).
 SymMatrix3x3(const double& a, const double& b, const double& c);

 //! Constructs a symmetric 3x3 matrix as a copy of an existing 3x3 matrix S.

 //! If 'S' is not symmetric, its upper triangular part is reflected to the lower.

 SymMatrix3x3(const Matrix3x3& S);

 bool operator==(const SymMatrix3x3& s) //!< True iff all entries are equal
  {return (M[0]==s.M[0] && M[1]==s.M[1] && M[2]==s.M[2] && M[3]==s.M[3] && M[4]==s.M[4] && M[5]==s.M[5]);}
 bool operator!=(const SymMatrix3x3& s) //!< True iff at least one different entry
  {return (M[0]!=s.M[0] || M[1]!=s.M[1] || M[2]!=s.M[2] || M[3]!=s.M[3] || M[4]!=s.M[4] || M[5]!=s.M[5]);}

 void operator+=(const SymMatrix3x3&); //!< Sum another matrix
 void operator-=(const SymMatrix3x3&); //!< Subtract another matrix
 void operator*=(const double&);       //!< Multiply by a scalar
 void operator/=(const double& d) {operator *=(1.0/d);}       //!< Divide by a scalar
 SymMatrix3x3 operator+(const SymMatrix3x3&) const; //!< Returns the sum of this and another matrix
 SymMatrix3x3 operator*(const double&) const; //!< Returns the product of this matrix with a scalar

 //! Initializes all elements to 'd'
 void operator=(const double& d) {M[0]=M[1]=M[2]=M[3]=M[4]=M[5]=d;}

 //! Returns the determinant
 double determinant() const {return (M[0]*M[2]*M[5])+(2.0*M[1]*M[3]*M[4])-(M[0]*M[4]*M[4])-(M[2]*M[3]*M[3])-(M[5]*M[1]*M[1]);}

 //! Returns TRUE iff the matrix is made of all zeroes.
 bool isNull() const {return (M[0]==0 && M[1]==0 && M[2]==0 && M[3]==0 && M[4]==0 && M[5]==0);}

 //! Returns Transpose(a,b,c)*M*(a,b,c)

 //! Returns the (scalar) result of multiplying the matrix on
 //! the left and on the right by the vector (a,b,c).
 double lrMultiply(const double& a, const double& b, const double& c) const;


 //! Returns the (scalar) result of v*M*w
 double lrMultiply(const double& v1, const double& v2, const double& v3,
                   const double& w1, const double& w2, const double& w3) const;

 bool invert(); //!< Inverts this matrix. Returns FALSE if not invertible, TRUE otherwise

 //! Returns the matrix trace
 double trace() const {return M[0]+M[2]+M[5];}

 //! Compute eigenvalues and eigenvectors of the matrix (Jacobi method).

 //! The calling function is responsible of verifying that the matrix
 //! is diagonalizable. Also, eigen_vals and eigen_vecs must be allocated
 //! prior to calling this method.\n
 //! eigen_vals are sorted in ascending order (i.e., eigen_vals[0] is the smallest one).\n
 //! eigen_vecs are sorted accordingly to the order of eigen_vals, that is,
 //! (eigen_vecs[0], eigen_vecs[1], eigen_vecs[2]) is the eigenvector corresponding to
 //! eigen_vals[0].
 void diagonalize(double eigen_vals[3], double eigen_vecs[9]) const;

 //! Compute the eigenvalues l1, l2 and l3.

 //! This method is much faster and precise than 'diagonalize', as it uses
 //! an analytical direct method instead of an iterative approach.
 void getEigenvalues(double *l1, double *l2, double *l3) const;

 //! Compute the eigenvector (a,b,c) corresponding to the minimum eigenvalue.

 //! This method is much faster and precise than 'diagonalize', as it uses
 //! an analytical direct method instead of an iterative approach.
 void getMinEigenvector(double *a, double *b, double *c) const;

 //! This method is much faster and precise than 'diagonalize', as it uses
 //! an analytical direct method instead of an iterative approach.
 void getMaxEigenvector(double *a, double *b, double *c) const;

 //! Prints the contents of the matrix to the specified FILE id.

 //! If no 'id' is specifyed, results are printed to stdout.

 void print(FILE *id =stdout) const;
};


//////////////////////////////////////////////////////////////////////////
//
// Symmetric 4x4 matrix
//
//////////////////////////////////////////////////////////////////////////

//! Symmetric 4x4 matrix.

//! Compact storage:	\n
//! a2 ab ac ad		\n
//! ab b2 bc bd		\n
//! ac bc c2 cd		\n
//! ad bd cd d2		\n

class SymMatrix4x4
{
 public:
 double a2,ab,ac,ad,b2,bc,bd,c2,cd,d2;	 	    //!< Actual matrix coeffs.

 SymMatrix4x4() {a2=ab=ac=ad=b2=bc=bd=c2=cd=d2=0;}  //!< Constructs a null matrix

 //! Extend a 3x3 symmetric matrix to homogeneous coordinates (ad=bd=cd=0 and d2=1).
 SymMatrix4x4(const SymMatrix3x3&);

 //! Quadric (a,b,c,d)*Transpose(a,b,c,d)
 SymMatrix4x4(const double& a, const double& b, const double& c, const double& d);

 bool operator==(const SymMatrix4x4&);	 	    //!< True iff equal
 bool operator!=(const SymMatrix4x4&);	 	    //!< True iff not equal

 void operator+=(const SymMatrix4x4&);	 	    //!< Sum another matrix
 SymMatrix4x4 operator+(const SymMatrix4x4&) const; //!< Returns the sum of this and another matrix
 SymMatrix4x4 operator*(const double&) const;	    //!< Returns the product of this matrix by a scalar

 //! Adds the quadric (a,b,c,d)*Transpose(a,b,c,d)
 void add(const double& a, const double& b, const double& c, const double& d);  

 //! Returns Transpose(a,b,c,d)*M*(a,b,c,d)

 //! Returns the (scalar) result of multiplying the matrix on
 //! the left and on the right by the vector (a,b,c,d).
 double lrMultiply(const double& a, const double& b, const double& c, const double& d) const;

 //! \brief Computes the vector (a,b,c) that minimizes the quantity lrMultiply(a,b,c,1).
 //! Returns FALSE if such a vector is not unique.
 bool getMinimizer(double *a, double *b, double *c) const;

 bool invert();	//!< Inverts this matrix. Returns FALSE if not invertible, TRUE otherwise
};


//////////////////////////////////////////////////////////////////////////
//
// Generic 4x4 matrix
//
//////////////////////////////////////////////////////////////////////////

//! Generic 4x4 matrix.

class Matrix4x4
{
 public:
 double matrix[4][4]; //!< Actual matrix coefficients

 //! Constructs an undefined matrix
 Matrix4x4();

 //! Constructs a diagonal matrix with 'd' values on the diagonal
 Matrix4x4(const double& d);

 //! Constructs a fully initialized matrix (parameters are in row dominant order M[0][0], M[0][1], ...).
 Matrix4x4(
	   const double&, const double&, const double&, const double&,
	   const double&, const double&, const double&, const double&,
	   const double&, const double&, const double&, const double&,
	   const double&, const double&, const double&, const double&
	  );

 //! Rotation matrix from a quaternion
 void setRotation(const double &, const double&, const double&, const double&);

 //! Translation matrix from a vector
 void setTranslation(const double &, const double&, const double&);

 Matrix4x4 operator*(const Matrix4x4&) const; //!< Returns the product of this and another matrix (rows by columns)

 void transform(double *, double *, double *); //! Transform the vector by left-multiplication with the matrix
};


#endif // MATRIX_H

