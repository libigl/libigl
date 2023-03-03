// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2021 Vladimir S. FONOV <vladimir.fonov@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.

/*
*  
*  C++ version based on the routines published in    
*  "Fast and Robust Triangle-Triangle Overlap Test      
*   Using Orientation Predicates"  P. Guigue - O. Devillers
*  
*  Works with Eigen data structures instead of plain C arrays
*  returns bool values
*  
*  Original notice: 
*
*  Triangle-Triangle Overlap Test Routines        
*  July, 2002                                                          
*  Updated December 2003                                                
*                                                                       
*  This file contains C implementation of algorithms for                
*  performing two and three-dimensional triangle-triangle intersection test 
*  The algorithms and underlying theory are described in                    
*                                                                           
* "Fast and Robust Triangle-Triangle Overlap Test 
*  Using Orientation Predicates"  P. Guigue - O. Devillers
*                                                 
*  Journal of Graphics Tools, 8(1), 2003                                    
*                                                                           
*  Several geometric predicates are defined.  Their parameters are all      
*  points.  Each point is an array of two or three double precision         
*  floating point numbers. The geometric predicates implemented in          
*  this file are:                                                            
*                                                                           
*    int tri_tri_overlap_test_3d(p1,q1,r1,p2,q2,r2)                         
*    int tri_tri_overlap_test_2d(p1,q1,r1,p2,q2,r2)                         
*                                                                           
*    int tri_tri_intersection_test_3d(p1,q1,r1,p2,q2,r2,
*                                     coplanar,source,target)               
*                                                                           
*       is a version that computes the segment of intersection when            
*       the triangles overlap (and are not coplanar)                        
*                                                                           
*    each function returns 1 if the triangles (including their              
*    boundary) intersect, otherwise 0                                       
*                                                                           
*                                                                           
*  Other information are available from the Web page                        
*  http://www.acm.org/jgt/papers/GuigueDevillers03/                         
*                                                                           
*/

#ifndef IGL_TRI_TRI_INTERSECT_CPP
#define IGL_TRI_TRI_INTERSECT_CPP

#include "Guigue2003_tri_tri_intersect.h"

// helper functions
namespace igl {
template <typename DerivedP1,typename DerivedQ1,typename DerivedR1,
typename DerivedP2,typename DerivedQ2,typename DerivedR2,
typename DerivedN1>
IGL_INLINE bool coplanar_tri_tri3d(
  const Eigen::MatrixBase<DerivedP1> &p1, const Eigen::MatrixBase<DerivedQ1> &q1, const Eigen::MatrixBase<DerivedR1> &r1,
  const Eigen::MatrixBase<DerivedP2> &p2, const Eigen::MatrixBase<DerivedQ2> &q2, const Eigen::MatrixBase<DerivedR2> &r2,
  const Eigen::MatrixBase<DerivedN1> &normal_1);

template <typename DerivedP1,typename DerivedQ1,typename DerivedR1,
typename DerivedP2,typename DerivedQ2,typename DerivedR2>
IGL_INLINE bool ccw_tri_tri_intersection_2d(
  const Eigen::MatrixBase<DerivedP1> &p1, const Eigen::MatrixBase<DerivedQ1> &q1, const Eigen::MatrixBase<DerivedR1> &r1,
  const Eigen::MatrixBase<DerivedP2> &p2, const Eigen::MatrixBase<DerivedQ2> &q2, const Eigen::MatrixBase<DerivedR2> &r2);  
}

/* some 3D macros */

#define _IGL_CROSS(dest,v1,v2)                       \
               dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
               dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
               dest[2]=v1[0]*v2[1]-v1[1]*v2[0];

#define _IGL_DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
 


#define _IGL_SUB(dest,v1,v2) dest[0]=v1[0]-v2[0]; \
                        dest[1]=v1[1]-v2[1]; \
                        dest[2]=v1[2]-v2[2]; 


#define _IGL_SCALAR(dest,alpha,v) dest[0] = alpha * v[0]; \
                             dest[1] = alpha * v[1]; \
                             dest[2] = alpha * v[2];



#define _IGL_CHECK_MIN_MAX(p1,q1,r1,p2,q2,r2) {\
  _IGL_SUB(v1,p2,q1)\
  _IGL_SUB(v2,p1,q1)\
  _IGL_CROSS(N1,v1,v2)\
  _IGL_SUB(v1,q2,q1)\
  if (_IGL_DOT(v1,N1) > 0.0f) return false;\
  _IGL_SUB(v1,p2,p1)\
  _IGL_SUB(v2,r1,p1)\
  _IGL_CROSS(N1,v1,v2)\
  _IGL_SUB(v1,r2,p1) \
  if (_IGL_DOT(v1,N1) > 0.0f) return false;\
  else return true; }



/* Permutation in a canonical form of T2's vertices */

#define _IGL_TRI_TRI_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2) { \
  if (dp2 > 0.0f) { \
     if (dq2 > 0.0f) _IGL_CHECK_MIN_MAX(p1,r1,q1,r2,p2,q2) \
     else if (dr2 > 0.0f) _IGL_CHECK_MIN_MAX(p1,r1,q1,q2,r2,p2)\
     else _IGL_CHECK_MIN_MAX(p1,q1,r1,p2,q2,r2) }\
  else if (dp2 < 0.0f) { \
    if (dq2 < 0.0f) _IGL_CHECK_MIN_MAX(p1,q1,r1,r2,p2,q2)\
    else if (dr2 < 0.0f) _IGL_CHECK_MIN_MAX(p1,q1,r1,q2,r2,p2)\
    else _IGL_CHECK_MIN_MAX(p1,r1,q1,p2,q2,r2)\
  } else { \
    if (dq2 < 0.0f) { \
      if (dr2 >= 0.0f)  _IGL_CHECK_MIN_MAX(p1,r1,q1,q2,r2,p2)\
      else _IGL_CHECK_MIN_MAX(p1,q1,r1,p2,q2,r2)\
    } \
    else if (dq2 > 0.0f) { \
      if (dr2 > 0.0f) _IGL_CHECK_MIN_MAX(p1,r1,q1,p2,q2,r2)\
      else  _IGL_CHECK_MIN_MAX(p1,q1,r1,q2,r2,p2)\
    } \
    else  { \
      if (dr2 > 0.0f) _IGL_CHECK_MIN_MAX(p1,q1,r1,r2,p2,q2)\
      else if (dr2 < 0.0f) _IGL_CHECK_MIN_MAX(p1,r1,q1,r2,p2,q2)\
      else return coplanar_tri_tri3d(p1,q1,r1,p2,q2,r2,N1);\
     }}}
  


/*
*
*  Three-dimensional Triangle-Triangle Overlap Test
*
*/

template <typename DerivedP1,typename DerivedQ1,typename DerivedR1,
typename DerivedP2,typename DerivedQ2,typename DerivedR2> 
IGL_INLINE bool igl::tri_tri_overlap_test_3d(
  const Eigen::MatrixBase<DerivedP1> &  p1, 
  const Eigen::MatrixBase<DerivedQ1> &  q1, 
  const Eigen::MatrixBase<DerivedR1> &  r1, 
  const Eigen::MatrixBase<DerivedP2> &  p2, 
  const Eigen::MatrixBase<DerivedQ2> &  q2, 
  const Eigen::MatrixBase<DerivedR2> &  r2)
{
  using Scalar    = typename DerivedP1::Scalar;
  using RowVector = typename Eigen::Matrix<Scalar,1,3>;

  Scalar dp1, dq1, dr1, dp2, dq2, dr2;
  RowVector v1, v2;
  RowVector N1, N2; 
  
  /* Compute distance signs  of p1, q1 and r1 to the plane of
     triangle(p2,q2,r2) */


  _IGL_SUB(v1,p2,r2)
  _IGL_SUB(v2,q2,r2)
  _IGL_CROSS(N2,v1,v2)

  _IGL_SUB(v1,p1,r2)
  dp1 = _IGL_DOT(v1,N2);
  _IGL_SUB(v1,q1,r2)
  dq1 = _IGL_DOT(v1,N2);
  _IGL_SUB(v1,r1,r2)
  dr1 = _IGL_DOT(v1,N2);
  
  if (((dp1 * dq1) > 0.0f) && ((dp1 * dr1) > 0.0f))  return false; 

  /* Compute distance signs  of p2, q2 and r2 to the plane of
     triangle(p1,q1,r1) */

  _IGL_SUB(v1,q1,p1)
  _IGL_SUB(v2,r1,p1)
  _IGL_CROSS(N1,v1,v2)

  _IGL_SUB(v1,p2,r1)
  dp2 = _IGL_DOT(v1,N1);
  _IGL_SUB(v1,q2,r1)
  dq2 = _IGL_DOT(v1,N1);
  _IGL_SUB(v1,r2,r1)
  dr2 = _IGL_DOT(v1,N1);
  
  if (((dp2 * dq2) > 0.0f) && ((dp2 * dr2) > 0.0f)) return false;

  /* Permutation in a canonical form of T1's vertices */


  if (dp1 > 0.0f) {
    if (dq1 > 0.0f) _IGL_TRI_TRI_3D(r1,p1,q1,p2,r2,q2,dp2,dr2,dq2)
    else if (dr1 > 0.0f) _IGL_TRI_TRI_3D(q1,r1,p1,p2,r2,q2,dp2,dr2,dq2)  
    else _IGL_TRI_TRI_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2)
  } else if (dp1 < 0.0f) {
    if (dq1 < 0.0f) _IGL_TRI_TRI_3D(r1,p1,q1,p2,q2,r2,dp2,dq2,dr2)
    else if (dr1 < 0.0f) _IGL_TRI_TRI_3D(q1,r1,p1,p2,q2,r2,dp2,dq2,dr2)
    else _IGL_TRI_TRI_3D(p1,q1,r1,p2,r2,q2,dp2,dr2,dq2)
  } else {
    if (dq1 < 0.0f) {
      if (dr1 >= 0.0f) _IGL_TRI_TRI_3D(q1,r1,p1,p2,r2,q2,dp2,dr2,dq2)
      else _IGL_TRI_TRI_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2)
    }
    else if (dq1 > 0.0f) {
      if (dr1 > 0.0f) _IGL_TRI_TRI_3D(p1,q1,r1,p2,r2,q2,dp2,dr2,dq2)
      else _IGL_TRI_TRI_3D(q1,r1,p1,p2,q2,r2,dp2,dq2,dr2)
    }
    else  {
      if (dr1 > 0.0f) _IGL_TRI_TRI_3D(r1,p1,q1,p2,q2,r2,dp2,dq2,dr2)
      else if (dr1 < 0.0f) _IGL_TRI_TRI_3D(r1,p1,q1,p2,r2,q2,dp2,dr2,dq2)
      else return coplanar_tri_tri3d(p1,q1,r1,p2,q2,r2,N1);
    }
  }
};



template <typename DerivedP1,typename DerivedQ1,typename DerivedR1,
typename DerivedP2,typename DerivedQ2,typename DerivedR2,
typename DerivedN1>
IGL_INLINE bool igl::coplanar_tri_tri3d(
  const Eigen::MatrixBase<DerivedP1> &p1, const Eigen::MatrixBase<DerivedQ1> &q1, const Eigen::MatrixBase<DerivedR1> &r1,
  const Eigen::MatrixBase<DerivedP2> &p2, const Eigen::MatrixBase<DerivedQ2> &q2, const Eigen::MatrixBase<DerivedR2> &r2,
  const Eigen::MatrixBase<DerivedN1> &normal_1)
{

  using Scalar= typename DerivedP1::Scalar;
  using RowVector2D = typename Eigen::Matrix<Scalar,1,2>;

  RowVector2D P1,Q1,R1;
  RowVector2D P2,Q2,R2;

  Scalar n_x, n_y, n_z;

  n_x = ((normal_1[0]<0.0)?-normal_1[0]:normal_1[0]);
  n_y = ((normal_1[1]<0.0)?-normal_1[1]:normal_1[1]);
  n_z = ((normal_1[2]<0.0)?-normal_1[2]:normal_1[2]);


  /* Projection of the triangles in 3D onto 2D such that the area of
     the projection is maximized. */


  if (( n_x > n_z ) && ( n_x >= n_y )) {
    // Project onto plane YZ

      P1[0] = q1[2]; P1[1] = q1[1];
      Q1[0] = p1[2]; Q1[1] = p1[1];
      R1[0] = r1[2]; R1[1] = r1[1]; 
    
      P2[0] = q2[2]; P2[1] = q2[1];
      Q2[0] = p2[2]; Q2[1] = p2[1];
      R2[0] = r2[2]; R2[1] = r2[1]; 

  } else if (( n_y > n_z ) && ( n_y >= n_x )) {
    // Project onto plane XZ

    P1[0] = q1[0]; P1[1] = q1[2];
    Q1[0] = p1[0]; Q1[1] = p1[2];
    R1[0] = r1[0]; R1[1] = r1[2]; 
 
    P2[0] = q2[0]; P2[1] = q2[2];
    Q2[0] = p2[0]; Q2[1] = p2[2];
    R2[0] = r2[0]; R2[1] = r2[2]; 
    
  } else {
    // Project onto plane XY

    P1[0] = p1[0]; P1[1] = p1[1]; 
    Q1[0] = q1[0]; Q1[1] = q1[1]; 
    R1[0] = r1[0]; R1[1] = r1[1]; 
    
    P2[0] = p2[0]; P2[1] = p2[1]; 
    Q2[0] = q2[0]; Q2[1] = q2[1]; 
    R2[0] = r2[0]; R2[1] = r2[1]; 
  }

  return tri_tri_overlap_test_2d(P1,Q1,R1,P2,Q2,R2);
    
};



/*
*                                                                
*  Three-dimensional Triangle-Triangle Intersection              
*
*/

/*
   This macro is called when the triangles surely intersect
   It constructs the segment of intersection of the two triangles
   if they are not coplanar.
*/

// NOTE: a faster, but possibly less precise, method of computing
// point B is described here: https://github.com/erich666/jgt-code/issues/5

#define _IGL_CONSTRUCT_INTERSECTION(p1,q1,r1,p2,q2,r2) { \
  _IGL_SUB(v1,q1,p1) \
  _IGL_SUB(v2,r2,p1) \
  _IGL_CROSS(N,v1,v2) \
  _IGL_SUB(v,p2,p1) \
  if (_IGL_DOT(v,N) > 0.0f) {\
    _IGL_SUB(v1,r1,p1) \
    _IGL_CROSS(N,v1,v2) \
    if (_IGL_DOT(v,N) <= 0.0f) { \
      _IGL_SUB(v2,q2,p1) \
      _IGL_CROSS(N,v1,v2) \
      if (_IGL_DOT(v,N) > 0.0f) { \
  _IGL_SUB(v1,p1,p2) \
  _IGL_SUB(v2,p1,r1) \
  alpha = _IGL_DOT(v1,N2) / _IGL_DOT(v2,N2); \
  _IGL_SCALAR(v1,alpha,v2) \
  _IGL_SUB(source,p1,v1) \
  _IGL_SUB(v1,p2,p1) \
  _IGL_SUB(v2,p2,r2) \
  alpha = _IGL_DOT(v1,N1) / _IGL_DOT(v2,N1); \
  _IGL_SCALAR(v1,alpha,v2) \
  _IGL_SUB(target,p2,v1) \
  return true; \
      } else { \
  _IGL_SUB(v1,p2,p1) \
  _IGL_SUB(v2,p2,q2) \
  alpha = _IGL_DOT(v1,N1) / _IGL_DOT(v2,N1); \
  _IGL_SCALAR(v1,alpha,v2) \
  _IGL_SUB(source,p2,v1) \
  _IGL_SUB(v1,p2,p1) \
  _IGL_SUB(v2,p2,r2) \
  alpha = _IGL_DOT(v1,N1) / _IGL_DOT(v2,N1); \
  _IGL_SCALAR(v1,alpha,v2) \
  _IGL_SUB(target,p2,v1) \
  return true; \
      } \
    } else { \
      return false; \
    } \
  } else { \
    _IGL_SUB(v2,q2,p1) \
    _IGL_CROSS(N,v1,v2) \
    if (_IGL_DOT(v,N) < 0.0f) { \
      return false; \
    } else { \
      _IGL_SUB(v1,r1,p1) \
      _IGL_CROSS(N,v1,v2) \
      if (_IGL_DOT(v,N) >= 0.0f) { \
  _IGL_SUB(v1,p1,p2) \
  _IGL_SUB(v2,p1,r1) \
  alpha = _IGL_DOT(v1,N2) / _IGL_DOT(v2,N2); \
  _IGL_SCALAR(v1,alpha,v2) \
  _IGL_SUB(source,p1,v1) \
  _IGL_SUB(v1,p1,p2) \
  _IGL_SUB(v2,p1,q1) \
  alpha = _IGL_DOT(v1,N2) / _IGL_DOT(v2,N2); \
  _IGL_SCALAR(v1,alpha,v2) \
  _IGL_SUB(target,p1,v1) \
  return true; \
      } else { \
  _IGL_SUB(v1,p2,p1) \
  _IGL_SUB(v2,p2,q2) \
  alpha = _IGL_DOT(v1,N1) / _IGL_DOT(v2,N1); \
  _IGL_SCALAR(v1,alpha,v2) \
  _IGL_SUB(source,p2,v1) \
  _IGL_SUB(v1,p1,p2) \
  _IGL_SUB(v2,p1,q1) \
  alpha = _IGL_DOT(v1,N2) / _IGL_DOT(v2,N2); \
  _IGL_SCALAR(v1,alpha,v2) \
  _IGL_SUB(target,p1,v1) \
  return true; \
      }}}} 

                

#define _IGL_TRI_TRI_INTER_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2) { \
  if (dp2 > 0.0f) { \
     if (dq2 > 0.0f) _IGL_CONSTRUCT_INTERSECTION(p1,r1,q1,r2,p2,q2) \
     else if (dr2 > 0.0f) _IGL_CONSTRUCT_INTERSECTION(p1,r1,q1,q2,r2,p2)\
     else _IGL_CONSTRUCT_INTERSECTION(p1,q1,r1,p2,q2,r2) }\
  else if (dp2 < 0.0f) { \
    if (dq2 < 0.0f) _IGL_CONSTRUCT_INTERSECTION(p1,q1,r1,r2,p2,q2)\
    else if (dr2 < 0.0f) _IGL_CONSTRUCT_INTERSECTION(p1,q1,r1,q2,r2,p2)\
    else _IGL_CONSTRUCT_INTERSECTION(p1,r1,q1,p2,q2,r2)\
  } else { \
    if (dq2 < 0.0f) { \
      if (dr2 >= 0.0f)  _IGL_CONSTRUCT_INTERSECTION(p1,r1,q1,q2,r2,p2)\
      else _IGL_CONSTRUCT_INTERSECTION(p1,q1,r1,p2,q2,r2)\
    } \
    else if (dq2 > 0.0f) { \
      if (dr2 > 0.0f) _IGL_CONSTRUCT_INTERSECTION(p1,r1,q1,p2,q2,r2)\
      else  _IGL_CONSTRUCT_INTERSECTION(p1,q1,r1,q2,r2,p2)\
    } \
    else  { \
      if (dr2 > 0.0f) _IGL_CONSTRUCT_INTERSECTION(p1,q1,r1,r2,p2,q2)\
      else if (dr2 < 0.0f) _IGL_CONSTRUCT_INTERSECTION(p1,r1,q1,r2,p2,q2)\
      else { \
        coplanar = true; \
  return coplanar_tri_tri3d(p1,q1,r1,p2,q2,r2,N1);\
     } \
  }} }
  

/*
   The following version computes the segment of intersection of the
   two triangles if it exists. 
   coplanar returns whether the triangles are coplanar
   source and target are the endpoints of the line segment of intersection 
*/

template <typename DerivedP1,typename DerivedQ1,typename DerivedR1,
typename DerivedP2,typename DerivedQ2,typename DerivedR2,
typename DerivedS,typename DerivedT>
IGL_INLINE bool igl::tri_tri_intersection_test_3d(
    const Eigen::MatrixBase<DerivedP1> & p1, const Eigen::MatrixBase<DerivedQ1> & q1, const Eigen::MatrixBase<DerivedR1> & r1, 
    const Eigen::MatrixBase<DerivedP2> & p2, const Eigen::MatrixBase<DerivedQ2> & q2, const Eigen::MatrixBase<DerivedR2> & r2,
    bool & coplanar, 
    Eigen::MatrixBase<DerivedS> & source, 
    Eigen::MatrixBase<DerivedT> & target )        
{
  using Scalar= typename DerivedP1::Scalar;
  using RowVector3D = typename Eigen::Matrix<Scalar,1,3>;

  Scalar dp1, dq1, dr1, dp2, dq2, dr2;
  RowVector3D v1, v2, v;
  RowVector3D N1, N2, N;
  Scalar alpha;
  // Compute distance signs  of p1, q1 and r1 
  // to the plane of triangle(p2,q2,r2)


  _IGL_SUB(v1,p2,r2)
  _IGL_SUB(v2,q2,r2)
  _IGL_CROSS(N2,v1,v2)

  _IGL_SUB(v1,p1,r2)
  dp1 = _IGL_DOT(v1,N2);
  _IGL_SUB(v1,q1,r2)
  dq1 = _IGL_DOT(v1,N2);
  _IGL_SUB(v1,r1,r2)
  dr1 = _IGL_DOT(v1,N2);
  
  coplanar = false;

  if (((dp1 * dq1) > 0.0f) && ((dp1 * dr1) > 0.0f))  return false; 

  // Compute distance signs  of p2, q2 and r2 
  // to the plane of triangle(p1,q1,r1)

  
  _IGL_SUB(v1,q1,p1)
  _IGL_SUB(v2,r1,p1)
  _IGL_CROSS(N1,v1,v2)

  _IGL_SUB(v1,p2,r1)
  dp2 = _IGL_DOT(v1,N1);
  _IGL_SUB(v1,q2,r1)
  dq2 = _IGL_DOT(v1,N1);
  _IGL_SUB(v1,r2,r1)
  dr2 = _IGL_DOT(v1,N1);
  
  if (((dp2 * dq2) > 0.0f) && ((dp2 * dr2) > 0.0f)) return false;

  // Permutation in a canonical form of T1's vertices


  if (dp1 > 0.0f) {
    if (dq1 > 0.0f) _IGL_TRI_TRI_INTER_3D(r1,p1,q1,p2,r2,q2,dp2,dr2,dq2)
    else if (dr1 > 0.0f) _IGL_TRI_TRI_INTER_3D(q1,r1,p1,p2,r2,q2,dp2,dr2,dq2)
  
    else _IGL_TRI_TRI_INTER_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2)
  } else if (dp1 < 0.0f) {
    if (dq1 < 0.0f) _IGL_TRI_TRI_INTER_3D(r1,p1,q1,p2,q2,r2,dp2,dq2,dr2)
    else if (dr1 < 0.0f) _IGL_TRI_TRI_INTER_3D(q1,r1,p1,p2,q2,r2,dp2,dq2,dr2)
    else _IGL_TRI_TRI_INTER_3D(p1,q1,r1,p2,r2,q2,dp2,dr2,dq2)
  } else {
    if (dq1 < 0.0f) {
      if (dr1 >= 0.0f) _IGL_TRI_TRI_INTER_3D(q1,r1,p1,p2,r2,q2,dp2,dr2,dq2)
      else _IGL_TRI_TRI_INTER_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2)
    }
    else if (dq1 > 0.0f) {
      if (dr1 > 0.0f) _IGL_TRI_TRI_INTER_3D(p1,q1,r1,p2,r2,q2,dp2,dr2,dq2)
      else _IGL_TRI_TRI_INTER_3D(q1,r1,p1,p2,q2,r2,dp2,dq2,dr2)
    }
    else  {
      if (dr1 > 0.0f) _IGL_TRI_TRI_INTER_3D(r1,p1,q1,p2,q2,r2,dp2,dq2,dr2)
      else if (dr1 < 0.0f) _IGL_TRI_TRI_INTER_3D(r1,p1,q1,p2,r2,q2,dp2,dr2,dq2)
      else {
  // triangles are co-planar

  coplanar = true;
  return coplanar_tri_tri3d(p1,q1,r1,p2,q2,r2,N1);
      }
    }
  }
};





/*
*
*  Two dimensional Triangle-Triangle Overlap Test    
*
*/


/* some 2D macros */

#define _IGL_ORIENT_2D(a, b, c)  ((a[0]-c[0])*(b[1]-c[1])-(a[1]-c[1])*(b[0]-c[0]))


#define INTERSECTION_TEST_VERTEX(P1, Q1, R1, P2, Q2, R2) {\
  if (_IGL_ORIENT_2D(R2,P2,Q1) >= 0.0f)\
    if (_IGL_ORIENT_2D(R2,Q2,Q1) <= 0.0f)\
      if (_IGL_ORIENT_2D(P1,P2,Q1) > 0.0f) {\
  if (_IGL_ORIENT_2D(P1,Q2,Q1) <= 0.0f) return true; \
  else return false;} else {\
  if (_IGL_ORIENT_2D(P1,P2,R1) >= 0.0f)\
    if (_IGL_ORIENT_2D(Q1,R1,P2) >= 0.0f) return true; \
    else return false;\
  else return false;}\
    else \
      if (_IGL_ORIENT_2D(P1,Q2,Q1) <= 0.0f)\
  if (_IGL_ORIENT_2D(R2,Q2,R1) <= 0.0f)\
    if (_IGL_ORIENT_2D(Q1,R1,Q2) >= 0.0f) return true; \
    else return false;\
  else return false;\
      else return false;\
  else\
    if (_IGL_ORIENT_2D(R2,P2,R1) >= 0.0f) \
      if (_IGL_ORIENT_2D(Q1,R1,R2) >= 0.0f)\
  if (_IGL_ORIENT_2D(P1,P2,R1) >= 0.0f) return true;\
  else return false;\
      else \
  if (_IGL_ORIENT_2D(Q1,R1,Q2) >= 0.0f) {\
    if (_IGL_ORIENT_2D(R2,R1,Q2) >= 0.0f) return true; \
    else return false; }\
  else return false; \
    else  return false; \
 };



#define _IGL_INTERSECTION_TEST_EDGE(P1, Q1, R1, P2, Q2, R2) { \
  if (_IGL_ORIENT_2D(R2,P2,Q1) >= 0.0f) {\
    if (_IGL_ORIENT_2D(P1,P2,Q1) >= 0.0f) { \
        if (_IGL_ORIENT_2D(P1,Q1,R2) >= 0.0f) return true; \
        else return false;} else { \
      if (_IGL_ORIENT_2D(Q1,R1,P2) >= 0.0f){ \
  if (_IGL_ORIENT_2D(R1,P1,P2) >= 0.0f) return true; else return false;} \
      else return false; } \
  } else {\
    if (_IGL_ORIENT_2D(R2,P2,R1) >= 0.0f) {\
      if (_IGL_ORIENT_2D(P1,P2,R1) >= 0.0f) {\
  if (_IGL_ORIENT_2D(P1,R1,R2) >= 0.0f) return true;  \
  else {\
    if (_IGL_ORIENT_2D(Q1,R1,R2) >= 0.0f) return true; else return false;}}\
      else  return false; }\
    else return false; }}




template <typename DerivedP1,typename DerivedQ1,typename DerivedR1,
typename DerivedP2,typename DerivedQ2,typename DerivedR2>
IGL_INLINE bool igl::ccw_tri_tri_intersection_2d(
  const Eigen::MatrixBase<DerivedP1> &p1, const Eigen::MatrixBase<DerivedQ1> &q1, const Eigen::MatrixBase<DerivedR1> &r1,
  const Eigen::MatrixBase<DerivedP2> &p2, const Eigen::MatrixBase<DerivedQ2> &q2, const Eigen::MatrixBase<DerivedR2> &r2) 
{
  if ( _IGL_ORIENT_2D(p2,q2,p1) >= 0.0f ) {
    if ( _IGL_ORIENT_2D(q2,r2,p1) >= 0.0f ) {
      if ( _IGL_ORIENT_2D(r2,p2,p1) >= 0.0f ) return true;
      else _IGL_INTERSECTION_TEST_EDGE(p1,q1,r1,p2,q2,r2)
    } else {  
      if ( _IGL_ORIENT_2D(r2,p2,p1) >= 0.0f ) 
  _IGL_INTERSECTION_TEST_EDGE(p1,q1,r1,r2,p2,q2)
      else INTERSECTION_TEST_VERTEX(p1,q1,r1,p2,q2,r2)}}
  else {
    if ( _IGL_ORIENT_2D(q2,r2,p1) >= 0.0f ) {
      if ( _IGL_ORIENT_2D(r2,p2,p1) >= 0.0f ) 
  _IGL_INTERSECTION_TEST_EDGE(p1,q1,r1,q2,r2,p2)
      else  INTERSECTION_TEST_VERTEX(p1,q1,r1,q2,r2,p2)}
    else INTERSECTION_TEST_VERTEX(p1,q1,r1,r2,p2,q2)}
};

template <typename DerivedP1,typename DerivedQ1,typename DerivedR1,
typename DerivedP2,typename DerivedQ2,typename DerivedR2>
IGL_INLINE bool igl::tri_tri_overlap_test_2d(
  const Eigen::MatrixBase<DerivedP1> &p1, const Eigen::MatrixBase<DerivedQ1> &q1, const Eigen::MatrixBase<DerivedR1> &r1,
  const Eigen::MatrixBase<DerivedP2> &p2, const Eigen::MatrixBase<DerivedQ2> &q2, const Eigen::MatrixBase<DerivedR2> &r2) 
{
  if ( _IGL_ORIENT_2D(p1,q1,r1) < 0.0f )
    if ( _IGL_ORIENT_2D(p2,q2,r2) < 0.0f )
      return ccw_tri_tri_intersection_2d(p1,r1,q1,p2,r2,q2);
    else
      return ccw_tri_tri_intersection_2d(p1,r1,q1,p2,q2,r2);
  else
    if ( _IGL_ORIENT_2D(p2,q2,r2) < 0.0f )
      return ccw_tri_tri_intersection_2d(p1,q1,r1,p2,r2,q2);
    else
      return ccw_tri_tri_intersection_2d(p1,q1,r1,p2,q2,r2);

};

//cleanup
#undef _IGL_ORIENT_2D
#undef _IGL_INTERSECTION_TEST_VERTEX
#undef _IGL_TRI_TRI_INTER_3D
#undef _IGL_SUB
#undef _IGL_DOT
#undef _IGL_CROSS
#undef _IGL_SCALAR
#undef _IGL_CHECK_MIN_MAX
#undef _IGL_INTERSECTION_TEST_EDGE


#endif //IGL_TRI_TRI_INTERSECT_CPP


#ifdef IGL_STATIC_LIBRARY
template bool igl::tri_tri_intersection_test_3d<Eigen::Block<Eigen::Matrix<double, -1, -1, 0, -1, -1> const, 1, -1, false>, Eigen::Block<Eigen::Matrix<double, -1, -1, 0, -1, -1> const, 1, -1, false>, Eigen::Block<Eigen::Matrix<double, -1, -1, 0, -1, -1> const, 1, -1, false>, Eigen::Block<Eigen::Matrix<double, -1, -1, 0, -1, -1> const, 1, -1, false>, Eigen::Block<Eigen::Matrix<double, -1, -1, 0, -1, -1> const, 1, -1, false>, Eigen::Block<Eigen::Matrix<double, -1, -1, 0, -1, -1> const, 1, -1, false>, Eigen::Matrix<double, 1, 3, 1, 1, 3>, Eigen::Matrix<double, 1, 3, 1, 1, 3> >(Eigen::MatrixBase<Eigen::Block<Eigen::Matrix<double, -1, -1, 0, -1, -1> const, 1, -1, false> > const&, Eigen::MatrixBase<Eigen::Block<Eigen::Matrix<double, -1, -1, 0, -1, -1> const, 1, -1, false> > const&, Eigen::MatrixBase<Eigen::Block<Eigen::Matrix<double, -1, -1, 0, -1, -1> const, 1, -1, false> > const&, Eigen::MatrixBase<Eigen::Block<Eigen::Matrix<double, -1, -1, 0, -1, -1> const, 1, -1, false> > const&, Eigen::MatrixBase<Eigen::Block<Eigen::Matrix<double, -1, -1, 0, -1, -1> const, 1, -1, false> > const&, Eigen::MatrixBase<Eigen::Block<Eigen::Matrix<double, -1, -1, 0, -1, -1> const, 1, -1, false> > const&, bool&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> >&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> >&);
#endif
