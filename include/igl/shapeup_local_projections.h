// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_SHAPEUP_LOCAL_PROJECTIONS_H
#define IGL_SHAPEUP_LOCAL_PROJECTIONS_H

#include <igl/igl_inline.h>
#include <igl/setdiff.h>
#include <igl/cat.h>
#include <Eigen/Core>
#include <vector>




namespace igl
{

	//Every function here defines a local projection for ShapeUp, and must have the following structure to qualify:
	//Input:
	//	P		#P by 3				the set of points, either the initial solution, or from previous iteration.
	//  SC		#Set by 1           cardinalities of sets in S
	//  S		#Sets by max(SC)    independent sets where the local projection applies. Values beyond column SC(i)-1 in row S(i,:) are "don't care"
	//Output:
	//	projP	#S by 3*max(SC) in format xyzxyzxyz,  where the projected points correspond to each set in S in the same order.
	typedef std::function<bool(const Eigen::PlainObjectBase<Eigen::MatrixXd>&, const Eigen::PlainObjectBase<Eigen::VectorXi>&, const Eigen::PlainObjectBase<Eigen::MatrixXi>&, Eigen::PlainObjectBase<Eigen::MatrixXd>&)> shapeup_projection_function;

  
  //This projection does nothing but render points into projP. Mostly used for "echoing" the global step
  IGL_INLINE bool shapeup_identity_projection(const Eigen::PlainObjectBase<Eigen::MatrixXd>& P, const Eigen::PlainObjectBase<Eigen::VectorXi>& SC, const Eigen::PlainObjectBase<Eigen::MatrixXi>& S,  Eigen::PlainObjectBase<Eigen::MatrixXd>& projP);
  
  //the projection assumes that the sets are vertices of polygons in cyclic order
  IGL_INLINE bool shapeup_regular_face_projection(const Eigen::PlainObjectBase<Eigen::MatrixXd>& P, const Eigen::PlainObjectBase<Eigen::VectorXi>& SC, const Eigen::PlainObjectBase<Eigen::MatrixXi>& S,  Eigen::PlainObjectBase<Eigen::MatrixXd>& projP);
  
 
}

#ifndef IGL_STATIC_LIBRARY
#include "shapeup_local_projections.cpp"
#endif

#endif
