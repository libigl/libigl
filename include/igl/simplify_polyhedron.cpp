// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2016 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "simplify_polyhedron.h"
#include "decimate.h"
#include "circulation.h"
#include "per_face_normals.h"
#include "infinite_cost_stopping_condition.h"
#include <functional>

IGL_INLINE void igl::simplify_polyhedron(
  const Eigen::MatrixXd & OV,
  const Eigen::MatrixXi & OF,
  Eigen::MatrixXd & V,
  Eigen::MatrixXi & F,
  Eigen::VectorXi & J)
{
  assert(false && "This is incorrect. In all but simple cases, this could introduce self-intersections");
  // TODO: to generalize to non-manifold and open meshes, the cost should be 0
  // if moving the vertex kept all incident faces in their original planes and
  // kept all incident boundary edges on their original lines.

  Eigen::MatrixXd N;
  // Function for computing cost of collapsing edge (0 if at least one
  // direction doesn't change pointset, inf otherwise) and placement (in lowest
  // cost direction).
  const auto & perfect= [&N](
    const int e,
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXi & E,
    const Eigen::VectorXi & EMAP,
    const Eigen::MatrixXi & EF,
    const Eigen::MatrixXi & EI,
    double & cost,
    Eigen::RowVectorXd & p)
  {
    // Function for ocmputing cost (0 or inf) of collapsing edge by placing
    // vertex at `positive` end of edge.
    const auto & perfect_directed = [&N](
      const int e,
      const bool positive,
      const Eigen::MatrixXd & V,
      const Eigen::MatrixXi & F,
      const Eigen::MatrixXi & E,
      const Eigen::VectorXi & EMAP,
      const Eigen::MatrixXi & EF,
      const Eigen::MatrixXi & EI,
      double & cost,
      Eigen::RowVectorXd & p)
    {
      const auto vi = E(e,positive);
      const auto vj = E(e,!positive);
      p = V.row(vj);
      std::vector<int> faces = igl::circulation(e,positive,F,E,EMAP,EF,EI);
      assert(faces.size() >= 2);
      const Eigen::RowVectorXd nfront = N.row(faces.front());
      const Eigen::RowVectorXd nback = N.row(faces.back());
      // loop around faces incident on vi, beginning by matching normals to
      // front, then allow one switch to back normal
      bool matching_front = true;
      cost = 0;
      for(auto f : faces)
      {
        const Eigen::RowVectorXd nf = N.row(f);
        const double epsilon = 1e-10;
        if(matching_front && (nf-nfront).norm()>epsilon)
        {
          matching_front = false;
        }
        if(!matching_front && (nf-nback).norm()>epsilon)
        {
          cost = std::numeric_limits<double>::infinity();
          break;
        }
      }
    }; 
    p.resize(3);
    double cost0, cost1;
    Eigen::RowVectorXd p0, p1;
    perfect_directed(e,false,V,F,E,EMAP,EF,EI,cost0,p0);
    perfect_directed(e,true,V,F,E,EMAP,EF,EI,cost1,p1);
    if(cost0 < cost1)
    {
      cost = cost0;
      p = p0;
    }else
    {
      cost = cost1;
      p = p1;
    }
  };
  igl::per_face_normals(OV,OF,N);
  igl::decimate(
    OV,OF,perfect,igl::infinite_cost_stopping_condition(perfect),V,F,J);
}

