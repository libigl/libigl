// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2019 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/iterative_closest_point.h>
#include <igl/random_dir.h>
#include <igl/PI.h>
#include <igl/AABB.h>
#include <igl/per_face_normals.h>
#include <igl/avg_edge_length.h>
#include <igl/per_vertex_normals.h>
#include <igl/barycenter.h>
#include <igl/slice.h>
#include <igl/slice_mask.h>
#include <iostream>

int main(int argc, char * argv[])
{
  Eigen::MatrixXd OVX,VX,VY;
  Eigen::MatrixXi FX,FY;
  igl::read_triangle_mesh( argc>1?argv[1]: TUTORIAL_SHARED_PATH "/decimated-max.obj",VY,FY);
  const double bbd = (VY.colwise().maxCoeff()-VY.colwise().minCoeff()).norm();
  FX = FY;
  {
    // sprinkle a noise so that we can see z-fighting when the match is perfect.
    const double h = igl::avg_edge_length(VY,FY);
    OVX = VY + 1e-2*h*Eigen::MatrixXd::Random(VY.rows(),VY.cols());
  }
  
  VX = OVX;

  igl::AABB<Eigen::MatrixXd,3> Ytree;
  Ytree.init(VY,FY);
  Eigen::MatrixXd NY;
  igl::per_face_normals(VY,FY,NY);

  igl::opengl::glfw::Viewer v;
  std::cout<<R"(
  [space]  conduct a single iterative closest point iteration
  R,r      reset to a random orientation and offset
)";

  const auto apply_random_rotation = [&]()
  {
    const Eigen::Matrix3d R = Eigen::AngleAxisd(
      2.*igl::PI*(double)rand()/RAND_MAX*0.3, igl::random_dir()).matrix();
    const Eigen::RowVector3d cen = 
      0.5*(VY.colwise().maxCoeff()+VY.colwise().minCoeff());
    VX = ((OVX*R).rowwise()+(cen-cen*R)).eval();
  };
  const auto single_iteration = [&]()
  {
    ////////////////////////////////////////////////////////////////////////
    // Perform single iteration of ICP method
    ////////////////////////////////////////////////////////////////////////
    Eigen::Matrix3d R;
    Eigen::RowVector3d t;
    igl::iterative_closest_point(VX,FX,VY,FY,Ytree,NY,1000,1,R,t);
    VX = ((VX*R).rowwise()+t).eval();
    v.data().set_mesh(VX,FX);
    v.data().compute_normals();
  };
  v.callback_pre_draw = [&](igl::opengl::glfw::Viewer &)->bool
  {
    if(v.core().is_animating)
    {
      single_iteration();
    }
    return false;
  };
  v.callback_key_pressed = 
    [&](igl::opengl::glfw::Viewer &,unsigned char key,int)->bool
  {
    switch(key)
    {
      case ' ':
      {
        v.core().is_animating = false;
        single_iteration();
        return true;
      }
      case 'R':
      case 'r':
        // Random rigid transformation
        apply_random_rotation();
        v.data().set_mesh(VX,FX);
        v.data().compute_normals();
        return true;
        break;
    }
    return false;
  };

  v.data().set_mesh(VY,FY);
  v.data().set_colors(Eigen::RowVector3d(1,1,1));
  v.data().show_lines = false;
  v.append_mesh();
  v.data().set_mesh(VX,FX);
  v.data().show_lines = false;
  v.launch();
}
