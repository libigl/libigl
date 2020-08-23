// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2020 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.

#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/blue_noise.h>
#include <igl/per_vertex_normals.h>
#include <igl/barycentric_interpolation.h>
#include <igl/blue_noise.h>
#include <igl/doublearea.h>
#include <igl/random_points_on_mesh.h>
#include <igl/PI.h>


int main(int argc, char *argv[])
{

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(
    argc>1?argv[1]: TUTORIAL_SHARED_PATH "/elephant.obj",V,F);
  Eigen::MatrixXd N;
  igl::per_vertex_normals(V,F,N);
  const double bbd = (V.colwise().maxCoeff()- V.colwise().minCoeff()).norm();

  Eigen::MatrixXd P_blue;
  Eigen::MatrixXd N_blue;
  Eigen::VectorXd A_blue;
  Eigen::MatrixXd C_blue;
  {
    const int n_desired = argc>2?atoi(argv[2]):50000;
    // Heuristic to  determine radius from desired number 
    const double r = [&V,&F](const int n)
    {
      Eigen::VectorXd A;
      igl::doublearea(V,F,A);
      return sqrt(((A.sum()*0.5/(n*0.6162910373))/igl::PI));
    }(n_desired);
    printf("blue noise radius: %g\n",r);
    Eigen::MatrixXd B;
    Eigen::VectorXi I;
    igl::blue_noise(V,F,r,B,I,P_blue);
    igl::barycentric_interpolation(N,F,B,I,N_blue);
    N_blue.rowwise().normalize();
  }
  Eigen::MatrixXd P_white;
  Eigen::MatrixXd N_white;
  Eigen::VectorXd A_white;
  Eigen::MatrixXd C_white;
  {
    Eigen::MatrixXd B;
    Eigen::VectorXi I;
    igl::random_points_on_mesh(P_blue.rows(),V,F,B,I,P_white);
    igl::barycentric_interpolation(N,F,B,I,N_white);
    N_white.rowwise().normalize();
  }

  // Color
  Eigen::RowVector3d C(1,1,1);

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V,F);
  viewer.data().show_lines = false;
  viewer.data().show_faces = false;
  viewer.data().point_size = 4;

  bool use_blue = true;
  const auto update = [&]()
  {
    const Eigen::RowVector3d orange(1.0,0.7,0.2);
    if(use_blue)
    {
      viewer.data().set_points(P_blue,Eigen::RowVector3d(0.3,0.4,1.0));
      viewer.data().set_edges_from_vector_field(P_blue,bbd*0.01*N_blue,orange);
    }else
    {
      viewer.data().set_points(P_white,Eigen::RowVector3d(1,1,1));
      viewer.data().set_edges_from_vector_field(P_white,bbd*0.01*N_white,orange);
    }
  };
  update();

  viewer.callback_key_pressed = 
    [&](igl::opengl::glfw::Viewer &,unsigned char key,int)->bool
  {
    switch(key)
    {
      case ' ':
      {
        use_blue = !use_blue;
        update();
        return true;
      }
    }
    return false;
  };
  std::cout<<R"(
[space]  toggle between blue and white noise samplings
)";

  viewer.launch();
}
