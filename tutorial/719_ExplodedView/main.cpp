// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2020 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.

#include <igl/opengl/glfw/Viewer.h>
#include <igl/readMESH.h>
#include <igl/barycenter.h>
#include <igl/volume.h>
#include <igl/exploded_view.h>
#include <igl/slice.h>
#include <igl/colormap.h>


int main(int argc, char *argv[])
{

  Eigen::MatrixXd V;
  Eigen::MatrixXi T,F;
  igl::readMESH(argc>1?argv[1]: TUTORIAL_SHARED_PATH "/octopus-low.mesh",V,T,F);
  // Some per-tet data
  Eigen::VectorXd D;
  {
    Eigen::MatrixXd BC;
    igl::barycenter(V,T,BC);
    Eigen::VectorXd vol;
    igl::volume(V,T,vol);
    const Eigen::RowVectorXd c = vol.transpose()*BC/vol.array().sum();
    D = (BC.rowwise()-c).rowwise().norm();
  }

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;

  double t = 1;
  double s = 1;
  const auto update = [&]()
  {
    Eigen::MatrixXd EV;
    Eigen::MatrixXi EF;
    Eigen::VectorXi I,J;
    igl::exploded_view(V,T,s,t,EV,EF,I,J);
    Eigen::VectorXd DJ;
    igl::slice(D,J,1,DJ);
    static bool first = true;
    if(first)
    {
      viewer.data().clear();
      viewer.data().set_mesh(EV,EF);
      first = false;
    }else
    {
      viewer.data().set_vertices(EV);
    }
    viewer.data().set_face_based(true);
    Eigen::MatrixXd C;
    igl::colormap(igl::COLOR_MAP_TYPE_VIRIDIS,DJ,true,C);
    viewer.data().set_colors(C);
  };
  int mod = 0;
  float prev_y;
  viewer.callback_mouse_move = [&](igl::opengl::glfw::Viewer &, int x, int y)->bool
  {
    if((mod & IGL_MOD_SHIFT)||(mod & IGL_MOD_ALT))
    {
      if(mod & IGL_MOD_SHIFT)
      {
        t = std::min(std::max(t+0.001*(y-prev_y),1.),2.);
      }
      if(mod & IGL_MOD_ALT)
      {
        s = std::min(std::max(s+0.001*(y-prev_y),0.),1.);
      }
      prev_y = y;
      update();
      return true;
    }
    return false;
  };
  // get modifier
  viewer.callback_key_down = 
    [&](igl::opengl::glfw::Viewer &, unsigned char key, int _mod)->bool
  {
    prev_y = viewer.current_mouse_y;
    mod = _mod;
    return false;
  };
  viewer.callback_key_up = 
    [&](igl::opengl::glfw::Viewer &, unsigned char key, int _mod)->bool
  {
    mod = _mod;
    return false;
  };
  update();
  viewer.launch();
}
