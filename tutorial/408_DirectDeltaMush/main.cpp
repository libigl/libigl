#include <igl/read_triangle_mesh.h>
#include <igl/readTGF.h>
#include <igl/readDMAT.h>
#include <igl/lbs_matrix.h>
#include <igl/deform_skeleton.h>
#include <igl/direct_delta_mush.h>
#include <igl/opengl/glfw/Viewer.h>
#include "tutorial_shared_path.h"
#include <Eigen/Geometry>
#include <vector>

int main(int argc, char * argv[]) 
{

  Eigen::MatrixXd V,U,C,W,T,M,Omega;
  Eigen::MatrixXi F,BE;
  igl::read_triangle_mesh(TUTORIAL_SHARED_PATH "/elephant.obj",V,F);
  igl::readTGF(           TUTORIAL_SHARED_PATH "/elephant.tgf",C,BE);
  igl::readDMAT(          TUTORIAL_SHARED_PATH "/elephant-weights.dmat",W);
  igl::readDMAT(          TUTORIAL_SHARED_PATH "/elephant-anim.dmat",T);
  // convert weight to piecewise-rigid weights to stress test DDM
  for (int i = 0; i < W.rows(); ++i)
  {
    int maxj;
    W.row(i).maxCoeff(&maxj);
    for (int j = 0; j < W.cols(); j++)
    {
      W(i, j) = double(maxj == j);
    }
  }

  igl::lbs_matrix(V,W,M);

  int p = 20;
  float lambda = 3; // 0 < lambda
  float kappa = 1; // 0 < kappa < lambda
  float alpha = 0.8; // 0 <= alpha < 1
  igl::direct_delta_mush_precomputation(V, F,W, p, lambda, kappa, alpha, Omega);

  igl::opengl::glfw::Viewer viewer;
  int frame = 0;
  const int pr_id = viewer.selected_data_index;
  viewer.append_mesh();
  const int ddm_id = viewer.selected_data_index;
  Eigen::RowVector3d offset(1.1*(V.col(0).maxCoeff()-V.col(0).minCoeff()),0,0);

  viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer &) -> bool
  {
    if(viewer.core().is_animating)
    {
      const Eigen::Map<Eigen::MatrixXd> Tf(
        T.data()+frame*T.rows(),4*W.cols(),3);
      U = (M*Tf).rowwise()-offset;

      Eigen::MatrixXd Cf;
      Eigen::MatrixXi BEf;
      igl::deform_skeleton(C,BE,Tf,Cf,BEf);
      viewer.data(pr_id).set_edges(Cf,BEf, Eigen::RowVector3d(1,1,1));

      viewer.data(pr_id).set_vertices(U);
      viewer.data(pr_id).compute_normals();

      {
        std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d>> 
          T_list(BE.rows());
        for (int e = 0; e < BE.rows(); e++)
        {
          T_list[e] = Eigen::Affine3d::Identity();
          T_list[e].matrix().block(0,0,3,4) = Tf.block(e*4,0,4,3).transpose();
        }
        igl::direct_delta_mush(V, T_list, Omega, U);
      }
      U.rowwise() += offset;
      viewer.data(ddm_id).set_vertices(U);
      viewer.data(ddm_id).compute_normals();

      frame++;
      if(frame == T.cols())
      {
        frame = 0;
        viewer.core().is_animating = false;
      }
    }
    return false;
  };
  viewer.callback_key_pressed = [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod)
  {
    switch(key)
    {
    case ' ':
      viewer.core().is_animating = !viewer.core().is_animating;
      return true;
    }
    return false;
  };



  for(auto & id : {pr_id,ddm_id})
  {
    if(id == pr_id)
    {
      viewer.data(id).set_mesh( (V.rowwise()-offset*1.0).eval() ,F);
      viewer.data(id).set_colors(Eigen::RowVector3d(214./255.,170./255.,148./255.));
      viewer.data(id).set_edges(C,BE, Eigen::RowVector3d(1,1,1));
    }else if(id == ddm_id){
      viewer.data(id).set_mesh( (V.rowwise()+offset*1.0).eval() ,F);
      viewer.data(id).set_colors(Eigen::RowVector3d(132./255.,214./255.,105./255.));
    }
    viewer.data(id).show_lines = false;
    viewer.data(id).set_face_based(true);
    viewer.data(id).show_overlay_depth = false;
  }
  viewer.core().is_animating = false;
  viewer.core().animation_max_fps = 24.;
  //viewer.data().set_colors(V,F);


  viewer.launch_init();
  viewer.core().align_camera_center(V);

  viewer.launch_rendering(true);
  viewer.launch_shut();
}
