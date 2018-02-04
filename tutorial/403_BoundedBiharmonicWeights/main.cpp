// Because of Mosek complications, we don't use static library if Mosek is used.
#ifdef LIBIGL_WITH_MOSEK
#ifdef IGL_STATIC_LIBRARY
#undef IGL_STATIC_LIBRARY
#endif
#endif

#include <igl/boundary_conditions.h>
#include <igl/colon.h>
#include <igl/column_to_quats.h>
#include <igl/directed_edge_parents.h>
#include <igl/forward_kinematics.h>
#include <igl/jet.h>
#include <igl/lbs_matrix.h>
#include <igl/deform_skeleton.h>
#include <igl/normalize_row_sums.h>
#include <igl/readDMAT.h>
#include <igl/readMESH.h>
#include <igl/readTGF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/bbw.h>
//#include <igl/embree/bone_heat.h>

#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <vector>
#include <algorithm>
#include <iostream>

#include "tutorial_shared_path.h"

typedef
  std::vector<Eigen::Quaterniond,Eigen::aligned_allocator<Eigen::Quaterniond> >
  RotationList;

const Eigen::RowVector3d sea_green(70./255.,252./255.,167./255.);
int selected = 0;
Eigen::MatrixXd V,W,U,C,M;
Eigen::MatrixXi T,F,BE;
Eigen::VectorXi P;
RotationList pose;
double anim_t = 1.0;
double anim_t_dir = -0.03;

bool pre_draw(igl::opengl::glfw::Viewer & viewer)
{
  using namespace Eigen;
  using namespace std;
  if(viewer.core.is_animating)
  {
    // Interpolate pose and identity
    RotationList anim_pose(pose.size());
    for(int e = 0;e<pose.size();e++)
    {
      anim_pose[e] = pose[e].slerp(anim_t,Quaterniond::Identity());
    }
    // Propagate relative rotations via FK to retrieve absolute transformations
    RotationList vQ;
    vector<Vector3d> vT;
    igl::forward_kinematics(C,BE,P,anim_pose,vQ,vT);
    const int dim = C.cols();
    MatrixXd T(BE.rows()*(dim+1),dim);
    for(int e = 0;e<BE.rows();e++)
    {
      Affine3d a = Affine3d::Identity();
      a.translate(vT[e]);
      a.rotate(vQ[e]);
      T.block(e*(dim+1),0,dim+1,dim) =
        a.matrix().transpose().block(0,0,dim+1,dim);
    }
    // Compute deformation via LBS as matrix multiplication
    U = M*T;

    // Also deform skeleton edges
    MatrixXd CT;
    MatrixXi BET;
    igl::deform_skeleton(C,BE,T,CT,BET);

    viewer.data().set_vertices(U);
    viewer.data().set_edges(CT,BET,sea_green);
    viewer.data().compute_normals();
    anim_t += anim_t_dir;
    anim_t_dir *= (anim_t>=1.0 || anim_t<=0.0?-1.0:1.0);
  }
  return false;
}

void set_color(igl::opengl::glfw::Viewer &viewer)
{
  Eigen::MatrixXd C;
  igl::jet(W.col(selected).eval(),true,C);
  viewer.data().set_colors(C);
}

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods)
{
  switch(key)
  {
    case ' ':
      viewer.core.is_animating = !viewer.core.is_animating;
      break;
    case '.':
      selected++;
      selected = std::min(std::max(selected,0),(int)W.cols()-1);
      set_color(viewer);
      break;
    case ',':
      selected--;
      selected = std::min(std::max(selected,0),(int)W.cols()-1);
      set_color(viewer);
      break;
  }
  return true;
}

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  igl::readMESH(TUTORIAL_SHARED_PATH "/hand.mesh",V,T,F);
  U=V;
  igl::readTGF(TUTORIAL_SHARED_PATH "/hand.tgf",C,BE);
  // retrieve parents for forward kinematics
  igl::directed_edge_parents(BE,P);

  // Read pose as matrix of quaternions per row
  MatrixXd Q;
  igl::readDMAT(TUTORIAL_SHARED_PATH "/hand-pose.dmat",Q);
  igl::column_to_quats(Q,pose);
  assert(pose.size() == BE.rows());

  // List of boundary indices (aka fixed value indices into VV)
  VectorXi b;
  // List of boundary conditions of each weight function
  MatrixXd bc;
  igl::boundary_conditions(V,T,C,VectorXi(),BE,MatrixXi(),b,bc);

  // compute BBW weights matrix
  igl::BBWData bbw_data;
  // only a few iterations for sake of demo
  bbw_data.active_set_params.max_iter = 8;
  bbw_data.verbosity = 2;
  if(!igl::bbw(V,T,b,bc,bbw_data,W))
  {
    return EXIT_FAILURE;
  }

  //MatrixXd Vsurf = V.topLeftCorner(F.maxCoeff()+1,V.cols());
  //MatrixXd Wsurf;
  //if(!igl::bone_heat(Vsurf,F,C,VectorXi(),BE,MatrixXi(),Wsurf))
  //{
  //  return false;
  //}
  //W.setConstant(V.rows(),Wsurf.cols(),1);
  //W.topLeftCorner(Wsurf.rows(),Wsurf.cols()) = Wsurf = Wsurf = Wsurf = Wsurf;

  // Normalize weights to sum to one
  igl::normalize_row_sums(W,W);
  // precompute linear blend skinning matrix
  igl::lbs_matrix(V,W,M);

  // Plot the mesh with pseudocolors
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(U, F);
  set_color(viewer);
  viewer.data().set_edges(C,BE,sea_green);
  viewer.data().show_lines = false;
  viewer.data().show_overlay_depth = false;
  viewer.data().line_width = 1;
  viewer.callback_pre_draw = &pre_draw;
  viewer.callback_key_down = &key_down;
  viewer.core.is_animating = false;
  viewer.core.animation_max_fps = 30.;
  cout<<
    "Press '.' to show next weight function."<<endl<<
    "Press ',' to show previous weight function."<<endl<<
    "Press [space] to toggle animation."<<endl;
  viewer.launch();
  return EXIT_SUCCESS;
}
