#include <igl/direct_delta_mush.h>
#include <igl/directed_edge_orientations.h>
#include <igl/directed_edge_parents.h>
#include <igl/forward_kinematics.h>
#include <igl/PI.h>
#include <igl/lbs_matrix.h>
#include <igl/deform_skeleton.h>
#include <igl/readDMAT.h>
#include <igl/readOBJ.h>
#include <igl/readTGF.h>
#include <igl/opengl/glfw/Viewer.h>

#include <vector>
#include <iostream>

#include "tutorial_shared_path.h"

#include <igl/writeDMAT.h>

const bool USE_SAVED_OMEGA = true;

typedef std::vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond>>
  RotationList;

typedef std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d>>
  TransformationList;

const Eigen::RowVector3d sea_green(70. / 255., 252. / 255., 167. / 255.);
Eigen::MatrixXd V, W, C, U, M, Omega;
Eigen::MatrixXi F, BE;
Eigen::VectorXi P;
Eigen::SparseMatrix<double> W_sparse;
std::vector<RotationList> poses;
double anim_t = 0.0;
double anim_t_dir = 0.015;
bool use_ddm = false;
bool recompute = true;
int p = 3;
float lambda = 0.75;
float kappa = 0.5;
float alpha = 0.5;

bool pre_draw(igl::opengl::glfw::Viewer & viewer)
{
  using namespace Eigen;
  using namespace std;
  if (recompute)
  {
    // Find pose interval
    const int begin = (int) floor(anim_t) % poses.size();
    const int end = (int) (floor(anim_t) + 1) % poses.size();
    const double t = anim_t - floor(anim_t);

    // Interpolate pose and identity
    RotationList anim_pose(poses[begin].size());
    for (int e = 0; e < poses[begin].size(); e++)
    {
      anim_pose[e] = poses[begin][e].slerp(t, poses[end][e]);
    }
    // Propagate relative rotations via FK to retrieve absolute transformations
    RotationList vQ;
    vector<Vector3d> vT;
    igl::forward_kinematics(C, BE, P, anim_pose, vQ, vT);
    const int dim = C.cols();
    MatrixXd T(BE.rows() * (dim + 1), dim);
    TransformationList T_list;
    for (int e = 0; e < BE.rows(); e++)
    {
      Affine3d a = Affine3d::Identity();
      a.translate(vT[e]);
      a.rotate(vQ[e]);
      T.block(e * (dim + 1), 0, dim + 1, dim) =
        a.matrix().transpose().block(0, 0, dim + 1, dim);
      T_list.push_back(a);
    }
    // Compute deformation via LBS as matrix multiplication
    if (use_ddm)
    {
      // igl::direct_delta_mush_pose_evaluation(T_list, Omega, U);
      igl::direct_delta_mush(V, F, C, BE, T_list, Omega, U);
    } else
    {
      U = M * T;
    }

    // Also deform skeleton edges
    MatrixXd CT;
    MatrixXi BET;
    igl::deform_skeleton(C, BE, T, CT, BET);

    viewer.data().set_vertices(U);
    viewer.data().set_edges(CT, BET, sea_green);
    viewer.data().compute_normals();
    if (viewer.core().is_animating)
    {
      anim_t += anim_t_dir;
    } else
    {
      recompute = false;
    }
  }
  return false;
}

bool key_down(igl::opengl::glfw::Viewer & viewer, unsigned char key, int mods)
{
  recompute = true;
  switch (key)
  {
    case 'D':
    case 'd':
      use_ddm = !use_ddm;
      return true;
    case ' ':
      viewer.core().is_animating = !viewer.core().is_animating;
      return true;
    default:
      return false;
  }
}

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  igl::readOBJ(TUTORIAL_SHARED_PATH "/arm.obj", V, F);
  U = V;
  igl::readTGF(TUTORIAL_SHARED_PATH "/arm.tgf", C, BE);
  // retrieve parents for forward kinematics
  igl::directed_edge_parents(BE, P);
  RotationList rest_pose;
  igl::directed_edge_orientations(C, BE, rest_pose);
  poses.resize(2, RotationList(4, Quaterniond::Identity()));
  const Quaterniond bend(AngleAxisd(-igl::PI * 0.7, Vector3d(0, 0, 1)));
  poses[1][2] = rest_pose[2] * bend * rest_pose[2].conjugate();

  igl::readDMAT(TUTORIAL_SHARED_PATH "/arm-weights.dmat", W);
  W_sparse = W.sparseView();
  igl::lbs_matrix(V, W, M);

  // Precomputation for Direct Delta Mush
  if (USE_SAVED_OMEGA) {
    igl::readDMAT(TUTORIAL_SHARED_PATH "/arm-weights-ddm-omega.dmat", Omega);
  }
  else {
    igl::direct_delta_mush_precomputation(V, F, C, BE, W_sparse, p, lambda, kappa, alpha, Omega);
    igl::writeDMAT(TUTORIAL_SHARED_PATH "/arm-weights-ddm-omega.dmat", Omega);
  }

  // Plot the mesh with pseudocolors
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(U, F);
  viewer.data().set_edges(C, BE, sea_green);
  viewer.data().show_lines = false;
  viewer.data().show_overlay_depth = false;
  viewer.data().line_width = 1;
  viewer.core().trackball_angle.normalize();
  viewer.callback_pre_draw = &pre_draw;
  viewer.callback_key_down = &key_down;
  viewer.core().is_animating = false;
  viewer.core().camera_zoom = 2.5;
  viewer.core().animation_max_fps = 30.;
  cout << "Press [d] to toggle between LBS and DDM" << endl
       << "Press [space] to toggle animation" << endl;
  viewer.launch();
}
