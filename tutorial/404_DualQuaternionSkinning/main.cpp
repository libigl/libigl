// Don't use static library for this example because of Mosek complications
//#define IGL_NO_MOSEK
#ifdef IGL_NO_MOSEK
#undef IGL_STATIC_LIBRARY
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
#include <igl/viewer/Viewer.h>
#include <igl/bbw/bbw.h>

#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <vector>
#include <algorithm>
#include <iostream>

typedef 
  std::vector<Eigen::Quaterniond,Eigen::aligned_allocator<Eigen::Quaterniond> >
  RotationList;

const Eigen::RowVector3d sea_green(70./255.,252./255.,167./255.);
Eigen::MatrixXd V,W;
Eigen::MatrixXi F,BE;
Eigen::VectorXi P;
RotationList pose;
double anim_t = 1.0;
double anim_t_dir = -0.03;

bool pre_draw(igl::Viewer & viewer)
{
  using namespace Eigen;
  using namespace std;
  if(viewer.options.is_animating)
  {
    // Interpolate pose and identity
    RotationList anim_pose(pose.size());
    for(int e = 0;e<pose.size();e++)
    {
      anim_pose[e] = pose[e].slerp(anim_t,Quaterniond::Identity());
    }
    // Propogate relative rotations via FK to retrieve absolute transformations
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
    
    viewer.set_vertices(U);
    viewer.set_edges(CT,BET,sea_green);
    viewer.compute_normals();
    anim_t += anim_t_dir;
    anim_t_dir *= (anim_t>=1.0 || anim_t<=0.0?-1.0:1.0);
  }
  return false;
}

bool key_down(igl::Viewer &viewer, unsigned char key, int mods)
{
  switch(key)
  {
    case ' ':
      viewer.options.is_animating = !viewer.options.is_animating;
      break;
  }
}

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  igl::readOBJ("../shared/arm.obj",V,F);
  U=V;
  igl::readTGF("../shared/arm.tgf",C,BE);
  // retrieve parents for forward kinematics
  directed_edge_parents(BE,P);
  igl::readDMAT("../shared/arm-weights.dmat",W);
  pose_0 identity
  pose_1 twist
  pose_2 bend

  // Plot the mesh with pseudocolors
  igl::Viewer viewer;
  viewer.set_mesh(U, F);
  viewer.set_edges(C,BE,sea_green);
  viewer.options.show_lines = false;
  viewer.options.show_overlay_depth = false;
  viewer.options.line_width = 1;
  viewer.options.trackball_angle.normalize();
  viewer.callback_pre_draw = &pre_draw;
  viewer.callback_key_down = &key_down;
  viewer.options.is_animating = false;
  viewer.options.animation_max_fps = 30.;
  viewer.launch();
}
