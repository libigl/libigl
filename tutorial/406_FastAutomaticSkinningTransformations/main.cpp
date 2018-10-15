#include <igl/colon.h>
#include <igl/directed_edge_orientations.h>
#include <igl/directed_edge_parents.h>
#include <igl/forward_kinematics.h>
#include <igl/PI.h>
#include <igl/partition.h>
#include <igl/mat_max.h>
#include <igl/lbs_matrix.h>
#include <igl/slice.h>
#include <igl/deform_skeleton.h>
#include <igl/dqs.h>
#include <igl/lbs_matrix.h>
#include <igl/columnize.h>
#include <igl/readDMAT.h>
#include <igl/readOBJ.h>
#include <igl/arap.h>
#include <igl/arap_dof.h>
#include <igl/opengl/glfw/Viewer.h>

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
Eigen::MatrixXd V,U,M;
Eigen::MatrixXi F;
Eigen::VectorXi S,b;
Eigen::MatrixXd L;
Eigen::RowVector3d mid;
double anim_t = 0.0;
double anim_t_dir = 0.03;
double bbd = 1.0;
bool resolve = true;
igl::ARAPData arap_data,arap_grouped_data;
igl::ArapDOFData<Eigen::MatrixXd,double> arap_dof_data;
Eigen::SparseMatrix<double> Aeq;

enum ModeType
{
  MODE_TYPE_ARAP = 0,
  MODE_TYPE_ARAP_GROUPED = 1,
  MODE_TYPE_ARAP_DOF = 2,
  NUM_MODE_TYPES = 4
} mode = MODE_TYPE_ARAP;

bool pre_draw(igl::opengl::glfw::Viewer & viewer)
{
  using namespace Eigen;
  using namespace std;
  if(resolve)
  {
    MatrixXd bc(b.size(),V.cols());
    VectorXd Beq(3*b.size());
    for(int i = 0;i<b.size();i++)
    {
      bc.row(i) = V.row(b(i));
      switch(i%4)
      {
        case 2:
          bc(i,0) += 0.15*bbd*sin(0.5*anim_t);
          bc(i,1) += 0.15*bbd*(1.-cos(0.5*anim_t));
          break;
        case 1:
          bc(i,1) += 0.10*bbd*sin(1.*anim_t*(i+1));
          bc(i,2) += 0.10*bbd*(1.-cos(1.*anim_t*(i+1)));
          break;
        case 0:
          bc(i,0) += 0.20*bbd*sin(2.*anim_t*(i+1));
          break;
      }
      Beq(3*i+0) = bc(i,0);
      Beq(3*i+1) = bc(i,1);
      Beq(3*i+2) = bc(i,2);
    }
    switch(mode)
    {
      default:
        assert("unknown mode");
      case MODE_TYPE_ARAP:
        igl::arap_solve(bc,arap_data,U);
        break;
      case MODE_TYPE_ARAP_GROUPED:
        igl::arap_solve(bc,arap_grouped_data,U);
        break;
      case MODE_TYPE_ARAP_DOF:
      {
        VectorXd L0 = L;
        arap_dof_update(arap_dof_data,Beq,L0,30,0,L);
        const auto & Ucol = M*L;
        U.col(0) = Ucol.block(0*U.rows(),0,U.rows(),1);
        U.col(1) = Ucol.block(1*U.rows(),0,U.rows(),1);
        U.col(2) = Ucol.block(2*U.rows(),0,U.rows(),1);
        break;
      }
    }
    viewer.data().set_vertices(U);
    viewer.data().set_points(bc,sea_green);
    viewer.data().compute_normals();
    if(viewer.core().is_animating)
    {
      anim_t += anim_t_dir;
    }else
    {
      resolve = false;
    }
  }
  return false;
}

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods)
{
  switch(key)
  {
    case '0':
      anim_t = 0;
      resolve = true;
      return true;
    case '.':
      mode = (ModeType)(((int)mode+1)%((int)NUM_MODE_TYPES-1));
      resolve = true;
      return true;
    case ',':
      mode = (ModeType)(((int)mode-1)%((int)NUM_MODE_TYPES-1));
      resolve = true;
      return true;
    case ' ':
      viewer.core().is_animating = !viewer.core().is_animating;
      if(viewer.core().is_animating)
      {
        resolve = true;
      }
      return true;
  }
  return false;
}

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  igl::readOBJ(TUTORIAL_SHARED_PATH "/armadillo.obj",V,F);
  U=V;
  MatrixXd W;
  igl::readDMAT(TUTORIAL_SHARED_PATH "/armadillo-weights.dmat",W);
  igl::lbs_matrix_column(V,W,M);

  // Cluster according to weights
  VectorXi G;
  {
    VectorXi S;
    VectorXd D;
    igl::partition(W,50,G,S,D);
  }

  // vertices corresponding to handles (those with maximum weight)
  {
    VectorXd maxW;
    igl::mat_max(W,1,maxW,b);
  }

  // Precomputation for FAST
  cout<<"Initializing Fast Automatic Skinning Transformations..."<<endl;
  // number of weights
  const int m = W.cols();
  Aeq.resize(m*3,m*3*(3+1));
  vector<Triplet<double> > ijv;
  for(int i = 0;i<m;i++)
  {
    RowVector4d homo;
    homo << V.row(b(i)),1.;
    for(int d = 0;d<3;d++)
    {
      for(int c = 0;c<(3+1);c++)
      {
        ijv.push_back(Triplet<double>(3*i + d,i + c*m*3 + d*m, homo(c)));
      }
    }
  }
  Aeq.setFromTriplets(ijv.begin(),ijv.end());
  igl::arap_dof_precomputation(V,F,M,G,arap_dof_data);
  igl::arap_dof_recomputation(VectorXi(),Aeq,arap_dof_data);
  // Initialize
  MatrixXd Istack = MatrixXd::Identity(3,3+1).replicate(1,m);
  igl::columnize(Istack,m,2,L);

  // Precomputation for ARAP
  cout<<"Initializing ARAP..."<<endl;
  arap_data.max_iter = 1;
  igl::arap_precomputation(V,F,V.cols(),b,arap_data);
  // Grouped arap
  cout<<"Initializing ARAP with grouped edge-sets..."<<endl;
  arap_grouped_data.max_iter = 2;
  arap_grouped_data.G = G;
  igl::arap_precomputation(V,F,V.cols(),b,arap_grouped_data);


  // bounding box diagonal
  bbd = (V.colwise().maxCoeff()- V.colwise().minCoeff()).norm();

  // Plot the mesh with pseudocolors
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(U, F);
  viewer.data().add_points(igl::slice(V,b,1),sea_green);
  viewer.data().show_lines = false;
  viewer.callback_pre_draw = &pre_draw;
  viewer.callback_key_down = &key_down;
  viewer.core().is_animating = false;
  viewer.core().animation_max_fps = 30.;
  cout<<
    "Press [space] to toggle animation."<<endl<<
    "Press '0' to reset pose."<<endl<<
    "Press '.' to switch to next deformation method."<<endl<<
    "Press ',' to switch to previous deformation method."<<endl;
  viewer.launch();
}
