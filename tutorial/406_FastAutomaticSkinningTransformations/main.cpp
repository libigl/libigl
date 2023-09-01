#include "precomputation.h"

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

  precomputation(V,F,W,M,b,L,arap_data,arap_grouped_data,arap_dof_data);

  // bounding box diagonal
  bbd = (V.colwise().maxCoeff()- V.colwise().minCoeff()).norm();

  // Plot the mesh with pseudocolors
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(U, F);
  viewer.data().add_points(V(b,Eigen::all),sea_green);
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
