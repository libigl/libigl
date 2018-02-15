#include <igl/colon.h>
#include <igl/directed_edge_orientations.h>
#include <igl/directed_edge_parents.h>
#include <igl/forward_kinematics.h>
#include <igl/PI.h>
#include <igl/lbs_matrix.h>
#include <igl/deform_skeleton.h>
#include <igl/dqs.h>
#include <igl/readDMAT.h>
#include <igl/readOFF.h>
#include <igl/arap.h>
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
Eigen::MatrixXd V,U;
Eigen::MatrixXi F;
Eigen::VectorXi S,b;
Eigen::RowVector3d mid;
double anim_t = 0.0;
double anim_t_dir = 0.03;
igl::ARAPData arap_data;

bool pre_draw(igl::opengl::glfw::Viewer & viewer)
{
  using namespace Eigen;
  using namespace std;
    MatrixXd bc(b.size(),V.cols());
    for(int i = 0;i<b.size();i++)
    {
      bc.row(i) = V.row(b(i));
      switch(S(b(i)))
      {
        case 0:
        {
          const double r = mid(0)*0.25;
          bc(i,0) += r*sin(0.5*anim_t*2.*igl::PI);
          bc(i,1) -= r+r*cos(igl::PI+0.5*anim_t*2.*igl::PI);
          break;
        }
        case 1:
        {
          const double r = mid(1)*0.15;
          bc(i,1) += r+r*cos(igl::PI+0.15*anim_t*2.*igl::PI);
          bc(i,2) -= r*sin(0.15*anim_t*2.*igl::PI);
          break;
        }
        case 2:
        {
          const double r = mid(1)*0.15;
          bc(i,2) += r+r*cos(igl::PI+0.35*anim_t*2.*igl::PI);
          bc(i,0) += r*sin(0.35*anim_t*2.*igl::PI);
          break;
        }
        default:
          break;
      }
    }
    igl::arap_solve(bc,arap_data,U);
    viewer.data().set_vertices(U);
    viewer.data().compute_normals();
  if(viewer.core.is_animating)
  {
    anim_t += anim_t_dir;
  }
  return false;
}

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods)
{
  switch(key)
  {
    case ' ':
      viewer.core.is_animating = !viewer.core.is_animating;
      return true;
  }
  return false;
}

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  igl::readOFF(TUTORIAL_SHARED_PATH "/decimated-knight.off",V,F);
  U=V;
  igl::readDMAT(TUTORIAL_SHARED_PATH "/decimated-knight-selection.dmat",S);

  // vertices in selection
  igl::colon<int>(0,V.rows()-1,b);
  b.conservativeResize(stable_partition( b.data(), b.data()+b.size(), 
   [](int i)->bool{return S(i)>=0;})-b.data());
  // Centroid
  mid = 0.5*(V.colwise().maxCoeff() + V.colwise().minCoeff());
  // Precomputation
  arap_data.max_iter = 100;
  igl::arap_precomputation(V,F,V.cols(),b,arap_data);

  // Set color based on selection
  MatrixXd C(F.rows(),3);
  RowVector3d purple(80.0/255.0,64.0/255.0,255.0/255.0);
  RowVector3d gold(255.0/255.0,228.0/255.0,58.0/255.0);
  for(int f = 0;f<F.rows();f++)
  {
    if( S(F(f,0))>=0 && S(F(f,1))>=0 && S(F(f,2))>=0)
    {
      C.row(f) = purple;
    }else
    {
      C.row(f) = gold;
    }
  }

  // Plot the mesh with pseudocolors
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(U, F);
  viewer.data().set_colors(C);
  viewer.callback_pre_draw = &pre_draw;
  viewer.callback_key_down = &key_down;
  viewer.core.is_animating = false;
  viewer.core.animation_max_fps = 30.;
  cout<<
    "Press [space] to toggle animation"<<endl;
  viewer.launch();
}
