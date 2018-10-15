#include <igl/colon.h>
#include <igl/harmonic.h>
#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include <algorithm>
#include <iostream>
#include "tutorial_shared_path.h"

double z_max = 1.0;
double z_dir = -0.03;
int k = 2;
bool resolve = true;
Eigen::MatrixXd V,U;
Eigen::VectorXd Z;
Eigen::MatrixXi F;
Eigen::VectorXi b;
Eigen::VectorXd bc;

bool pre_draw(igl::opengl::glfw::Viewer & viewer)
{
  using namespace Eigen;
  if(resolve)
  {
    igl::harmonic(V,F,b,bc,k,Z);
    resolve = false;
  }
  U.col(2) = z_max*Z;
  viewer.data().set_vertices(U);
  viewer.data().compute_normals();
  if(viewer.core().is_animating)
  {
    z_max += z_dir;
    z_dir *= (z_max>=1.0 || z_max<=0.0?-1.0:1.0);
  }
  return false;
}

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods)
{
  switch(key)
  {
    case ' ':
      viewer.core().is_animating = !viewer.core().is_animating;
      break;
    case '.':
      k++;
      k = (k>4?4:k);
      resolve = true;
      break;
    case ',':
      k--;
      k = (k<1?1:k);
      resolve = true;
      break;
  }
  return true;
}

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  igl::readOBJ(TUTORIAL_SHARED_PATH "/bump-domain.obj",V,F);
  U=V;
  // Find boundary vertices outside annulus
  typedef Matrix<bool,Dynamic,1> VectorXb;
  VectorXb is_outer = (V.rowwise().norm().array()-1.0)>-1e-15;
  VectorXb is_inner = (V.rowwise().norm().array()-0.15)<1e-15;
  VectorXb in_b = is_outer.array() || is_inner.array();
  igl::colon<int>(0,V.rows()-1,b);
  b.conservativeResize(stable_partition( b.data(), b.data()+b.size(),
   [&in_b](int i)->bool{return in_b(i);})-b.data());
  bc.resize(b.size(),1);
  for(int bi = 0;bi<b.size();bi++)
  {
    bc(bi) = (is_outer(b(bi))?0.0:1.0);
  }


  // Pseudo-color based on selection
  MatrixXd C(F.rows(),3);
  RowVector3d purple(80.0/255.0,64.0/255.0,255.0/255.0);
  RowVector3d gold(255.0/255.0,228.0/255.0,58.0/255.0);
  for(int f = 0;f<F.rows();f++)
  {
    if( in_b(F(f,0)) && in_b(F(f,1)) && in_b(F(f,2)))
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
  viewer.data().show_lines = false;
  viewer.data().set_colors(C);
  viewer.core().trackball_angle = Eigen::Quaternionf(0.81,-0.58,-0.03,-0.03);
  viewer.core().trackball_angle.normalize();
  viewer.callback_pre_draw = &pre_draw;
  viewer.callback_key_down = &key_down;
  viewer.core().is_animating = true;
  viewer.core().animation_max_fps = 30.;
  cout<<
    "Press [space] to toggle animation."<<endl<<
    "Press '.' to increase k."<<endl<<
    "Press ',' to decrease k."<<endl;
  viewer.launch();
}
