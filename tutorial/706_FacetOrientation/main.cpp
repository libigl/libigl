#include <igl/read_triangle_mesh.h>
#include <igl/randperm.h>
#include <igl/orientable_patches.h>
#include <igl/slice.h>
#include <igl/hsv_to_rgb.h>
#include <igl/embree/reorient_facets_raycast.h>
#include <igl/opengl/glfw/Viewer.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

igl::opengl::glfw::Viewer viewer;
Eigen::MatrixXd V;
std::vector<Eigen::VectorXi> C(2);
std::vector<Eigen::MatrixXd> RGBcolors(2);
Eigen::MatrixXi F;
std::vector<Eigen::MatrixXi> FF(2);
bool is_showing_reoriented = false;
bool facetwise = false;

int main(int argc, char * argv[])
{
  using namespace std;
  cout<<R"(
Usage:

[space]  Toggle between original and reoriented faces
F,f      Toggle between patchwise and facetwise reorientation
S,s      Scramble colors
)";
  igl::read_triangle_mesh(TUTORIAL_SHARED_PATH "/truck.obj",V,F);

  const auto & scramble_colors = []()
  {
    for(int pass = 0;pass<2;pass++)
    {
      Eigen::MatrixXi R;
      igl::randperm(C[pass].maxCoeff()+1,R);
      C[pass] = igl::slice(R,Eigen::MatrixXi(C[pass]));
      Eigen::MatrixXd HSV(C[pass].rows(),3);
      HSV.col(0) = 
        360.*C[pass].array().cast<double>()/(double)C[pass].maxCoeff();
      HSV.rightCols(2).setConstant(1.0);
      igl::hsv_to_rgb(HSV,RGBcolors[pass]);
    }
    viewer.data().set_colors(RGBcolors[facetwise]);
  };

  viewer.callback_key_pressed = 
    [&scramble_colors]
    (igl::opengl::glfw::Viewer& /*viewer*/, unsigned int key, int mod)->bool
  {
    switch(key)
    {
    default:
      return false;
    case 'F':
    case 'f':
    {
      facetwise = !facetwise;
      break;
    }
    case 'S':
    case 's':
    {
      scramble_colors();
      return true;
    }
    case ' ':
    {
      is_showing_reoriented = !is_showing_reoriented;
      break;
    }
    }
    viewer.data().clear();
    viewer.data().set_mesh(V,is_showing_reoriented?FF[facetwise]:F);
    viewer.data().set_colors(RGBcolors[facetwise]);
    return true;
  };


  // Compute patches
  for(int pass = 0;pass<2;pass++)
  {
    Eigen::VectorXi I;
    igl::embree::reorient_facets_raycast(
      V,F,F.rows()*100,10,pass==1,false,false,I,C[pass]);
    // apply reorientation
    FF[pass].conservativeResize(F.rows(),F.cols());
    for(int i = 0;i<I.rows();i++)
    {
      if(I(i))
      {
        FF[pass].row(i) = (F.row(i).reverse()).eval();
      }else
      {
        FF[pass].row(i) = F.row(i);
      }
    }
  }

  viewer.data().set_mesh(V,is_showing_reoriented?FF[facetwise]:F);
  viewer.data().set_face_based(true);
  scramble_colors();
  viewer.launch();
}


