#include <igl/opengl/gl.h>
#include <igl/arap.h>
#include <igl/biharmonic_coordinates.h>
#include <igl/cat.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/matrix_to_list.h>
#include <igl/parula.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/readDMAT.h>
#include <igl/readMESH.h>
#include <igl/remove_unreferenced.h>
#include <igl/slice.h>
#include <igl/writeDMAT.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Sparse>
#include <iostream>
#include <queue>

#include "tutorial_shared_path.h"

struct Mesh
{
  Eigen::MatrixXd V,U;
  Eigen::MatrixXi T,F;
} low,high,scene;

Eigen::MatrixXd W;
igl::ARAPData arap_data;

int main(int argc, char * argv[])
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;
  if(!readMESH(TUTORIAL_SHARED_PATH "/octopus-low.mesh",low.V,low.T,low.F))
  {
    cout<<"failed to load mesh"<<endl;
  }
  if(!readMESH(TUTORIAL_SHARED_PATH "/octopus-high.mesh",high.V,high.T,high.F))
  {
    cout<<"failed to load mesh"<<endl;
  }

  // Precomputation
  {
    Eigen::VectorXi b;
    {
      Eigen::VectorXi J = Eigen::VectorXi::LinSpaced(high.V.rows(),0,high.V.rows()-1);
      Eigen::VectorXd sqrD;
      Eigen::MatrixXd _2;
      cout<<"Finding closest points..."<<endl;
      igl::point_mesh_squared_distance(low.V,high.V,J,sqrD,b,_2);
      assert(sqrD.minCoeff() < 1e-7 && "low.V should exist in high.V");
    }
    // force perfect positioning, rather have popping in low-res than high-res.
    // The correct/elaborate thing to do is express original low.V in terms of
    // linear interpolation (or extrapolation) via elements in (high.V,high.F)
    igl::slice(high.V,b,1,low.V);
    // list of points --> list of singleton lists
    std::vector<std::vector<int> > S;
    igl::matrix_to_list(b,S);
    cout<<"Computing weights for "<<b.size()<<
      " handles at "<<high.V.rows()<<" vertices..."<<endl;
    // Technically k should equal 3 for smooth interpolation in 3d, but 2 is
    // faster and looks OK
    const int k = 2;
    igl::biharmonic_coordinates(high.V,high.T,S,k,W);
    cout<<"Reindexing..."<<endl;
    // Throw away interior tet-vertices, keep weights and indices of boundary
    VectorXi I,J;
    igl::remove_unreferenced(high.V.rows(),high.F,I,J);
    for_each(high.F.data(),high.F.data()+high.F.size(),[&I](int & a){a=I(a);});
    for_each(b.data(),b.data()+b.size(),[&I](int & a){a=I(a);});
    igl::slice(MatrixXd(high.V),J,1,high.V);
    igl::slice(MatrixXd(W),J,1,W);
  }

  // Resize low res (high res will also be resized by affine precision of W)
  low.V.rowwise() -= low.V.colwise().mean();
  low.V /= (low.V.maxCoeff()-low.V.minCoeff());
  low.V.rowwise() += RowVector3d(0,1,0);
  low.U = low.V;
  high.U = high.V;

  arap_data.with_dynamics = true;
  arap_data.max_iter = 10;
  arap_data.energy = ARAP_ENERGY_TYPE_DEFAULT;
  arap_data.h = 0.01;
  arap_data.ym = 0.001;
  if(!arap_precomputation(low.V,low.T,3,VectorXi(),arap_data))
  {
    cerr<<"arap_precomputation failed."<<endl;
    return EXIT_FAILURE;
  }
  // Constant gravitational force
  Eigen::SparseMatrix<double> M;
  igl::massmatrix(low.V,low.T,igl::MASSMATRIX_TYPE_DEFAULT,M);
  const size_t n = low.V.rows();
  arap_data.f_ext =  M * RowVector3d(0,-9.8,0).replicate(n,1);
  // Random initial velocities to wiggle things
  arap_data.vel = MatrixXd::Random(n,3);
  
  igl::opengl::glfw::Viewer viewer;
  // Create one huge mesh containing both meshes
  igl::cat(1,low.U,high.U,scene.U);
  igl::cat(1,low.F,MatrixXi(high.F.array()+low.V.rows()),scene.F);
  // Color each mesh
  viewer.data().set_mesh(scene.U,scene.F);
  MatrixXd C(scene.F.rows(),3);
  C<<
    RowVector3d(0.8,0.5,0.2).replicate(low.F.rows(),1),
    RowVector3d(0.3,0.4,1.0).replicate(high.F.rows(),1);
  viewer.data().set_colors(C);

  viewer.callback_key_pressed = 
    [&](igl::opengl::glfw::Viewer & viewer,unsigned int key,int mods)->bool
  {
    switch(key)
    {
      default: 
        return false;
      case ' ':
        viewer.core().is_animating = !viewer.core().is_animating;
        return true;
      case 'r':
        low.U = low.V;
        return true;
    }
  };
  viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer & viewer)->bool
  {
    glEnable(GL_CULL_FACE);
    if(viewer.core().is_animating)
    {
      arap_solve(MatrixXd(0,3),arap_data,low.U);
      for(int v = 0;v<low.U.rows();v++)
      {
        // collide with y=0 plane
        const int y = 1;
        if(low.U(v,y) < 0)
        {
          low.U(v,y) = -low.U(v,y);
          // ~ coefficient of restitution
          const double cr = 1.1;
          arap_data.vel(v,y) = - arap_data.vel(v,y) / cr;
        }
      }

      scene.U.block(0,0,low.U.rows(),low.U.cols()) = low.U;
      high.U = W * (low.U.rowwise() + RowVector3d(1,0,0));
      scene.U.block(low.U.rows(),0,high.U.rows(),high.U.cols()) = high.U;

      viewer.data().set_vertices(scene.U);
      viewer.data().compute_normals();
    }
    return false;
  };
  viewer.data().show_lines = false;
  viewer.core().is_animating = true;
  viewer.core().animation_max_fps = 30.;
  viewer.data().set_face_based(true);
  cout<<R"(
[space] to toggle animation
'r'     to reset positions 
      )";
  viewer.core().rotation_type = 
    igl::opengl::ViewerCore::ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP;
  viewer.launch();
}
