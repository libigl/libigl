#include <igl/circulation.h>
#include <igl/collapse_edge.h>
#include <igl/edge_flaps.h>
#include <igl/shortest_edge_and_midpoint.h>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <iostream>
#include <set>

#include "tutorial_shared_path.h"

int main(int argc, char * argv[])
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  cout<<"Usage: ./703_Decimation_bin [filename.(off|obj|ply)]"<<endl;
  cout<<"  [space]  toggle animation."<<endl;
  cout<<"  'r'  reset."<<endl;
  // Load a closed manifold mesh
  string filename(TUTORIAL_SHARED_PATH "/fertility.off");
  if(argc>=2)
  {
    filename = argv[1];
  }
  MatrixXd V,OV;
  MatrixXi F,OF;
  read_triangle_mesh(filename,OV,OF);

  igl::opengl::glfw::Viewer viewer;

  // Prepare array-based edge data structures and priority queue
  VectorXi EMAP;
  MatrixXi E,EF,EI;
  typedef std::set<std::pair<double,int> > PriorityQueue;
  PriorityQueue Q;
  std::vector<PriorityQueue::iterator > Qit;
  // If an edge were collapsed, we'd collapse it to these points:
  MatrixXd C;
  int num_collapsed;

  // Function to reset original mesh and data structures
  const auto & reset = [&]()
  {
    F = OF;
    V = OV;
    edge_flaps(F,E,EMAP,EF,EI);
    Qit.resize(E.rows());

    C.resize(E.rows(),V.cols());
    VectorXd costs(E.rows());
    Q.clear();
    for(int e = 0;e<E.rows();e++)
    {
      double cost = e;
      RowVectorXd p(1,3);
      shortest_edge_and_midpoint(e,V,F,E,EMAP,EF,EI,cost,p);
      C.row(e) = p;
      Qit[e] = Q.insert(std::pair<double,int>(cost,e)).first;
    }
    num_collapsed = 0;
    viewer.data().clear();
    viewer.data().set_mesh(V,F);
    viewer.data().set_face_based(true);
  };

  const auto &pre_draw = [&](igl::opengl::glfw::Viewer & viewer)->bool
  {
    // If animating then collapse 10% of edges
    if(viewer.core.is_animating && !Q.empty())
    {
      bool something_collapsed = false;
      // collapse edge
      const int max_iter = std::ceil(0.01*Q.size());
      for(int j = 0;j<max_iter;j++)
      {
        if(!collapse_edge(
          shortest_edge_and_midpoint, V,F,E,EMAP,EF,EI,Q,Qit,C))
        {
          break;
        }
        something_collapsed = true;
        num_collapsed++;
      }

      if(something_collapsed)
      {
        viewer.data().clear();
        viewer.data().set_mesh(V,F);
        viewer.data().set_face_based(true);
      }
    }
    return false;
  };

  const auto &key_down =
    [&](igl::opengl::glfw::Viewer &viewer,unsigned char key,int mod)->bool
  {
    switch(key)
    {
      case ' ':
        viewer.core.is_animating ^= 1;
        break;
      case 'R':
      case 'r':
        reset();
        break;
      default:
        return false;
    }
    return true;
  };

  reset();
  viewer.core.background_color.setConstant(1);
  viewer.core.is_animating = true;
  viewer.callback_key_down = key_down;
  viewer.callback_pre_draw = pre_draw;
  return viewer.launch();
}
