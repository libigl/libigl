#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/jet.h>
#include <igl/shapeup.h>
#include <igl/shapeup_local_projections.h>
#include <igl/quad_planarity.h>
#include <igl/readDMAT.h>
#include <igl/readOFF.h>
#include <igl/slice.h>
#include <igl/viewer/Viewer.h>
#include <vector>
#include <cstdlib>

#include "tutorial_shared_path.h"

// Quad mesh loaded
Eigen::MatrixXd VQC;
Eigen::MatrixXi FQC;
Eigen::MatrixXi E;
Eigen::MatrixXi FQCtri;
Eigen::MatrixXd PQC0, PQC1, PQC2, PQC3;
// Euclidean-regular quad mesh
Eigen::MatrixXd VQCregular;
Eigen::MatrixXi FQCtriregular;
Eigen::MatrixXd PQC0regular, PQC1regular, PQC2regular, PQC3regular;

igl::ShapeupData su_data;



// Scale for visualizing the fields
double global_scale; //TODO: not used


bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier)
{
  using namespace std;
  using namespace Eigen;

  // Plot the original quad mesh
  if (key == '1')
  {
      cout<<"before setting mesh 1"<<endl;
    // Draw the triangulated quad mesh
    viewer.data.set_mesh(VQC, FQCtri);

    // Assign a color to each quad that corresponds to its planarity
    VectorXd planarity;
    igl::quad_planarity( VQC, FQC, planarity);
    MatrixXd Ct;
    igl::jet(planarity, 0, 0.01, Ct);
    MatrixXd C(FQCtri.rows(),3);
    C << Ct, Ct;
    viewer.data.set_colors(C);

    // Plot a line for each edge of the quad mesh
    viewer.data.add_edges(PQC0, PQC1, Eigen::RowVector3d(0,0,0));
    viewer.data.add_edges(PQC1, PQC2, Eigen::RowVector3d(0,0,0));
    viewer.data.add_edges(PQC2, PQC3, Eigen::RowVector3d(0,0,0));
    viewer.data.add_edges(PQC3, PQC0, Eigen::RowVector3d(0,0,0));
  }

  // Plot the planarized quad mesh
  if (key == '2')
  {
    // Draw the triangulated quad mesh
      cout<<"before setting mesh 2"<<endl;
    viewer.data.set_mesh(VQCregular, FQCtri);

    // Assign a color to each quad that corresponds to its planarity
    VectorXd planarity;
    igl::quad_planarity( VQCregular, FQC, planarity);
    MatrixXd Ct;
    igl::jet(planarity, 0, 0.01, Ct);
    MatrixXd C(FQCtri.rows(),3);
    C << Ct, Ct;
    viewer.data.set_colors(C);

    // Plot a line for each edge of the quad mesh
    viewer.data.add_edges(PQC0regular, PQC1regular, Eigen::RowVector3d(0,0,0));
    viewer.data.add_edges(PQC1regular, PQC2regular, Eigen::RowVector3d(0,0,0));
    viewer.data.add_edges(PQC2regular, PQC3regular, Eigen::RowVector3d(0,0,0));
    viewer.data.add_edges(PQC3regular, PQC0regular, Eigen::RowVector3d(0,0,0));
  }

  return false;
}

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;

  // Load a quad mesh
  igl::readOFF(TUTORIAL_SHARED_PATH "/halftunnel.off", VQC, FQC);

  // Convert it in a triangle mesh
  FQCtri.resize(2*FQC.rows(), 3);
  FQCtri <<  FQC.col(0),FQC.col(1),FQC.col(2),
             FQC.col(2),FQC.col(3),FQC.col(0);
  igl::slice( VQC, FQC.col(0).eval(), 1, PQC0);
  igl::slice( VQC, FQC.col(1).eval(), 1, PQC1);
  igl::slice( VQC, FQC.col(2).eval(), 1, PQC2);
  igl::slice( VQC, FQC.col(3).eval(), 1, PQC3);

  // Create a planar version with ShapeUp
  //igl::planarize_quad_mesh(VQC, FQC, 100, 0.005, VQCregular);
    
  E.resize(FQC.size(),2);
  E.col(0)<<FQC.col(0),FQC.col(1),FQC.col(2),FQC.col(3);
  E.col(1)<<FQC.col(1),FQC.col(2),FQC.col(3),FQC.col(0);
    
  VectorXi b;
  VectorXd w(FQC.rows());
  MatrixXd bc;
    
  VectorXi array_of_fours=VectorXi::Constant(FQC.rows(),4);
    cout<<"before pre-computation"<<endl;
    std::function<bool(const MatrixXd&, const VectorXi&, const MatrixXi&, MatrixXd&)> localFunction=std::function<bool(const MatrixXd&, const VectorXi&, const MatrixXi&, MatrixXd&)>(igl::shapeup_identity_projection);
    //shapeup_precomputation(VQC, array_of_fours,FQC,E,b,w, localFunction,su_data);
    cout<<"after pre-computation"<<endl;
    shapeup_solve(bc,VQC,su_data,VQCregular);
    cout<<"after computation"<<endl;
    

  // Convert the planarized mesh to triangles
  igl::slice( VQCregular, FQC.col(0).eval(), 1, PQC0regular);
  igl::slice( VQCregular, FQC.col(1).eval(), 1, PQC1regular);
  igl::slice( VQCregular, FQC.col(2).eval(), 1, PQC2regular);
  igl::slice( VQCregular, FQC.col(3).eval(), 1, PQC3regular);

  // Launch the viewer
  igl::viewer::Viewer viewer;
  key_down(viewer,'1',0);
  viewer.core.invert_normals = true;
  viewer.core.show_lines = false;
  viewer.callback_key_down = &key_down;
  viewer.launch();
}
