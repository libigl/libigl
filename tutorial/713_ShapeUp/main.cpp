#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/jet.h>
#include <igl/shapeup.h>
#include <igl/quad_planarity.h>
#include <igl/readDMAT.h>
#include <igl/readOFF.h>
#include <igl/slice.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/PI.h>
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

void quadAngleRegularity(const Eigen::MatrixXd& V, const Eigen::MatrixXi& Q, Eigen::VectorXd& angleRegularity)
{
  angleRegularity.conservativeResize(Q.rows());
  angleRegularity.setZero();
  for (int i=0;i<Q.rows();i++){
    for (int j=0;j<4;j++){
      Eigen::RowVectorXd v21=(V.row(Q(i,j))-V.row(Q(i,(j+1)%4))).normalized();
      Eigen::RowVectorXd v23=(V.row(Q(i,(j+2)%4))-V.row(Q(i,(j+1)%4))).normalized();
  
      angleRegularity(i)+=(abs(acos(v21.dot(v23))-igl::PI/2.0)/(igl::PI/2.0))/4.0;
    }
  }
}


bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
  using namespace std;
  using namespace Eigen;

  // Plot the original quad mesh
  
  if (key == '1')
  {
    viewer.data().clear();
    // Draw the triangulated quad mesh
    viewer.data().set_mesh(VQC, FQCtri);

    // Assign a color to each quad that corresponds to the average deviation of each angle from pi/2
    VectorXd angleRegularity(FQC.rows());
    quadAngleRegularity( VQC, FQC, angleRegularity);
    MatrixXd Ct;
    igl::jet(angleRegularity, 0.0, 0.05, Ct);
    MatrixXd C(FQCtri.rows(),3);
    C << Ct, Ct;
    viewer.data().set_colors(C);

    // Plot a line for each edge of the quad mesh
    viewer.data().add_edges(PQC0, PQC1, Eigen::RowVector3d(0,0,0));
    viewer.data().add_edges(PQC1, PQC2, Eigen::RowVector3d(0,0,0));
    viewer.data().add_edges(PQC2, PQC3, Eigen::RowVector3d(0,0,0));
    viewer.data().add_edges(PQC3, PQC0, Eigen::RowVector3d(0,0,0));
  }

  // Plot the planarized quad mesh
  if (key == '2')
  {
    viewer.data().clear();
    // Draw the triangulated quad mesh
    viewer.data().set_mesh(VQCregular, FQCtri);

    // Assign a color to each quad that corresponds to its planarity
    VectorXd angleRegularity(FQC.rows());
    quadAngleRegularity( VQCregular, FQC, angleRegularity);
    MatrixXd Ct;
    igl::jet(angleRegularity, 0, 0.05, Ct);
    MatrixXd C(FQCtri.rows(),3);
    C << Ct, Ct;
    viewer.data().set_colors(C);

    // Plot a line for each edge of the quad mesh
    viewer.data().add_edges(PQC0regular, PQC1regular, Eigen::RowVector3d(0,0,0));
    viewer.data().add_edges(PQC1regular, PQC2regular, Eigen::RowVector3d(0,0,0));
    viewer.data().add_edges(PQC2regular, PQC3regular, Eigen::RowVector3d(0,0,0));
    viewer.data().add_edges(PQC3regular, PQC0regular, Eigen::RowVector3d(0,0,0));
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
  
  VectorXi b(1); b(0)=0;  //setting the first vertex to be the same.
  
  VectorXd wShape=VectorXd::Constant(FQC.rows(),1.0);
  VectorXd wSmooth=VectorXd::Constant(E.rows(),1.0);
  MatrixXd bc(1,3); bc<<VQC.row(0);
  
  VectorXi array_of_fours=VectorXi::Constant(FQC.rows(),4);
  igl::shapeup_projection_function localFunction(igl::shapeup_regular_face_projection);
  
  su_data.maxIterations=200;
  shapeup_precomputation(VQC, array_of_fours,FQC,E,b,wShape, wSmooth,su_data);
  shapeup_solve(bc,localFunction, VQC,su_data, false,VQCregular);
  

  // Convert the planarized mesh to triangles
  igl::slice( VQCregular, FQC.col(0).eval(), 1, PQC0regular);
  igl::slice( VQCregular, FQC.col(1).eval(), 1, PQC1regular);
  igl::slice( VQCregular, FQC.col(2).eval(), 1, PQC2regular);
  igl::slice( VQCregular, FQC.col(3).eval(), 1, PQC3regular);

  // Launch the viewer
  igl::opengl::glfw::Viewer viewer;
  key_down(viewer,'1',0);
  viewer.data().invert_normals = true;
  viewer.data().show_lines = false;
  viewer.callback_key_down = &key_down;
  viewer.launch();
}
