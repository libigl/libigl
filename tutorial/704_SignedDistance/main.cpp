#include <igl/cat.h>
#include <igl/edge_lengths.h>
#include <igl/parula.h>
#include <igl/per_edge_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/readMESH.h>
#include <igl/signed_distance.h>
#include <igl/slice_mask.h>
#include <igl/marching_tets.h>
#include <igl/upsample.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/writeOBJ.h>
#include <Eigen/Sparse>
#include <iostream>

#include "tutorial_shared_path.h"

Eigen::MatrixXd V;
Eigen::MatrixXi T,F;
igl::AABB<Eigen::MatrixXd,3> tree;
Eigen::MatrixXd FN,VN,EN;
Eigen::MatrixXi E;
Eigen::VectorXi EMAP;
double max_distance = 1;

double slice_z = 0.5;
bool overlay = false;

void update_visualization(igl::opengl::glfw::Viewer & viewer)
{
  using namespace Eigen;
  using namespace std;
  Eigen::Vector4d plane(
    0,0,1,-((1-slice_z)*V.col(2).minCoeff()+slice_z*V.col(2).maxCoeff()));
  MatrixXd V_vis;
  MatrixXi F_vis;
  // Extract triangle mesh slice through volume mesh and subdivide nasty
  // triangles
  {
    VectorXi J;
    SparseMatrix<double> bary;
    {
      // Value of plane's implicit function at all vertices
      const VectorXd IV = 
        (V.col(0)*plane(0) + 
         V.col(1)*plane(1) + 
         V.col(2)*plane(2)).array()
        + plane(3);
      igl::marching_tets(V,T,IV,V_vis,F_vis,J,bary);
      igl::writeOBJ("vis.obj",V_vis,F_vis);
    }
    while(true)
    {
      MatrixXd l;
      igl::edge_lengths(V_vis,F_vis,l);
      l /= (V_vis.colwise().maxCoeff() - V_vis.colwise().minCoeff()).norm();
      const double max_l = 0.03;
      if(l.maxCoeff()<max_l)
      {
        break;
      }
      Array<bool,Dynamic,1> bad = l.array().rowwise().maxCoeff() > max_l;
      MatrixXi F_vis_bad, F_vis_good;
      igl::slice_mask(F_vis,bad,1,F_vis_bad);
      igl::slice_mask(F_vis,(bad!=true).eval(),1,F_vis_good);
      igl::upsample(V_vis,F_vis_bad);
      F_vis = igl::cat(1,F_vis_bad,F_vis_good);
    }
  }

  // Compute signed distance
  VectorXd S_vis;
  {
    VectorXi I;
    MatrixXd N,C;
    // Bunny is a watertight mesh so use pseudonormal for signing
    signed_distance_pseudonormal(V_vis,V,F,tree,FN,VN,EN,EMAP,S_vis,I,C,N);
  }
  // push to [0,1] range
  S_vis.array() = 0.5*(S_vis.array()/max_distance)+0.5;
  MatrixXd C_vis;
  // color without normalizing
  igl::parula(S_vis,false,C_vis);


  const auto & append_mesh = [&C_vis,&F_vis,&V_vis](
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const RowVector3d & color)
  {
    F_vis.conservativeResize(F_vis.rows()+F.rows(),3);
    F_vis.bottomRows(F.rows()) = F.array()+V_vis.rows();
    V_vis.conservativeResize(V_vis.rows()+V.rows(),3);
    V_vis.bottomRows(V.rows()) = V;
    C_vis.conservativeResize(C_vis.rows()+V.rows(),3);
    C_vis.bottomRows(V.rows()).rowwise() = color;
  };
  if(overlay)
  {
    append_mesh(V,F,RowVector3d(0.8,0.8,0.8));
  }
  viewer.data().clear();
  viewer.data().set_mesh(V_vis,F_vis);
  viewer.data().set_colors(C_vis);
  viewer.core.lighting_factor = overlay;
}

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int mod)
{
  switch(key)
  {
    default:
      return false;
    case ' ':
      overlay ^= true;
      break;
    case '.':
      slice_z = std::min(slice_z+0.01,0.99);
      break;
    case ',':
      slice_z = std::max(slice_z-0.01,0.01);
      break;
  }
  update_visualization(viewer);
  return true;
}

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;

  cout<<"Usage:"<<endl;
  cout<<"[space]  toggle showing surface."<<endl;
  cout<<"'.'/','  push back/pull forward slicing plane."<<endl;
  cout<<endl;

  // Load mesh: (V,T) tet-mesh of convex hull, F contains original surface
  // triangles
  igl::readMESH(TUTORIAL_SHARED_PATH "/bunny.mesh",V,T,F);


  // Encapsulated call to point_mesh_squared_distance to determine bounds
  {
    VectorXd sqrD;
    VectorXi I;
    MatrixXd C;
    igl::point_mesh_squared_distance(V,V,F,sqrD,I,C);
    max_distance = sqrt(sqrD.maxCoeff());
  }

  // Precompute signed distance AABB tree
  tree.init(V,F);
  // Precompute vertex,edge and face normals
  igl::per_face_normals(V,F,FN);
  igl::per_vertex_normals(
    V,F,igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE,FN,VN);
  igl::per_edge_normals(
    V,F,igl::PER_EDGE_NORMALS_WEIGHTING_TYPE_UNIFORM,FN,EN,E,EMAP);

  // Plot the generated mesh
  igl::opengl::glfw::Viewer viewer;
  update_visualization(viewer);
  viewer.callback_key_down = &key_down;
  viewer.data().show_lines = false;
  viewer.launch();
}
