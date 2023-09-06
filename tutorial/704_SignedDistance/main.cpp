#include <igl/cat.h>
#include <igl/edge_lengths.h>
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


Eigen::MatrixXd V;
Eigen::MatrixXi T,F;

igl::AABB<Eigen::MatrixXd,3> tree;
igl::FastWindingNumberBVH fwn_bvh;

Eigen::MatrixXd FN,VN,EN;
Eigen::MatrixXi E;
Eigen::VectorXi EMAP;
double max_distance = 1;

double slice_z = 0.5;
bool overlay = false;

bool useFastWindingNumber = false;


const Eigen::MatrixXd CM = 
  (Eigen::MatrixXd(50,3)<<
  242,242,242,
  247,251,253,
  228,234,238,
  233,243,249,
  214,227,234,
  217,234,244,
  199,218,230,
  203,226,240,
  186,211,226,
  187,217,236,
  171,203,222,
  173,209,232,
  157,195,218,
  158,201,228,
  142,187,214,
  143,193,223,
  129,179,210,
  128,185,219,
  114,171,206,
  112,176,215,
  100,163,202,
  98,168,211,
  86,156,198,
  82,159,207,
  71,148,194,
  255,247,223,
  242,230,204,
  255,235,206,
  242,219,189,
  255,225,191,
  242,209,175,
  255,214,176,
  242,198,159,
  255,203,160,
  242,188,145,
  255,192,145,
  242,177,129,
  255,181,128,
  242,167,115,
  255,170,113,
  242,157,101,
  255,159,97,
  242,146,85,
  255,148,82,
  242,136,71,
  255,137,65,
  242,125,55,
  255,126,50,
  242,115,41,
  255,116,36).finished()/255.0;

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

  if (!useFastWindingNumber)
  {
    VectorXi I;
    MatrixXd N,C;
    // Bunny is a watertight mesh so use pseudonormal for signing
    signed_distance_pseudonormal(V_vis,V,F,tree,FN,VN,EN,EMAP,S_vis,I,C,N);
  } else {
    signed_distance_fast_winding_number(V_vis, V, F, tree, fwn_bvh, S_vis);
  }    

  const auto & append_mesh = [&F_vis,&V_vis](
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const RowVector3d & color)
  {
    F_vis.conservativeResize(F_vis.rows()+F.rows(),3);
    F_vis.bottomRows(F.rows()) = F.array()+V_vis.rows();
    V_vis.conservativeResize(V_vis.rows()+V.rows(),3);
    V_vis.bottomRows(V.rows()) = V;
  };
  if(overlay)
  {
    append_mesh(V,F,RowVector3d(0.8,0.8,0.8));
  }
  viewer.data().clear();
  viewer.data().set_mesh(V_vis,F_vis);
  viewer.data().set_colormap(CM);
  viewer.data().set_data(S_vis);
  viewer.core().lighting_factor = overlay;
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
    case '1':
      useFastWindingNumber = true;
      break;
    case '2':
      useFastWindingNumber = false;
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
  cout<< "1/2 toggle between fast winding number (1) and pseudonormal (2) signing. \n";
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

  // Fast winding and Pseudo normal depend on differnt AABB trees... We initialize both here.

  // Pseudonormal setup...
  // Precompute signed distance AABB tree
  tree.init(V,F);
  // Precompute vertex,edge and face normals
  igl::per_face_normals(V,F,FN);
  igl::per_vertex_normals(
    V,F,igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE,FN,VN);
  igl::per_edge_normals(
    V,F,igl::PER_EDGE_NORMALS_WEIGHTING_TYPE_UNIFORM,FN,EN,E,EMAP);

  // fast winding number setup (just init fwn bvh)
  igl::fast_winding_number(V, F, 2, fwn_bvh);

  // Plot the generated mesh
  igl::opengl::glfw::Viewer viewer;
  update_visualization(viewer);
  viewer.callback_key_down = &key_down;
  viewer.data().show_lines = false;
  viewer.launch();
}
