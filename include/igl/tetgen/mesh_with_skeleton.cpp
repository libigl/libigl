#include "mesh_with_skeleton.h"

#include <igl/sample_edges.h>
#include <igl/cat.h>
#include <igl/tetgen/tetrahedralize.h>
#include <igl/writeOFF.h>

#include <iostream>

bool igl::mesh_with_skeleton(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & C,
  const Eigen::VectorXi & /*P*/,
  const Eigen::MatrixXi & BE,
  const Eigen::MatrixXi & CE,
  const int samples_per_bone,
  Eigen::MatrixXd & VV,
  Eigen::MatrixXi & TT,
  Eigen::MatrixXi & FF)
{
  using namespace Eigen;
  using namespace igl;
  using namespace std;
  // Collect all edges that need samples:
  MatrixXi BECE = cat(1,BE,CE);
  MatrixXd S;
  // Sample each edge with 10 samples. (Choice of 10 doesn't seem to matter so
  // much, but could under some circumstances)
  sample_edges(C,BECE,samples_per_bone,S);
  // Vertices we'll constrain tet mesh to meet
  MatrixXd VS = cat(1,V,S);
  // Boundary faces
  MatrixXi BF;
  // Use tetgen to mesh the interior of surface, this assumes surface:
  //   * has no holes
  //   * has no non-manifold edges or vertices
  //   * has consistent orientation
  //   * has no self-intersections
  //   * has no 0-volume pieces
  // Default settings pq100 tell tetgen to mesh interior of triangle mesh and
  // to produce a graded tet mesh
  writeOFF("mesh_with_skeleton.off",VS,F);
  cerr<<"tetgen begin()"<<endl;
  int status = tetrahedralize( VS,F,"pq100Y",VV,TT,FF);
  cerr<<"tetgen end()"<<endl;
  if(FF.rows() != F.rows())
  {
    // Issue a warning if the surface has changed
    cerr<<"mesh_with_skeleton: Warning: boundary faces != input faces"<<endl;
  }
  if(status != 0)
  {
    cerr<<
      "***************************************************************"<<endl<<
      "***************************************************************"<<endl<<
      "***************************************************************"<<endl<<
      "***************************************************************"<<endl<<
      "* mesh_with_skeleton: tetgen failed. Just meshing convex hull *"<<endl<<
      "***************************************************************"<<endl<<
      "***************************************************************"<<endl<<
      "***************************************************************"<<endl<<
      "***************************************************************"<<endl;
    // If meshing convex hull then use more regular mesh
    status = tetrahedralize(VS,F,"q1.414",VV,TT,FF);
    // I suppose this will fail if the skeleton is outside the mesh
    assert(FF.maxCoeff() < VV.rows());
    if(status != 0)
    {
      cerr<<"mesh_with_skeleton: tetgen failed again."<<endl;
      return false;
    }
  }

  return true;
}

