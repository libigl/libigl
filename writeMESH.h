#ifndef IGL_WRITEMESH_H
#define IGL_WRITEMESH_H

#include <string>
#include <vector>

namespace igl
{
  // save a tetrahedral volume mesh to a .mesh file
  //
  // Templates:
  //   Scalar  type for positions and vectors (will be cast as double)
  //   Index  type for indices (will be cast to int)
  // Input:
  //   mesh_file_name  path of .mesh file
  // Outputs:
  //   V  double matrix of vertex positions  #V by 3
  //   T  #T list of tet indices into vertex positions
  //   F  #F list of face indices into vertex positions
  template <typename Scalar, typename Index>
  inline bool writeMESH(
    const std::string mesh_file_name,
    std::vector<std::vector<Scalar > > & V,
    std::vector<std::vector<Index > > & T,
    std::vector<std::vector<Index > > & F);

  // Input:
  //   mesh_file_name  path of .mesh file
  // Outputs:
  //   V  eigen double matrix #V by 3
  //   T  eigen int matrix #T by 4
  //   F  eigen int matrix #F by 3
  inline bool writeMESH(
    const std::string str,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& T,
    Eigen::MatrixXi& F);
}

// Implementation
#include <cstdio>
#include "verbose.h"

template <typename Scalar, typename Index>
inline bool igl::writeMESH(
  const std::string mesh_file_name,
  std::vector<std::vector<Scalar > > & V,
  std::vector<std::vector<Index > > & T,
  std::vector<std::vector<Index > > & F)
{
  // not implemented but should be
  assert(false);
  return false;
}

#include <Eigen/Core>

inline bool igl::writeMESH(
  const std::string str,
  Eigen::MatrixXd& V,
  Eigen::MatrixXi& T,
  Eigen::MatrixXi& F)
{
  using namespace std;
  using namespace igl;
  using namespace Eigen;

  FILE * mesh_file = fopen(str.c_str(),"w");
  if(NULL==mesh_file)
  {
    fprintf(stderr,"IOError: %s could not be opened...",str.c_str());
    return false;
  }
  // print header
  fprintf(mesh_file,"MeshVersionFormatted 1\n");
  fprintf(mesh_file,"Dimension 3\n");
  // print tet vertices
  fprintf(mesh_file,"Vertices\n");
  // print number of tet vertices
  int number_of_tet_vertices = V.rows();
  fprintf(mesh_file,"%d\n",number_of_tet_vertices);
  // loop over tet vertices
  for(int i = 0;i<number_of_tet_vertices;i++)
  {
    // print position of ith tet vertex
    fprintf(mesh_file,"%lg %lg %lg 1\n",
      (double)V(i,0),
      (double)V(i,1),
      (double)V(i,2));
  }
  verbose("WARNING: save_mesh() assumes that vertices have"
      " same indices in surface as volume...\n");
  // print faces
  fprintf(mesh_file,"Triangles\n");
  // print number of triangles
  int number_of_triangles = F.rows();
  fprintf(mesh_file,"%d\n",number_of_triangles);
  // loop over faces
  for(int i = 0;i<number_of_triangles;i++)
  {
    // loop over vertices in face
    fprintf(mesh_file,"%d %d %d 1\n", 
      (int)F(i,0)+1, 
      (int)F(i,1)+1, 
      (int)F(i,2)+1);
  }
  // print tetrahedra
  fprintf(mesh_file,"Tetrahedra\n");
  int number_of_tetrahedra = T.rows();
  // print number of tetrahedra
  fprintf(mesh_file,"%d\n",number_of_tetrahedra);
  // loop over tetrahedra
  for(int i = 0; i < number_of_tetrahedra;i++)
  {
    // mesh standard uses 1-based indexing
    fprintf(mesh_file, "%d %d %d %d 1\n",
      (int)T(i,0)+1,
      (int)T(i,1)+1,
      (int)T(i,2)+1,
      (int)T(i,3)+1);
  }
  fclose(mesh_file);
  return true;
}

#endif
