// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "writePLY.h"
#include <vector>

#include <igl/ply.h>
#include <vector>

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedN,
  typename DerivedUV>
IGL_INLINE bool igl::writePLY(
  const std::string & filename,
  const Eigen::PlainObjectBase<DerivedV> & V,
  const Eigen::PlainObjectBase<DerivedF> & F,
  const Eigen::PlainObjectBase<DerivedN> & N,
  const Eigen::PlainObjectBase<DerivedUV> & UV,
  const bool ascii)
{
  // Largely based on obj2ply.c

  typedef struct Vertex
  {
    double x,y,z,w;          /* position */
    double nx,ny,nz;         /* surface normal */
    double s,t;              /* texture coordinates */
  } Vertex;

  typedef struct Face
  {
    unsigned char nverts;    /* number of vertex indices in list */
    int *verts;              /* vertex index list */
  } Face;

  PlyProperty vert_props[] =
  { /* list of property information for a vertex */
    {"x", PLY_DOUBLE, PLY_DOUBLE, offsetof(Vertex,x), 0, 0, 0, 0},
    {"y", PLY_DOUBLE, PLY_DOUBLE, offsetof(Vertex,y), 0, 0, 0, 0},
    {"z", PLY_DOUBLE, PLY_DOUBLE, offsetof(Vertex,z), 0, 0, 0, 0},
    {"nx", PLY_DOUBLE, PLY_DOUBLE, offsetof(Vertex,nx), 0, 0, 0, 0},
    {"ny", PLY_DOUBLE, PLY_DOUBLE, offsetof(Vertex,ny), 0, 0, 0, 0},
    {"nz", PLY_DOUBLE, PLY_DOUBLE, offsetof(Vertex,nz), 0, 0, 0, 0},
    {"s", PLY_DOUBLE, PLY_DOUBLE, offsetof(Vertex,s), 0, 0, 0, 0},
    {"t", PLY_DOUBLE, PLY_DOUBLE, offsetof(Vertex,t), 0, 0, 0, 0},
  };

  PlyProperty face_props[] =
  { /* list of property information for a face */
    {"vertex_indices", PLY_INT, PLY_INT, offsetof(Face,verts),
      1, PLY_UCHAR, PLY_UCHAR, offsetof(Face,nverts)},
  };
  const bool has_normals = N.rows() > 0;
  const bool has_texture_coords = UV.rows() > 0;
  std::vector<Vertex> vlist(V.rows());
  std::vector<Face> flist(F.rows());
  for(size_t i = 0;i<(size_t)V.rows();i++)
  {
    vlist[i].x = V(i,0);
    vlist[i].y = V(i,1);
    vlist[i].z = V(i,2);
    if(has_normals)
    {
      vlist[i].nx = N(i,0);
      vlist[i].ny = N(i,1);
      vlist[i].nz = N(i,2);
    }
    if(has_texture_coords)
    {
      vlist[i].s = UV(i,0);
      vlist[i].t = UV(i,1);
    }
  }
  for(size_t i = 0;i<(size_t)F.rows();i++)
  {
    flist[i].nverts = F.cols();
    flist[i].verts = new int[F.cols()];
    for(size_t c = 0;c<(size_t)F.cols();c++)
    {
      flist[i].verts[c] = F(i,c);
    }
  }

  const char * elem_names[] = {"vertex","face"};
  FILE * fp = fopen(filename.c_str(),"w");
  if(fp==NULL)
  {
    return false;
  }
  PlyFile * ply = ply_write(fp, 2,elem_names,
      (ascii ? PLY_ASCII : PLY_BINARY_LE));
  if(ply==NULL)
  {
    return false;
  }

  std::vector<PlyProperty> plist;
  plist.push_back(vert_props[0]);
  plist.push_back(vert_props[1]);
  plist.push_back(vert_props[2]);
  if (has_normals)
  {
    plist.push_back(vert_props[3]);
    plist.push_back(vert_props[4]);
    plist.push_back(vert_props[5]);
  }
  if (has_texture_coords)
  {
    plist.push_back(vert_props[6]);
    plist.push_back(vert_props[7]);
  }
  ply_describe_element(ply, "vertex", V.rows(),plist.size(),
    &plist[0]);

  ply_describe_element(ply, "face", F.rows(),1,&face_props[0]);
  ply_header_complete(ply);
  ply_put_element_setup(ply, "vertex");
  for(const auto v : vlist)
  {
    ply_put_element(ply, (void *) &v);
  }
  ply_put_element_setup(ply, "face");
  for(const auto f : flist)
  {
    ply_put_element(ply, (void *) &f);
  }

  ply_close(ply);
  fclose(fp);
  for(size_t i = 0;i<(size_t)F.rows();i++)
  {
    delete[] flist[i].verts;
  }
  return true;
}

template <
  typename DerivedV,
  typename DerivedF>
IGL_INLINE bool igl::writePLY(
  const std::string & filename,
  const Eigen::PlainObjectBase<DerivedV> & V,
  const Eigen::PlainObjectBase<DerivedF> & F,
  const bool ascii)
{
  Eigen::MatrixXd N,UV;
  return writePLY(filename,V,F,N,UV,ascii);
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template bool igl::writePLY<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(std::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, bool);
template bool igl::writePLY<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3> >(std::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, bool);
#endif
