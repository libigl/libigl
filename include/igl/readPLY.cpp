// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "readPLY.h"
#include "list_to_matrix.h"
#include "tinyply.h"
#include <fstream>
#include <iostream>

template <
  typename Vtype,
  typename Ftype,
  typename Ntype,
  typename UVtype>
IGL_INLINE bool igl::readPLY(
  const std::string & filename,
  std::vector<std::vector<Vtype> > & V,
  std::vector<std::vector<Ftype> > & F,
  std::vector<std::vector<Ntype> > & N,
  std::vector<std::vector<UVtype> >  & UV)
{
  try {
    std::ifstream ply_stream(filename, std::ios::binary);
    if (ply_stream.fail())
    {
          std::cerr << "ReadPLY: Error opening file " << filename << std::endl;
          return false;
    }
    return readPLY(ply_stream,V,F,N,UV);
  } catch (const std::exception& e) {
    std::cerr << "ReadPLY error: " << filename << e.what() << std::endl;
    return false;
  }
}

template <
  typename Vtype,
  typename Ftype,
  typename Ntype,
  typename UVtype>
IGL_INLINE bool igl::readPLY(
  std::istream & ply_stream,
  std::vector<std::vector<Vtype> > & V,
  std::vector<std::vector<Ftype> > & F,
  std::vector<std::vector<Ntype> > & N,
  std::vector<std::vector<UVtype> > & UV)
{   
  igl::tinyply::PlyFile file;
  file.parse_header(ply_stream);

  /*
  std::cout << "igl::readPLY Header..........................................................\n";
  for (auto c : file.get_comments())
      std::cout << "Comment: " << c << std::endl;
  for (auto e : file.get_elements())
  {
      std::cout << "element - " << e.name << " (" << e.size << ")" << std::endl;
      for (auto p : e.properties)
          std::cout << "\tproperty - " << p.name << " (" << igl::tinyply::PropertyTable[p.propertyType].str << ")" << std::endl;
  }
  std::cout << "igl::readPLY Data..........................................................\n";
  //*/

  // Tinyply treats parsed data as untyped byte buffers.
  std::shared_ptr<igl::tinyply::PlyData> vertices, normals, faces, texcoords;

  // The header information can be used to programmatically extract properties on elements
  // known to exist in the header prior to reading the data. For brevity of this sample, properties
  // like vertex position are hard-coded:  
  try {
    vertices = file.request_properties_from_element("vertex", { "x", "y", "z" });
  }
  catch (const std::exception & e) { }

  try {
    normals = file.request_properties_from_element("vertex", { "nx", "ny", "nz" });
  }
  catch (const std::exception & e) { }

  //Try texture coordinates with several names
  try {
    //texture_u texture_v are the names used by meshlab to store textures
    texcoords = file.request_properties_from_element("vertex", { "texture_u", "texture_v" });
  }
  catch (const std::exception & e) { }
  if (!texcoords)
  {
      try {
        //u v are the naive names
        texcoords = file.request_properties_from_element("vertex", { "u", "v" });
      }
      catch (const std::exception & e) { }

  }
  if (!texcoords)
  {
      try {
        //s t were the names used by blender and the previous libigl PLY reader.
        texcoords = file.request_properties_from_element("vertex", { "s", "t" });
      }
      catch (const std::exception & e) { }

  }

  // Providing a list size hint (the last argument) is a 2x performance improvement. If you have
  // arbitrary ply files, it is best to leave this 0.
  try {
    faces = file.request_properties_from_element("face", { "vertex_indices" }, 3);
  }
  catch (const std::exception & e) { }

  // Parse the geometry data
  file.read(ply_stream);

  if (vertices) {
    //std::cout << "\tRead " << vertices->count << " total vertices "<< std::endl;
    V.resize(vertices->count, std::vector<Vtype>(3));

    if (vertices->t == igl::tinyply::Type::FLOAT32)
    {
      const float *v = reinterpret_cast<const float *>(vertices->buffer.get());
      for (size_t i=0; i< vertices->count; ++i)
      {
        V[i][0]=Vtype(v[3*i]);
        V[i][1]=Vtype(v[3*i+1]);
        V[i][2]=Vtype(v[3*i+2]);
      }
    }
    else if (vertices->t == igl::tinyply::Type::FLOAT64)
    {
      const double *v = reinterpret_cast<const double *>(vertices->buffer.get());
      for (size_t i=0; i< vertices->count; ++i)
      {
        V[i][0]=Vtype(v[3*i]);
        V[i][1]=Vtype(v[3*i+1]);
        V[i][2]=Vtype(v[3*i+2]);
      }
    }
  }
  else {
    V.resize(0);
  }

  if (normals) {
    //std::cout << "\tRead " << normals->count << " total vertex normals " << std::endl;
    N.resize(normals->count, std::vector<Ntype>(3));

    if (normals->t == igl::tinyply::Type::FLOAT32)
    {
      const float *v = reinterpret_cast<const float *>(normals->buffer.get());
      for (size_t i=0; i< normals->count; ++i)
      {
        N[i][0]=Ntype(v[3*i]);
        N[i][1]=Ntype(v[3*i+1]);
        N[i][2]=Ntype(v[3*i+2]);
      }
    }
    else if (normals->t == igl::tinyply::Type::FLOAT64)
    {
      const double *v = reinterpret_cast<const double *>(normals->buffer.get());
      for (size_t i=0; i< normals->count; ++i)
      {
        N[i][0]=Ntype(v[3*i]);
        N[i][1]=Ntype(v[3*i+1]);
        N[i][2]=Ntype(v[3*i+2]);
      }
    }
  }
  else {
      N.resize(0);
  }

  if (texcoords) {
    //std::cout << "\tRead " << texcoords->count << " total vertex texcoords " << std::endl;
    UV.resize(texcoords->count, std::vector<UVtype>(2));

    if (texcoords->t == igl::tinyply::Type::FLOAT32)
    {
      const float *uv = reinterpret_cast<const float *>(texcoords->buffer.get());
      for (size_t i=0; i< texcoords->count; ++i)
      {
        UV[i][0]=UVtype(uv[2*i]);
        UV[i][1]=UVtype(uv[2*i+1]);
      }
    }
    else if (texcoords->t == igl::tinyply::Type::FLOAT64)
    {
      const double *uv = reinterpret_cast<const double *>(texcoords->buffer.get());
      for (size_t i=0; i< texcoords->count; ++i)
      {
        UV[i][0]=UVtype(uv[2*i]);
        UV[i][1]=UVtype(uv[2*i+1]);
      }
    }
  }
  else {
    UV.resize(0);
  }

  if (faces) {
    //std::cout << "\tRead " << faces->count << " total faces (triangles) " << std::endl;
    F.resize(faces->count, std::vector<Ftype>(3));

    //Set triangles to this mesh
    unsigned *f = reinterpret_cast<unsigned *>(faces->buffer.get());
    for (size_t i=0; i< faces->count; ++i)
    {
        F[i][0]=Ftype(f[3*i]);
        F[i][1]=Ftype(f[3*i+1]);
        F[i][2]=Ftype(f[3*i+2]);
    }
  }
  else {
    F.resize(0);
  }

  return true;
}

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedN,
  typename DerivedUV>
IGL_INLINE bool igl::readPLY(
  const std::string & filename,
  Eigen::PlainObjectBase<DerivedV> & V,
  Eigen::PlainObjectBase<DerivedF> & F,
  Eigen::PlainObjectBase<DerivedN> & N,
  Eigen::PlainObjectBase<DerivedUV> & UV)
{
  std::vector<std::vector<typename DerivedV::Scalar> > vV;
  std::vector<std::vector<typename DerivedF::Scalar> > vF;
  std::vector<std::vector<typename DerivedN::Scalar> > vN;
  std::vector<std::vector<typename DerivedUV::Scalar> > vUV;
  if(!readPLY(filename,vV,vF,vN,vUV))
  {
    return false;
  }
  return 
    list_to_matrix(vV,V) &&
    list_to_matrix(vF,F) &&
    list_to_matrix(vN,N) &&
    list_to_matrix(vUV,UV);
}

template <
  typename DerivedV,
  typename DerivedF>
IGL_INLINE bool igl::readPLY(
  const std::string & filename,
  Eigen::PlainObjectBase<DerivedV> & V,
  Eigen::PlainObjectBase<DerivedF> & F)
{
  Eigen::MatrixXd N,UV;
  return readPLY(filename,V,F,N,UV);
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template bool igl::readPLY<double, int, double, double>(std::string const &, std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >&, std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > >&, std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >&, std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >&);
template bool igl::readPLY<Eigen::Matrix<float, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3> >(std::string const&, Eigen::PlainObjectBase<Eigen::Matrix<float, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&);
template bool igl::readPLY<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3> >(std::string const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&);
template bool igl::readPLY<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >(std::string const &, Eigen::PlainObjectBase<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > &, Eigen::PlainObjectBase<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > &, Eigen::PlainObjectBase<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > &, Eigen::PlainObjectBase<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > &);
template bool igl::readPLY<Eigen::Matrix<double, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 3, 1, -1, 3> >(std::string const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 1, -1, 3> >&);
template bool igl::readPLY<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> >(std::string const &, Eigen::PlainObjectBase<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > &, Eigen::PlainObjectBase<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> > &);
template bool igl::readPLY<Eigen::Matrix<float, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 3, 1, -1, 3> >(const std::string &, Eigen::PlainObjectBase<Eigen::Matrix<float, -1, 3, 1, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 1, -1, 3> >&);
template bool igl::readPLY<Eigen::Matrix<float, -1, 3, 1, -1, 3>, Eigen::Matrix<unsigned int, -1, 3, 1, -1, 3> >(std::string const&, Eigen::PlainObjectBase<Eigen::Matrix<float, -1, 3, 1, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<unsigned int, -1, 3, 1, -1, 3> >&);
template bool igl::readPLY<Eigen::Matrix<double, -1, -1, 1, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(std::string const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 1, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
#endif
