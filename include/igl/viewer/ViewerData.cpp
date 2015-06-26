// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "ViewerData.h"

#include <iostream>

#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>

#ifdef ENABLE_SERIALIZATION
#include <igl/serialize.h>
namespace igl {
  namespace serialization {

    IGL_INLINE void serialization(bool s,igl::viewer::ViewerData& obj,std::vector<char>& buffer)
    {
      SERIALIZE_MEMBER(V);
      SERIALIZE_MEMBER(F);

      SERIALIZE_MEMBER(F_normals);
      SERIALIZE_MEMBER(F_material_ambient);
      SERIALIZE_MEMBER(F_material_diffuse);
      SERIALIZE_MEMBER(F_material_specular);

      SERIALIZE_MEMBER(V_normals);
      SERIALIZE_MEMBER(V_material_ambient);
      SERIALIZE_MEMBER(V_material_diffuse);
      SERIALIZE_MEMBER(V_material_specular);

      SERIALIZE_MEMBER(V_uv);
      SERIALIZE_MEMBER(F_uv);

      SERIALIZE_MEMBER(texture_R);
      SERIALIZE_MEMBER(texture_G);
      SERIALIZE_MEMBER(texture_B);

      SERIALIZE_MEMBER(lines);
      SERIALIZE_MEMBER(points);

      SERIALIZE_MEMBER(labels_positions);
      SERIALIZE_MEMBER(labels_strings);

      SERIALIZE_MEMBER(face_based);
    }

    IGL_INLINE void serialize(const igl::viewer::ViewerData& obj,std::vector<char>& buffer)
    {
      serialization(true,const_cast<igl::viewer::ViewerData&>(obj),buffer);
    }

    IGL_INLINE void deserialize(igl::viewer::ViewerData& obj,const std::vector<char>& buffer)
    {
      serialization(false,obj,const_cast<std::vector<char>&>(buffer));
      obj.dirty = igl::viewer::ViewerData::DIRTY_ALL;
    }
  }
}
#endif

IGL_INLINE igl::viewer::ViewerData::ViewerData()
: dirty(DIRTY_ALL)
{
  clear();
};

IGL_INLINE void igl::viewer::ViewerData::set_face_based(bool newvalue)
{
  if (face_based != newvalue)
  {
    face_based = newvalue;
    dirty = DIRTY_ALL;
  }
}

// Helpers that draws the most common meshes
IGL_INLINE void igl::viewer::ViewerData::set_mesh(const Eigen::MatrixXd& _V, const Eigen::MatrixXi& _F)
{
  using namespace std;

  Eigen::MatrixXd V_temp;

  // If V only has two columns, pad with a column of zeros
  if (_V.cols() == 2)
  {
    V_temp = Eigen::MatrixXd::Zero(_V.rows(),3);
    V_temp.block(0,0,_V.rows(),2) = _V;
  }
  else
    V_temp = _V;

  if (V.rows() == 0 && F.rows() == 0)
  {
    V = V_temp;
    F = _F;

    compute_normals();
    uniform_colors(Eigen::Vector3d(51.0/255.0,43.0/255.0,33.3/255.0),
                   Eigen::Vector3d(255.0/255.0,228.0/255.0,58.0/255.0),
                   Eigen::Vector3d(255.0/255.0,235.0/255.0,80.0/255.0));

    grid_texture();
  }
  else
  {
    if (_V.rows() == V.rows() && _F.rows() == F.rows())
    {
      V = V_temp;
      F = _F;
    }
    else
      cerr << "ERROR (set_mesh): The new mesh has a different number of vertices/faces. Please clear the mesh before plotting.";
  }
  dirty |= DIRTY_FACE | DIRTY_POSITION;
}

IGL_INLINE void igl::viewer::ViewerData::set_vertices(const Eigen::MatrixXd& _V)
{
  V = _V;
  assert(F.size() == 0 || F.maxCoeff() < V.rows());
  dirty |= DIRTY_POSITION;
}

IGL_INLINE void igl::viewer::ViewerData::set_normals(const Eigen::MatrixXd& N)
{
  using namespace std;
  if (N.rows() == V.rows())
  {
    set_face_based(false);
    V_normals = N;
  }
  else if (N.rows() == F.rows() || N.rows() == F.rows()*3)
  {
    set_face_based(true);
    F_normals = N;
  }
  else
    cerr << "ERROR (set_normals): Please provide a normal per face, per corner or per vertex.";
  dirty |= DIRTY_NORMAL;
}

IGL_INLINE void igl::viewer::ViewerData::set_colors(const Eigen::MatrixXd &C)
{
  using namespace std;
  using namespace Eigen;
  // Ambient color should be darker color
  const auto ambient = [](const MatrixXd & C)->MatrixXd
  {
    return 0.1*C;
  };
  // Specular color should be a less saturated and darker color: dampened
  // highlights
  const auto specular = [](const MatrixXd & C)->MatrixXd
  {
    const double grey = 0.3;
    return grey+0.1*(C.array()-grey);
  };
  if (C.rows() == 1)
  {
    for (unsigned i=0;i<V_material_diffuse.rows();++i)
    {
      V_material_diffuse.row(i) = C.row(0);
    }
    V_material_ambient = ambient(V_material_diffuse);
    V_material_specular = specular(V_material_diffuse);

    for (unsigned i=0;i<F_material_diffuse.rows();++i)
    {
      F_material_diffuse.row(i) = C.row(0);
    }
    F_material_ambient = ambient(F_material_diffuse);
    F_material_specular = specular(F_material_diffuse);
  }
  else if (C.rows() == V.rows())
  {
    set_face_based(false);
    V_material_diffuse = C;
    V_material_ambient = ambient(V_material_diffuse);
    V_material_specular = specular(V_material_diffuse);
  }
  else if (C.rows() == F.rows())
  {
    set_face_based(true);
    F_material_diffuse = C;
    F_material_ambient = ambient(F_material_diffuse);
    F_material_specular = specular(F_material_diffuse);
  }
  else
    cerr << "ERROR (set_colors): Please provide a single color, or a color per face or per vertex.";
  dirty |= DIRTY_DIFFUSE;

}

IGL_INLINE void igl::viewer::ViewerData::set_uv(const Eigen::MatrixXd& UV)
{
  using namespace std;
  if (UV.rows() == V.rows())
  {
    set_face_based(false);
    V_uv = UV;
  }
  else
    cerr << "ERROR (set_UV): Please provide uv per vertex.";
  dirty |= DIRTY_UV;
}

IGL_INLINE void igl::viewer::ViewerData::set_uv(const Eigen::MatrixXd& UV_V, const Eigen::MatrixXi& UV_F)
{
  set_face_based(true);
  V_uv = UV_V.block(0,0,UV_V.rows(),2);
  F_uv = UV_F;
  dirty |= DIRTY_UV;
}


IGL_INLINE void igl::viewer::ViewerData::set_texture(
  const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& R,
  const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& G,
  const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& B)
{
  texture_R = R;
  texture_G = G;
  texture_B = B;
  dirty |= DIRTY_TEXTURE;
}

IGL_INLINE void igl::viewer::ViewerData::set_points(
  const Eigen::MatrixXd& P,
  const Eigen::MatrixXd& C)
{
  // clear existing points
  points.resize(0,0);
  add_points(P,C);
}

IGL_INLINE void igl::viewer::ViewerData::add_points(const Eigen::MatrixXd& P,  const Eigen::MatrixXd& C)
{
  Eigen::MatrixXd P_temp;

  // If P only has two columns, pad with a column of zeros
  if (P.cols() == 2)
  {
    P_temp = Eigen::MatrixXd::Zero(P.rows(),3);
    P_temp.block(0,0,P.rows(),2) = P;
  }
  else
    P_temp = P;

  int lastid = points.rows();
  points.conservativeResize(points.rows() + P_temp.rows(),6);
  for (unsigned i=0; i<P_temp.rows(); ++i)
    points.row(lastid+i) << P_temp.row(i), i<C.rows() ? C.row(i) : C.row(C.rows()-1);

  dirty |= DIRTY_OVERLAY_POINTS;
}

IGL_INLINE void igl::viewer::ViewerData::set_edges(
  const Eigen::MatrixXd& P,
  const Eigen::MatrixXi& E,
  const Eigen::MatrixXd& C)
{
  using namespace Eigen;
  lines.resize(E.rows(),9);
  assert(C.cols() == 3);
  for(int e = 0;e<E.rows();e++)
  {
    RowVector3d color;
    if(C.size() == 3)
    {
      color<<C;
    }else if(C.rows() == E.rows())
    {
      color<<C.row(e);
    }
    lines.row(e)<< P.row(E(e,0)), P.row(E(e,1)), color;
  }
  dirty |= DIRTY_OVERLAY_LINES;
}

IGL_INLINE void igl::viewer::ViewerData::add_edges(const Eigen::MatrixXd& P1, const Eigen::MatrixXd& P2, const Eigen::MatrixXd& C)
{
  Eigen::MatrixXd P1_temp,P2_temp;

  // If P1 only has two columns, pad with a column of zeros
  if (P1.cols() == 2)
  {
    P1_temp = Eigen::MatrixXd::Zero(P1.rows(),3);
    P1_temp.block(0,0,P1.rows(),2) = P1;
    P2_temp = Eigen::MatrixXd::Zero(P2.rows(),3);
    P2_temp.block(0,0,P2.rows(),2) = P2;
  }
  else
  {
    P1_temp = P1;
    P2_temp = P2;
  }

  int lastid = lines.rows();
  lines.conservativeResize(lines.rows() + P1_temp.rows(),9);
  for (unsigned i=0; i<P1_temp.rows(); ++i)
    lines.row(lastid+i) << P1_temp.row(i), P2_temp.row(i), i<C.rows() ? C.row(i) : C.row(C.rows()-1);

  dirty |= DIRTY_OVERLAY_LINES;
}

IGL_INLINE void igl::viewer::ViewerData::add_label(const Eigen::VectorXd& P,  const std::string& str)
{
  Eigen::RowVectorXd P_temp;

  // If P only has two columns, pad with a column of zeros
  if (P.size() == 2)
  {
    P_temp = Eigen::RowVectorXd::Zero(3);
    P_temp << P, 0;
  }
  else
    P_temp = P;

  int lastid = labels_positions.rows();
  labels_positions.conservativeResize(lastid+1, 3);
  labels_positions.row(lastid) = P_temp;
  labels_strings.push_back(str);
}

IGL_INLINE void igl::viewer::ViewerData::clear()
{
  V                       = Eigen::MatrixXd (0,3);
  F                       = Eigen::MatrixXi (0,3);

  F_material_ambient      = Eigen::MatrixXd (0,3);
  F_material_diffuse      = Eigen::MatrixXd (0,3);
  F_material_specular     = Eigen::MatrixXd (0,3);

  V_material_ambient      = Eigen::MatrixXd (0,3);
  V_material_diffuse      = Eigen::MatrixXd (0,3);
  V_material_specular     = Eigen::MatrixXd (0,3);

  F_normals               = Eigen::MatrixXd (0,3);
  V_normals               = Eigen::MatrixXd (0,3);

  V_uv                    = Eigen::MatrixXd (0,2);
  F_uv                    = Eigen::MatrixXi (0,3);

  lines                   = Eigen::MatrixXd (0,9);
  points                  = Eigen::MatrixXd (0,6);
  labels_positions        = Eigen::MatrixXd (0,3);
  labels_strings.clear();

  face_based = false;
}

IGL_INLINE void igl::viewer::ViewerData::compute_normals()
{
  igl::per_face_normals(V, F, F_normals);
  igl::per_vertex_normals(V, F, F_normals, V_normals);
  dirty |= DIRTY_NORMAL;
}

IGL_INLINE void igl::viewer::ViewerData::uniform_colors(Eigen::Vector3d ambient, Eigen::Vector3d diffuse, Eigen::Vector3d specular)
{
  V_material_ambient.resize(V.rows(),3);
  V_material_diffuse.resize(V.rows(),3);
  V_material_specular.resize(V.rows(),3);

  for (unsigned i=0; i<V.rows();++i)
  {
    V_material_ambient.row(i) = ambient;
    V_material_diffuse.row(i) = diffuse;
    V_material_specular.row(i) = specular;
  }

  F_material_ambient.resize(F.rows(),3);
  F_material_diffuse.resize(F.rows(),3);
  F_material_specular.resize(F.rows(),3);

  for (unsigned i=0; i<F.rows();++i)
  {
    F_material_ambient.row(i) = ambient;
    F_material_diffuse.row(i) = diffuse;
    F_material_specular.row(i) = specular;
  }
  dirty |= DIRTY_SPECULAR | DIRTY_DIFFUSE | DIRTY_AMBIENT;
}

IGL_INLINE void igl::viewer::ViewerData::grid_texture()
{
  if (V_uv.rows() == 0)
  {
    V_uv = V.block(0, 0, V.rows(), 2);
    V_uv.col(0) = V_uv.col(0).array() - V_uv.col(0).minCoeff();
    V_uv.col(0) = V_uv.col(0).array() / V_uv.col(0).maxCoeff();
    V_uv.col(1) = V_uv.col(1).array() - V_uv.col(1).minCoeff();
    V_uv.col(1) = V_uv.col(1).array() / V_uv.col(1).maxCoeff();
    V_uv = V_uv.array() * 10;
    dirty |= DIRTY_TEXTURE;
  }

  unsigned size = 128;
  unsigned size2 = size/2;
  texture_R.resize(size, size);
  for (unsigned i=0; i<size; ++i)
  {
    for (unsigned j=0; j<size; ++j)
    {
      texture_R(i,j) = 0;
      if ((i<size2 && j<size2) || (i>=size2 && j>=size2))
        texture_R(i,j) = 255;
    }
  }

  texture_G = texture_R;
  texture_B = texture_R;
  dirty |= DIRTY_TEXTURE;
}
