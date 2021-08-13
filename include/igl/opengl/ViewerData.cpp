// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "ViewerData.h"
#include "ViewerCore.h"

#include "../per_face_normals.h"
#include "../material_colors.h"
#include "../per_vertex_normals.h"

// Really? Just for GL_NEAREST?
#include "gl.h"

#include <iostream>


IGL_INLINE igl::opengl::ViewerData::ViewerData()
: dirty(MeshGL::DIRTY_ALL),
  show_faces        (~unsigned(0)),
  show_lines        (~unsigned(0)),
  face_based        (false),
  double_sided      (false),
  invert_normals    (false),
  show_overlay      (~unsigned(0)),
  show_overlay_depth(~unsigned(0)),
  show_vertex_labels(0),
  show_face_labels  (0),
  show_custom_labels(0),
  show_texture      (false),
  use_matcap        (false),
  point_size(30),
  line_width(0.5f),
  label_size(1),
  line_color(0,0,0,1),
  label_color(0,0,0.04,1),
  shininess(35.0f),
  id(-1),
  is_visible        (~unsigned(0))
{
  clear();
};

IGL_INLINE void igl::opengl::ViewerData::set_face_based(bool newvalue)
{
  if (face_based != newvalue)
  {
    face_based = newvalue;
    dirty = MeshGL::DIRTY_ALL;
  }
}

// Helpers that draws the most common meshes
IGL_INLINE void igl::opengl::ViewerData::set_mesh(
    const Eigen::MatrixXd& _V, const Eigen::MatrixXi& _F)
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
    uniform_colors(
      Eigen::Vector3d(GOLD_AMBIENT[0], GOLD_AMBIENT[1], GOLD_AMBIENT[2]),
      Eigen::Vector3d(GOLD_DIFFUSE[0], GOLD_DIFFUSE[1], GOLD_DIFFUSE[2]),
      Eigen::Vector3d(GOLD_SPECULAR[0], GOLD_SPECULAR[1], GOLD_SPECULAR[2]));

    // Generates a checkerboard texture
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
      cerr << "ERROR (set_mesh): The new mesh has a different number of vertices/faces. Please clear the mesh before plotting."<<endl;
  }
  dirty |= MeshGL::DIRTY_FACE | MeshGL::DIRTY_POSITION;
}

IGL_INLINE void igl::opengl::ViewerData::set_vertices(const Eigen::MatrixXd& _V)
{
  V = _V;
  assert(F.size() == 0 || F.maxCoeff() < V.rows());
  dirty |= MeshGL::DIRTY_POSITION;
}

IGL_INLINE void igl::opengl::ViewerData::set_normals(const Eigen::MatrixXd& N)
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
    cerr << "ERROR (set_normals): Please provide a normal per face, per corner or per vertex."<<endl;
  dirty |= MeshGL::DIRTY_NORMAL;
}

IGL_INLINE void igl::opengl::ViewerData::set_visible(bool value, unsigned int core_id /*= 1*/)
{
  if (value)
    is_visible |= core_id;
  else
  is_visible &= ~core_id;
}

IGL_INLINE void igl::opengl::ViewerData::copy_options(const ViewerCore &from, const ViewerCore &to)
{
  to.set(show_overlay      , from.is_set(show_overlay)      );
  to.set(show_overlay_depth, from.is_set(show_overlay_depth));
  to.set(show_texture      , from.is_set(show_texture)      );
  to.set(use_matcap        , from.is_set(use_matcap)        );
  to.set(show_faces        , from.is_set(show_faces)        );
  to.set(show_lines        , from.is_set(show_lines)        );
}

IGL_INLINE void igl::opengl::ViewerData::set_colors(const Eigen::MatrixXd &C)
{
  using namespace std;
  using namespace Eigen;
  // This Gouraud coloring should be deprecated in favor of Phong coloring in
  // set-data
  if(C.rows()>0 && C.cols() == 1)
  {
    assert(false && "deprecated: call set_data directly instead");
    return set_data(C);
  }
  // Ambient color should be darker color
  const auto ambient = [](const MatrixXd & C)->MatrixXd
  {
    MatrixXd T = 0.1*C;
    T.col(3) = C.col(3);
    return T;
  };
  // Specular color should be a less saturated and darker color: dampened
  // highlights
  const auto specular = [](const MatrixXd & C)->MatrixXd
  {
    const double grey = 0.3;
    MatrixXd T = grey+0.1*(C.array()-grey);
    T.col(3) = C.col(3);
    return T;
  };
  if (C.rows() == 1)
  {
    for (unsigned i=0;i<V_material_diffuse.rows();++i)
    {
      if (C.cols() == 3)
        V_material_diffuse.row(i) << C.row(0),1;
      else if (C.cols() == 4)
        V_material_diffuse.row(i) << C.row(0);
    }
    V_material_ambient = ambient(V_material_diffuse);
    V_material_specular = specular(V_material_diffuse);

    for (unsigned i=0;i<F_material_diffuse.rows();++i)
    {
      if (C.cols() == 3)
        F_material_diffuse.row(i) << C.row(0),1;
      else if (C.cols() == 4)
        F_material_diffuse.row(i) << C.row(0);
    }
    F_material_ambient = ambient(F_material_diffuse);
    F_material_specular = specular(F_material_diffuse);
  }
  else if(C.rows() == V.rows() || C.rows() == F.rows())
  {
    // face based colors?
    if((C.rows()==F.rows()) && (C.rows() != V.rows() || face_based))
    {
      set_face_based(true);
      for (unsigned i=0;i<F_material_diffuse.rows();++i)
      {
        if (C.cols() == 3)
          F_material_diffuse.row(i) << C.row(i), 1;
        else if (C.cols() == 4)
          F_material_diffuse.row(i) << C.row(i);
      }
      F_material_ambient = ambient(F_material_diffuse);
      F_material_specular = specular(F_material_diffuse);
    }
    else/*(C.rows() == V.rows())*/
    {
      set_face_based(false);
      for (unsigned i=0;i<V_material_diffuse.rows();++i)
      {
        if (C.cols() == 3)
          V_material_diffuse.row(i) << C.row(i), 1;
        else if (C.cols() == 4)
          V_material_diffuse.row(i) << C.row(i);
      }
      V_material_ambient = ambient(V_material_diffuse);
      V_material_specular = specular(V_material_diffuse);
    }
  }
  else
    cerr << "ERROR (set_colors): Please provide a single color, or a color per face or per vertex."<<endl;
  dirty |= MeshGL::DIRTY_DIFFUSE | MeshGL::DIRTY_SPECULAR | MeshGL::DIRTY_AMBIENT;

}

IGL_INLINE void igl::opengl::ViewerData::set_uv(const Eigen::MatrixXd& UV)
{
  using namespace std;
  if (UV.rows() == V.rows())
  {
    set_face_based(false);
    V_uv = UV;
  }
  else
    cerr << "ERROR (set_UV): Please provide uv per vertex."<<endl;;
  dirty |= MeshGL::DIRTY_UV;
}

IGL_INLINE void igl::opengl::ViewerData::set_uv(const Eigen::MatrixXd& UV_V, const Eigen::MatrixXi& UV_F)
{
  set_face_based(true);
  V_uv = UV_V.block(0,0,UV_V.rows(),2);
  F_uv = UV_F;
  dirty |= MeshGL::DIRTY_UV;
}

IGL_INLINE void igl::opengl::ViewerData::set_texture(
  const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& R,
  const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& G,
  const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& B)
{
  texture_R = R;
  texture_G = G;
  texture_B = B;
  texture_A = Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>::Constant(R.rows(),R.cols(),255);
  dirty |= MeshGL::DIRTY_TEXTURE;
}

IGL_INLINE void igl::opengl::ViewerData::set_texture(
  const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& R,
  const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& G,
  const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& B,
  const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& A)
{
  texture_R = R;
  texture_G = G;
  texture_B = B;
  texture_A = A;
  dirty |= MeshGL::DIRTY_TEXTURE;
}

IGL_INLINE void igl::opengl::ViewerData::set_data(
  const Eigen::VectorXd & D,
  double caxis_min,
  double caxis_max,
  igl::ColorMapType cmap,
  int num_steps)
{
  if(!show_texture)
  {
    Eigen::MatrixXd CM;
    igl::colormap(cmap,Eigen::VectorXd::LinSpaced(num_steps,0,1).eval(),0,1,CM);
    set_colormap(CM);
  }
  set_uv(((D.array()-caxis_min)/(caxis_max-caxis_min)).replicate(1,2));
}

IGL_INLINE void igl::opengl::ViewerData::set_data(const Eigen::VectorXd & D, igl::ColorMapType cmap, int num_steps)
{
  const double caxis_min = D.minCoeff();
  const double caxis_max = D.maxCoeff();
  return set_data(D,caxis_min,caxis_max,cmap,num_steps);
}

IGL_INLINE void igl::opengl::ViewerData::set_colormap(const Eigen::MatrixXd & CM)
{
  assert(CM.cols() == 3 && "colormap CM should have 3 columns");
  // Convert to R,G,B textures
  const Eigen::Matrix<unsigned char,Eigen::Dynamic, Eigen::Dynamic> R =
    (CM.col(0)*255.0).cast<unsigned char>();
  const Eigen::Matrix<unsigned char,Eigen::Dynamic, Eigen::Dynamic> G =
    (CM.col(1)*255.0).cast<unsigned char>();
  const Eigen::Matrix<unsigned char,Eigen::Dynamic, Eigen::Dynamic> B =
    (CM.col(2)*255.0).cast<unsigned char>();
  set_colors(Eigen::RowVector3d(1,1,1));
  set_texture(R,G,B);
  show_texture = ~unsigned(0);
  meshgl.tex_filter = GL_NEAREST;
  meshgl.tex_wrap = GL_CLAMP_TO_EDGE;
}

IGL_INLINE void igl::opengl::ViewerData::set_points(
  const Eigen::MatrixXd& P,
  const Eigen::MatrixXd& C)
{
  // clear existing points
  points.resize(0,0);
  add_points(P,C);
}

IGL_INLINE void igl::opengl::ViewerData::add_points(const Eigen::MatrixXd& P,  const Eigen::MatrixXd& C)
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

  dirty |= MeshGL::DIRTY_OVERLAY_POINTS;
}

IGL_INLINE void igl::opengl::ViewerData::clear_points()
{
  points.resize(0, 6);
}

IGL_INLINE void igl::opengl::ViewerData::set_edges(
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
  dirty |= MeshGL::DIRTY_OVERLAY_LINES;
}

IGL_INLINE void igl::opengl::ViewerData::set_edges_from_vector_field(
  const Eigen::MatrixXd& P, 
  const Eigen::MatrixXd& V, 
  const Eigen::MatrixXd& C)
{
  assert(P.rows() == V.rows());
  Eigen::MatrixXi E(P.rows(),2);
  const Eigen::MatrixXd PV = 
    (Eigen::MatrixXd(P.rows()+V.rows(),3)<<P,P+V).finished();
  for(int i = 0;i<P.rows();i++)
  {
    E(i,0) = i;
    E(i,1) = i+P.rows();
  }
  const Eigen::MatrixXd CC = C.replicate<2,1>();
  set_edges(PV,E, C.rows() == 1?C:C.replicate<2,1>());
}

IGL_INLINE void igl::opengl::ViewerData::add_edges(const Eigen::MatrixXd& P1, const Eigen::MatrixXd& P2, const Eigen::MatrixXd& C)
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

  dirty |= MeshGL::DIRTY_OVERLAY_LINES;
}

IGL_INLINE void igl::opengl::ViewerData::clear_edges()
{
  lines.resize(0, 9);
}

IGL_INLINE void igl::opengl::ViewerData::add_label(const Eigen::VectorXd& P,  const std::string& str)
{
  Eigen::RowVectorXd P_temp;

  // If P only has two columns, pad with a column of zeros
  if (P.size() == 2)
  {
    P_temp = Eigen::RowVectorXd::Zero(3);
    P_temp << P.transpose(), 0;
  }
  else
    P_temp = P;

  int lastid = labels_positions.rows();
  labels_positions.conservativeResize(lastid+1, 3);
  labels_positions.row(lastid) = P_temp;
  labels_strings.push_back(str);

  dirty |= MeshGL::DIRTY_CUSTOM_LABELS;
}

IGL_INLINE void igl::opengl::ViewerData::set_labels(const Eigen::MatrixXd& P, const std::vector<std::string>& str)
{
  assert(P.rows() == str.size() && "position # and label # do not match!");
  assert(P.cols() == 3 && "dimension of label positions incorrect!");
  labels_positions = P;
  labels_strings = str;

  dirty |= MeshGL::DIRTY_CUSTOM_LABELS;
}

IGL_INLINE void igl::opengl::ViewerData::clear_labels()
{
  labels_positions.resize(0,3);
  labels_strings.clear();
}

IGL_INLINE void igl::opengl::ViewerData::clear()
{
  V                       = Eigen::MatrixXd (0,3);
  F                       = Eigen::MatrixXi (0,3);

  F_material_ambient      = Eigen::MatrixXd (0,4);
  F_material_diffuse      = Eigen::MatrixXd (0,4);
  F_material_specular     = Eigen::MatrixXd (0,4);

  V_material_ambient      = Eigen::MatrixXd (0,4);
  V_material_diffuse      = Eigen::MatrixXd (0,4);
  V_material_specular     = Eigen::MatrixXd (0,4);

  F_normals               = Eigen::MatrixXd (0,3);
  V_normals               = Eigen::MatrixXd (0,3);

  V_uv                    = Eigen::MatrixXd (0,2);
  F_uv                    = Eigen::MatrixXi (0,3);

  lines                   = Eigen::MatrixXd (0,9);
  points                  = Eigen::MatrixXd (0,6);

  vertex_labels_positions = Eigen::MatrixXd (0,3);
  face_labels_positions   = Eigen::MatrixXd (0,3);
  labels_positions        = Eigen::MatrixXd (0,3);
  vertex_labels_strings.clear();
  face_labels_strings.clear();
  labels_strings.clear();

  face_based = false;
  double_sided = false;
  invert_normals = false;
  show_texture = false;
  use_matcap = false;
}

IGL_INLINE void igl::opengl::ViewerData::compute_normals()
{
  if(V.cols() == 2)
  {
    F_normals = Eigen::RowVector3d(0,0,1).replicate(F.rows(),1);
    V_normals = Eigen::RowVector3d(0,0,1).replicate(V.rows(),1);
  }else
  {
    assert(V.cols() == 3);
    igl::per_face_normals(V, F, F_normals);
    igl::per_vertex_normals(V, F, F_normals, V_normals);
  }
  dirty |= MeshGL::DIRTY_NORMAL;
}

IGL_INLINE void igl::opengl::ViewerData::uniform_colors(
  const Eigen::Vector3d& ambient,
  const Eigen::Vector3d& diffuse,
  const Eigen::Vector3d& specular)
{
  Eigen::Vector4d ambient4;
  Eigen::Vector4d diffuse4;
  Eigen::Vector4d specular4;

  ambient4 << ambient, 1;
  diffuse4 << diffuse, 1;
  specular4 << specular, 1;

  uniform_colors(ambient4,diffuse4,specular4);
}

IGL_INLINE void igl::opengl::ViewerData::uniform_colors(
  const Eigen::Vector4d& ambient,
  const Eigen::Vector4d& diffuse,
  const Eigen::Vector4d& specular)
{
  V_material_ambient.resize(V.rows(),4);
  V_material_diffuse.resize(V.rows(),4);
  V_material_specular.resize(V.rows(),4);

  for (unsigned i=0; i<V.rows();++i)
  {
    V_material_ambient.row(i) = ambient;
    V_material_diffuse.row(i) = diffuse;
    V_material_specular.row(i) = specular;
  }

  F_material_ambient.resize(F.rows(),4);
  F_material_diffuse.resize(F.rows(),4);
  F_material_specular.resize(F.rows(),4);

  for (unsigned i=0; i<F.rows();++i)
  {
    F_material_ambient.row(i) = ambient;
    F_material_diffuse.row(i) = diffuse;
    F_material_specular.row(i) = specular;
  }
  dirty |= MeshGL::DIRTY_SPECULAR | MeshGL::DIRTY_DIFFUSE | MeshGL::DIRTY_AMBIENT;
}

IGL_INLINE void igl::opengl::ViewerData::normal_matcap()
{
  const int size = 512;
  texture_R.resize(size, size);
  texture_G.resize(size, size);
  texture_B.resize(size, size);
  const Eigen::Vector3d navy(0.3,0.3,0.5);
  static const auto clamp = [](double t){ return std::max(std::min(t,1.0),0.0);};
  for(int i = 0;i<size;i++)
  {
    const double x = (double(i)/double(size-1)*2.-1.);
    for(int j = 0;j<size;j++)
    {
      const double y = (double(j)/double(size-1)*2.-1.);
      const double z = sqrt(1.0-std::min(x*x+y*y,1.0));
      Eigen::Vector3d C = Eigen::Vector3d(x*0.5+0.5,y*0.5+0.5,z);
      texture_R(i,j) = clamp(C(0))*255;
      texture_G(i,j) = clamp(C(1))*255;
      texture_B(i,j) = clamp(C(2))*255;
    }
  }
  texture_A.setConstant(texture_R.rows(),texture_R.cols(),255);
  dirty |= MeshGL::DIRTY_TEXTURE;
}

IGL_INLINE void igl::opengl::ViewerData::grid_texture()
{
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
  texture_A = Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>::Constant(texture_R.rows(),texture_R.cols(),255);
  dirty |= MeshGL::DIRTY_TEXTURE;
}

// Populate VBOs of a particular label stype (Vert, Face, Custom)
IGL_INLINE void igl::opengl::ViewerData::update_labels(
  igl::opengl::MeshGL& meshgl,
  igl::opengl::MeshGL::TextGL& GL_labels,
  const Eigen::MatrixXd& positions,
  const std::vector<std::string>& strings
){
  if (positions.rows()>0)
  {
    int numCharsToRender = 0;
    for(size_t p=0; p<positions.rows(); p++)
    {
      numCharsToRender += strings.at(p).length();
    }
    GL_labels.label_pos_vbo.resize(numCharsToRender, 3);
    GL_labels.label_char_vbo.resize(numCharsToRender, 1);
    GL_labels.label_offset_vbo.resize(numCharsToRender, 1);
    GL_labels.label_indices_vbo.resize(numCharsToRender, 1);
    int idx=0;
    assert(strings.size() == positions.rows());
    for(size_t s=0; s<strings.size(); s++)
    {
      const auto & label = strings.at(s);
      for(size_t c=0; c<label.length(); c++)
      {
        GL_labels.label_pos_vbo.row(idx) = positions.row(s).cast<float>();
        GL_labels.label_char_vbo(idx) = (float)(label.at(c));
        GL_labels.label_offset_vbo(idx) = c;
        GL_labels.label_indices_vbo(idx) = idx;
        idx++;
      }
    }
  }
}

IGL_INLINE void igl::opengl::ViewerData::updateGL(
  const igl::opengl::ViewerData& data,
  const bool invert_normals,
  igl::opengl::MeshGL& meshgl
  )
{
  if (!meshgl.is_initialized)
  {
    meshgl.init();
  }

  bool per_corner_uv = (data.F_uv.rows() == data.F.rows());
  bool per_corner_normals = (data.F_normals.rows() == 3 * data.F.rows());

  meshgl.dirty |= data.dirty;

  // Input:
  //   X  #F by dim quantity
  // Output:
  //   X_vbo  #F*3 by dim scattering per corner
  const auto per_face = [&data](
      const Eigen::MatrixXd & X,
      MeshGL::RowMatrixXf & X_vbo)
  {
    assert(X.cols() == 4);
    X_vbo.resize(data.F.rows()*3,4);
    for (unsigned i=0; i<data.F.rows();++i)
      for (unsigned j=0;j<3;++j)
        X_vbo.row(i*3+j) = X.row(i).cast<float>();
  };

  // Input:
  //   X  #V by dim quantity
  // Output:
  //   X_vbo  #F*3 by dim scattering per corner
  const auto per_corner = [&data](
      const Eigen::MatrixXd & X,
      MeshGL::RowMatrixXf & X_vbo)
  {
    X_vbo.resize(data.F.rows()*3,X.cols());
    for (unsigned i=0; i<data.F.rows();++i)
      for (unsigned j=0;j<3;++j)
        X_vbo.row(i*3+j) = X.row(data.F(i,j)).cast<float>();
  };

  if (!data.face_based)
  {
    if (!(per_corner_uv || per_corner_normals))
    {
      // Vertex positions
      if (meshgl.dirty & MeshGL::DIRTY_POSITION)
        meshgl.V_vbo = data.V.cast<float>();

      // Vertex normals
      if (meshgl.dirty & MeshGL::DIRTY_NORMAL)
      {
        meshgl.V_normals_vbo = data.V_normals.cast<float>();
        if (invert_normals)
          meshgl.V_normals_vbo = -meshgl.V_normals_vbo;
      }

      // Per-vertex material settings
      if (meshgl.dirty & MeshGL::DIRTY_AMBIENT)
        meshgl.V_ambient_vbo = data.V_material_ambient.cast<float>();
      if (meshgl.dirty & MeshGL::DIRTY_DIFFUSE)
        meshgl.V_diffuse_vbo = data.V_material_diffuse.cast<float>();
      if (meshgl.dirty & MeshGL::DIRTY_SPECULAR)
        meshgl.V_specular_vbo = data.V_material_specular.cast<float>();

      // Face indices
      if (meshgl.dirty & MeshGL::DIRTY_FACE)
        meshgl.F_vbo = data.F.cast<unsigned>();

      // Texture coordinates
      if (meshgl.dirty & MeshGL::DIRTY_UV)
      {
        meshgl.V_uv_vbo = data.V_uv.cast<float>();
      }
    }
    else
    {

      // Per vertex properties with per corner UVs
      if (meshgl.dirty & MeshGL::DIRTY_POSITION)
      {
        per_corner(data.V,meshgl.V_vbo);
      }

      if (meshgl.dirty & MeshGL::DIRTY_AMBIENT)
      {
        meshgl.V_ambient_vbo.resize(data.F.rows()*3,4);
        per_corner(data.V_material_ambient,meshgl.V_ambient_vbo);
      }
      if (meshgl.dirty & MeshGL::DIRTY_DIFFUSE)
      {
        meshgl.V_diffuse_vbo.resize(data.F.rows()*3,4);
        per_corner(data.V_material_diffuse,meshgl.V_diffuse_vbo);
      }
      if (meshgl.dirty & MeshGL::DIRTY_SPECULAR)
      {
        meshgl.V_specular_vbo.resize(data.F.rows()*3,4);
        per_corner(data.V_material_specular,meshgl.V_specular_vbo);
      }

      if (meshgl.dirty & MeshGL::DIRTY_NORMAL)
      {
        meshgl.V_normals_vbo.resize(data.F.rows()*3,3);
        per_corner(data.V_normals,meshgl.V_normals_vbo);

        if (invert_normals)
          meshgl.V_normals_vbo = -meshgl.V_normals_vbo;
      }

      if (meshgl.dirty & MeshGL::DIRTY_FACE)
      {
        meshgl.F_vbo.resize(data.F.rows(),3);
        for (unsigned i=0; i<data.F.rows();++i)
          meshgl.F_vbo.row(i) << i*3+0, i*3+1, i*3+2;
      }

      if ( (meshgl.dirty & MeshGL::DIRTY_UV) && data.V_uv.rows()>0)
      {
        meshgl.V_uv_vbo.resize(data.F.rows()*3,2);
        for (unsigned i=0; i<data.F.rows();++i)
          for (unsigned j=0;j<3;++j)
            meshgl.V_uv_vbo.row(i*3+j) =
              data.V_uv.row(per_corner_uv ?
                data.F_uv(i,j) : data.F(i,j)).cast<float>();
      }
    }
  } else
  {
    if (meshgl.dirty & MeshGL::DIRTY_POSITION)
    {
      per_corner(data.V,meshgl.V_vbo);
    }
    if (meshgl.dirty & MeshGL::DIRTY_AMBIENT)
    {
      per_face(data.F_material_ambient,meshgl.V_ambient_vbo);
    }
    if (meshgl.dirty & MeshGL::DIRTY_DIFFUSE)
    {
      per_face(data.F_material_diffuse,meshgl.V_diffuse_vbo);
    }
    if (meshgl.dirty & MeshGL::DIRTY_SPECULAR)
    {
      per_face(data.F_material_specular,meshgl.V_specular_vbo);
    }

    if (meshgl.dirty & MeshGL::DIRTY_NORMAL)
    {
      meshgl.V_normals_vbo.resize(data.F.rows()*3,3);
      for (unsigned i=0; i<data.F.rows();++i)
        for (unsigned j=0;j<3;++j)
          meshgl.V_normals_vbo.row(i*3+j) =
             per_corner_normals ?
               data.F_normals.row(i*3+j).cast<float>() :
               data.F_normals.row(i).cast<float>();

      if (invert_normals)
        meshgl.V_normals_vbo = -meshgl.V_normals_vbo;
    }

    if (meshgl.dirty & MeshGL::DIRTY_FACE)
    {
      meshgl.F_vbo.resize(data.F.rows(),3);
      for (unsigned i=0; i<data.F.rows();++i)
        meshgl.F_vbo.row(i) << i*3+0, i*3+1, i*3+2;
    }

    if( (meshgl.dirty & MeshGL::DIRTY_UV) && data.V_uv.rows()>0)
    {
        meshgl.V_uv_vbo.resize(data.F.rows()*3,2);
        for (unsigned i=0; i<data.F.rows();++i)
          for (unsigned j=0;j<3;++j)
            meshgl.V_uv_vbo.row(i*3+j) = data.V_uv.row(per_corner_uv ? data.F_uv(i,j) : data.F(i,j)).cast<float>();
    }
  }

  if (meshgl.dirty & MeshGL::DIRTY_TEXTURE)
  {
    meshgl.tex_u = data.texture_R.rows();
    meshgl.tex_v = data.texture_R.cols();
    meshgl.tex.resize(data.texture_R.size()*4);
    for (unsigned i=0;i<data.texture_R.size();++i)
    {
      meshgl.tex(i*4+0) = data.texture_R(i);
      meshgl.tex(i*4+1) = data.texture_G(i);
      meshgl.tex(i*4+2) = data.texture_B(i);
      meshgl.tex(i*4+3) = data.texture_A(i);
    }
  }

  if (meshgl.dirty & MeshGL::DIRTY_OVERLAY_LINES)
  {
    meshgl.lines_V_vbo.resize(data.lines.rows()*2,3);
    meshgl.lines_V_colors_vbo.resize(data.lines.rows()*2,3);
    meshgl.lines_F_vbo.resize(data.lines.rows()*2,1);
    for (unsigned i=0; i<data.lines.rows();++i)
    {
      meshgl.lines_V_vbo.row(2*i+0) = data.lines.block<1, 3>(i, 0).cast<float>();
      meshgl.lines_V_vbo.row(2*i+1) = data.lines.block<1, 3>(i, 3).cast<float>();
      meshgl.lines_V_colors_vbo.row(2*i+0) = data.lines.block<1, 3>(i, 6).cast<float>();
      meshgl.lines_V_colors_vbo.row(2*i+1) = data.lines.block<1, 3>(i, 6).cast<float>();
      meshgl.lines_F_vbo(2*i+0) = 2*i+0;
      meshgl.lines_F_vbo(2*i+1) = 2*i+1;
    }
  }

  if (meshgl.dirty & MeshGL::DIRTY_OVERLAY_POINTS)
  {
    meshgl.points_V_vbo.resize(data.points.rows(),3);
    meshgl.points_V_colors_vbo.resize(data.points.rows(),3);
    meshgl.points_F_vbo.resize(data.points.rows(),1);
    for (unsigned i=0; i<data.points.rows();++i)
    {
      meshgl.points_V_vbo.row(i) = data.points.block<1, 3>(i, 0).cast<float>();
      meshgl.points_V_colors_vbo.row(i) = data.points.block<1, 3>(i, 3).cast<float>();
      meshgl.points_F_vbo(i) = i;
    }
  }

  if (meshgl.dirty & MeshGL::DIRTY_FACE_LABELS)
  {
    if(face_labels_positions.rows()==0)
    {
      face_labels_positions.conservativeResize(F.rows(), 3);
      Eigen::MatrixXd faceNormals = F_normals.normalized();
      for (int f=0; f<F.rows();++f)
      {
        std::string faceName = std::to_string(f);
        face_labels_positions.row(f) = V.row(F.row(f)(0));
        face_labels_positions.row(f) += V.row(F.row(f)(1));
        face_labels_positions.row(f) += V.row(F.row(f)(2));
        face_labels_positions.row(f) /= 3.;
        face_labels_positions.row(f) = (faceNormals*0.05).row(f) + face_labels_positions.row(f);
        face_labels_strings.push_back(faceName);
      }
    }
    update_labels(
      meshgl,
      meshgl.face_labels,
      face_labels_positions,
      face_labels_strings
    );
  }

  if (meshgl.dirty & MeshGL::DIRTY_VERTEX_LABELS)
  {
    if(vertex_labels_positions.rows()==0)
    {
      vertex_labels_positions.conservativeResize(V.rows(), 3);
      Eigen::MatrixXd normalized = V_normals.normalized();
      for (int v=0; v<V.rows();++v)
      {
        std::string vertName = std::to_string(v);
        vertex_labels_positions.row(v) = (normalized*0.1).row(v) + V.row(v);
        vertex_labels_strings.push_back(vertName);
      }
    }
    update_labels(
      meshgl,
      meshgl.vertex_labels,
      vertex_labels_positions,
      vertex_labels_strings
    );
  }

  if (meshgl.dirty & MeshGL::DIRTY_CUSTOM_LABELS)
  {
    update_labels(
      meshgl,
      meshgl.custom_labels,
      labels_positions,
      labels_strings
    );
  }
}
