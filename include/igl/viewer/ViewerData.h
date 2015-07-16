// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_VIEWER_VIEWER_DATA_H
#define IGL_VIEWER_VIEWER_DATA_H

#include <cstdint>
#include <vector>

#include <Eigen/Core>

#include <igl/igl_inline.h>

namespace igl
{
namespace viewer
{

// TODO: write documentation

class ViewerData
{
public:
  ViewerData();

  enum DirtyFlags
  {
    DIRTY_NONE           = 0x0000,
    DIRTY_POSITION       = 0x0001,
    DIRTY_UV             = 0x0002,
    DIRTY_NORMAL         = 0x0004,
    DIRTY_AMBIENT        = 0x0008,
    DIRTY_DIFFUSE        = 0x0010,
    DIRTY_SPECULAR       = 0x0020,
    DIRTY_TEXTURE        = 0x0040,
    DIRTY_FACE           = 0x0080,
    DIRTY_MESH           = 0x00FF,
    DIRTY_OVERLAY_LINES  = 0x0100,
    DIRTY_OVERLAY_POINTS = 0x0200,
    DIRTY_ALL            = 0x03FF
  };

  // Empy all fields
  IGL_INLINE void clear();

  // Change the visualization mode, invalidating the cache if necessary
  IGL_INLINE void set_face_based(bool newvalue);

  // Helpers that can draw the most common meshes
  IGL_INLINE void set_mesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);
  IGL_INLINE void set_vertices(const Eigen::MatrixXd& V);
  IGL_INLINE void set_normals(const Eigen::MatrixXd& N);

  // Set the color of the mesh
  //
  // Inputs:
  //   C  #V|#F|1 by 3 list of colors
  IGL_INLINE void set_colors(const Eigen::MatrixXd &C);
  IGL_INLINE void set_uv(const Eigen::MatrixXd& UV);
  IGL_INLINE void set_uv(const Eigen::MatrixXd& UV_V, const Eigen::MatrixXi& UV_F);
  IGL_INLINE void set_texture(
                    const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& R,
                    const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& G,
                    const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& B);

  // Sets points given a list of point vertices. In constrast to `set_points`
  // this will (purposefully) clober existing points.
  //
  // Inputs:
  //   P  #P by 3 list of vertex positions
  //   C  #P|1 by 3 color(s)
  IGL_INLINE void set_points(
    const Eigen::MatrixXd& P,  
    const Eigen::MatrixXd& C);
  IGL_INLINE void add_points(const Eigen::MatrixXd& P,  const Eigen::MatrixXd& C);
  // Sets edges given a list of edge vertices and edge indices. In constrast
  // to `add_edges` this will (purposefully) clober existing edges.
  //
  // Inputs:
  //   P  #P by 3 list of vertex positions
  //   E  #E by 2 list of edge indices into P
  //   C  #E|1 by 3 color(s)
  IGL_INLINE void set_edges (const Eigen::MatrixXd& P, const Eigen::MatrixXi& E, const Eigen::MatrixXd& C);
  IGL_INLINE void add_edges (const Eigen::MatrixXd& P1, const Eigen::MatrixXd& P2, const Eigen::MatrixXd& C);
  IGL_INLINE void add_label (const Eigen::VectorXd& P,  const std::string& str);

  // Computes the normals of the mesh
  IGL_INLINE void compute_normals();

  // Assigns uniform colors to all faces/vertices
  IGL_INLINE void uniform_colors(Eigen::Vector3d ambient, Eigen::Vector3d diffuse, Eigen::Vector3d specular);

  // Generates a default grid texture
  IGL_INLINE void grid_texture();

  Eigen::MatrixXd V; // Vertices of the current mesh (#V x 3)
  Eigen::MatrixXi F; // Faces of the mesh (#F x 3)

  // Per face attributes
  Eigen::MatrixXd F_normals; // One normal per face

  Eigen::MatrixXd F_material_ambient; // Per face ambient color
  Eigen::MatrixXd F_material_diffuse; // Per face diffuse color
  Eigen::MatrixXd F_material_specular; // Per face specular color

  // Per vertex attributes
  Eigen::MatrixXd V_normals; // One normal per vertex

  Eigen::MatrixXd V_material_ambient; // Per vertex ambient color
  Eigen::MatrixXd V_material_diffuse; // Per vertex diffuse color
  Eigen::MatrixXd V_material_specular; // Per vertex specular color

  // UV parametrization
  Eigen::MatrixXd V_uv; // UV vertices
  Eigen::MatrixXi F_uv; // optional faces for UVs

  // Texture
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> texture_R;
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> texture_G;
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> texture_B;

  // Overlays

  // Lines plotted over the scene
  // (Every row contains 9 doubles in the following format S_x, S_y, S_z, T_x, T_y, T_z, C_r, C_g, C_b),
  // with S and T the coordinates of the two vertices of the line in global coordinates, and C the color in floating point rgb format
  Eigen::MatrixXd lines;

  // Points plotted over the scene
  // (Every row contains 6 doubles in the following format P_x, P_y, P_z, C_r, C_g, C_b),
  // with P the position in global coordinates of the center of the point, and C the color in floating point rgb format
  Eigen::MatrixXd points;

  // Text labels plotted over the scene
  // Textp contains, in the i-th row, the position in global coordinates where the i-th label should be anchored
  // Texts contains in the i-th position the text of the i-th label
  Eigen::MatrixXd           labels_positions;
  std::vector<std::string>  labels_strings;

  // Marks dirty buffers that need to be uploaded to OpenGL
  uint32_t dirty;

  // Enable per-face or per-vertex properties
  bool face_based;
  /*********************************/
};

}
}

#ifndef IGL_STATIC_LIBRARY
#  include "ViewerData.cpp"
#endif

#endif
