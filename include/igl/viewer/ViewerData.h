#ifndef IGL_VIEWER_DATA_H
#define IGL_VIEWER_DATA_H

#include <igl/igl_inline.h>
#include <Eigen/Core>

namespace igl
{

class ViewerData
#ifdef ENABLE_XML_SERIALIZATION
: public ::igl::XMLSerialization
#endif
{
public:
  ViewerData()
  #ifdef ENABLE_XML_SERIALIZATION
  : XMLSerialization("Data"), dirty(DIRTY_ALL)
  #endif
  {};

  IGL_INLINE void InitSerialization();

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
  Eigen::Matrix<char,Eigen::Dynamic,Eigen::Dynamic> texture_R;
  Eigen::Matrix<char,Eigen::Dynamic,Eigen::Dynamic> texture_G;
  Eigen::Matrix<char,Eigen::Dynamic,Eigen::Dynamic> texture_B;

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
  std::vector<std::string > labels_strings;

  // Marks dirty buffers that need to be uploaded to OpenGL
  uint32_t dirty;

  // Caches the two-norm between the min/max point of the bounding box
  float object_scale;
  /*********************************/
};


}

#ifndef IGL_STATIC_LIBRARY
#  include "ViewerData.cpp"
#endif

#endif
