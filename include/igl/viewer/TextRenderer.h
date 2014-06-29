/* This class extends the font rendering code in AntTweakBar
   so that it can be used to render text at arbitrary 3D positions */

#ifndef IGL_TEXT_RENDERER_H
#define IGL_TEXT_RENDERER_H

#include <igl/igl_inline.h>


namespace igl
{

class TextRenderer : public CTwGraphOpenGLCore
{
public:
  IGL_INLINE TextRenderer();

  IGL_INLINE virtual int Init();
  IGL_INLINE virtual int Shut();

  IGL_INLINE void BeginDraw(const Eigen::Matrix4f &view, const Eigen::Matrix4f &proj,
    const Eigen::Vector4f &_viewport, float _object_scale);

  IGL_INLINE void EndDraw();

  IGL_INLINE void DrawText(Eigen::Vector3d pos, Eigen::Vector3d normal, const std::string &text);

protected:
  igl::Viewer::OpenGL_shader m_shader;
  std::map<std::string, void *> m_textObjects;
  GLuint m_shaderHandleBackup;
  GLuint m_TriTexUniLocationDepth;
  Eigen::Matrix4f view_matrix, proj_matrix;
  Eigen::Vector4f viewport;
  float object_scale;
};

}

#ifndef IGL_STATIC_LIBRARY
#  include "TextRenderer.cpp"
#endif

#endif
