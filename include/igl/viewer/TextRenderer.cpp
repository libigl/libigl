// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Wenzel Jacob <wenzel@inf.ethz.ch>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "TextRenderer.h"
#include <igl/project.h>

  IGL_INLINE igl::TextRenderer::TextRenderer() : m_shaderHandleBackup(0) { }

  IGL_INLINE int igl::TextRenderer::Init()
  {
    int retval = CTwGraphOpenGLCore::Init();
    if (retval == 1)
    {
      std::string vertexShader =
          "#version 150\n"
          "uniform vec2 offset;"
          "uniform vec2 wndSize;"
          "uniform vec4 color;"
          "uniform float depth;"
          "in vec2 vertex;"
          "in vec2 uv;"
          "out vec4 fcolor;"
          "out vec2 fuv;"
          "void main() {"
          "  gl_Position = vec4(2.0*(vertex.x+offset.x-0.5)/wndSize.x - 1.0,"
          "                     1.0 - 2.0*(vertex.y+offset.y-0.5)/wndSize.y,"
          "                     depth, 1);"
          " fuv = uv;"
          " fcolor = color;"
          "}";

      std::string fragmentShader =
        "#version 150\n"
        "uniform sampler2D tex;"
        "in vec2 fuv;"
        "in vec4 fcolor;"
        "out vec4 outColor;"
        "void main() { outColor.rgb = fcolor.bgr; outColor.a = fcolor.a * texture(tex, fuv).r; }";

      if (!m_shader.init(vertexShader, fragmentShader, "outColor"))
        return 0;

      /* Adjust location bindings */
      glBindAttribLocation(m_shader.program_shader, 0, "vertex");
      glBindAttribLocation(m_shader.program_shader, 1, "uv");
      glBindAttribLocation(m_shader.program_shader, 2, "color");
      glLinkProgram(m_shader.program_shader);

      m_shaderHandleBackup = m_TriTexUniProgram;
      m_TriTexUniProgram = m_shader.program_shader;
      m_TriTexUniLocationOffset = m_shader.uniform("offset");
      m_TriTexUniLocationWndSize = m_shader.uniform("wndSize");
      m_TriTexUniLocationColor = m_shader.uniform("color");
      m_TriTexUniLocationTexture = m_shader.uniform("tex");
      m_TriTexUniLocationDepth = m_shader.uniform("depth");
    }
    return retval;
  }

  IGL_INLINE int igl::TextRenderer::Shut()
  {
    for (auto kv : m_textObjects)
      DeleteTextObj(kv.second);
    m_shader.free();
    m_TriTexUniProgram = m_shaderHandleBackup;
    return CTwGraphOpenGLCore::Shut();
  }

  IGL_INLINE void igl::TextRenderer::BeginDraw(const Eigen::Matrix4f &view, const Eigen::Matrix4f &proj,
    const Eigen::Vector4f &_viewport, float _object_scale)
  {
    viewport = _viewport;
    proj_matrix = proj;
    view_matrix = view;
    CTwGraphOpenGLCore::BeginDraw(viewport[2], viewport[3]);
    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_FALSE);
    object_scale = _object_scale;
  }

  IGL_INLINE void igl::TextRenderer::EndDraw()
  {
    /* Limit the number of cached text objects */
    for (auto it = m_textObjects.cbegin(); it != m_textObjects.cend(); )
    {
      if (m_textObjects.size() < 1000000)
        break;
      DeleteTextObj(it->second);
      m_textObjects.erase(it++);
    }

    glDepthMask(GL_TRUE);
    CTwGraphOpenGLCore::EndDraw();
  }

  IGL_INLINE void igl::TextRenderer::DrawText(Eigen::Vector3d pos, Eigen::Vector3d normal, const std::string &text)
  {
    pos += normal * 0.005f * object_scale;
    Eigen::Vector3f coord = igl::project(Eigen::Vector3f(pos(0), pos(1), pos(2)),
        view_matrix, proj_matrix, viewport);
    auto it = m_textObjects.find(text);
    void *text_obj = nullptr;
    if (it == m_textObjects.end())
    {
      text_obj = NewTextObj();
      BuildText(text_obj, &text, NULL, NULL, 1, g_DefaultNormalFont, 0, 0);
      m_textObjects[text] = text_obj;
    } else {
      text_obj = it->second;
    }
    m_shader.bind();
    glUniform1f(m_TriTexUniLocationDepth, 2*(coord(2)-0.5f));
    CTwGraphOpenGLCore::DrawText(text_obj, coord[0], viewport[3] - coord[1], COLOR32_BLUE, 0);
  }
