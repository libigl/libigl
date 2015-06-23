// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Wenzel Jacob <wenzel@inf.ethz.ch>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "TextRenderer.h"
#include <igl/project.h>
#include <nanovg.h>
#define NANOVG_GL3
#include <nanovg_gl.h>

using namespace std;

  IGL_INLINE igl::TextRenderer::TextRenderer() { }

  IGL_INLINE int igl::TextRenderer::Init()
  {
    cerr << "Init TextRenderer" << endl;
    #ifdef NDEBUG
      ctx = nvgCreateGL3(NVG_STENCIL_STROKES | NVG_ANTIALIAS);
    #else
      ctx = nvgCreateGL3(NVG_STENCIL_STROKES | NVG_ANTIALIAS | NVG_DEBUG);
    #endif

    nvgFontSize(ctx, 16);
    nvgFontFace(ctx, "sans");
    nvgTextAlign(ctx, NVG_ALIGN_LEFT | NVG_ALIGN_MIDDLE);
    NVGcolor c;
    c.r = 0.2;
    c.g = 0.2;
    c.b = 255;
    c.a = 255;

    nvgFillColor(ctx, c);
    nvgStrokeColor(ctx, c);

  }

  IGL_INLINE int igl::TextRenderer::Shut()
  {
    cerr << "Shut TextRenderer" << endl;

    nvgDeleteGL3(ctx);
  }

  IGL_INLINE void igl::TextRenderer::BeginDraw(const Eigen::Matrix4f &view, const Eigen::Matrix4f &proj,
    const Eigen::Vector4f &_viewport, float _object_scale)
  {
    cerr << "BeginDraw TextRenderer" << endl;


    viewport = _viewport;
    proj_matrix = proj;
    view_matrix = view;
    // CTwGraphOpenGLCore::BeginDraw(viewport[2], viewport[3]);
    // glEnable(GL_DEPTH_TEST);
    // glDepthMask(GL_FALSE);
    object_scale = _object_scale;


    Eigen::Vector2i mFBSize;
    Eigen::Vector2i mSize;

    GLFWwindow* mGLFWWindow = glfwGetCurrentContext();
    glfwGetFramebufferSize(mGLFWWindow,&mFBSize[0],&mFBSize[1]);
    glfwGetWindowSize(mGLFWWindow,&mSize[0],&mSize[1]);
    glViewport(0,0,mFBSize[0],mFBSize[1]);

    /* Calculate pixel ratio for hi-dpi devices. */
    float mPixelRatio = (float)mFBSize[0] / (float)mSize[0];
    nvgBeginFrame(ctx,mSize[0],mSize[1],mPixelRatio);

  }

  IGL_INLINE void igl::TextRenderer::EndDraw()
  {
    // /* Limit the number of cached text objects */
    // for (auto it = m_textObjects.cbegin(); it != m_textObjects.cend(); )
    // {
    //   if (m_textObjects.size() < 1000000)
    //     break;
    //   DeleteTextObj(it->second);
    //   m_textObjects.erase(it++);
    // }

    // glDepthMask(GL_TRUE);
    // CTwGraphOpenGLCore::EndDraw();
    nvgEndFrame(ctx);
  }

  IGL_INLINE void igl::TextRenderer::DrawText(Eigen::Vector3d pos, Eigen::Vector3d normal, const std::string &text)
  {
    pos += normal * 0.005f * object_scale;
    Eigen::Vector3f coord = igl::project(Eigen::Vector3f(pos(0), pos(1), pos(2)),
        view_matrix, proj_matrix, viewport);

    // auto it = m_textObjects.find(text);
    // void *text_obj = nullptr;
    // if (it == m_textObjects.end())
    // {
    //   text_obj = NewTextObj();
    //   BuildText(text_obj, &text, NULL, NULL, 1, g_DefaultNormalFont, 0, 0);
    //   m_textObjects[text] = text_obj;
    // } else {
    //   text_obj = it->second;
    // }
    // m_shader.bind();
    // glUniform1f(m_TriTexUniLocationDepth, 2*(coord(2)-0.5f));
    //CTwGraphOpenGLCore::DrawText(text_obj, coord[0], viewport[3] - coord[1], COLOR32_BLUE, 0);

    nvgText(ctx, coord[0], viewport[3] - coord[1], text.c_str(), nullptr);


    nvgText(ctx, 10, 10, "Ciao", nullptr);

    cerr << "Draw TextRenderer " << coord[0] << " " << viewport[3] - coord[1] << " " << text << endl;


  }
