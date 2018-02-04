// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Wenzel Jacob <wenzel@inf.ethz.ch>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifdef IGL_VIEWER_WITH_NANOGUI
#include "TextRenderer.h"
#include "TextRenderer_fonts.h"
#include <igl/project.h>

#include <nanogui/opengl.h>
#include <nanovg.h>

#include <Eigen/Dense>

#define NANOVG_GL3
#include <nanovg_gl.h>


IGL_INLINE igl::opengl::glfw::TextRenderer::TextRenderer(): ctx(nullptr) {}

IGL_INLINE int igl::opengl::glfw::TextRenderer::Init()
{
  using namespace std;
  #ifdef NDEBUG
    ctx = nvgCreateGL3(NVG_STENCIL_STROKES | NVG_ANTIALIAS);
  #else
    ctx = nvgCreateGL3(NVG_STENCIL_STROKES | NVG_ANTIALIAS | NVG_DEBUG);
  #endif

  nvgCreateFontMem(ctx, "sans", igl_roboto_regular_ttf,
                             igl_roboto_regular_ttf_size, 0);

  return 0;
}

IGL_INLINE int igl::opengl::glfw::TextRenderer::Shut()
{
  using namespace std;
  if(ctx)
    nvgDeleteGL3(ctx);
  return 0;
}

IGL_INLINE void igl::opengl::glfw::TextRenderer::BeginDraw(
  const Eigen::Matrix4f &view,
  const Eigen::Matrix4f &proj,
  const Eigen::Vector4f &_viewport,
  float _object_scale)
{
  using namespace std;
  viewport = _viewport;
  proj_matrix = proj;
  view_matrix = view;
  object_scale = _object_scale;

  Eigen::Vector2i mFBSize;
  Eigen::Vector2i mSize;

  GLFWwindow* mGLFWWindow = glfwGetCurrentContext();
  glfwGetFramebufferSize(mGLFWWindow,&mFBSize[0],&mFBSize[1]);
  glfwGetWindowSize(mGLFWWindow,&mSize[0],&mSize[1]);
  glViewport(0,0,mFBSize[0],mFBSize[1]);

  glClear(GL_STENCIL_BUFFER_BIT);

  /* Calculate pixel ratio for hi-dpi devices. */
  mPixelRatio = (float)mFBSize[0] / (float)mSize[0];
  nvgBeginFrame(ctx,mSize[0],mSize[1],mPixelRatio);
}

IGL_INLINE void igl::opengl::glfw::TextRenderer::EndDraw()
{
  using namespace std;
  nvgEndFrame(ctx);
}

IGL_INLINE void igl::opengl::glfw::TextRenderer::DrawText(
  Eigen::Vector3d pos, Eigen::Vector3d normal, const std::string &text)
{
  using namespace std;
  pos += normal * 0.005f * object_scale;
  Eigen::Vector3f coord = igl::project(Eigen::Vector3f(pos(0), pos(1), pos(2)),
      view_matrix, proj_matrix, viewport);

  nvgFontSize(ctx, 16*mPixelRatio);
  nvgFontFace(ctx, "sans");
  nvgTextAlign(ctx, NVG_ALIGN_LEFT | NVG_ALIGN_MIDDLE);
  nvgFillColor(ctx, nvgRGBA(10,10,250,255));
  nvgText(ctx, coord[0]/mPixelRatio, (viewport[3] - coord[1])/mPixelRatio, text.c_str(), NULL);
}
#endif
