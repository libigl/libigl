// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "init_render_to_texture.h"
#include <cassert>

IGL_INLINE void igl::opengl::init_render_to_texture(
  const size_t width,
  const size_t height,
  GLuint & tex_id,
  GLuint & fbo_id,
  GLuint & dfbo_id)
{
  using namespace std;
  // Delete if already exists
  glDeleteTextures(1,&tex_id);
  glDeleteFramebuffersEXT(1,&fbo_id);
  glDeleteFramebuffersEXT(1,&dfbo_id);
  // http://www.opengl.org/wiki/Framebuffer_Object_Examples#Quick_example.2C_render_to_texture_.282D.29
  glGenTextures(1, &tex_id);
  glBindTexture(GL_TEXTURE_2D, tex_id);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  //NULL means reserve texture memory, but texels are undefined
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, width, height, 0, GL_BGRA, GL_FLOAT, NULL);
  glBindTexture(GL_TEXTURE_2D, 0);
  glGenFramebuffersEXT(1, &fbo_id);
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo_id);
  //Attach 2D texture to this FBO
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, tex_id, 0);

  glGenRenderbuffersEXT(1, &dfbo_id);
  glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, dfbo_id);
  glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT24, width, height);
  //Attach depth buffer to FBO (for this example it's not really needed, but if
  //drawing a 3D scene it would be necessary to attach something)
  glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, dfbo_id);

  //Does the GPU support current FBO configuration?
  GLenum status;
  status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
  assert(status == GL_FRAMEBUFFER_COMPLETE_EXT);
  // Unbind to clean up
  glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, 0);
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
}
