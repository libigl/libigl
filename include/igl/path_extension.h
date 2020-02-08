// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2019 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_PATH_EXTENSION_H
#define IGL_PATH_EXTENSION_H
#include "igl_inline.h"
#include <string>

namespace igl
{
  // Input:
  //  path  string containing input path
  // Return the extension determined by igl::pathinfo
  IGL_INLINE std::string path_extension(
    const std::string & path);

}

#ifndef IGL_STATIC_LIBRARY
#  include "path_extension.cpp"
#endif

#endif
