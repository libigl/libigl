// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
//               2014 Christian Schüller <schuellchr@gmail.com> 
//
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_HIT_H
#define IGL_HIT_H

namespace igl
{
  // Reimplementation of the embree::Hit struct from embree1.0
  // 
  // TODO: template on floating point type
  struct Hit
  {
    int id; // primitive id
    int gid; // geometry id (not used)
    // barycentric coordinates so that 
    //   pos = V.row(F(id,0))*(1-u-v)+V.row(F(id,1))*u+V.row(F(id,2))*v;
    float u,v; 
    // parametric distance so that
    //   pos = origin + t * dir
    float t; 
  };
}
#endif 
