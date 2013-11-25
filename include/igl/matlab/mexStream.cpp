// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "mexStream.h"
#include "mex.h"

std::streamsize igl::mexStream::xsputn(const char *s, std::streamsize n) 
{
  mexPrintf("%.*s",n,s);
  mexEvalString("drawnow;"); // to dump string.
  return n;
}

int igl::mexStream::overflow(int c) 
{
    if (c != EOF) {
      mexPrintf("%.1s",&c);
      mexEvalString("drawnow;"); // to dump string.
    }
    return 1;
}
