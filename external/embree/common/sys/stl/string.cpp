// ======================================================================== //
// Copyright 2009-2014 Intel Corporation                                    //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#include <algorithm>
#include <ctype.h>
#include "string.h"

namespace std
{
  char to_lower(char c) { return char(tolower(int(c))); }
  char to_upper(char c) { return char(toupper(int(c))); }
  string strlwr(const string& s) { string dst(s); std::transform(dst.begin(), dst.end(), dst.begin(), to_lower); return dst; }
  string strupr(const string& s) { string dst(s); std::transform(dst.begin(), dst.end(), dst.begin(), to_upper); return dst; }
}
