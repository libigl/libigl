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

#pragma once

namespace embree
{
  struct DiscreteTessellationPattern
  {
    DiscreteTessellationPattern (float ftess, const int sublevel = 0)
    : sublevel(sublevel)
    {
      for (size_t i=0; i<sublevel; i++) ftess *= 0.5f;
      int tess = ceil(ftess);
      for (size_t i=0; i<sublevel; i++) tess *= 2;
      rcp_tess = 1.0f/float(tess);
      N = tess;
    }
    
    /* returns number of intervals (points-1) */
    __forceinline int size() const {
      return N;
    }
    
    __forceinline float operator() (int i) const {
      return min(float(i)*rcp_tess,1.0f);
    }
    
  private:
    int sublevel;
    float rcp_tess;
    int   N;
  };
}
