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
  struct FractionalTessellationPattern
  {
    FractionalTessellationPattern (float tess, const bool sublevel = false)
    : sublevel(sublevel)
    {
      rcp_tess = 1.0f/tess;
      if (sublevel) tess *= 0.5f;
      Nl = floor(tess);
      Nh = Nl+1;
      th = tess-float(Nh);
      tl = float(Nl)-tess;
      N = Nl+1+(Nl % 2);
      if (sublevel) N *= 2;
    }
    
    /* returns number of intervals (points-1) */
    __forceinline int size() const {
      return N;
    }
    
    __forceinline float operator() (int i) const
    {
      float ofs = 0.0f;
      if (sublevel && i>=N/2) {
        ofs = 0.5f;
        i-=N/2;
      }
      
      const int Nl2 = Nl/2;
      if (Nl % 2 == 0) {
        const float f0 = float(i > Nl2) * th;
        return min((float(i) + f0)*rcp_tess+ofs,1.0f);
      } 
      else {
        const float f0 = (i > (Nl2+0)) * th;
        const float f1 = (i > (Nl2+1)) * tl;
        const float f2 = (i > (Nl2+2)) * th;
        return min((float(i) + f0 + f1 + f2)*rcp_tess+ofs,1.0f);
      }
    }
    
  private:
    bool sublevel;
    float rcp_tess;
    int   Nl, Nh, N;
    float tl, th;
  };
}
