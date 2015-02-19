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

#include "common/default.h"
#include "common/geometry.h"

namespace embree
{
  namespace isa
  {
    class BicubicBezierPatch
    {
    public:
      Vec3fa cp[4][4];
      
      
      // ==================================================
      // ======== evaluate using the full 4x4 matrix ======
      // ==================================================
      
      __forceinline Vec3fa lerp(const Vec3fa &u, const Vec3fa &v, const Vec3fa &x, const Vec3fa &y) const
      {
        return u * x + v * y;
      }
      
      __forceinline Vec3fa eval(const float uu,
                                const float vv) const
      {
        const Vec3fa u = Vec3fa(uu);
        const Vec3fa v = Vec3fa(vv);
        const Vec3fa one = Vec3fa(1.0f);
        const Vec3fa one_minus_u = one - u;
        const Vec3fa one_minus_v = one - v;
        
        Vec3fa u_ptr[4]; 
        
        for (size_t i=0;i<4;i++)
	{
	  const Vec3fa v00 = one_minus_v * cp[0][i] + v * cp[1][i];
	  const Vec3fa v01 = one_minus_v * cp[1][i] + v * cp[2][i];
	  const Vec3fa v02 = one_minus_v * cp[2][i] + v * cp[3][i];
	  const Vec3fa v10 = one_minus_v * v00 + v * v01;
	  const Vec3fa v11 = one_minus_v * v01 + v * v02;
	  const Vec3fa v20 = one_minus_v * v10 + v * v11;
	  u_ptr[i] = v20;
	}
        
        const Vec3fa u00 = one_minus_u * u_ptr[0] + u * u_ptr[1];
        const Vec3fa u01 = one_minus_u * u_ptr[1] + u * u_ptr[2];
        const Vec3fa u02 = one_minus_u * u_ptr[2] + u * u_ptr[3];
        const Vec3fa u10 = one_minus_u * u00 + u * u01;
        const Vec3fa u11 = one_minus_u * u01 + u * u02;
        const Vec3fa u20 = one_minus_u * u10 + u * u11;
        return u20;
      }
      
      
      // ==================================================
      // ======== evaluate assuming a 4x3 matrix (U) ======
      // ==================================================
      
      __forceinline Vec3fa evalU_4x3(const float uu,
                                     const float vv) const
      {
        const Vec3fa u = Vec3fa(uu);
        const Vec3fa v = Vec3fa(vv);
        const Vec3fa one = Vec3fa(1.0f);
        const Vec3fa one_minus_u = one - u;
        const Vec3fa one_minus_v = one - v;
        
        Vec3fa u_ptr[4]; 
        
        for (size_t i=0;i<4;i++)
	{
	  const Vec3fa v00 = one_minus_v * cp[0][i] + v * cp[1][i];
	  const Vec3fa v01 = one_minus_v * cp[1][i] + v * cp[2][i];
	  const Vec3fa v02 = one_minus_v * cp[2][i] + v * cp[3][i];
	  const Vec3fa v10 = one_minus_v * v00 + v * v01;
	  const Vec3fa v11 = one_minus_v * v01 + v * v02;
	  const Vec3fa v20 = one_minus_v * v10 + v * v11;
	  u_ptr[i] = v20;
	}
        
        const Vec3fa u00 =  one_minus_u * u_ptr[0] + u * u_ptr[1];
        const Vec3fa u01 =  one_minus_u * u_ptr[1] + u * u_ptr[2];
        const Vec3fa u10 =  one_minus_u * u00 + u * u01;
        return u10;
      }
      
      // ==================================================
      // ======== evaluate assuming a 3x4 matrix (V) ======
      // ==================================================
      
      __forceinline Vec3fa evalV_3x4(const float uu,
                                     const float vv) const
      {
        const Vec3fa u = Vec3fa(uu);
        const Vec3fa v = Vec3fa(vv);
        const Vec3fa one = Vec3fa(1.0f);
        const Vec3fa one_minus_u = one - u;
        const Vec3fa one_minus_v = one - v;
        Vec3fa u_ptr[4]; 
        
        for (size_t i=0;i<4;i++)
	{
	  const Vec3fa v00 = one_minus_v * cp[0][i] + v * cp[1][i];
	  const Vec3fa v01 = one_minus_v * cp[1][i] + v * cp[2][i];
	  const Vec3fa v10 = one_minus_v * v00 + v * v01;
	  u_ptr[i] = v10;
	}
        
        const Vec3fa u00 = one_minus_u * u_ptr[0] + u * u_ptr[1];
        const Vec3fa u01 = one_minus_u * u_ptr[1] + u * u_ptr[2];
        const Vec3fa u02 = one_minus_u * u_ptr[2] + u * u_ptr[3];
        const Vec3fa u10 = one_minus_u * u00 + u * u01;
        const Vec3fa u11 = one_minus_u * u01 + u * u02;
        const Vec3fa u20 = one_minus_u * u10 + u * u11;
        return u20;
      }
      
      
      
      
      __forceinline void init(const Vec3fa cc_patch[4][4])
      {
        const Vec3fa b11 = ((cc_patch[1][1]) * 4.0f + ((cc_patch[1][2]) + (cc_patch[2][1])) * 2.0f + (cc_patch[2][2])) * 1.0f / 9.0f;
        const Vec3fa b12 = ((cc_patch[1][2]) * 4.0f + ((cc_patch[1][1]) + (cc_patch[2][2])) * 2.0f + (cc_patch[2][1])) * 1.0f / 9.0f;
        const Vec3fa b22 = ((cc_patch[2][2]) * 4.0f + ((cc_patch[2][1]) + (cc_patch[1][2])) * 2.0f + (cc_patch[1][1])) * 1.0f / 9.0f;
        const Vec3fa b21 = ((cc_patch[2][1]) * 4.0f + ((cc_patch[1][1]) + (cc_patch[2][2])) * 2.0f + (cc_patch[1][2])) * 1.0f / 9.0f;
        
        cp[1][1] = b11;
        cp[1][2] = b12;
        cp[2][2] = b22;
        cp[2][1] = b21;
        
        // --- edges ---
        
        const Vec3fa b01 = ((cc_patch[1][1]) * 8.0f + (cc_patch[1][2]) * 4.0f + ((cc_patch[0][1]) + (cc_patch[2][1])) * 2.0f + (cc_patch[0][2]) + (cc_patch[2][2])) * 1.0f / 18.0f;
        const Vec3fa b02 = ((cc_patch[1][2]) * 8.0f + (cc_patch[1][1]) * 4.0f + ((cc_patch[0][2]) + (cc_patch[2][2])) * 2.0f + (cc_patch[0][1]) + (cc_patch[2][1])) * 1.0f / 18.0f;
        
        cp[0][1] = b01;
        cp[0][2] = b02;
        
        const Vec3fa b13 = ((cc_patch[1][2]) * 8.0f + (cc_patch[2][2]) * 4.0f + ((cc_patch[1][1]) + (cc_patch[1][3])) * 2.0f + (cc_patch[2][1]) + (cc_patch[2][3])) * 1.0f / 18.0f;
        const Vec3fa b23 = ((cc_patch[2][2]) * 8.0f + (cc_patch[1][2]) * 4.0f + ((cc_patch[2][1]) + (cc_patch[2][3])) * 2.0f + (cc_patch[1][1]) + (cc_patch[1][3])) * 1.0f / 18.0f;
        
        cp[1][3] = b13;
        cp[2][3] = b23;
        
        const Vec3fa b32 = ((cc_patch[2][2]) * 8.0f + (cc_patch[2][1]) * 4.0f + ((cc_patch[1][2]) + (cc_patch[3][2])) * 2.0f + (cc_patch[1][1]) + (cc_patch[3][1])) * 1.0f / 18.0f;
        const Vec3fa b31 = ((cc_patch[2][1]) * 8.0f + (cc_patch[2][2]) * 4.0f + ((cc_patch[1][1]) + (cc_patch[3][1])) * 2.0f + (cc_patch[1][2]) + (cc_patch[3][2])) * 1.0f / 18.0f;
        
        cp[3][2] = b32;
        cp[3][1] = b31;
        
        const Vec3fa b20 = ((cc_patch[2][1]) * 8.0f + (cc_patch[1][1]) * 4.0f + ((cc_patch[2][0]) + (cc_patch[2][2])) * 2.0f + (cc_patch[1][0]) + (cc_patch[1][2])) * 1.0f / 18.0f;
        const Vec3fa b10 = ((cc_patch[1][1]) * 8.0f + (cc_patch[2][1]) * 4.0f + ((cc_patch[1][0]) + (cc_patch[1][2])) * 2.0f + (cc_patch[2][0]) + (cc_patch[2][2])) * 1.0f / 18.0f;
        
        cp[2][0] = b20;
        cp[1][0] = b10;
        
        // --- corner ----
        
        const Vec3fa b00 = ((cc_patch[1][1]) * 16.0f + ((cc_patch[0][1]) + (cc_patch[1][2]) + (cc_patch[2][1]) + (cc_patch[1][0])) * 4.0f + 
                            ((cc_patch[0][0]) + (cc_patch[0][2]) + (cc_patch[2][2]) + (cc_patch[2][0]))) * 1.0f / 36.0f;
        const Vec3fa b03 = ((cc_patch[1][2]) * 16.0f + ((cc_patch[0][2]) + (cc_patch[1][3]) + (cc_patch[2][2]) + (cc_patch[1][1])) * 4.0f + 
                            ((cc_patch[0][1]) + (cc_patch[0][3]) + (cc_patch[2][3]) + (cc_patch[2][1]))) * 1.0f / 36.0f;
        const Vec3fa b33 = ((cc_patch[2][2]) * 16.0f + ((cc_patch[1][2]) + (cc_patch[2][3]) + (cc_patch[3][2]) + (cc_patch[2][1])) * 4.0f + 
                            ((cc_patch[1][1]) + (cc_patch[1][3]) + (cc_patch[3][3]) + (cc_patch[3][1]))) * 1.0f / 36.0f;
        const Vec3fa b30 = ((cc_patch[2][1]) * 16.0f + ((cc_patch[1][1]) + (cc_patch[2][2]) + (cc_patch[3][1]) + (cc_patch[2][0])) * 4.0f + 
                            ((cc_patch[1][0]) + (cc_patch[1][2]) + (cc_patch[3][2]) + (cc_patch[3][0]))) * 1.0f / 36.0f;
        
        cp[0][0] = b00;
        cp[0][3] = b03;
        cp[3][3] = b33;
        cp[3][0] = b30;    
      }
      
      
      
    };
    
    __forceinline void createTangentPatchU4x3(const BicubicBezierPatch &source, BicubicBezierPatch &dest)
    {
      for (int y=0;y<4;y++) 
        for (int x=0;x<3;x++) 
          dest.cp[y][x] = ((source.cp[y][x+1]) - (source.cp[y][x])) * 3.0f;
    }
    
    __forceinline void createTangentPatchV3x4(const BicubicBezierPatch &source, BicubicBezierPatch &dest)
    {
      for (int y=0;y<3;y++) 
        for (int x=0;x<4;x++) 
          dest.cp[y][x] = ((source.cp[y+1][x]) - (source.cp[y][x])) * 3.0f;
    }
    
    
    __forceinline Vec3fa getNormalFromTangentPatches(const BicubicBezierPatch &tangentU, 
                                                     const BicubicBezierPatch &tangentV,
                                                     const float u,
                                                     const float v) 
    {
      const Vec3fa tU = tangentU.evalU_4x3(u,v);
      const Vec3fa tV = tangentV.evalV_3x4(u,v);
      return cross(tU,tV);
    }
    
    
    
    __forceinline std::ostream &operator<<(std::ostream &o, const BicubicBezierPatch &patch)
    {
      for (int y=0;y<4;y++) 
	for (int x=0;x<4;x++) 
	  o << "y = " << y << " x = " << x << " : " << patch.cp[y][x] << std::endl;
      return o;
    } 
    
    
  };
}
