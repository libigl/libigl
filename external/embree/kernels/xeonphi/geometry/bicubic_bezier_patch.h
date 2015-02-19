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
#if defined(__MIC__)

  __forceinline mic_f bc(const Vec3fa &v)
  {
    return broadcast4to16f(&v);
  }

  __forceinline mic3f vec3fa_to_mic3f(const Vec3fa &v) 
  {
    return mic3f(mic_f(v.x),mic_f(v.y),mic_f(v.z));
  }

  class BicubicBezierPatch
  {
  public:
    Vec3fa cp[4][4];

    BicubicBezierPatch() {}

    BicubicBezierPatch(Vec3fa matrix[4][4]) 
      {
	for (size_t y=0;y<4;y++)
	  for (size_t x=0;x<4;x++)
	    cp[y][x] = matrix[y][x];
	
      }

    // ==================================================
    // ======== evaluate using the full 4x4 matrix ======
    // ==================================================

    __forceinline mic_f lerp(const mic_f &u, const mic_f &v, const mic_f &x, const mic_f &y) const
    {
      return u * x + v * y;
    }

    __forceinline mic_f eval(const float uu,
			     const float vv) const
    {
      const mic_f u = mic_f(uu);
      const mic_f v = mic_f(vv);
      const mic_f one = mic_f::one();
      const mic_f one_minus_u = one - u;
      const mic_f one_minus_v = one - v;
      const mic_f v00 = one_minus_v * load16f((float*)&cp[0][0]) + v * load16f((float*)&cp[1][0]);
      const mic_f v01 = one_minus_v * load16f((float*)&cp[1][0]) + v * load16f((float*)&cp[2][0]);
      const mic_f v02 = one_minus_v * load16f((float*)&cp[2][0]) + v * load16f((float*)&cp[3][0]);
      const mic_f v10 = one_minus_v * v00 + v * v01;
      const mic_f v11 = one_minus_v * v01 + v * v02;
      const mic_f v20 = one_minus_v * v10 + v * v11;

      const Vec3fa *const u_ptr = (Vec3fa*)&v20; // FIXME
      const mic_f u00 = one_minus_u * bc(u_ptr[0]) + u * bc(u_ptr[1]);
      const mic_f u01 = one_minus_u * bc(u_ptr[1]) + u * bc(u_ptr[2]);
      const mic_f u02 = one_minus_u * bc(u_ptr[2]) + u * bc(u_ptr[3]);
      const mic_f u10 = one_minus_u * u00 + u * u01;
      const mic_f u11 = one_minus_u * u01 + u * u02;
      const mic_f u20 = one_minus_u * u10 + u * u11;
      return u20;
    }


    // ==================================================
    // ======== evaluate assuming a 4x3 matrix (U) ======
    // ==================================================

    __forceinline mic_f evalU_4x3(const float uu,
				  const float vv) const
    {
      const mic_f u = mic_f(uu);
      const mic_f v = mic_f(vv);
      const mic_f one = mic_f::one();
      const mic_f one_minus_u = one - u;
      const mic_f one_minus_v = one - v;

      const mic_f v00 = one_minus_v * load16f((float*)&cp[0][0]) + v * load16f((float*)&cp[1][0]);
      const mic_f v01 = one_minus_v * load16f((float*)&cp[1][0]) + v * load16f((float*)&cp[2][0]);
      const mic_f v02 = one_minus_v * load16f((float*)&cp[2][0]) + v * load16f((float*)&cp[3][0]);
      const mic_f v10 = one_minus_v * v00 + v * v01;
      const mic_f v11 = one_minus_v * v01 + v * v02;
      const mic_f v20 = one_minus_v * v10 + v * v11;

      const Vec3fa *const u_ptr = (Vec3fa*)&v20; //FIXME
      const mic_f u00 =  one_minus_u * bc(u_ptr[0]) + u * bc(u_ptr[1]);
      const mic_f u01 =  one_minus_u * bc(u_ptr[1]) + u * bc(u_ptr[2]);
      const mic_f u10 =  one_minus_u * u00 + u * u01;

      return u10;
    }

    __forceinline mic3f evalU_4x3(const mic_f u,
				  const mic_f v) const
    {
      const mic_f one = mic_f::one();
      const mic_f one_minus_u = one - u;
      const mic_f one_minus_v = one - v;

      const mic3f u00 = one_minus_u * vec3fa_to_mic3f(cp[0][0]) + u * vec3fa_to_mic3f(cp[0][1]);
      const mic3f u01 = one_minus_u * vec3fa_to_mic3f(cp[0][1]) + u * vec3fa_to_mic3f(cp[0][2]);
      const mic3f u0  = one_minus_u * u00+ u * u01;
      const mic3f u10 = one_minus_u * vec3fa_to_mic3f(cp[1][0]) + u * vec3fa_to_mic3f(cp[1][1]);
      const mic3f u11 = one_minus_u * vec3fa_to_mic3f(cp[1][1]) + u * vec3fa_to_mic3f(cp[1][2]);
      const mic3f u1  = one_minus_u * u10+ u * u11;
      const mic3f u20 = one_minus_u * vec3fa_to_mic3f(cp[2][0]) + u * vec3fa_to_mic3f(cp[2][1]);
      const mic3f u21 = one_minus_u * vec3fa_to_mic3f(cp[2][1]) + u * vec3fa_to_mic3f(cp[2][2]);
      const mic3f u2  = one_minus_u * u20+ u * u21;

      const mic3f u30 = one_minus_u * vec3fa_to_mic3f(cp[3][0]) + u * vec3fa_to_mic3f(cp[3][1]);
      const mic3f u31 = one_minus_u * vec3fa_to_mic3f(cp[3][1]) + u * vec3fa_to_mic3f(cp[3][2]);
      const mic3f u3  = one_minus_u * u30 + u * u31;

      const mic3f v00 = one_minus_v * u0 + v * u1;
      const mic3f v01 = one_minus_v * u1 + v * u2;
      const mic3f v02 = one_minus_v * u2 + v * u3;

      const mic3f v10 = one_minus_v * v00 + v * v01;
      const mic3f v11 = one_minus_v * v01 + v * v02;
      const mic3f v20 = one_minus_v * v10 + v * v11;
      return v20;
    }

    // ==================================================
    // ======== evaluate assuming a 3x4 matrix (V) ======
    // ==================================================

    __forceinline mic_f evalV_3x4(const float uu,
				  const float vv) const
    {
      const mic_f u = mic_f(uu);
      const mic_f v = mic_f(vv);
      const mic_f one = mic_f::one();
      const mic_f one_minus_u = one - u;
      const mic_f one_minus_v = one - v;
      const mic_f v00 = one_minus_v * load16f((float*)&cp[0][0]) + v * load16f((float*)&cp[1][0]);
      const mic_f v01 = one_minus_v * load16f((float*)&cp[1][0]) + v * load16f((float*)&cp[2][0]);
      const mic_f v10 = one_minus_v * v00 + v * v01;

      const Vec3fa *const u_ptr = (Vec3fa*)&v10;
      const mic_f u00 = one_minus_u * bc(u_ptr[0]) + u * bc(u_ptr[1]);
      const mic_f u01 = one_minus_u * bc(u_ptr[1]) + u * bc(u_ptr[2]);
      const mic_f u02 = one_minus_u * bc(u_ptr[2]) + u * bc(u_ptr[3]);
      const mic_f u10 = one_minus_u * u00 + u * u01;
      const mic_f u11 = one_minus_u * u01 + u * u02;
      const mic_f u20 = one_minus_u * u10 + u * u11;
      return u20;
    }

    __forceinline mic3f evalV_3x4(const mic_f u,
				  const mic_f v) const
    {
      const mic_f one = mic_f::one();
      const mic_f one_minus_u = one - u;
      const mic_f one_minus_v = one - v;      
      const mic3f v00 = one_minus_v * vec3fa_to_mic3f(cp[0][0]) + v * vec3fa_to_mic3f(cp[1][0]);
      const mic3f v01 = one_minus_v * vec3fa_to_mic3f(cp[1][0]) + v * vec3fa_to_mic3f(cp[2][0]);
      const mic3f v0  = one_minus_v * v00 + v * v01;
      const mic3f v10 = one_minus_v * vec3fa_to_mic3f(cp[0][1]) + v * vec3fa_to_mic3f(cp[1][1]);
      const mic3f v11 = one_minus_v * vec3fa_to_mic3f(cp[1][1]) + v * vec3fa_to_mic3f(cp[2][1]);
      const mic3f v1  = one_minus_v * v10 + v * v11;
      const mic3f v20 = one_minus_v * vec3fa_to_mic3f(cp[0][2]) + v * vec3fa_to_mic3f(cp[1][2]);
      const mic3f v21 = one_minus_v * vec3fa_to_mic3f(cp[1][2]) + v * vec3fa_to_mic3f(cp[2][2]);
      const mic3f v2  = one_minus_v * v20 + v * v21;
      const mic3f v30 = one_minus_v * vec3fa_to_mic3f(cp[0][3]) + v * vec3fa_to_mic3f(cp[1][3]);
      const mic3f v31 = one_minus_v * vec3fa_to_mic3f(cp[1][3]) + v * vec3fa_to_mic3f(cp[2][3]);
      const mic3f v3  = one_minus_v * v30 + v * v31;

      const mic3f u00 = one_minus_u * v0 + u * v1;
      const mic3f u10 = one_minus_u * v1 + u * v2;
      const mic3f u20 = one_minus_u * v2 + u * v3;

      const mic3f u01 = one_minus_u * u00 + u * u10;
      const mic3f u11 = one_minus_u * u10 + u * u20;
      const mic3f u02 = one_minus_u * u01 + u * u11;
      return u02;
    }

    // ==============================================================
    // ======== fast evaluation code using four samples at once ======
    // ==============================================================
  
    __forceinline mic_f eval4(const mic_f u,
			      const mic_f v) const
    {
      const mic_f one = mic_f::one();
      const mic_f one_minus_u = one - u;
      const mic_f one_minus_v = one - v;      

      mic_f dest[4];
#pragma unroll(4)
      for (int i=0;i<4;i++)
	{
	  const mic_f v00 = lerp(one_minus_v,v,broadcast4to16f((const float*)&cp[0][i]),broadcast4to16f((const float*)&cp[1][i]));
	  const mic_f v01 = lerp(one_minus_v,v,broadcast4to16f((const float*)&cp[1][i]),broadcast4to16f((const float*)&cp[2][i]));
	  const mic_f v02 = lerp(one_minus_v,v,broadcast4to16f((const float*)&cp[2][i]),broadcast4to16f((const float*)&cp[3][i]));

	  const mic_f v10 = lerp(one_minus_v,v,v00,v01);
	  const mic_f v11 = lerp(one_minus_v,v,v01,v02);
	  const mic_f v20 = lerp(one_minus_v,v,v10,v11);
	  dest[i] = v20;
	}

      const mic_f u00 = lerp(one_minus_u,u,dest[0],dest[1]);
      const mic_f u01 = lerp(one_minus_u,u,dest[1],dest[2]);
      const mic_f u02 = lerp(one_minus_u,u,dest[2],dest[3]);

      const mic_f u10 = lerp(one_minus_u,u,u00,u01);
      const mic_f u11 = lerp(one_minus_u,u,u01,u02);
      const mic_f u20 = lerp(one_minus_u,u,u10,u11);

      return u20;
    }

    // ===========================================
    // ======== evaluate 16 samples at once ======
    // ===========================================


    __forceinline mic3f eval(const mic_f u,
			     const mic_f v) const
    {
      const mic_f one = mic_f::one();
      const mic_f one_minus_u = one - u;
      const mic_f one_minus_v = one - v;      

      const mic_f B0_u = one_minus_u * one_minus_u * one_minus_u;
      const mic_f B0_v = one_minus_v * one_minus_v * one_minus_v;
      const mic_f B1_u = 3.0f * one_minus_u * one_minus_u * u;
      const mic_f B1_v = 3.0f * one_minus_v * one_minus_v * v;
      const mic_f B2_u = 3.0f * one_minus_u * u * u;
      const mic_f B2_v = 3.0f * one_minus_v * v * v;
      const mic_f B3_u = u * u * u;
      const mic_f B3_v = v * v * v;

      const mic_f x = \
	(B0_u * cp[0][0].x + B1_u * cp[0][1].x + B2_u * cp[0][2].x + B3_u * cp[0][3].x) * B3_v + 
	(B0_u * cp[1][0].x + B1_u * cp[1][1].x + B2_u * cp[1][2].x + B3_u * cp[1][3].x) * B2_v + 
	(B0_u * cp[2][0].x + B1_u * cp[2][1].x + B2_u * cp[2][2].x + B3_u * cp[2][3].x) * B1_v + 
	(B0_u * cp[3][0].x + B1_u * cp[3][1].x + B2_u * cp[3][2].x + B3_u * cp[3][3].x) * B0_v; 

      const mic_f y = \
	(B0_u * cp[0][0].y + B1_u * cp[0][1].y + B2_u * cp[0][2].y + B3_u * cp[0][3].y) * B3_v + 
	(B0_u * cp[1][0].y + B1_u * cp[1][1].y + B2_u * cp[1][2].y + B3_u * cp[1][3].y) * B2_v + 
	(B0_u * cp[2][0].y + B1_u * cp[2][1].y + B2_u * cp[2][2].y + B3_u * cp[2][3].y) * B1_v + 
	(B0_u * cp[3][0].y + B1_u * cp[3][1].y + B2_u * cp[3][2].y + B3_u * cp[3][3].y) * B0_v; 

      const mic_f z = \
	(B0_u * cp[0][0].z + B1_u * cp[0][1].z + B2_u * cp[0][2].z + B3_u * cp[0][3].z) * B3_v + 
	(B0_u * cp[1][0].z + B1_u * cp[1][1].z + B2_u * cp[1][2].z + B3_u * cp[1][3].z) * B2_v + 
	(B0_u * cp[2][0].z + B1_u * cp[2][1].z + B2_u * cp[2][2].z + B3_u * cp[2][3].z) * B1_v + 
	(B0_u * cp[3][0].z + B1_u * cp[3][1].z + B2_u * cp[3][2].z + B3_u * cp[3][3].z) * B0_v; 
      return mic3f(x,y,z);
    }


    __forceinline void init(const Vec3fa cc_patch[4][4]) // from a b-spline patch
    {
      const mic_f b11 = (bc(cc_patch[1][1]) * 4.0f + (bc(cc_patch[1][2]) + bc(cc_patch[2][1])) * 2.0f + bc(cc_patch[2][2])) * 1.0f / 9.0f;
      const mic_f b12 = (bc(cc_patch[1][2]) * 4.0f + (bc(cc_patch[1][1]) + bc(cc_patch[2][2])) * 2.0f + bc(cc_patch[2][1])) * 1.0f / 9.0f;
      const mic_f b22 = (bc(cc_patch[2][2]) * 4.0f + (bc(cc_patch[2][1]) + bc(cc_patch[1][2])) * 2.0f + bc(cc_patch[1][1])) * 1.0f / 9.0f;
      const mic_f b21 = (bc(cc_patch[2][1]) * 4.0f + (bc(cc_patch[1][1]) + bc(cc_patch[2][2])) * 2.0f + bc(cc_patch[1][2])) * 1.0f / 9.0f;

      store4f((float*)&cp[1][1],b11);
      store4f((float*)&cp[1][2],b12);
      store4f((float*)&cp[2][2],b22);
      store4f((float*)&cp[2][1],b21);

      // --- edges ---

      const mic_f b01 = (bc(cc_patch[1][1]) * 8.0f + bc(cc_patch[1][2]) * 4.0f + (bc(cc_patch[0][1]) + bc(cc_patch[2][1])) * 2.0f + bc(cc_patch[0][2]) + bc(cc_patch[2][2])) * 1.0f / 18.0f;
      const mic_f b02 = (bc(cc_patch[1][2]) * 8.0f + bc(cc_patch[1][1]) * 4.0f + (bc(cc_patch[0][2]) + bc(cc_patch[2][2])) * 2.0f + bc(cc_patch[0][1]) + bc(cc_patch[2][1])) * 1.0f / 18.0f;

      store4f((float*)&cp[0][1],b01);
      store4f((float*)&cp[0][2],b02);

      const mic_f b13 = (bc(cc_patch[1][2]) * 8.0f + bc(cc_patch[2][2]) * 4.0f + (bc(cc_patch[1][1]) + bc(cc_patch[1][3])) * 2.0f + bc(cc_patch[2][1]) + bc(cc_patch[2][3])) * 1.0f / 18.0f;
      const mic_f b23 = (bc(cc_patch[2][2]) * 8.0f + bc(cc_patch[1][2]) * 4.0f + (bc(cc_patch[2][1]) + bc(cc_patch[2][3])) * 2.0f + bc(cc_patch[1][1]) + bc(cc_patch[1][3])) * 1.0f / 18.0f;

      store4f((float*)&cp[1][3],b13);
      store4f((float*)&cp[2][3],b23);

      const mic_f b32 = (bc(cc_patch[2][2]) * 8.0f + bc(cc_patch[2][1]) * 4.0f + (bc(cc_patch[1][2]) + bc(cc_patch[3][2])) * 2.0f + bc(cc_patch[1][1]) + bc(cc_patch[3][1])) * 1.0f / 18.0f;
      const mic_f b31 = (bc(cc_patch[2][1]) * 8.0f + bc(cc_patch[2][2]) * 4.0f + (bc(cc_patch[1][1]) + bc(cc_patch[3][1])) * 2.0f + bc(cc_patch[1][2]) + bc(cc_patch[3][2])) * 1.0f / 18.0f;

      store4f((float*)&cp[3][2],b32);
      store4f((float*)&cp[3][1],b31);

      const mic_f b20 = (bc(cc_patch[2][1]) * 8.0f + bc(cc_patch[1][1]) * 4.0f + (bc(cc_patch[2][0]) + bc(cc_patch[2][2])) * 2.0f + bc(cc_patch[1][0]) + bc(cc_patch[1][2])) * 1.0f / 18.0f;
      const mic_f b10 = (bc(cc_patch[1][1]) * 8.0f + bc(cc_patch[2][1]) * 4.0f + (bc(cc_patch[1][0]) + bc(cc_patch[1][2])) * 2.0f + bc(cc_patch[2][0]) + bc(cc_patch[2][2])) * 1.0f / 18.0f;

      store4f((float*)&cp[2][0],b20);
      store4f((float*)&cp[1][0],b10);

      // --- corner ----

      const mic_f b00 = (bc(cc_patch[1][1]) * 16.0f + (bc(cc_patch[0][1]) + bc(cc_patch[1][2]) + bc(cc_patch[2][1]) + bc(cc_patch[1][0])) * 4.0f + 
			 (bc(cc_patch[0][0]) + bc(cc_patch[0][2]) + bc(cc_patch[2][2]) + bc(cc_patch[2][0]))) * 1.0f / 36.0f;
      const mic_f b03 = (bc(cc_patch[1][2]) * 16.0f + (bc(cc_patch[0][2]) + bc(cc_patch[1][3]) + bc(cc_patch[2][2]) + bc(cc_patch[1][1])) * 4.0f + 
			 (bc(cc_patch[0][1]) + bc(cc_patch[0][3]) + bc(cc_patch[2][3]) + bc(cc_patch[2][1]))) * 1.0f / 36.0f;
      const mic_f b33 = (bc(cc_patch[2][2]) * 16.0f + (bc(cc_patch[1][2]) + bc(cc_patch[2][3]) + bc(cc_patch[3][2]) + bc(cc_patch[2][1])) * 4.0f + 
			 (bc(cc_patch[1][1]) + bc(cc_patch[1][3]) + bc(cc_patch[3][3]) + bc(cc_patch[3][1]))) * 1.0f / 36.0f;
      const mic_f b30 = (bc(cc_patch[2][1]) * 16.0f + (bc(cc_patch[1][1]) + bc(cc_patch[2][2]) + bc(cc_patch[3][1]) + bc(cc_patch[2][0])) * 4.0f + 
			 (bc(cc_patch[1][0]) + bc(cc_patch[1][2]) + bc(cc_patch[3][2]) + bc(cc_patch[3][0]))) * 1.0f / 36.0f;

      store4f((float*)&cp[0][0],b00);
      store4f((float*)&cp[0][3],b03);
      store4f((float*)&cp[3][3],b33);
      store4f((float*)&cp[3][0],b30);    
    }


 
  };

  __forceinline void createTangentPatchU4x3(const BicubicBezierPatch &source, BicubicBezierPatch &dest)
  {
    for (int y=0;y<4;y++) 
      for (int x=0;x<3;x++) 
	store4f(&dest.cp[y][x],(bc(source.cp[y][x+1]) - bc(source.cp[y][x])) * 3.0f);
  }

  __forceinline void createTangentPatchV3x4(const BicubicBezierPatch &source, BicubicBezierPatch &dest)
  {
    for (int y=0;y<3;y++) 
      for (int x=0;x<4;x++) 
	store4f(&dest.cp[y][x],(bc(source.cp[y+1][x]) - bc(source.cp[y][x])) * 3.0f);
  }


  __forceinline mic_f getNormalFromTangentPatches(const BicubicBezierPatch &tangentU, 
						  const BicubicBezierPatch &tangentV,
						  const float u,
						  const float v) 
  {
    const mic_f tU = tangentU.evalU_4x3(u,v);
    const mic_f tV = tangentV.evalV_3x4(u,v);
    return lcross_xyz(tU,tV);
  }

  __forceinline mic3f getNormalFromTangentPatches(const BicubicBezierPatch &tangentU, 
						  const BicubicBezierPatch &tangentV,
						  const mic_f u,
						  const mic_f v,
						  const mic_f one_minus_u,
						  const mic_f one_minus_v) 
  {
    const mic3f tV = tangentV.evalV_3x4(u,v);
    const mic3f tU = tangentU.evalU_4x3(u,v);
    const mic3f normal = cross(tU,tV); 
    return normal;
  }


  __forceinline std::ostream &operator<<(std::ostream &o, const BicubicBezierPatch &patch)
    {
      for (int y=0;y<4;y++) 
	for (int x=0;x<4;x++) 
	  o << "y = " << y << " x = " << x << " : " << patch.cp[y][x] << std::endl;
      return o;
    } 


#endif

};
