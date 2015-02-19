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

#include "catmullclark_patch.h"

namespace embree
{
  template<typename Tessellator>
  struct FeatureAdaptiveSubdivisionBSpline
  {
    Tessellator& tessellator;

    __forceinline FeatureAdaptiveSubdivisionBSpline (int primID, const SubdivMesh::HalfEdge* h, const BufferT<Vec3fa>& vertices, Tessellator& tessellator)
      : tessellator(tessellator)
    {
#if 0
      const Vec2f uv[4] = { Vec2f(0.0f,0.0f),Vec2f(0.0f,1.0f),Vec2f(1.0f,1.0f),Vec2f(1.0f,0.0f) };
      int neighborSubdiv[4];
      const CatmullClarkPatch patch(h,vertices);
      for (size_t i=0; i<4; i++) {
	neighborSubdiv[i] = h->hasOpposite() ? h->opposite()->noRegularFace() : 0; h = h->next();
      }
      subdivide(patch,0,uv,neighborSubdiv);
#else
      int neighborSubdiv[GeneralCatmullClarkPatch::SIZE];
      const GeneralCatmullClarkPatch patch(h,vertices);
      for (size_t i=0; i<patch.size(); i++) {
	neighborSubdiv[i] = h->hasOpposite() ? h->opposite()->noRegularFace() : 0; h = h->next();
      }
      subdivide(patch,0,neighborSubdiv);
#endif
    }

    void subdivide(const GeneralCatmullClarkPatch& patch, int depth, int neighborSubdiv[GeneralCatmullClarkPatch::SIZE])
    {
      /* convert into standard quad patch if possible */
      if (likely(patch.isQuadPatch())) 
      {
	const Vec2f uv[4] = { Vec2f(0.0f,0.0f), Vec2f(0.0f,1.0f), Vec2f(1.0f,1.0f), Vec2f(1.0f,0.0f) };
	CatmullClarkPatch qpatch; patch.init(qpatch);
	subdivide(qpatch,depth,uv,neighborSubdiv);
	return;
      }

      /* subdivide patch */
      size_t N;
      CatmullClarkPatch patches[GeneralCatmullClarkPatch::SIZE]; 
      patch.subdivide(patches,N);

      /* check if subpatches need further subdivision */
      //const bool noleaf = depth > 1;
      bool childSubdiv[GeneralCatmullClarkPatch::SIZE];
      for (size_t i=0; i<N; i++)
        childSubdiv[i] = /*noleaf &&*/ !patches[i].isRegularOrFinal(depth);

      /* parametrization for triangles */
      if (N == 3) {
	const Vec2f uv_0(0.0f,0.0f);
	const Vec2f uv01(0.5f,0.0f);
	const Vec2f uv_1(1.0f,0.0f);
	const Vec2f uv12(0.5f,0.5f);
	const Vec2f uv_2(0.0f,1.0f);
	const Vec2f uv20(0.0f,0.5f);
	const Vec2f uvcc(1.0f/3.0f);
	const Vec2f uv0[4] = { uv_0,uv01,uvcc,uv20 };
	const Vec2f uv1[4] = { uv_1,uv12,uvcc,uv01 };
	const Vec2f uv2[4] = { uv_2,uv20,uvcc,uv12 };
	const int neighborSubdiv0[4] = { max(0,childSubdiv[0]-1),childSubdiv[1],childSubdiv[2],max(0,childSubdiv[2]-1) };
	const int neighborSubdiv1[4] = { max(0,childSubdiv[1]-1),childSubdiv[2],childSubdiv[0],max(0,childSubdiv[0]-1) };
	const int neighborSubdiv2[4] = { max(0,childSubdiv[2]-1),childSubdiv[0],childSubdiv[1],max(0,childSubdiv[1]-1) };
	subdivide(patches[0],depth+1, uv0, neighborSubdiv0);
	subdivide(patches[1],depth+1, uv1, neighborSubdiv1);
	subdivide(patches[2],depth+1, uv2, neighborSubdiv2);
      } 

      /* parametrization for quads */
      else if (N == 4) {
	const Vec2f uv_0(0.0f,0.0f);
	const Vec2f uv01(0.5f,0.0f);
	const Vec2f uv_1(1.0f,0.0f);
	const Vec2f uv12(1.0f,0.5f);
	const Vec2f uv_2(1.0f,1.0f);
	const Vec2f uv23(0.5f,1.0f);
	const Vec2f uv_3(0.0f,1.0f);
	const Vec2f uv30(0.0f,0.5f);
	const Vec2f uvcc(0.5f,0.5f);
	const Vec2f uv0[4] = { uv_0,uv01,uvcc,uv30 };
	const Vec2f uv1[4] = { uv_1,uv12,uvcc,uv01 };
	const Vec2f uv2[4] = { uv_2,uv23,uvcc,uv12 };
	const Vec2f uv3[4] = { uv_3,uv30,uvcc,uv23 };
	const int neighborSubdiv0[4] = { max(0,childSubdiv[0]-1),childSubdiv[1],childSubdiv[3],max(0,childSubdiv[3]-1) };
	const int neighborSubdiv1[4] = { max(0,childSubdiv[1]-1),childSubdiv[2],childSubdiv[0],max(0,childSubdiv[0]-1) };
	const int neighborSubdiv2[4] = { max(0,childSubdiv[2]-1),childSubdiv[3],childSubdiv[1],max(0,childSubdiv[1]-1) };
	const int neighborSubdiv3[4] = { max(0,childSubdiv[3]-1),childSubdiv[0],childSubdiv[2],max(0,childSubdiv[2]-1) };
	subdivide(patches[0],depth+1, uv0, neighborSubdiv0);
	subdivide(patches[1],depth+1, uv1, neighborSubdiv1);
	subdivide(patches[2],depth+1, uv2, neighborSubdiv2);
	subdivide(patches[3],depth+1, uv3, neighborSubdiv3);
      } 

      /* parametrization for arbitrary polygons */
      else 
      {
	for (size_t i=0; i<N; i++) 
	{
	  const Vec2f uv[4] = { Vec2f(float(i)+0.0f,0.0f),Vec2f(float(i)+0.0f,1.0f),Vec2f(float(i)+1.0f,1.0f),Vec2f(float(i)+1.0f,0.0f) };
	  const int neighborSubdiv[4] = { false,childSubdiv[(i+1)%N],childSubdiv[(i-1)%N],false };
	  subdivide(patches[i],depth+1,uv,neighborSubdiv);
	  
	}
      }
    }

    __forceinline void tessellate(const CatmullClarkPatch& patch, int depth, const Vec2f uv[4], const int neighborSubdiv_i[4]) {
      tessellator(patch,uv,neighborSubdiv_i);
    }

    void subdivide(const CatmullClarkPatch& patch, int depth, const Vec2f uv[4], const int neighborSubdiv_i[4])
    {
      if (depth == 0)
	if (patch.isRegularOrFinal(depth))
	  return tessellator(patch,uv,neighborSubdiv_i);

      CatmullClarkPatch patches[4]; 
      patch.subdivide(patches);

      int neighborSubdiv[4];
      for (size_t i=0; i<4; i++) 
	neighborSubdiv[i] = max(0,neighborSubdiv_i[i]-1);

      const int childSubdiv0 = !patches[0].isRegularOrFinal(depth);
      const int childSubdiv1 = !patches[1].isRegularOrFinal(depth);
      const int childSubdiv2 = !patches[2].isRegularOrFinal(depth);
      const int childSubdiv3 = !patches[3].isRegularOrFinal(depth);

      const Vec2f uv01 = 0.5f*(uv[0]+uv[1]);
      const Vec2f uv12 = 0.5f*(uv[1]+uv[2]);
      const Vec2f uv23 = 0.5f*(uv[2]+uv[3]);
      const Vec2f uv30 = 0.5f*(uv[3]+uv[0]);
      const Vec2f uvcc = 0.25f*(uv[0]+uv[1]+uv[2]+uv[3]);

      const Vec2f uv0[4] = { uv[0],uv01,uvcc,uv30 };
      const Vec2f uv1[4] = { uv01,uv[1],uv12,uvcc };
      const Vec2f uv2[4] = { uvcc,uv12,uv[2],uv23 };
      const Vec2f uv3[4] = { uv30,uvcc,uv23,uv[3] };

      const int neighborSubdiv0[4] = { neighborSubdiv[0],childSubdiv1,childSubdiv3,neighborSubdiv[3] };
      const int neighborSubdiv1[4] = { neighborSubdiv[0],neighborSubdiv[1],childSubdiv2,childSubdiv0 };
      const int neighborSubdiv2[4] = { childSubdiv1,neighborSubdiv[1],neighborSubdiv[2],childSubdiv3 };
      const int neighborSubdiv3[4] = { childSubdiv0,childSubdiv2,neighborSubdiv[2],neighborSubdiv[3] };

      if (childSubdiv0) subdivide  (patches[0], depth+1, uv0, neighborSubdiv0);
      else              tessellate (patches[0], depth+1, uv0, neighborSubdiv0);
      
      if (childSubdiv1) subdivide  (patches[1], depth+1, uv1, neighborSubdiv1);
      else              tessellate (patches[1], depth+1, uv1, neighborSubdiv1);
      
      if (childSubdiv2) subdivide  (patches[2], depth+1, uv2, neighborSubdiv2);
      else              tessellate (patches[2], depth+1, uv2, neighborSubdiv2);
      
      if (childSubdiv3) subdivide  (patches[3], depth+1, uv3, neighborSubdiv3);
      else              tessellate (patches[3], depth+1, uv3, neighborSubdiv3);
    }
  };

   template<typename Tessellator>
     inline void feature_adaptive_subdivision_bspline (int primID, const SubdivMesh::HalfEdge* h, const BufferT<Vec3fa>& vertices, Tessellator tessellator)
   {
     FeatureAdaptiveSubdivisionBSpline<Tessellator>(primID,h,vertices,tessellator);
   }
}
