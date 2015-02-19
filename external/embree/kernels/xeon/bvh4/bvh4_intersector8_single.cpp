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

#include "bvh4_intersector8_single.h"
#include "bvh4_intersector1.h"

#include "geometry/bezier1v_intersector8.h"
#include "geometry/bezier1i_intersector8.h"
#include "geometry/subdivpatch1_intersector1.h"
#include "geometry/subdivpatch1cached_intersector1.h"
#include "geometry/grid_intersector1.h"

namespace embree
{
  namespace isa
  {
    template<int types, bool robust, typename PrimitiveIntersector8>
    void BVH4Intersector8Single<types,robust,PrimitiveIntersector8>::intersect(avxb* valid_i, BVH4* bvh, Ray8& ray)
    {
      /* load ray */
      const avxb valid0 = *valid_i;
      avx3f ray_org = ray.org;
      avx3f ray_dir = ray.dir;
      avxf ray_tnear = ray.tnear, ray_tfar  = ray.tfar;
      const avx3f rdir = rcp_safe(ray_dir);
      const avx3f org(ray_org), org_rdir = org * rdir;
      ray_tnear = select(valid0,ray_tnear,avxf(pos_inf));
      ray_tfar  = select(valid0,ray_tfar ,avxf(neg_inf));
      const avxf inf = avxf(pos_inf);
      Precalculations pre(valid0,ray);

      /* compute near/far per ray */
      avx3i nearXYZ;
      nearXYZ.x = select(rdir.x >= 0.0f,avxi(0*(int)sizeof(ssef)),avxi(1*(int)sizeof(ssef)));
      nearXYZ.y = select(rdir.y >= 0.0f,avxi(2*(int)sizeof(ssef)),avxi(3*(int)sizeof(ssef)));
      nearXYZ.z = select(rdir.z >= 0.0f,avxi(4*(int)sizeof(ssef)),avxi(5*(int)sizeof(ssef)));

      /* we have no packet implementation for OBB nodes yet */
      size_t bits = movemask(valid0);
      for (size_t i=__bsf(bits); bits!=0; bits=__btc(bits,i), i=__bsf(bits)) {
	intersect1(bvh, bvh->root, i, pre, ray, ray_org, ray_dir, rdir, ray_tnear, ray_tfar, nearXYZ);
      }
      AVX_ZERO_UPPER();
    }
    
    template<int types, bool robust, typename PrimitiveIntersector8>
    void BVH4Intersector8Single<types,robust, PrimitiveIntersector8>::occluded(avxb* valid_i, BVH4* bvh, Ray8& ray)
    {
      /* load ray */
      const avxb valid = *valid_i;
      avxb terminated = !valid;
      avx3f ray_org = ray.org, ray_dir = ray.dir;
      avxf ray_tnear = ray.tnear, ray_tfar  = ray.tfar;
      const avx3f rdir = rcp_safe(ray_dir);
      const avx3f org(ray_org), org_rdir = org * rdir;
      ray_tnear = select(valid,ray_tnear,avxf(pos_inf));
      ray_tfar  = select(valid,ray_tfar ,avxf(neg_inf));
      const avxf inf = avxf(pos_inf);
      Precalculations pre(valid,ray);

      /* compute near/far per ray */
      avx3i nearXYZ;
      nearXYZ.x = select(rdir.x >= 0.0f,avxi(0*(int)sizeof(ssef)),avxi(1*(int)sizeof(ssef)));
      nearXYZ.y = select(rdir.y >= 0.0f,avxi(2*(int)sizeof(ssef)),avxi(3*(int)sizeof(ssef)));
      nearXYZ.z = select(rdir.z >= 0.0f,avxi(4*(int)sizeof(ssef)),avxi(5*(int)sizeof(ssef)));

      /* we have no packet implementation for OBB nodes yet */
      size_t bits = movemask(valid);
      for (size_t i=__bsf(bits); bits!=0; bits=__btc(bits,i), i=__bsf(bits)) {
	if (occluded1(bvh,bvh->root,i,pre,ray,ray_org,ray_dir,rdir,ray_tnear,ray_tfar,nearXYZ))
	  terminated[i] = -1;
      }
      store8i(valid & terminated,&ray.geomID,0);
      AVX_ZERO_UPPER();
    }

    template<typename Intersector1>
    void BVH4Intersector8FromIntersector1<Intersector1>::intersect(avxb* valid_i, BVH4* bvh, Ray8& ray)
    {
      Ray rays[8];
      ray.get(rays);
      size_t bits = movemask(*valid_i);
      for (size_t i=__bsf(bits); bits!=0; bits=__btc(bits,i), i=__bsf(bits)) {
	Intersector1::intersect(bvh,rays[i]);
      }
      ray.set(rays);
      AVX_ZERO_UPPER();
    }
    
    template<typename Intersector1>
    void BVH4Intersector8FromIntersector1<Intersector1>::occluded(avxb* valid_i, BVH4* bvh, Ray8& ray)
    {
      Ray rays[8];
      ray.get(rays);
      size_t bits = movemask(*valid_i);
      for (size_t i=__bsf(bits); bits!=0; bits=__btc(bits,i), i=__bsf(bits)) {
	Intersector1::intersect(bvh,rays[i]);
      }
      ray.set(rays);
      AVX_ZERO_UPPER();
    }
    
    DEFINE_INTERSECTOR8(BVH4Bezier1vIntersector8Single_OBB, BVH4Intersector8Single<0x101 COMMA false COMMA LeafIterator8_1<Bezier1vIntersector8<LeafMode> > >);
    DEFINE_INTERSECTOR8(BVH4Bezier1iIntersector8Single_OBB, BVH4Intersector8Single<0x101 COMMA false COMMA LeafIterator8_1<Bezier1iIntersector8<LeafMode> > >);
    DEFINE_INTERSECTOR8(BVH4Bezier1iMBIntersector8Single_OBB,BVH4Intersector8Single<0x1010 COMMA false COMMA LeafIterator8_1<Bezier1iIntersector8MB<LeafMode> > >);

    DEFINE_INTERSECTOR8(BVH4Subdivpatch1Intersector8, BVH4Intersector8FromIntersector1<BVH4Intersector1<0x1 COMMA false COMMA LeafIterator1<SubdivPatch1Intersector1 > > >);
    DEFINE_INTERSECTOR8(BVH4Subdivpatch1CachedIntersector8,BVH4Intersector8FromIntersector1<BVH4Intersector1<0x1 COMMA false COMMA SubdivPatch1CachedIntersector1> >);

    DEFINE_INTERSECTOR8(BVH4GridIntersector8, BVH4Intersector8FromIntersector1<BVH4Intersector1<0x1 COMMA false COMMA GridIntersector1> >);
    DEFINE_INTERSECTOR8(BVH4GridLazyIntersector8, BVH4Intersector8FromIntersector1<BVH4Intersector1<0x1 COMMA false COMMA Switch2Intersector1<GridIntersector1 COMMA GridLazyIntersector1> > >);
  }
}
