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

#include "primitive.h"
#include "discrete_tessellation.h"
#include "common/subdiv/feature_adaptive_eval.h"

#define GRID_COMPRESS_BOUNDS 1

namespace embree
{
  struct Grid
  {
    struct CompressedBounds16;
    struct UncompressedBounds16;

#if GRID_COMPRESS_BOUNDS
    typedef CompressedBounds16 Bounds16;
#else
    typedef UncompressedBounds16 Bounds16;
#endif

    struct UncompressedBounds16
    {
      static const size_t N = 4;

      /*! Sets bounding box of child. */
      __forceinline void clear() 
      {
	for (size_t i=0; i<4; i++) {
	  lower_x[i] = lower_y[i] = lower_z[i] = pos_inf;
	  upper_x[i] = upper_y[i] = upper_z[i] = neg_inf;
	}
      }

      /*! set bounds of all children (only required for compression) */
      __forceinline void set(const BBox3fa& bounds) {
      }
      
      /*! Sets bounding box of child. */
      __forceinline void set(size_t i, const BBox3fa& bounds) 
      {
        assert(i < 4*N);
        ((float*)&lower_x)[i] = bounds.lower.x; 
	((float*)&lower_y)[i] = bounds.lower.y; 
	((float*)&lower_z)[i] = bounds.lower.z;
        ((float*)&upper_x)[i] = bounds.upper.x; 
	((float*)&upper_y)[i] = bounds.upper.y; 
	((float*)&upper_z)[i] = bounds.upper.z;
      }
      
      /*! intersection with single rays */
      template<bool robust>
      __forceinline size_t intersect(size_t i, size_t _nearX, size_t _nearY, size_t _nearZ,
                                     const sse3f& org, const sse3f& rdir, const sse3f& org_rdir, const ssef& tnear, const ssef& tfar) const
      {
        const size_t nearX = 4*_nearX, nearY = 4*_nearY, nearZ = 4*_nearZ; 
        const size_t farX  = nearX ^ (4*sizeof(ssef)), farY  = nearY ^ (4*sizeof(ssef)), farZ  = nearZ ^ (4*sizeof(ssef));
#if defined (__AVX2__)
        const ssef tNearX = msub(load4f((const char*)&lower_x[i]+nearX), rdir.x, org_rdir.x);
        const ssef tNearY = msub(load4f((const char*)&lower_x[i]+nearY), rdir.y, org_rdir.y);
        const ssef tNearZ = msub(load4f((const char*)&lower_x[i]+nearZ), rdir.z, org_rdir.z);
        const ssef tFarX  = msub(load4f((const char*)&lower_x[i]+farX ), rdir.x, org_rdir.x);
        const ssef tFarY  = msub(load4f((const char*)&lower_x[i]+farY ), rdir.y, org_rdir.y);
        const ssef tFarZ  = msub(load4f((const char*)&lower_x[i]+farZ ), rdir.z, org_rdir.z);
#else
        const ssef tNearX = (load4f((const char*)&lower_x[i]+nearX) - org.x) * rdir.x;
        const ssef tNearY = (load4f((const char*)&lower_x[i]+nearY) - org.y) * rdir.y;
        const ssef tNearZ = (load4f((const char*)&lower_x[i]+nearZ) - org.z) * rdir.z;
        const ssef tFarX  = (load4f((const char*)&lower_x[i]+farX ) - org.x) * rdir.x;
        const ssef tFarY  = (load4f((const char*)&lower_x[i]+farY ) - org.y) * rdir.y;
        const ssef tFarZ  = (load4f((const char*)&lower_x[i]+farZ ) - org.z) * rdir.z;
#endif

        if (robust) {
          const float round_down = 1.0f-2.0f*float(ulp);
          const float round_up   = 1.0f+2.0f*float(ulp);
          const ssef tNear = max(tNearX,tNearY,tNearZ,tnear);
          const ssef tFar  = min(tFarX ,tFarY ,tFarZ ,tfar);
          const sseb vmask = round_down*tNear <= round_up*tFar;
          const size_t mask = movemask(vmask);
          return mask;
        }
        
#if defined(__SSE4_1__)
        const ssef tNear = maxi(maxi(tNearX,tNearY),maxi(tNearZ,tnear));
        const ssef tFar  = mini(mini(tFarX ,tFarY ),mini(tFarZ ,tfar ));
        const sseb vmask = cast(tNear) > cast(tFar);
        const size_t mask = movemask(vmask)^0xf;
#else
        const ssef tNear = max(tNearX,tNearY,tNearZ,tnear);
        const ssef tFar  = min(tFarX ,tFarY ,tFarZ ,tfar);
        const sseb vmask = tNear <= tFar;
        const size_t mask = movemask(vmask);
#endif
        return mask;
      }

      template<bool robust>
      __forceinline size_t intersect(size_t nearX, size_t nearY, size_t nearZ,
                                     const sse3f& org, const sse3f& rdir, const sse3f& org_rdir, const ssef& tnear, const ssef& tfar) const
      {
        const size_t mask0 = intersect<robust>(0, nearX, nearY, nearZ, org, rdir, org_rdir, tnear, tfar);
        const size_t mask1 = intersect<robust>(1, nearX, nearY, nearZ, org, rdir, org_rdir, tnear, tfar);
        const size_t mask2 = intersect<robust>(2, nearX, nearY, nearZ, org, rdir, org_rdir, tnear, tfar);
        const size_t mask3 = intersect<robust>(3, nearX, nearY, nearZ, org, rdir, org_rdir, tnear, tfar);
        return mask0 | (mask1 << 4) | (mask2 << 8) | (mask3 << 12);
      }

#if defined (__AVX__)

      /*! intersection with single rays */
      template<bool robust>
      __forceinline size_t intersect(size_t i, size_t _nearX, size_t _nearY, size_t _nearZ,
                                     const avx3f& org, const avx3f& rdir, const avx3f& org_rdir, const avxf& tnear, const avxf& tfar) const
      {
        const size_t nearX = 4*_nearX, nearY = 4*_nearY, nearZ = 4*_nearZ; 
        const size_t farX  = nearX ^ (4*sizeof(ssef)), farY  = nearY ^ (4*sizeof(ssef)), farZ  = nearZ ^ (4*sizeof(ssef));
#if defined (__AVX2__)
        const avxf tNearX = msub(load8f((const char*)&lower_x[i]+nearX), rdir.x, org_rdir.x);
        const avxf tNearY = msub(load8f((const char*)&lower_x[i]+nearY), rdir.y, org_rdir.y);
        const avxf tNearZ = msub(load8f((const char*)&lower_x[i]+nearZ), rdir.z, org_rdir.z);
        const avxf tFarX  = msub(load8f((const char*)&lower_x[i]+farX ), rdir.x, org_rdir.x);
        const avxf tFarY  = msub(load8f((const char*)&lower_x[i]+farY ), rdir.y, org_rdir.y);
        const avxf tFarZ  = msub(load8f((const char*)&lower_x[i]+farZ ), rdir.z, org_rdir.z);
#else
        const avxf tNearX = (load8f((const char*)&lower_x[i]+nearX) - org.x) * rdir.x;
        const avxf tNearY = (load8f((const char*)&lower_x[i]+nearY) - org.y) * rdir.y;
        const avxf tNearZ = (load8f((const char*)&lower_x[i]+nearZ) - org.z) * rdir.z;
        const avxf tFarX  = (load8f((const char*)&lower_x[i]+farX ) - org.x) * rdir.x;
        const avxf tFarY  = (load8f((const char*)&lower_x[i]+farY ) - org.y) * rdir.y;
        const avxf tFarZ  = (load8f((const char*)&lower_x[i]+farZ ) - org.z) * rdir.z;
#endif

        if (robust) {
          const float round_down = 1.0f-2.0f*float(ulp);
          const float round_up   = 1.0f+2.0f*float(ulp);
          const avxf tNear = max(tNearX,tNearY,tNearZ,tnear);
          const avxf tFar  = min(tFarX ,tFarY ,tFarZ ,tfar);
          const avxb vmask = round_down*tNear <= round_up*tFar;
          const size_t mask = movemask(vmask);
          return mask;
        }
        
/*#if defined(__AVX2__) // FIXME: not working for cube
        const avxf tNear = maxi(maxi(tNearX,tNearY),maxi(tNearZ,tnear));
        const avxf tFar  = mini(mini(tFarX ,tFarY ),mini(tFarZ ,tfar ));
        const avxb vmask = cast(tNear) > cast(tFar);
        const size_t mask = movemask(vmask)^0xf;
	#else*/
        const avxf tNear = max(tNearX,tNearY,tNearZ,tnear);
        const avxf tFar  = min(tFarX ,tFarY ,tFarZ ,tfar);
        const avxb vmask = tNear <= tFar;
        const size_t mask = movemask(vmask);
//#endif
        return mask;
      }

      template<bool robust>
      __forceinline size_t intersect(size_t nearX, size_t nearY, size_t nearZ,
                                     const avx3f& org, const avx3f& rdir, const avx3f& org_rdir, const avxf& tnear, const avxf& tfar) const
      {
        const size_t mask01 = intersect<robust>(0, nearX, nearY, nearZ, org, rdir, org_rdir, tnear, tfar);
        const size_t mask23 = intersect<robust>(2, nearX, nearY, nearZ, org, rdir, org_rdir, tnear, tfar);
        return mask01 | (mask23 << 8);
      }
      
#endif

    public:
      ssef lower_x[4];           //!< X dimension of lower bounds of all 4 children.
      ssef upper_x[4];           //!< X dimension of upper bounds of all 4 children.
      ssef lower_y[4];           //!< Y dimension of lower bounds of all 4 children.
      ssef upper_y[4];           //!< Y dimension of upper bounds of all 4 children.
      ssef lower_z[4];           //!< Z dimension of lower bounds of all 4 children.
      ssef upper_z[4];           //!< Z dimension of upper bounds of all 4 children.
    };

    struct CompressedBounds16
    {
      static const size_t N = 4;

      /*! Sets bounding box of child. */
      __forceinline void clear() 
      {
	for (size_t i=0; i<16; i++) {
	  lower_x[i] = lower_y[i] = lower_z[i] = 255;
	  upper_x[i] = upper_y[i] = upper_z[i] = 0;
	}
      }

      __forceinline void set(const BBox3fa& bounds) {
        offset = bounds.lower;
        scale = max(Vec3fa(1E-20f),bounds.size()/255.0f);
      }

      /*! Sets bounding box of child. */
      __forceinline void set(size_t i, const BBox3fa& bounds) 
      {
	assert(i < 4*N);
        const Vec3fa lower = clamp(floor((bounds.lower-offset)/scale),Vec3fa(0.0f),Vec3fa(255.0f));
        const Vec3fa upper = clamp(ceil ((bounds.upper-offset)/scale),Vec3fa(0.0f),Vec3fa(255.0f));
        lower_x[i] = (unsigned char) lower.x;
        lower_y[i] = (unsigned char) lower.y;
        lower_z[i] = (unsigned char) lower.z;
        upper_x[i] = (unsigned char) upper.x;
        upper_y[i] = (unsigned char) upper.y;
        upper_z[i] = (unsigned char) upper.z;
      }
      
      /*! intersection with single rays */
      template<bool robust>
      __forceinline size_t intersect(size_t i, size_t nearX, size_t nearY, size_t nearZ,
                                     const sse3f& org, const sse3f& rdir, const sse3f& org_rdir, const ssef& tnear, const ssef& tfar) const
      {
        const size_t farX  = nearX ^ 16, farY  = nearY ^ 16, farZ  = nearZ ^ 16;

        const sse3f vscale(scale), voffset(offset);
        const ssef near_x = madd(ssef::load(&this->lower_x[i]+nearX),vscale.x,voffset.x);
        const ssef near_y = madd(ssef::load(&this->lower_x[i]+nearY),vscale.y,voffset.y);
        const ssef near_z = madd(ssef::load(&this->lower_x[i]+nearZ),vscale.z,voffset.z);
        const ssef far_x  = madd(ssef::load(&this->lower_x[i]+farX),vscale.x,voffset.x);
        const ssef far_y  = madd(ssef::load(&this->lower_x[i]+farY),vscale.y,voffset.y);
        const ssef far_z  = madd(ssef::load(&this->lower_x[i]+farZ),vscale.z,voffset.z);

#if defined (__AVX2__)
        const ssef tNearX = msub(near_x, rdir.x, org_rdir.x);
        const ssef tNearY = msub(near_y, rdir.y, org_rdir.y);
        const ssef tNearZ = msub(near_z, rdir.z, org_rdir.z);
        const ssef tFarX  = msub(far_x , rdir.x, org_rdir.x);
        const ssef tFarY  = msub(far_y , rdir.y, org_rdir.y);
        const ssef tFarZ  = msub(far_z , rdir.z, org_rdir.z);
#else
        const ssef tNearX = (near_x - org.x) * rdir.x;
        const ssef tNearY = (near_y - org.y) * rdir.y;
        const ssef tNearZ = (near_z - org.z) * rdir.z;
        const ssef tFarX  = (far_x  - org.x) * rdir.x;
        const ssef tFarY  = (far_y  - org.y) * rdir.y;
        const ssef tFarZ  = (far_z  - org.z) * rdir.z;
#endif

        if (robust) {
          const float round_down = 1.0f-2.0f*float(ulp);
          const float round_up   = 1.0f+2.0f*float(ulp);
          const ssef tNear = max(tNearX,tNearY,tNearZ,tnear);
          const ssef tFar  = min(tFarX ,tFarY ,tFarZ ,tfar);
          const sseb vmask = round_down*tNear <= round_up*tFar;
          const size_t mask = movemask(vmask);
          return mask;
        }
        
#if defined(__SSE4_1__)
        const ssef tNear = maxi(maxi(tNearX,tNearY),maxi(tNearZ,tnear));
        const ssef tFar  = mini(mini(tFarX ,tFarY ),mini(tFarZ ,tfar ));
        const sseb vmask = cast(tNear) > cast(tFar);
        const size_t mask = movemask(vmask)^0xf;
#else
        const ssef tNear = max(tNearX,tNearY,tNearZ,tnear);
        const ssef tFar  = min(tFarX ,tFarY ,tFarZ ,tfar);
        const sseb vmask = tNear <= tFar;
        const size_t mask = movemask(vmask);
#endif
        return mask;
      }

      template<bool robust>
      __forceinline size_t intersect(size_t nearX, size_t nearY, size_t nearZ,
                                     const sse3f& org, const sse3f& rdir, const sse3f& org_rdir, const ssef& tnear, const ssef& tfar) const
      {
        const size_t mask0 = intersect<robust>( 0, nearX, nearY, nearZ, org, rdir, org_rdir, tnear, tfar);
        const size_t mask1 = intersect<robust>( 4, nearX, nearY, nearZ, org, rdir, org_rdir, tnear, tfar);
        const size_t mask2 = intersect<robust>( 8, nearX, nearY, nearZ, org, rdir, org_rdir, tnear, tfar);
        const size_t mask3 = intersect<robust>(12, nearX, nearY, nearZ, org, rdir, org_rdir, tnear, tfar);
        return mask0 | (mask1 << 4) | (mask2 << 8) | (mask3 << 12);
      }

#if defined (__AVX__)

      /*! intersection with single rays */
      template<bool robust>
      __forceinline size_t intersect(size_t i, size_t nearX, size_t nearY, size_t nearZ,
                                     const avx3f& org, const avx3f& rdir, const avx3f& org_rdir, const avxf& tnear, const avxf& tfar) const
      {
        const size_t farX  = nearX ^ 16, farY  = nearY ^ 16, farZ  = nearZ ^ 16;

        const avx3f vscale(scale), voffset(offset);
        const avxf near_x = madd(avxf::load(&this->lower_x[i]+nearX),vscale.x,voffset.x);
        const avxf near_y = madd(avxf::load(&this->lower_x[i]+nearY),vscale.y,voffset.y);
        const avxf near_z = madd(avxf::load(&this->lower_x[i]+nearZ),vscale.z,voffset.z);
        const avxf far_x  = madd(avxf::load(&this->lower_x[i]+farX),vscale.x,voffset.x);
        const avxf far_y  = madd(avxf::load(&this->lower_x[i]+farY),vscale.y,voffset.y);
        const avxf far_z  = madd(avxf::load(&this->lower_x[i]+farZ),vscale.z,voffset.z);

#if defined (__AVX2__)

        const avxf tNearX = msub(near_x, rdir.x, org_rdir.x);
        const avxf tNearY = msub(near_y, rdir.y, org_rdir.y);
        const avxf tNearZ = msub(near_z, rdir.z, org_rdir.z);
        const avxf tFarX  = msub(far_x , rdir.x, org_rdir.x);
        const avxf tFarY  = msub(far_y , rdir.y, org_rdir.y);
        const avxf tFarZ  = msub(far_z , rdir.z, org_rdir.z);
#else

        const avxf tNearX = (near_x - org.x) * rdir.x;
        const avxf tNearY = (near_y - org.y) * rdir.y;
        const avxf tNearZ = (near_z - org.z) * rdir.z;
        const avxf tFarX  = (far_x  - org.x) * rdir.x;
        const avxf tFarY  = (far_y  - org.y) * rdir.y;
        const avxf tFarZ  = (far_z  - org.z) * rdir.z;
#endif

        if (robust) {
          const float round_down = 1.0f-2.0f*float(ulp);
          const float round_up   = 1.0f+2.0f*float(ulp);
          const avxf tNear = max(tNearX,tNearY,tNearZ,tnear);
          const avxf tFar  = min(tFarX ,tFarY ,tFarZ ,tfar);
          const avxb vmask = round_down*tNear <= round_up*tFar;
          const size_t mask = movemask(vmask);
          return mask;
        }
        
/*#if defined(__AVX2__) // FIXME: not working for cube
        const avxf tNear = maxi(maxi(tNearX,tNearY),maxi(tNearZ,tnear));
        const avxf tFar  = mini(mini(tFarX ,tFarY ),mini(tFarZ ,tfar ));
        const avxb vmask = cast(tNear) > cast(tFar);
        const size_t mask = movemask(vmask)^0xf;
	#else*/
        const avxf tNear = max(tNearX,tNearY,tNearZ,tnear);
        const avxf tFar  = min(tFarX ,tFarY ,tFarZ ,tfar);
        const avxb vmask = tNear <= tFar;
        const size_t mask = movemask(vmask);
//#endif
        return mask;
      }

      template<bool robust>
      __forceinline size_t intersect(size_t nearX, size_t nearY, size_t nearZ,
                                     const avx3f& org, const avx3f& rdir, const avx3f& org_rdir, const avxf& tnear, const avxf& tfar) const
      {
        const size_t mask01 = intersect<robust>(0, nearX, nearY, nearZ, org, rdir, org_rdir, tnear, tfar);
        const size_t mask23 = intersect<robust>(8, nearX, nearY, nearZ, org, rdir, org_rdir, tnear, tfar);
        return mask01 | (mask23 << 8);
      }
      
#endif

    public:
      Vec3fa offset;               //!< offset to decompress bounds
      Vec3fa scale;                //!< scale  to decompress bounds
      unsigned char lower_x[16]; 
      unsigned char upper_x[16]; 
      unsigned char lower_y[16]; 
      unsigned char upper_y[16]; 
      unsigned char lower_z[16]; 
      unsigned char upper_z[16]; 
    };

//#endif

    struct EagerLeaf
    {
      struct Quads
      {
	enum { QUAD2X2 = 0, QUAD1X2 = 1, QUAD2X1 = 2, QUAD1X1 = 3, NONE = 4 };
	
	__forceinline Quads () : type(NONE), ofs(0) {}
	__forceinline Quads (unsigned char type, unsigned char ofs) : type(type), ofs(ofs) {}

	friend std::ostream& operator<<(std::ostream& cout, const Quads& a) {
	  return cout << "{ type = " << (int) a.type << ", ofs = " << (int) a.ofs << " }";
	}

	unsigned char type;
	unsigned char ofs;
      };

      __forceinline const BBox3fa getBounds(const size_t x0, const size_t x1, const size_t y0, const size_t y1) const 
      {
	BBox3fa bounds = empty;
	for (size_t y=y0; y<=y1; y++)
	  for (size_t x=x0; x<=x1; x++) 
	    bounds.extend(grid.point(x,y));

	return bounds;
      }

      __forceinline EagerLeaf (const Grid& grid) 
	: grid(grid) {}

      __forceinline const BBox3fa init (size_t x0, size_t x1, size_t y0, size_t y1) 
      {
	BBox3fa box_list[16];

	BBox3fa box = empty;
	bounds.clear();
	size_t i=0;
	for (size_t y=y0; y<y1; y+=2)
	{
	  for (size_t x=x0; x<x1; x+=2)
	  {
            assert(i < 16);
	    const bool right = x+1 == x1;
	    const bool bottom = y+1 == y1;
	    const size_t kx0 = x, kx1 = min(x+2,x1);
	    const size_t ky0 = y, ky1 = min(y+2,y1);
	    new (&quads[i]) Quads(2*bottom+right,y*grid.width+x);
	    const BBox3fa b = getBounds(kx0,kx1,ky0,ky1);
	    box_list[i++] = b;
	    box.extend(b);
	  }
	}

	bounds.set(box);
	for (size_t j=0; j<i; j++)
	  bounds.set(j,box_list[j]);

	return box;
      }

      friend std::ostream& operator<<(std::ostream& cout, const EagerLeaf& a) {
	cout << "{ " << std::endl;
	//cout << "  bounds = " << a.bounds << std::endl;
	for (size_t i=0; i<16; i++) cout << "  quads[" << i << "] = " << a.quads[i] << ", " << std::endl;
	cout << "  grid = " << &a.grid << std::endl;
	return cout << "}";
      }

      Bounds16 bounds;
      Quads quads[16];
      const Grid& grid;
    };

  public:

    __forceinline Grid(unsigned width, unsigned height, unsigned geomID, unsigned primID)
      : width(width), height(height), geomID(geomID), primID(primID) 
    {
      assert(width <= 17);
      assert(height <= 17);
    }

    static __forceinline Grid* create(FastAllocator::Thread& alloc, const size_t width, const size_t height, const unsigned geomID, const unsigned primID) {
      return new (alloc.malloc(sizeof(Grid)-17*17*sizeof(Vec3fa)+width*height*sizeof(Vec3fa))) Grid(width,height,geomID,primID);
    }
    
    __forceinline       Vec3fa& point(const size_t x, const size_t y)       { assert(y*width+x < width*height); return P[y*width+x]; }
    __forceinline const Vec3fa& point(const size_t x, const size_t y) const { assert(y*width+x < width*height); return P[y*width+x]; }

    __forceinline BBox3fa bounds() const
    {
      BBox3fa bounds = empty;
      for (size_t y=0; y<height; y++)
	for (size_t x=0; x<width; x++)
	  bounds.extend(point(x,y));
      return bounds;
    }

    static size_t getNumEagerLeaves(size_t width, size_t height) {
      const size_t w = (((width +1)/2)+3)/4;
      const size_t h = (((height+1)/2)+3)/4;
      return w*h;
    }

    static size_t getNumLazyLeaves(size_t width, size_t height) {
      const size_t w = (width +15)/16;
      const size_t h = (height+15)/16;
      return w*h;
    }

    int stitch(const int x, const int fine, const int coarse) {
      return (2*x+1)*coarse/(2*fine);
    }
    
    __forceinline void displace(Scene* scene, RTCDisplacementFunc func, void* userPtr, const Vec2f* luv, const Vec2f* uv, const Vec3fa* Ng)
    {
      /* calculate uv coordinates */
      __aligned(64) float qu[17*17+16], qv[17*17+16];
      __aligned(64) float qx[17*17+16], qy[17*17+16], qz[17*17+16];
      __aligned(64) float nx[17*17+16], ny[17*17+16], nz[17*17+16];
      for (size_t y=0; y<height; y++) 
      {
        for (size_t x=0; x<width; x++) 
        {
          qu[y*width+x] = uv[y*width+x].x;
          qv[y*width+x] = uv[y*width+x].y;
          qx[y*width+x] = P[y*width+x].x;
          qy[y*width+x] = P[y*width+x].y;
          qz[y*width+x] = P[y*width+x].z;
          nx[y*width+x] = Ng[y*width+x].x;
          ny[y*width+x] = Ng[y*width+x].y;
          nz[y*width+x] = Ng[y*width+x].z;
        }
      }
      
      /* call displacement shader */
      func(userPtr,geomID,primID,
	   (float*)qu,(float*)qv,
	   (float*)nx,(float*)ny,(float*)nz,
	   (float*)qx,(float*)qy,(float*)qz,
	   width*height);
      
      /* add displacements */
      for (size_t y=0; y<height; y++) {
        for (size_t x=0; x<width; x++) {
          const Vec3fa P0 = P[y*width+x];
          const Vec3fa P1 = Vec3fa(qx[y*width+x],qy[y*width+x],qz[y*width+x]);
#if defined(DEBUG)
	  SubdivMesh* mesh = (SubdivMesh*) scene->get(geomID);
          if (!mesh->displBounds.empty() && !inside(mesh->displBounds,P1-P0))
            THROW_RUNTIME_ERROR("displacement out of bounds");
#endif
          P[y*width+x] = P1;
        }
      }
    }

    __forceinline bool stitch_x(const CatmullClarkPatch& patch, const size_t y0, const size_t y_ofs, const size_t border, const size_t x0, const size_t x1,
				const DiscreteTessellationPattern& fine, const DiscreteTessellationPattern& coarse, 
				Vec2f luv[17*17], Vec3fa Ng[17*17])
    {
      if (unlikely(y0 != border || fine.size() == coarse.size()))
	return false;

      const size_t x0s = stitch(x0,fine.size(),coarse.size());
      const size_t x1s = stitch(x1,fine.size(),coarse.size());
      assert(x1s-x0s < 17);
      
      Vec3fa p_y0[17], Ng_y0[17];
      feature_adaptive_eval (patch, y0!=0,y0!=0,x0s,x1s,2,coarse.size()+1, p_y0,Ng_y0,1,17);

      Vec2f luv_y0[17];
      int y = y0-y_ofs;
      for (int x=0; x<=x1-x0; x++) {
	const size_t xs = stitch(x0+x,fine.size(),coarse.size());
	const float fx = coarse(xs);
	luv[x*width+y].y = fx;
	P  [x*width+y] = p_y0 [xs-x0s];
	Ng [x*width+y] = Ng_y0[xs-x0s];
      }
      return true;
    }

    __forceinline bool stitch_y(const CatmullClarkPatch& patch, const size_t y0, const size_t y_ofs, const size_t border, const size_t x0, const size_t x1,
				const DiscreteTessellationPattern& fine, const DiscreteTessellationPattern& coarse, 
				Vec2f luv[17*17], Vec3fa Ng[17*17])
    {
      if (unlikely(y0 != border || fine.size() == coarse.size()))
	return false;

      const size_t x0s = stitch(x0,fine.size(),coarse.size());
      const size_t x1s = stitch(x1,fine.size(),coarse.size());
      assert(x1s-x0s < 17);
      
      Vec3fa p_y0[17], Ng_y0[17];
      feature_adaptive_eval (patch, x0s,x1s,y0!=0,y0!=0,coarse.size()+1,2, p_y0,Ng_y0,17,1);
      
      Vec2f luv_y0[17];
      int y = y0-y_ofs;
      for (int x=0; x<=x1-x0; x++) {
	const size_t xs = stitch(x0+x,fine.size(),coarse.size());
	const float fx = coarse(xs);
	luv[y*width+x].x = fx;
	P  [y*width+x] = p_y0 [xs-x0s];
	Ng [y*width+x] = Ng_y0[xs-x0s];
      }
      return true;
    }

    __forceinline void calculateLocalUVs(const size_t x0, const size_t x1,
					 const size_t y0, const size_t y1,
					 const DiscreteTessellationPattern& pattern_x,
					 const DiscreteTessellationPattern& pattern_y,
					 Vec2f luv[17*17])
    {
      for (int y=0; y<height; y++) {
        const float fy = pattern_y(y0+y);
        for (int x=0; x<width; x++) {
	  const float fx = pattern_x(x0+x);
          luv[y*width+x] = Vec2f(fx,fy);
        }
      }
    }

    __forceinline void stitchLocalUVs(const size_t x0, const size_t x1,
				      const size_t y0, const size_t y1,
				      const DiscreteTessellationPattern& pattern0, 
				      const DiscreteTessellationPattern& pattern1, 
				      const DiscreteTessellationPattern& pattern2, 
				      const DiscreteTessellationPattern& pattern3, 
				      const DiscreteTessellationPattern& pattern_x,
				      const DiscreteTessellationPattern& pattern_y,
				      Vec2f luv[17*17])
    {
      if (unlikely(y0 == 0 && pattern_x.size() != pattern0.size())) {
        const float fy = pattern_y(y0);
        for (int x=x0; x<=x1; x++) {
          const float fx = pattern0(stitch(x,pattern_x.size(),pattern0.size()));
	  assert((y0-y0)*width+(x-x0) < 17*17);
          luv[(y0-y0)*width+(x-x0)] = Vec2f(fx,fy);
        }
      }

      if (unlikely(y1 == pattern_y.size() && pattern_x.size() != pattern2.size())) {
	const float fy = pattern_y(y1);
	for (int x=x0; x<=x1; x++) {
	  const float fx = pattern2(stitch(x,pattern_x.size(),pattern2.size()));
	  assert((y1-y0)*width+(x-x0) < 17*17);
	  luv[(y1-y0)*width+(x-x0)] = Vec2f(fx,fy);
	}
      }

      if (unlikely(x0 == 0 && pattern_y.size() != pattern3.size())) {
        const float fx = pattern_x(x0);
        for (int y=y0; y<=y1; y++) {
          const float fy = pattern3(stitch(y,pattern_y.size(),pattern3.size()));
	  assert((y-y0)*width+(x0-x0) < 17*17);
          luv[(y-y0)*width+(x0-x0)] = Vec2f(fx,fy);
        }
      }

      if (unlikely(x1 == pattern_x.size() && pattern_y.size() != pattern1.size())) {
	const float fx = pattern_x(x1);
	for (int y=y0; y<=y1; y++) {
	  const float fy = pattern1(stitch(y,pattern_y.size(),pattern1.size()));
	  assert((y-y0)*width+(x1-x0) < 17*17);
	  luv[(y-y0)*width+(x1-x0)] = Vec2f(fx,fy);
	}
      }
    }

    __forceinline void calculateGlobalUVs(const Vec2f& uv0, const Vec2f& uv1, const Vec2f& uv2, const Vec2f& uv3, 
					  Vec2f luv[17*17], Vec2f guv[17*17])
    {
      for (int y=0; y<height; y++) {
        for (int x=0; x<width; x++) {
	  const float fx =  luv[y*width+x].x;
	  const float fy =  luv[y*width+x].y;
	  const Vec2f uv01 = (1.0f-fx) * uv0  + fx * uv1;
	  const Vec2f uv32 = (1.0f-fx) * uv3  + fx * uv2;
	  const Vec2f uvxy = (1.0f-fy) * uv01 + fy * uv32;
	  guv[y*width+x] = uvxy;
	  const int iu = clamp(uvxy.x * 0xFFFF, 0.0f, float(0xFFFF));
	  const int iv = clamp(uvxy.y * 0xFFFF, 0.0f, float(0xFFFF));
	  point(x,y).a = (iv << 16) | iu;
        }
      }
    }

    template<typename Patch>
    __forceinline BBox3fa calculatePositionAndNormal(const Patch& patch, Vec2f luv[17*17], Vec3fa Ng[17*17])
    {
      BBox3fa bounds = empty;
      for (int y=0; y<height; y++) {
        for (int x=0; x<width; x++) {
	  const Vec2f& uv = luv[y*width+x];
	  const Vec3fa P = patch.eval(uv.x,uv.y);
	  bounds.extend(P);
	  point(x,y)    = P;
	  Ng[y*width+x] = normalize_safe(patch.normal(luv[y*width+x].x, luv[y*width+x].y)); // FIXME: enable only for displacement mapping
        }
      }
      return bounds;
    }

    template<typename Patch>
    __forceinline void calculatePositionAndNormal(const Patch& patch, 
						  const size_t x0, const size_t x1,
						  const size_t y0, const size_t y1,
						  const DiscreteTessellationPattern& pattern0, 
						  const DiscreteTessellationPattern& pattern1, 
						  const DiscreteTessellationPattern& pattern2, 
						  const DiscreteTessellationPattern& pattern3, 
						  const DiscreteTessellationPattern& pattern_x,
						  const DiscreteTessellationPattern& pattern_y,
						  Vec2f luv[17*17], Vec3fa Ng[17*17])
    {
      /* evaluate position and normal */
      size_t swidth  = pattern_x.size()+1;
      size_t sheight = pattern_y.size()+1;
#if 0
      feature_adaptive_eval (patch, x0,x1,y0,y1, swidth,sheight, P,Ng,width,height);
#else
      const bool st = stitch_y(patch,y0,y0,0        ,x0,x1,pattern_x,pattern0,luv,Ng);
      const bool sr = stitch_x(patch,x1,x0,swidth-1 ,y0,y1,pattern_y,pattern1,luv,Ng);
      const bool sb = stitch_y(patch,y1,y0,sheight-1,x0,x1,pattern_x,pattern2,luv,Ng);
      const bool sl = stitch_x(patch,x0,x0,0        ,y0,y1,pattern_y,pattern3,luv,Ng);
      feature_adaptive_eval (patch, x0+sl,x1-sr,y0+st,y1-sb, swidth,sheight, P+st*width+sl,Ng+st*width+sl,width,height);
#endif
    }

   __forceinline size_t createEagerPrims(FastAllocator::Thread& alloc, PrimRef* prims, 
					 const size_t x0, const size_t x1,
					 const size_t y0, const size_t y1)
	{
	  size_t i=0;
	  for (size_t y=y0; y<y1; y+=8) {
	    for (size_t x=x0; x<x1; x+=8) {
	      const size_t rx0 = x-x0, rx1 = min(x+8,x1)-x0;
	      const size_t ry0 = y-y0, ry1 = min(y+8,y1)-y0;
	      EagerLeaf* leaf = new (alloc.malloc(sizeof(EagerLeaf))) EagerLeaf(*this);
	      const BBox3fa bounds = leaf->init(rx0,rx1,ry0,ry1);
	      prims[i++] = PrimRef(bounds,BVH4::encodeTypedLeaf(leaf,0));
	    }
	  }
	  return i;
	}


    template<typename Patch>
    __forceinline void build(Scene* scene, const Patch& patch,
			     const size_t x0, const size_t x1,
			     const size_t y0, const size_t y1,
			     const Vec2f& uv0, const Vec2f& uv1, const Vec2f& uv2, const Vec2f& uv3,
			     const DiscreteTessellationPattern& pattern0, 
			     const DiscreteTessellationPattern& pattern1, 
			     const DiscreteTessellationPattern& pattern2, 
			     const DiscreteTessellationPattern& pattern3, 
			     const DiscreteTessellationPattern& pattern_x,
			     const DiscreteTessellationPattern& pattern_y)
    {
      /* calculate local UVs */
      Vec2f luv[17*17]; 
      calculateLocalUVs(x0,x1,y0,y1,pattern_x,pattern_y,luv);

      /* stitch local UVs */
      stitchLocalUVs(x0,x1,y0,y1,pattern0,pattern1,pattern2,pattern3,pattern_x,pattern_y,luv);
      
      /* evaluate position and normal */
      Vec3fa Ng[17*17];
      calculatePositionAndNormal(patch,luv,Ng);

      /* calculate global UVs */
      Vec2f guv[17*17]; 
      calculateGlobalUVs(uv0,uv1,uv2,uv3,luv,guv);

      /* perform displacement */
      SubdivMesh* mesh = (SubdivMesh*) scene->get(geomID);
      if (mesh->displFunc) 
	displace(scene,mesh->displFunc,mesh->userPtr,luv,guv,Ng);
    }

    __forceinline void build(Scene* scene, const CatmullClarkPatch& patch,
			     const size_t x0, const size_t x1,
			     const size_t y0, const size_t y1,
			     const Vec2f& uv0, const Vec2f& uv1, const Vec2f& uv2, const Vec2f& uv3,
			     const DiscreteTessellationPattern& pattern0, 
			     const DiscreteTessellationPattern& pattern1, 
			     const DiscreteTessellationPattern& pattern2, 
			     const DiscreteTessellationPattern& pattern3, 
			     const DiscreteTessellationPattern& pattern_x,
			     const DiscreteTessellationPattern& pattern_y)
    {
      /* calculate local UVs */
      Vec2f luv[17*17]; 
      calculateLocalUVs(x0,x1,y0,y1,pattern_x,pattern_y,luv);

      /* evaluate position and normal */
      Vec3fa Ng[17*17];
      calculatePositionAndNormal(patch,x0,x1,y0,y1,pattern0,pattern1,pattern2,pattern3,pattern_x,pattern_y,luv,Ng);

      /* calculate global UVs */
      Vec2f guv[17*17]; 
      calculateGlobalUVs(uv0,uv1,uv2,uv3,luv,guv);

      /* perform displacement */
      SubdivMesh* mesh = (SubdivMesh*) scene->get(geomID);
      if (mesh->displFunc) 
	displace(scene,mesh->displFunc,mesh->userPtr,luv,guv,Ng);
    }

    __forceinline BBox3fa buildLazyBounds(Scene* scene, const CatmullClarkPatch& patch,
					  const size_t x0, const size_t x1,
					  const size_t y0, const size_t y1,
					  const Vec2f& uv0, const Vec2f& uv1, const Vec2f& uv2, const Vec2f& uv3,
					  const DiscreteTessellationPattern& pattern0, 
					  const DiscreteTessellationPattern& pattern1, 
					  const DiscreteTessellationPattern& pattern2, 
					  const DiscreteTessellationPattern& pattern3, 
					  const DiscreteTessellationPattern& pattern_x,
					  const DiscreteTessellationPattern& pattern_y)
    {
      /* calculate local UVs */
      Vec2f luv[17*17]; 
      calculateLocalUVs(x0,x1,y0,y1,pattern_x,pattern_y,luv);

      /* evaluate position and normal */
      Vec3fa Ng[17*17];
      calculatePositionAndNormal(patch,x0,x1,y0,y1,pattern0,pattern1,pattern2,pattern3,pattern_x,pattern_y,luv,Ng);

      /* calculate global UVs */
      Vec2f guv[17*17]; 
      calculateGlobalUVs(uv0,uv1,uv2,uv3,luv,guv);

      /* try to approximate bounding box */
      SubdivMesh* mesh = (SubdivMesh*) scene->get(geomID);
      const BBox3fa dbounds = mesh->displBounds;
      if (!dbounds.empty()) {
	const BBox3fa gbounds = bounds();
	if (all(gt_mask(8.0f*gbounds.size(),dbounds.size()))) {
	  return gbounds+dbounds;
	}
      }

      /* perform displacement */
      if (mesh->displFunc) 
	displace(scene,mesh->displFunc,mesh->userPtr,luv,guv,Ng);

      return bounds();
    }
    
    template<typename Patch>
    static size_t createEager(unsigned geomID, unsigned primID, 
			      Scene* scene, const Patch& patch,
			      FastAllocator::Thread& alloc, PrimRef* prims,
			      const size_t x0, const size_t x1,
			      const size_t y0, const size_t y1,
			      const Vec2f uv[4], 
			      const DiscreteTessellationPattern& pattern0, 
			      const DiscreteTessellationPattern& pattern1, 
			      const DiscreteTessellationPattern& pattern2, 
			      const DiscreteTessellationPattern& pattern3, 
			      const DiscreteTessellationPattern& pattern_x,
			      const DiscreteTessellationPattern& pattern_y)
    {
      size_t N = 0;
      for (size_t y=y0; y<y1; y+=16)
      {
	for (size_t x=x0; x<x1; x+=16) 
	{
	  const size_t lx0 = x, lx1 = min(lx0+16,x1);
	  const size_t ly0 = y, ly1 = min(ly0+16,y1);
	  Grid* leaf = Grid::create(alloc,lx1-lx0+1,ly1-ly0+1,geomID,primID);
	  leaf->build(scene,patch,lx0,lx1,ly0,ly1,uv[0],uv[1],uv[2],uv[3],pattern0,pattern1,pattern2,pattern3,pattern_x,pattern_y);
	  size_t n = leaf->createEagerPrims(alloc,prims,lx0,lx1,ly0,ly1);
	  prims += n;
	  N += n;
	}
      }
      return N;
    }

    std::pair<BBox3fa,BVH4::NodeRef> createLazyPrims(FastAllocator::Thread& alloc,
						     const size_t x0, const size_t x1,
						     const size_t y0, const size_t y1)
      {
	if (x1-x0 <= 8 && y1-y0 <= 8) {
	  EagerLeaf* leaf = new (alloc.malloc(sizeof(EagerLeaf))) EagerLeaf(*this);
	  const BBox3fa bounds = leaf->init(x0,x1,y0,y1);
	  return std::pair<BBox3fa,BVH4::NodeRef>(bounds,BVH4::encodeTypedLeaf(leaf,0));
	}
	
	struct Range2D { 
	  __forceinline Range2D () {}
	  __forceinline Range2D (size_t x0, size_t x1, size_t y0, size_t y1) : x0(x0), x1(x1), y0(y0), y1(y1) {}
	  __forceinline size_t width () const { return x1-x0; }
	  __forceinline size_t height() const { return y1-y0; }
	  size_t x0,x1,y0,y1; 
	};
	
	size_t N = 0;
	Range2D ranges[4];
	ranges[N++] = Range2D(x0,x1,y0,y1);
	for (size_t j=0; j<2; j++)
	{
	  for (size_t i=0; i<N && N<4; i++) 
	  {
	    /* ignore leaves */
	    Range2D& r = ranges[i];
	    if (r.width() <= 8 && r.height() <= 8) 
	      continue;
	    
	    /* split range along largest dimension */
	    if (r.width() > r.height()) {
	      const size_t x0 = r.x0, xc = (r.x0+r.x1)/2, x1 = r.x1;
	      const size_t y0 = r.y0,                     y1 = r.y1;
	      new (&ranges[i]  ) Range2D(x0,xc,y0,y1);
	      new (&ranges[N++]) Range2D(xc,x1,y0,y1);
	    } else {
	      const size_t x0 = r.x0, x1 = r.x1;
	      const size_t y0 = r.y0, yc = (r.y0+r.y1)/2, y1 = r.y1;
	      new (&ranges[i]  ) Range2D(x0,x1,y0,yc);
	      new (&ranges[N++]) Range2D(x0,x1,yc,y1);
	    }
	  }
	}

	/* recurse */
	BBox3fa bounds = empty;
	BVH4::Node* node = (BVH4::Node*) alloc.malloc(sizeof(BVH4::Node),16); node->clear();
	for (size_t i=0; i<N; i++) {
	  std::pair<BBox3fa,BVH4::NodeRef> child = createLazyPrims(alloc,ranges[i].x0,ranges[i].x1,ranges[i].y0,ranges[i].y1);
	  node->set(i,child.first,child.second);
	  bounds.extend(child.first);
	}
	return std::pair<BBox3fa,BVH4::NodeRef>(bounds,BVH4::encodeNode2(node));
      }
    
    struct LazyLeaf
    {
      /*! Construction from vertices and IDs. */
      __forceinline LazyLeaf (BVH4* bvh,
			      BVH4::NodeRef* parent,
			      const unsigned int geomID, 
			      const unsigned int primID, 
			      const unsigned int quadID,
			      const unsigned int x0,
			      const unsigned int x1,
			      const unsigned int y0,
			      const unsigned int y1)
	
	: initializing(0), initialized(0), bvh(bvh), parent(parent), geomID(geomID), primID(primID), quadID(quadID), x0(x0), x1(x1), y0(y0), y1(y1) {}
    

      __forceinline BBox3fa bounds()
      {
	BBox3fa box = empty;

	/* build sub-BVH */
	Scene* scene = bvh->scene;
	SubdivMesh* mesh = scene->getSubdivMesh(geomID);
	BVH4::NodeRef node = BVH4::emptyNode;

	feature_adaptive_subdivision_eval(mesh->getHalfEdge(primID),mesh->getVertexBuffer(), // FIXME: only recurse into one sub-quad
					  [&](const CatmullClarkPatch& patch, const Vec2f uv[4], const int subdiv[4], const int id)
	{
	  if (id != quadID) return;

	  const float l0 = patch.ring[0].edge_level;
	  const float l1 = patch.ring[1].edge_level;
	  const float l2 = patch.ring[2].edge_level;
	  const float l3 = patch.ring[3].edge_level;
	  const DiscreteTessellationPattern pattern0(l0,subdiv[0]);
	  const DiscreteTessellationPattern pattern1(l1,subdiv[1]);
	  const DiscreteTessellationPattern pattern2(l2,subdiv[2]);
	  const DiscreteTessellationPattern pattern3(l3,subdiv[3]);
	  const DiscreteTessellationPattern pattern_x = pattern0.size() > pattern2.size() ? pattern0 : pattern2;
	  const DiscreteTessellationPattern pattern_y = pattern1.size() > pattern3.size() ? pattern1 : pattern3;
	  const int nx = pattern_x.size();
	  const int ny = pattern_y.size();
	  const size_t width  = x1-x0+1;
	  const size_t height = y1-y0+1;
	  assert(width <= 17 && height <= 17);
	  Grid grid(width,height,geomID,primID);
	  box = grid.buildLazyBounds(scene,patch,x0,x1,y0,y1,uv[0],uv[1],uv[2],uv[3],pattern0,pattern1,pattern2,pattern3,pattern_x,pattern_y);
	});

	return box;
      }

      __forceinline size_t initialize()
      {
	/* let only one thread lazily build this object */
	if (atomic_add(&initializing,1) != 0) {
	  while (!initialized) __pause_cpu();
	  return (size_t)*parent;
	}

	/* build sub-BVH */
	Scene* scene = bvh->scene;
	SubdivMesh* mesh = scene->getSubdivMesh(geomID);
	BVH4::NodeRef node = BVH4::emptyNode;

	FastAllocator::Thread& alloc = *bvh->alloc2.instance();
	
	feature_adaptive_subdivision_eval(mesh->getHalfEdge(primID),mesh->getVertexBuffer(), // FIXME: only recurse into one sub-quad
					  [&](const CatmullClarkPatch& patch, const Vec2f uv[4], const int subdiv[4], const int id)
	{
	  if (id != quadID || node != BVH4::emptyNode) return;
	  
	  const float l0 = patch.ring[0].edge_level;
	  const float l1 = patch.ring[1].edge_level;
	  const float l2 = patch.ring[2].edge_level;
	  const float l3 = patch.ring[3].edge_level;
	  const DiscreteTessellationPattern pattern0(l0,subdiv[0]);
	  const DiscreteTessellationPattern pattern1(l1,subdiv[1]);
	  const DiscreteTessellationPattern pattern2(l2,subdiv[2]);
	  const DiscreteTessellationPattern pattern3(l3,subdiv[3]);
	  const DiscreteTessellationPattern pattern_x = pattern0.size() > pattern2.size() ? pattern0 : pattern2;
	  const DiscreteTessellationPattern pattern_y = pattern1.size() > pattern3.size() ? pattern1 : pattern3;
	  const int nx = pattern_x.size();
	  const int ny = pattern_y.size();
	  Grid* leaf = Grid::create(alloc,x1-x0+1,y1-y0+1,geomID,primID);
	  leaf->build(scene,patch,x0,x1,y0,y1,uv[0],uv[1],uv[2],uv[3],pattern0,pattern1,pattern2,pattern3,pattern_x,pattern_y);
	  node = leaf->createLazyPrims(alloc,0,x1-x0,0,y1-y0).second;
	});
	
	/* link to sub-BVH */
	*parent = node;
	__memory_barrier();
	initialized = 1;
	return (size_t)node;
      }
      
    public:
      volatile atomic_t initializing; // FIXME: reduce size
      volatile char initialized; 
            
    public:
      BVH4* bvh;
      BVH4::NodeRef* parent;
      unsigned int geomID;
      unsigned int primID;
      unsigned short quadID;
      unsigned short x0;
      unsigned short x1;
      unsigned short y0;
      unsigned short y1;
    };
    
    static size_t createLazy(BVH4* bvh, BVH4::NodeRef* parent, unsigned geomID, unsigned primID, unsigned quadID, unsigned width, unsigned height, 
			     FastAllocator::Thread& alloc, PrimRef* prims)
    {
      size_t N = 0;
      for (size_t y=0; y<height; y+=16)
      {
	for (size_t x=0; x<width; x+=16) 
	{
	  const size_t lx0 = x, lx1 = min(lx0+16,width);
	  const size_t ly0 = y, ly1 = min(ly0+16,height);
	  LazyLeaf* leaf = new (alloc.malloc(sizeof(LazyLeaf))) LazyLeaf(bvh,parent,geomID,primID,quadID,lx0,lx1,ly0,ly1);
	  const BBox3fa bounds = leaf->bounds();
	  prims[N++] = PrimRef(bounds,BVH4::encodeTypedLeaf(leaf,1));
	}
      }
      return N;
    }
    
  public:
    unsigned width;
    unsigned height;
    unsigned primID;
    unsigned geomID;
    Vec3fa P[17*17];
  };
}
