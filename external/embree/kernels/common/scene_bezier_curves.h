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
#include "common/primref.h"
#include "common/buffer.h"

namespace embree
{
    struct BezierCurves : public Geometry
    {
      struct Vertex {
        float x,y,z,r;
      };

      static const GeometryTy geom_type = BEZIER_CURVES;

    public:
      BezierCurves (Scene* parent, RTCGeometryFlags flags, size_t numCurves, size_t numVertices, size_t numTimeSteps); 
    
      void write(std::ofstream& file);

    public:
      void enabling();
      void disabling();
      void setMask (unsigned mask);
      void setBuffer(RTCBufferType type, void* ptr, size_t offset, size_t stride);
      void* map(RTCBufferType type);
      void unmap(RTCBufferType type);
      void setUserData (void* ptr, bool ispc);
      void immutable ();
      bool verify ();

    public:

      /*! returns number of bezier curves */
      __forceinline size_t size() const {
	return numCurves;
      }

      /*! returns the i'th curve */
      __forceinline const int& curve(size_t i) const {
        assert(i < numCurves);
        return curves[i];
      }

      /*! returns i'th vertex of j'th timestep */
      __forceinline const Vec3fa& vertex(size_t i, size_t j = 0) const {
        assert(i < numVertices);
        assert(j < numTimeSteps);
        return (Vec3fa&)vertices[j][i];
      }

      /*! returns i'th radius of j'th timestep */
      __forceinline float radius(size_t i, size_t j = 0) const {
        assert(i < numVertices);
        assert(j < numTimeSteps);
        return vertices[j][i].r;
      }

      __forceinline unsigned int maxSubdivisionSteps(const float eps, const Vec3fa& p0, const Vec3fa& p1, const Vec3fa& p2, const Vec3fa& p3) const {
	const Vec3fa d0 = abs(p0 - 2.0f * p1 + p2);
	const Vec3fa d1 = abs(p1 - 2.0f * p2 + p3);
	const float d0_max = max(d0.x,d0.y,d0.z);
	const float d1_max = max(d1.x,d1.y,d1.z);
	const float L0 = max(d0_max,d1_max);
	const float r0 = (sqrtf(2.0f) * 4 * (4-1) * L0) / (8.0f * eps);      
	return (unsigned int)logf(r0);
      }

      /*! check if the i'th primitive is valid */
      __forceinline bool valid(size_t i, BBox3fa* bbox = NULL) const 
      {
	const int index = curve(i);
	if (index+3 >= numVertices) return false;

	for (size_t j=0; j<numTimeSteps; j++) 
	{
	  const float r0 = radius(index+0,j);
	  const float r1 = radius(index+1,j);
	  const float r2 = radius(index+2,j);
	  const float r3 = radius(index+3,j);
	  if (!inFloatRange(r0) || !inFloatRange(r1) || !inFloatRange(r2) || !inFloatRange(r3))
	    return false;

	  const Vec3fa& v0 = vertex(index+0,j);
	  const Vec3fa& v1 = vertex(index+1,j);
	  const Vec3fa& v2 = vertex(index+2,j);
	  const Vec3fa& v3 = vertex(index+3,j);
	  if (!inFloatRange(v0) || !inFloatRange(v1) || !inFloatRange(v2) || !inFloatRange(v3))
	    return false;
	}

	if (bbox) *bbox = bounds(i);
	return true;
      }

      /*! calculates bounding box of i'th bezier curve */
      __forceinline BBox3fa bounds(size_t i, size_t j = 0) const 
      {
        const int index = curve(i);
        const float r0 = radius(index+0,j);
        const float r1 = radius(index+1,j);
        const float r2 = radius(index+2,j);
        const float r3 = radius(index+3,j);
        const Vec3fa& v0 = vertex(index+0,j);
        const Vec3fa& v1 = vertex(index+1,j);
        const Vec3fa& v2 = vertex(index+2,j);
        const Vec3fa& v3 = vertex(index+3,j);
        const BBox3fa b = merge(BBox3fa(v0),BBox3fa(v1),BBox3fa(v2),BBox3fa(v3));
        return enlarge(b,Vec3fa(max(r0,r1,r2,r3)));
      }

      /*! calculates bounding box of i'th bezier curve */
      __forceinline BBox3fa bounds(const AffineSpace3fa& space, size_t i, size_t j = 0) const 
      {
        const int index = curve(i);
        const float r0 = radius(index+0,j);
        const float r1 = radius(index+1,j);
        const float r2 = radius(index+2,j);
        const float r3 = radius(index+3,j);
        const Vec3fa& v0 = xfmPoint(space,vertex(index+0,j));
        const Vec3fa& v1 = xfmPoint(space,vertex(index+1,j));
        const Vec3fa& v2 = xfmPoint(space,vertex(index+2,j));
        const Vec3fa& v3 = xfmPoint(space,vertex(index+3,j));
        const BBox3fa b = merge(BBox3fa(v0),BBox3fa(v1),BBox3fa(v2),BBox3fa(v3));
        return enlarge(b,Vec3fa(max(r0,r1,r2,r3)));
      }

      /*! calculates residual bounding box of i'th bezier curve */
      __forceinline BBox3fa bounds(const AffineSpace3fa& space0, const AffineSpace3fa& space1, size_t i) const 
      {
#if 1
        const int index = curve(i);
        const float r0 = radius(index+0,0) + radius(index+0,1); // FIXME: can one use max here?
        const float r1 = radius(index+1,0) + radius(index+1,1);
        const float r2 = radius(index+2,0) + radius(index+2,1);
        const float r3 = radius(index+3,0) + radius(index+3,1);

        const Vec3fa a0 = xfmPoint(space0,vertex(index+0,1));
        const Vec3fa a1 = xfmPoint(space0,vertex(index+1,1));
        const Vec3fa a2 = xfmPoint(space0,vertex(index+2,1));
        const Vec3fa a3 = xfmPoint(space0,vertex(index+3,1));

        const Vec3fa b0 = xfmPoint(space1,vertex(index+0,0));
        const Vec3fa b1 = xfmPoint(space1,vertex(index+1,0));
        const Vec3fa b2 = xfmPoint(space1,vertex(index+2,0));
        const Vec3fa b3 = xfmPoint(space1,vertex(index+3,0));

        const Vec3fa v0 = a0 + b0;
        const Vec3fa v1 = a1 + b1;
        const Vec3fa v2 = a2 + b2;
        const Vec3fa v3 = a3 + b3;

        const BBox3fa b = merge(BBox3fa(v0),BBox3fa(v1),BBox3fa(v2),BBox3fa(v3));
        return enlarge(b,Vec3fa(max(r0,r1,r2,r3))); // FIXME: is radius used properly?
#else

        const int index = curve(i);
        const float r0 = 0.5*radius(index+0,0) + 0.5*radius(index+0,1);
        const float r1 = 0.5*radius(index+1,0) + 0.5*radius(index+1,1);
        const float r2 = 0.5*radius(index+2,0) + 0.5*radius(index+2,1);
        const float r3 = 0.5*radius(index+3,0) + 0.5*radius(index+3,1);

        const Vec3fa p0 = 0.5f*vertex(index+0,0) + 0.5f*vertex(index+0,1);
        const Vec3fa p1 = 0.5f*vertex(index+1,0) + 0.5f*vertex(index+1,1);
        const Vec3fa p2 = 0.5f*vertex(index+2,0) + 0.5f*vertex(index+2,1);
        const Vec3fa p3 = 0.5f*vertex(index+3,0) + 0.5f*vertex(index+3,1);

        const AffineSpace3fa space = 0.5f*space0 + 0.5f*space1;
        //const AffineSpace3fa space = frame(normalize(0.5f*space0.row2() + 0.5f*space1.row2())).transposed();
        const Vec3fa v0 = xfmPoint(space,p0);
        const Vec3fa v1 = xfmPoint(space,p1);
        const Vec3fa v2 = xfmPoint(space,p2);
        const Vec3fa v3 = xfmPoint(space,p3);

        const BBox3fa b = merge(BBox3fa(v0),BBox3fa(v1),BBox3fa(v2),BBox3fa(v3));
        return enlarge(b,Vec3fa(max(r0,r1,r2,r3)));

#endif
      }

      __forceinline const Vec3fa *fristVertexPtr(size_t i) const { // FIXME: remove, use buffer to access vertices instead!
        return &vertex(curve(i));
      }

#if defined(__MIC__)

      __forceinline mic2f bounds_mic2f(size_t i) const 
      {
        const int index = curve(i);
        const Vec3fa& cp0 = vertex(index+0);
        const Vec3fa& cp1 = vertex(index+1);
        const Vec3fa& cp2 = vertex(index+2);
        const Vec3fa& cp3 = vertex(index+3);
	
#if 0
	const mic_f v0 = broadcast4to16f((float*)&cp0);
	const mic_f v1 = broadcast4to16f((float*)&cp1);
	const mic_f v2 = broadcast4to16f((float*)&cp2);
	const mic_f v3 = broadcast4to16f((float*)&cp3);
#else
	const mic_m m_4f = 0xf;
	const mic_f v0 = permute<0,0,0,0>(uload16f(m_4f,(float*)&cp0));
	const mic_f v1 = permute<0,0,0,0>(uload16f(m_4f,(float*)&cp1));
	const mic_f v2 = permute<0,0,0,0>(uload16f(m_4f,(float*)&cp2));
	const mic_f v3 = permute<0,0,0,0>(uload16f(m_4f,(float*)&cp3));
#endif

	const mic_f b_min = min(min(v0,v1),min(v2,v3));
	const mic_f b_max = max(max(v0,v1),max(v2,v3));

	const mic_f b_min_r = b_min - swDDDD(b_max);
	const mic_f b_max_r = b_max + swDDDD(b_max);

        return mic2f(b_min_r,b_max_r);
      }
      
#endif

      __forceinline BBox3fa subBounds(size_t curveID, size_t segmentID) const 
      {
	assert(curveID < numCurves);
	assert(segmentID < 8);
        const int index = curve(curveID);
        const float r0 = radius(index+0);
        const float r1 = radius(index+1);
        const float r2 = radius(index+2);
        const float r3 = radius(index+3);
        const Vec3fa& v0 = vertex(index+0);
        const Vec3fa& v1 = vertex(index+1);
        const Vec3fa& v2 = vertex(index+2);
        const Vec3fa& v3 = vertex(index+3);

	BBox3fa b(empty);
	{
	  unsigned int step = segmentID;
	  float t1 = (float)step / 8.0f;
	  float t0 = 1.0f - t1;
	  const float coeff0 = t0 * t0 * t0;
	  const float coeff1 = 3.0f * t1* t0 * t0;
	  const float coeff2 = 3.0f * t1* t1 * t0;
	  const float coeff3 = t1 * t1 * t1;
	  const Vec3fa p = coeff0 * v0 + coeff1 * v1 + coeff2 * v2 + coeff3 * v3; 
	  b.extend(p);
	}

	{
	  unsigned int step = segmentID+1;
	  float t1 = (float)step / 8.0f;
	  float t0 = 1.0f - t1;
	  const float coeff0 = t0 * t0 * t0;
	  const float coeff1 = 3.0f * t1* t0 * t0;
	  const float coeff2 = 3.0f * t1* t1 * t0;
	  const float coeff3 = t1 * t1 * t1;
	  const Vec3fa p = coeff0 * v0 + coeff1 * v1 + coeff2 * v2 + coeff3 * v3; 
	  b.extend(p);
	}
        return enlarge(b,Vec3fa(max(r0,r1,r2,r3)));
      }

    public:
      unsigned int mask;                //!< for masking out geometry
      unsigned char numTimeSteps;       //!< number of time steps (1 or 2)

      BufferT<int> curves;              //!< array of curve indices
      size_t numCurves;                 //!< number of triangles

      BufferT<Vertex> vertices[2];      //!< vertex array
      size_t numVertices;               //!< number of vertices
    };
}
