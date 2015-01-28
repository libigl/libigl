// ======================================================================== //
// Copyright 2009-2013 Intel Corporation                                    //
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

#ifndef __EMBREE_QUADRATIC_BEZIER_CURVES_SCENE_H__
#define __EMBREE_QUADRATIC_BEZIER_CURVES_SCENE_H__

#include "common/default.h"
#include "common/geometry.h"
#include "common/buildsource.h"
#include "common/buffer.h"

namespace embree
{
  namespace QuadraticBezierCurvesScene
  {

    struct QuadraticBezierCurves : public Geometry, public BuildSource
    {
      struct Vertex {
        float x,y,z,r;
      };

    public:
      QuadraticBezierCurves (Scene* parent, RTCGeometryFlags flags, size_t numCurves, size_t numVertices, size_t numTimeSteps); 
      
    public:
      void setMask (unsigned mask);
      void enable ();
      void update ();
      void disable ();
      void erase ();
      void immutable ();
      bool verify ();
      void setBuffer(RTCBufferType type, void* ptr, size_t offset, size_t stride);
      void* map(RTCBufferType type);
      void unmap(RTCBufferType type);

      void enabling();
      void disabling();

    public:

      bool isEmpty () const { 
        return numCurves == 0;
      }
      
      size_t groups () const { 
        return 1;
      }
      
      size_t prims (size_t group, size_t* pnumVertices) const {
        if (pnumVertices) *pnumVertices = numVertices*numTimeSteps;
        return numCurves;
      }

      const BBox3f bounds(size_t group, size_t prim) const {
        return bounds(prim);
      }

    public:

      __forceinline const int& curve(size_t i) const {
        assert(i < numCurves);
        return curves[i];
      }

      __forceinline const Vec3fa& vertex(size_t i, size_t j = 0) const {
        assert(i < numVertices);
        assert(j < 2);
        return (Vec3fa&)vertices[j][i];
      }

      __forceinline float radius(size_t i, size_t j = 0) const {
        assert(i < numVertices);
        assert(j < 2);
        return vertices[j][i].r;
      }

      __forceinline BBox3f bounds(size_t i) const 
      {
        const int index = curve(i);
        const float r0 = radius(index+0);
        const float r1 = radius(index+1);
        const float r2 = radius(index+2);
        const float r3 = radius(index+3);
        const Vec3fa& v0 = vertex(index+0);
        const Vec3fa& v1 = vertex(index+1);
        const Vec3fa& v2 = vertex(index+2);
        const Vec3fa& v3 = vertex(index+3);
        const BBox3f b = merge(BBox3f(v0),BBox3f(v1),BBox3f(v2),BBox3f(v3));
        return enlarge(b,Vec3fa(max(r0,r1,r2,r3)));
      }

      __forceinline bool anyMappedBuffers() const {
        return curves.isMapped() || vertices[0].isMapped() || vertices[1].isMapped();
      }

    public:
      unsigned mask;                    //!< for masking out geometry
      bool built;                       //!< geometry got built
      unsigned char numTimeSteps;       //!< number of time steps (1 or 2)

      BufferT<int> curves;              //!< array of curve indices
      bool needCurves;                  //!< set if curve indices required by acceleration structure
      size_t numCurves;                 //!< number of triangles

      BufferT<Vertex> vertices[2];      //!< vertex array
      bool needVertices;                //!< set if vertex array required by acceleration structure
      size_t numVertices;               //!< number of vertices
    };

  }
}

#endif
