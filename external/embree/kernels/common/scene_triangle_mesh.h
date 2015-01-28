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

#ifndef __EMBREE_TRIANGLE_MESH_SCENE_H__
#define __EMBREE_TRIANGLE_MESH_SCENE_H__

#include "common/default.h"
#include "common/geometry.h"
#include "common/buildsource.h"
#include "common/buffer.h"

namespace embree
{
  namespace TriangleMeshScene
  {

    /*! Triangle Mesh */
    struct TriangleMesh : public Geometry, public BuildSource
    {
      struct Triangle {
        unsigned int v[3];
      };

    public:
      TriangleMesh (Scene* parent, RTCGeometryFlags flags, size_t numTriangles, size_t numVertices, size_t numTimeSteps); 
      
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
      void setUserData (void* ptr, bool ispc);

      void enabling();
      void disabling();

    public:

      bool isEmpty () const { 
        return numTriangles == 0;
      }
      
      size_t groups () const { 
        return 1;
      }
      
      size_t prims (size_t group, size_t* pnumVertices) const {
        if (pnumVertices) *pnumVertices = numVertices*numTimeSteps;
        return numTriangles;
      }

      const BBox3f bounds(size_t group, size_t prim) const {
        return bounds(prim);
      }

      void bounds(size_t group, size_t begin, size_t end, BBox3f* bounds_o) const 
      {
        BBox3f b = empty;
        for (size_t i=begin; i<end; i++) b.extend(bounds(i));
        *bounds_o = b;
      }

      void split (const PrimRef& prim, int dim, float pos, PrimRef& left_o, PrimRef& right_o) const;

    public:

      __forceinline const Triangle& triangle(size_t i) const {
        assert(i < numTriangles);
        return triangles[i];
      }

      __forceinline const Vec3fa& vertex(size_t i, size_t j = 0) const {
        assert(i < numVertices);
        assert(j < 2);
        return vertices[j][i];
      }

      __forceinline BBox3f bounds(size_t index) const 
      {
        const Triangle& tri = triangle(index);
        const Vec3fa& v0 = vertex(tri.v[0]);
        const Vec3fa& v1 = vertex(tri.v[1]);
        const Vec3fa& v2 = vertex(tri.v[2]);
	return BBox3f( min(min(v0,v1),v2), max(max(v0,v1),v2) );
      }

      __forceinline bool anyMappedBuffers() const {
        return triangles.isMapped() || vertices[0].isMapped() || vertices[1].isMapped();
      }

    public:
      unsigned mask;                    //!< for masking out geometry
      bool built;                       //!< geometry got built
      unsigned char numTimeSteps;       //!< number of time steps (1 or 2)

      BufferT<Triangle> triangles;      //!< array of triangles
      bool needTriangles;               //!< set if triangle array required by acceleration structure
      size_t numTriangles;              //!< number of triangles

      BufferT<Vec3fa> vertices[2];      //!< vertex array
      bool needVertices;                //!< set if vertex array required by acceleration structure
      size_t numVertices;               //!< number of vertices
    };
  }
}

#endif
