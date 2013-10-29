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

#ifndef __EMBREE_BUILD_TRIANGLE_MESH_H__
#define __EMBREE_BUILD_TRIANGLE_MESH_H__

#include "geometry.h"
#include "../common/registry_accel.h"
#include "../common/registry_builder.h"
#include "../common/registry_intersector.h"

namespace embree
{
  /*! Triangle interface structure to the builder. The builders get an
   *  indexed face set as input, consisting of an array of vertices
   *  and triangles. If the topmost bit of id0 is set, the vertex IDs
   *  of this triangle are assumed to refer to the vertex position and
   *  a motion vector for motion blur. */
  struct TriangleMesh : public RTCGeometry
  {
    typedef RTCVertex Vertex;
    typedef RTCTriangle Triangle;

    /*! acceleration structure registry */
    static AccelRegistry accels; 

    /*! builder registry */
    static BuilderRegistry builders; 

    /*! intersector registrys */
    static IntersectorRegistry<Intersector1> intersectors1;
    
#if defined(__SSE__)
    static IntersectorRegistry<Intersector4> intersectors4;
#endif
    
#if defined(__AVX__)
    static IntersectorRegistry<Intersector8> intersectors8;
#endif
    
#if defined(__MIC__)
    static IntersectorRegistry<Intersector16> intersectors16;
#endif
    
    /*! Constructor */
    TriangleMesh (size_t numTriangles, size_t numVertices, const BBox3f& bounds, const char* accelTy);

    /*! Destruction. */
    ~TriangleMesh ();

    /*! clear registries */
    static void clearRegistry ();

    /*! returns number of triangles */
    size_t size() const;
    
    /*! returns bounds of ith triangle */
    BBox3f bounds(size_t i) const;

    /*! builds acceleration structure */
    void build(TaskScheduler::Event* event, std::string builderName);

    /*! frees unused data */
    void freeze();
    
    /*! returns intersector */
    Intersector1* intersector1(std::string travName) const;

#if defined(__SSE__)
    Intersector4* intersector4(std::string travName) const;
#endif

#if defined(__AVX__)
    Intersector8* intersector8(std::string travName) const;
#endif

#if defined(__MIC__)
    Intersector16* intersector16(std::string travName) const;
#endif

    /*! returns the ith triangle. */
    __forceinline const RTCTriangle& triangle(size_t i) const {
      return triangles[i];
    }

    /*! returns the ith vertex. */
    __forceinline const Vec3fa& vertex(size_t i) const {
      return ((Vec3fa*)vertices)[i];
    }

    /*! returns pointer to vertex array */
    void* getVertices() { 
      return vertices; 
    }

    /*! returns number of vertices */
    size_t getNumVertices() const { 
      return numVertices; 
    }

  public:
    RTCTriangle* triangles;  //!< array of triangles
    size_t numTriangles;       //!< number of triangles in array
    RTCVertex* vertices;     //!< array of vertices, has to be aligned to 16 bytes
    size_t numVertices;        //!< number of vertices in array
  };
}

#endif
