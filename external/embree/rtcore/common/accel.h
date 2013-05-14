// ======================================================================== //
// Copyright 2009-2012 Intel Corporation                                    //
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

#ifndef __EMBREE_ACCEL_H__
#define __EMBREE_ACCEL_H__

#include "default.h"

namespace embree
{
  /*! Malloc to allocate arrays for rtcCreateAccel function. */
  void* rtcMalloc(size_t bytes);

  /*! Frees all unused internal memory. */
  void rtcFreeMemory();

  /*! Triangle interface structure to the builder. The builders get an
   *  indexed face set as input, consisting of an array of vertices
   *  and triangles. If the topmost bit of id0 is set, the vertex IDs
   *  of this triangle are assumed to refer to the vertex position and
   *  a motion vector for motion blur. */
  struct BuildTriangle
  {
    /*! Constructs a builder triangle. */
    __forceinline BuildTriangle(int v0,      //!< reference to 1st vertex of the triangle
                                int v1,      //!< reference to 2nd vertex of the triangle
                                int v2,      //!< reference to 3rd vertex of the triangle
                                int id0 = 0, //!< 1st optional user ID
                                int id1 = 0) //!< 2nd optional user ID
      : v0(v0), v1(v1), v2(v2), id0(id0), id1(id1) { }
    
  public:
    int v0, v1, v2;       //!< references to 1st, 2nd, and 3rd triangle vertex
    int id0, id1;         //!< two user IDs
  };

  /*! Vertex type */
  struct BuildVertex {
  public:
    BuildVertex (float x, float y, float z) : x(x), y(y), z(z), align(0) {}
  public:
    float x,y,z; 
  private:
    int align;
  };

  /*! Abstract acceleration structure interface. */
  class Accel : public RefCount {
  public:

    /*! Construction with intersector type. */
    Accel (const std::string& intTy) : intTy(intTy) {}

    /*! A virtual destructor is required. */
    virtual ~Accel () {}

    /*! type safe interface query */
    template<typename Interface> Ref<Interface> queryInterface() {
      return dynamic_cast<Interface*>(query(Interface::name).ptr);
    }

  private:

    /*! Query interface to the acceleration structure. */
    virtual Ref<RefCount> query(const char* name) = 0;

  protected:

    /*! intersector type */
    const std::string intTy;
  };

  /*! Creates acceleration structure of specified type. Input is a
   *  mesh represented as an indexed face set. The triangle and vertex
   *  arrays have to be allocated with the rtcMalloc function and
   *  are deleted by the rtcCreateAccel function when no longer needed. */
  Ref<Accel> rtcCreateAccel(const char* accelTy,             //!< type of acceleration structure to use
                            const char* triTy,               //!< type of triangle representation to use
                            const BuildTriangle* triangles,  //!< array of triangles
                            size_t numTriangles,             //!< number of triangles in array
                            const BuildVertex* vertices,     //!< array of vertices, has to be aligned to 16 bytes
                            size_t numVertices,              //!< number of vertices in array
                            const BBox3f& bounds = empty,    //!< optional approximate bounding box of the geometry
                            bool freeArrays = true);         //!< if true, triangle and vertex arrays are freed when no longer needed
}

#endif
