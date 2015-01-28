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

#include "scene_triangle_mesh.h"
#include "scene.h"

namespace embree
{
  TriangleMeshScene::TriangleMesh::TriangleMesh (Scene* parent, RTCGeometryFlags flags, size_t numTriangles, size_t numVertices, size_t numTimeSteps)
    : Geometry(parent,TRIANGLE_MESH,numTriangles,flags), 
      mask(-1), built(false), numTimeSteps(numTimeSteps),
      numTriangles(numTriangles), needTriangles(false),
      numVertices(numVertices), needVertices(false)
  {
    triangles.init(numTriangles,sizeof(Triangle));
    for (size_t i=0; i<numTimeSteps; i++) {
      vertices[i].init(numVertices,sizeof(Vec3fa));
    }
    enabling();
  }
  
  void TriangleMeshScene::TriangleMesh::enabling() 
  { 
    if (numTimeSteps == 1) atomic_add(&parent->numTriangleMeshes ,1); 
    else                   atomic_add(&parent->numTriangleMeshes2,1); 
  }
  
  void TriangleMeshScene::TriangleMesh::disabling() 
  { 
    if (numTimeSteps == 1) atomic_add(&parent->numTriangleMeshes ,-1); 
    else                   atomic_add(&parent->numTriangleMeshes2,-1); 
  }

  void TriangleMeshScene::TriangleMesh::split (const PrimRef& prim, int dim, float pos, PrimRef& left_o, PrimRef& right_o) const
  {
    const TriangleMeshScene::TriangleMesh::Triangle& tri = triangle(prim.primID());
    const Vec3fa& v0 = vertex(tri.v[0]);
    const Vec3fa& v1 = vertex(tri.v[1]);
    const Vec3fa& v2 = vertex(tri.v[2]);
    splitTriangle(prim,dim,pos,v0,v1,v2,left_o,right_o);
  }
  
  void TriangleMeshScene::TriangleMesh::setMask (unsigned mask) 
  {
    if (parent->isStatic() && parent->isBuild()) {
      recordError(RTC_INVALID_OPERATION);
      return;
    }
    this->mask = mask; 
  }

  void TriangleMeshScene::TriangleMesh::enable () 
  {
    if (parent->isStatic() || anyMappedBuffers()) {
      recordError(RTC_INVALID_OPERATION);
      return;
    }
    Geometry::enable();
  }

  void TriangleMeshScene::TriangleMesh::update () 
  {
    if (parent->isStatic() || anyMappedBuffers()) {
      recordError(RTC_INVALID_OPERATION);
      return;
    }
    Geometry::update();
  }

  void TriangleMeshScene::TriangleMesh::disable () 
  {
    if (parent->isStatic() || anyMappedBuffers()) {
      recordError(RTC_INVALID_OPERATION);
      return;
    }
    Geometry::disable();
  }

  void TriangleMeshScene::TriangleMesh::erase () 
  {
    if (parent->isStatic() || anyMappedBuffers()) {
      recordError(RTC_INVALID_OPERATION);
      return;
    }
    Geometry::erase();
  }

  void TriangleMeshScene::TriangleMesh::setBuffer(RTCBufferType type, void* ptr, size_t offset, size_t stride) 
  { 
    if (parent->isStatic() && parent->isBuild()) {
      recordError(RTC_INVALID_OPERATION);
      return;
    }

    /* verify that all accesses are 4 bytes aligned */
    if (((size_t(ptr) + offset) & 0x3) || (stride & 0x3)) {
      recordError(RTC_INVALID_OPERATION);
      return;
    }

    /* verify that all vertex accesses are 16 bytes aligned */
#if defined(__MIC__)
    if (type == RTC_VERTEX_BUFFER0 || type == RTC_VERTEX_BUFFER1) {
      if (((size_t(ptr) + offset) & 0xF) || (stride & 0xF)) {
        recordError(RTC_INVALID_OPERATION);
        return;
      }
    }
#endif

    switch (type) {
    case RTC_INDEX_BUFFER  : 
      triangles.set(ptr,offset,stride); 
      break;
    case RTC_VERTEX_BUFFER0: 
      vertices[0].set(ptr,offset,stride); 
      if (numVertices) {
        /* test if array is properly padded */
        volatile int w = *((int*)&vertices[0][numVertices-1]+3); // FIXME: is failing hard avoidable?
      }
      break;
    case RTC_VERTEX_BUFFER1: 
      vertices[1].set(ptr,offset,stride); 
      if (numVertices) {
        /* test if array is properly padded */
        volatile int w = *((int*)&vertices[1][numVertices-1]+3); // FIXME: is failing hard avoidable?
      }
      break;
    default: 
      recordError(RTC_INVALID_ARGUMENT); break;
    }
  }

  void* TriangleMeshScene::TriangleMesh::map(RTCBufferType type) 
  {
    if (parent->isStatic() && parent->isBuild()) {
      recordError(RTC_INVALID_OPERATION);
      return NULL;
    }

    switch (type) {
    case RTC_INDEX_BUFFER  : return triangles  .map(parent->numMappedBuffers);
    case RTC_VERTEX_BUFFER0: return vertices[0].map(parent->numMappedBuffers);
    case RTC_VERTEX_BUFFER1: return vertices[1].map(parent->numMappedBuffers);
    default: 
      recordError(RTC_INVALID_ARGUMENT); 
      return NULL;
    }
  }

  void TriangleMeshScene::TriangleMesh::unmap(RTCBufferType type) 
  {
    if (parent->isStatic() && parent->isBuild()) {
      recordError(RTC_INVALID_OPERATION);
      return;
    }

    switch (type) {
    case RTC_INDEX_BUFFER  : triangles  .unmap(parent->numMappedBuffers); break;
    case RTC_VERTEX_BUFFER0: vertices[0].unmap(parent->numMappedBuffers); break;
    case RTC_VERTEX_BUFFER1: vertices[1].unmap(parent->numMappedBuffers); break;
    default                : recordError(RTC_INVALID_ARGUMENT); break;
    }
  }

  void TriangleMeshScene::TriangleMesh::setUserData (void* ptr, bool ispc) {
    userPtr = ptr;
  }

  void TriangleMeshScene::TriangleMesh::immutable () 
  {
    built = true;
    bool freeTriangles = !(needTriangles || parent->needTriangles);
    bool freeVertices  = !(needVertices  || parent->needVertices);
    if (freeTriangles) triangles.free();
    if (freeVertices ) vertices[0].free();
    if (freeVertices ) vertices[1].free();
  }

  bool TriangleMeshScene::TriangleMesh::verify () 
  {
    float range = sqrtf(0.5f*FLT_MAX);
    for (size_t i=0; i<numTriangles; i++) {
      if (triangles[i].v[0] >= numVertices) return false;
      if (triangles[i].v[1] >= numVertices) return false;
      if (triangles[i].v[2] >= numVertices) return false;
    }
    for (size_t j=0; j<numTimeSteps; j++) {
      BufferT<Vec3fa>& verts = vertices[j];
      for (size_t i=0; i<numVertices; i++) {
        if (verts[i].x < -range || verts[i].x > range) return false;
        if (verts[i].y < -range || verts[i].y > range) return false;
        if (verts[i].z < -range || verts[i].z > range) return false;
      }
    }
    return true;
  }
}
