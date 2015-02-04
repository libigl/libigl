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

#include "scene_triangle_mesh.h"
#include "scene.h"

namespace embree
{
  TriangleMesh::TriangleMesh (Scene* parent, RTCGeometryFlags flags, size_t numTriangles, size_t numVertices, size_t numTimeSteps)
    : Geometry(parent,TRIANGLE_MESH,numTriangles,flags), 
      mask(-1), numTimeSteps(numTimeSteps),
      numTriangles(numTriangles), numVertices(numVertices)
  {
    triangles.init(numTriangles,sizeof(Triangle));
    for (size_t i=0; i<numTimeSteps; i++) {
      vertices[i].init(numVertices,sizeof(Vec3fa));
    }
    enabling();
  }
  
  void TriangleMesh::enabling() 
  { 
    if (numTimeSteps == 1) { atomic_add(&parent->numTriangles ,numTriangles); }
    else                   { atomic_add(&parent->numTriangles2,numTriangles); }
  }
  
  void TriangleMesh::disabling() 
  { 
    if (numTimeSteps == 1) { atomic_add(&parent->numTriangles ,-(ssize_t)numTriangles); }
    else                   { atomic_add(&parent->numTriangles2,-(ssize_t)numTriangles); }
  }

  void TriangleMesh::setMask (unsigned mask) 
  {
    if (parent->isStatic() && parent->isBuild()) {
      process_error(RTC_INVALID_OPERATION,"static geometries cannot get modified");
      return;
    }
    this->mask = mask; 
  }

  void TriangleMesh::setBuffer(RTCBufferType type, void* ptr, size_t offset, size_t stride) 
  { 
    if (parent->isStatic() && parent->isBuild()) {
      process_error(RTC_INVALID_OPERATION,"static geometries cannot get modified");
      return;
    }

    /* verify that all accesses are 4 bytes aligned */
    if (((size_t(ptr) + offset) & 0x3) || (stride & 0x3)) {
      process_error(RTC_INVALID_OPERATION,"data must be 4 bytes aligned");
      return;
    }

    /* verify that all vertex accesses are 16 bytes aligned */
#if defined(__MIC__)
    // if (type == RTC_VERTEX_BUFFER0 || type == RTC_VERTEX_BUFFER1) {
    //   if (((size_t(ptr) + offset) & 0xF) || (stride & 0xF)) {
    //     process_error(RTC_INVALID_OPERATION,"data must be 16 bytes aligned");
    //     return;
    //   }
    // }
#endif

    switch (type) {
    case RTC_INDEX_BUFFER  : 
      triangles.set(ptr,offset,stride); 
      break;
    case RTC_VERTEX_BUFFER0: 
      vertices[0].set(ptr,offset,stride); 
      if (numVertices) {
        /* test if array is properly padded */
        volatile int w = *((int*)vertices[0].getPtr(numVertices-1)+3); // FIXME: is failing hard avoidable?
      }
      break;
    case RTC_VERTEX_BUFFER1: 
      vertices[1].set(ptr,offset,stride); 
      if (numVertices) {
        /* test if array is properly padded */
        volatile int w = *((int*)vertices[1].getPtr(numVertices-1)+3); // FIXME: is failing hard avoidable?
      }
      break;
    default: 
      process_error(RTC_INVALID_ARGUMENT,"unknown buffer type");
      break;
    }
  }

  void* TriangleMesh::map(RTCBufferType type) 
  {
    if (parent->isStatic() && parent->isBuild()) {
      process_error(RTC_INVALID_OPERATION,"static geometries cannot get modified");
      return NULL;
    }

    switch (type) {
    case RTC_INDEX_BUFFER  : return triangles  .map(parent->numMappedBuffers);
    case RTC_VERTEX_BUFFER0: return vertices[0].map(parent->numMappedBuffers);
    case RTC_VERTEX_BUFFER1: return vertices[1].map(parent->numMappedBuffers);
    default                : process_error(RTC_INVALID_ARGUMENT,"unknown buffer type"); return NULL;
    }
  }

  void TriangleMesh::unmap(RTCBufferType type) 
  {
    if (parent->isStatic() && parent->isBuild()) {
      process_error(RTC_INVALID_OPERATION,"static geometries cannot get modified");
      return;
    }

    switch (type) {
    case RTC_INDEX_BUFFER  : triangles  .unmap(parent->numMappedBuffers); break;
    case RTC_VERTEX_BUFFER0: vertices[0].unmap(parent->numMappedBuffers); break;
    case RTC_VERTEX_BUFFER1: vertices[1].unmap(parent->numMappedBuffers); break;
    default                : process_error(RTC_INVALID_ARGUMENT,"unknown buffer type"); break;
    }
  }

  void TriangleMesh::setUserData (void* ptr, bool ispc) {
    userPtr = ptr;
  }

  void TriangleMesh::immutable () 
  {
    bool freeTriangles = !parent->needTriangles;
    bool freeVertices  = !parent->needVertices;
    if (freeTriangles) triangles.free();
    if (freeVertices ) vertices[0].free();
    if (freeVertices ) vertices[1].free();
  }

  bool TriangleMesh::verify () 
  {
    for (size_t i=0; i<numTriangles; i++) {     
      if (triangles[i].v[0] >= numVertices) return false; 
      if (triangles[i].v[1] >= numVertices) return false; 
      if (triangles[i].v[2] >= numVertices) return false; 
    }
    for (size_t j=0; j<numTimeSteps; j++) {
      BufferT<Vec3fa>& verts = vertices[j];
      for (size_t i=0; i<numVertices; i++) {
	if (!inFloatRange(verts[i])) 
	  return false;
      }
    }
    return true;
  }

  void TriangleMesh::write(std::ofstream& file)
  {
    int type = TRIANGLE_MESH;
    file.write((char*)&type,sizeof(int));
    file.write((char*)&numTimeSteps,sizeof(int));
    file.write((char*)&numVertices,sizeof(int));
    file.write((char*)&numTriangles,sizeof(int));

    for (size_t j=0; j<numTimeSteps; j++) {
      while ((file.tellp() % 16) != 0) { char c = 0; file.write(&c,1); }
      for (size_t i=0; i<numVertices; i++) file.write((char*)vertexPtr(i,j),sizeof(Vec3fa));  
    }

    while ((file.tellp() % 16) != 0) { char c = 0; file.write(&c,1); }
    for (size_t i=0; i<numTriangles; i++) file.write((char*)&triangle(i),sizeof(Triangle));  

    while ((file.tellp() % 16) != 0) { char c = 0; file.write(&c,1); }
    for (size_t i=0; i<numTriangles; i++) file.write((char*)&triangle(i),sizeof(Triangle));  
  }
}
