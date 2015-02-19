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

#include "scene_bezier_curves.h"
#include "scene.h"

namespace embree
{
  BezierCurves::BezierCurves (Scene* parent, RTCGeometryFlags flags, size_t numCurves, size_t numVertices, size_t numTimeSteps) 
    : Geometry(parent,BEZIER_CURVES,numCurves,flags), 
      mask(-1), numTimeSteps(numTimeSteps), numCurves(numCurves), numVertices(numVertices)
  {
    curves.init(numCurves,sizeof(int));
    for (size_t i=0; i<numTimeSteps; i++) {
      vertices[i].init(numVertices,sizeof(Vertex));
    }
    enabling();
  }

  void BezierCurves::enabling() 
  { 
    if (numTimeSteps == 1) atomic_add(&parent->numBezierCurves ,numCurves); 
    else                   atomic_add(&parent->numBezierCurves2,numCurves); 
  }
  
  void BezierCurves::disabling() 
  { 
    if (numTimeSteps == 1) atomic_add(&parent->numBezierCurves ,-(ssize_t)numCurves); 
	else                   atomic_add(&parent->numBezierCurves2,-(ssize_t)numCurves);
  }
  
  void BezierCurves::setMask (unsigned mask) 
  {
    if (parent->isStatic() && parent->isBuild()) {
      process_error(RTC_INVALID_OPERATION,"static geometries cannot get modified");
      return;
    }
    this->mask = mask; 
  }

  void BezierCurves::setBuffer(RTCBufferType type, void* ptr, size_t offset, size_t stride) 
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
    if (type == RTC_VERTEX_BUFFER0 || type == RTC_VERTEX_BUFFER1) {
      if (((size_t(ptr) + offset) & 0xF) || (stride & 0xF)) {
        process_error(RTC_INVALID_OPERATION,"data must be 16 bytes aligned");
        return;
      }
    }
#endif

    switch (type) {
    case RTC_INDEX_BUFFER  : 
      curves.set(ptr,offset,stride); 
      break;
    case RTC_VERTEX_BUFFER0: 
      vertices[0].set(ptr,offset,stride); 
      break;
    case RTC_VERTEX_BUFFER1: 
      vertices[1].set(ptr,offset,stride); 
      break;
    default: 
      process_error(RTC_INVALID_ARGUMENT,"unknown buffer type");
      break;
    }
  }

  void* BezierCurves::map(RTCBufferType type) 
  {
    if (parent->isStatic() && parent->isBuild()) {
      process_error(RTC_INVALID_OPERATION,"static geometries cannot get modified");
      return NULL;
    }

    switch (type) {
    case RTC_INDEX_BUFFER  : return curves.map(parent->numMappedBuffers);
    case RTC_VERTEX_BUFFER0: return vertices[0].map(parent->numMappedBuffers);
    case RTC_VERTEX_BUFFER1: return vertices[1].map(parent->numMappedBuffers);
    default                : process_error(RTC_INVALID_ARGUMENT,"unknown buffer type"); return NULL;
    }
  }

  void BezierCurves::unmap(RTCBufferType type) 
  {
    if (parent->isStatic() && parent->isBuild()) {
      process_error(RTC_INVALID_OPERATION,"static geometries cannot get modified");
      return;
    }

    switch (type) {
    case RTC_INDEX_BUFFER  : curves.unmap(parent->numMappedBuffers); break;
    case RTC_VERTEX_BUFFER0: vertices[0].unmap(parent->numMappedBuffers); break;
    case RTC_VERTEX_BUFFER1: vertices[1].unmap(parent->numMappedBuffers); break;
    default                : process_error(RTC_INVALID_ARGUMENT,"unknown buffer type"); break;
    }
  }

  void BezierCurves::setUserData (void* ptr, bool ispc) {
    userPtr = ptr;
  }

  void BezierCurves::immutable () 
  {
    bool freeCurves    = true; //!parent->needCurves;
    bool freeVertices  = !parent->needVertices;
    if (freeCurves   ) curves.free();
    if (freeVertices ) vertices[0].free();
    if (freeVertices ) vertices[1].free();
  }

  bool BezierCurves::verify () 
  {
    for (size_t i=0; i<numCurves; i++) {
      if (curves[i]+3 >= numVertices) return false;
    }
    for (size_t j=0; j<numTimeSteps; j++) {
      BufferT<Vertex>& verts = vertices[j];
      for (size_t i=0; i<numVertices; i++) {
        if (!inFloatRange(verts[i].x)) return false;
	if (!inFloatRange(verts[i].y)) return false;
	if (!inFloatRange(verts[i].z)) return false;
	if (!inFloatRange(verts[i].r)) return false;
      }
    }
    return true;
  }

  void BezierCurves::write(std::ofstream& file)
  {
    int type = BEZIER_CURVES;
    file.write((char*)&type,sizeof(int));
    file.write((char*)&numTimeSteps,sizeof(int));
    file.write((char*)&numVertices,sizeof(int));
    file.write((char*)&numCurves,sizeof(int));

    for (size_t j=0; j<numTimeSteps; j++) {
      while ((file.tellp() % 16) != 0) { char c = 0; file.write(&c,1); }
      for (size_t i=0; i<numVertices; i++) file.write((char*)&vertex(i,j),sizeof(Vec3fa));  
    }

    while ((file.tellp() % 16) != 0) { char c = 0; file.write(&c,1); }
    for (size_t i=0; i<numCurves; i++) file.write((char*)&curve(i),sizeof(int));  
  }
}
