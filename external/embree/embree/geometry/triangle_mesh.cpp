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

#include "triangle_mesh.h"

namespace embree
{
  AccelRegistry TriangleMesh::accels; 

  BuilderRegistry TriangleMesh::builders; 

  IntersectorRegistry<Intersector1> TriangleMesh::intersectors1;
    
#if defined(__SSE__)
  IntersectorRegistry<Intersector4> TriangleMesh::intersectors4;
#endif
    
#if defined(__AVX__)
  IntersectorRegistry<Intersector8> TriangleMesh::intersectors8;
#endif
    
#if defined(__MIC__)
  IntersectorRegistry<Intersector16> TriangleMesh::intersectors16;
#endif
    
  TriangleMesh::TriangleMesh (size_t numTriangles, size_t numVertices, const BBox3f& bounds, const char* accelTy) 
    : RTCGeometry(bounds), numTriangles(numTriangles), numVertices(numVertices)
  {
    triangles = (RTCTriangle*) alignedMalloc(numTriangles*sizeof(RTCTriangle));
    vertices  = (RTCVertex  *) alignedMalloc(numVertices*sizeof(RTCVertex   ));
    accel = accels.create(accelTy,this);
  }
  
  TriangleMesh::~TriangleMesh ()
  {
    if (triangles) alignedFree(triangles);
    if (vertices ) alignedFree(vertices);
  }

  void TriangleMesh::clearRegistry ()
  {
    accels.clear();
    builders.clear();
    intersectors1.clear();
#if defined(__SSE__)
    intersectors4.clear();
#endif
#if defined(__AVX__)
    intersectors8.clear();
#endif
#if defined(__MIC__)
    intersectors16.clear();
#endif
  }

  size_t TriangleMesh::size() const {
    return numTriangles;
  }

  BBox3f TriangleMesh::bounds(size_t i) const 
  {
    BBox3f b = empty;
    RTCTriangle& tri = triangles[i];
    b.grow(vertex(tri.v0));
    b.grow(vertex(tri.v1));
    b.grow(vertex(tri.v2));
    return b;
  }

  void TriangleMesh::build(TaskScheduler::Event* event, std::string builderName) {
    builders.build(builderName,event,accel);
  }

  void TriangleMesh::freeze()
  {
    if (triangles) alignedFree(triangles); triangles = NULL;
    if (accel && accel->trity.needVertices) return;
    if (vertices) alignedFree(vertices); vertices  = NULL;
  }
  
  Intersector1* TriangleMesh::intersector1(std::string travName) const{
    assert(accel);
    return intersectors1.get(accel->name(),travName,accel);
  }

#if defined(__SSE__)
  Intersector4* TriangleMesh::intersector4(std::string travName) const {
    assert(accel);
    return intersectors4.get(accel->name(),travName,accel);
  }
#endif

#if defined(__AVX__)
  Intersector8* TriangleMesh::intersector8(std::string travName) const {
    assert(accel);
    return intersectors8.get(accel->name(),travName,accel);
  }
#endif

#if defined(__MIC__)
  Intersector16* TriangleMesh::intersector16(std::string travName) const {
    assert(accel);
    return intersectors16.get(accel->name(),travName,accel);
  }
#endif
}
