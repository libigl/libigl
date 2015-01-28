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

#include "triangle1.h"
#include "common/scene.h"

namespace embree
{
  SceneTriangle1 SceneTriangle1::type;
  TriangleMeshTriangle1 TriangleMeshTriangle1::type;

  Triangle1Type::Triangle1Type () 
      : PrimitiveType("triangle1",sizeof(Triangle1),1,false,1) {} 
  
  size_t Triangle1Type::blocks(size_t x) const {
    return x;
  }
    
  size_t Triangle1Type::size(const char* This) const {
    return 1;
  }

  void SceneTriangle1::pack(char* dst, atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, void* geom) const 
  {
    Scene* scene = (Scene*) geom;
    const PrimRef& prim = *prims;
    const unsigned geomID = prim.geomID();
    const unsigned primID = prim.primID();
    const TriangleMeshScene::TriangleMesh* mesh = scene->getTriangleMesh(geomID);
    const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(primID);
    const Vec3fa& p0 = mesh->vertex(tri.v[0]);
    const Vec3fa& p1 = mesh->vertex(tri.v[1]);
    const Vec3fa& p2 = mesh->vertex(tri.v[2]);
    new (dst) Triangle1(p0,p1,p2,mesh->id,primID,mesh->mask);
    prims++;
  }
  
  void SceneTriangle1::pack(char* dst, const PrimRef* prims, size_t num, void* geom) const 
  {
    Scene* scene = (Scene*) geom;
    const PrimRef& prim = *prims;
    const unsigned geomID = prim.geomID();
    const unsigned primID = prim.primID();
    const TriangleMeshScene::TriangleMesh* mesh = scene->getTriangleMesh(geomID);
    const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(primID);
    const Vec3fa& p0 = mesh->vertex(tri.v[0]);
    const Vec3fa& p1 = mesh->vertex(tri.v[1]);
    const Vec3fa& p2 = mesh->vertex(tri.v[2]);
    new (dst) Triangle1(p0,p1,p2,mesh->id,primID,mesh->mask);
    prims++;
  }
    
  BBox3f SceneTriangle1::update(char* prim, size_t num, void* geom) const 
  {
    BBox3f bounds = empty;
    Scene* scene = (Scene*) geom;
    
    for (size_t j=0; j<num; j++) 
    {
      Triangle1& dst = ((Triangle1*) prim)[j];
      const unsigned geomID = dst.geomID();
      const unsigned primID = dst.primID();
      const TriangleMeshScene::TriangleMesh* mesh = scene->getTriangleMesh(geomID);
      const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(primID);
      const Vec3fa v0 = mesh->vertex(tri.v[0]);
      const Vec3fa v1 = mesh->vertex(tri.v[1]);
      const Vec3fa v2 = mesh->vertex(tri.v[2]);
      new (&dst) Triangle1(v0,v1,v2,geomID,primID,mesh->mask);
      bounds.extend(merge(BBox3f(v0),BBox3f(v1),BBox3f(v2)));
    }
    return bounds; 
  }

  void TriangleMeshTriangle1::pack(char* dst, atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, void* geom) const 
  {
    const PrimRef& prim = *prims;
    const unsigned geomID = prim.geomID();
    const unsigned primID = prim.primID();
    const TriangleMeshScene::TriangleMesh* mesh = (TriangleMeshScene::TriangleMesh*) geom;
    const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(primID);
    const Vec3fa& p0 = mesh->vertex(tri.v[0]);
    const Vec3fa& p1 = mesh->vertex(tri.v[1]);
    const Vec3fa& p2 = mesh->vertex(tri.v[2]);
    new (dst) Triangle1(p0,p1,p2,mesh->id,primID,mesh->mask);
    prims++;
  }
  
  void TriangleMeshTriangle1::pack(char* dst, const PrimRef* prims, size_t num, void* geom) const 
  {
    const PrimRef& prim = *prims;
    const unsigned geomID = prim.geomID();
    const unsigned primID = prim.primID();
    const TriangleMeshScene::TriangleMesh* mesh = (TriangleMeshScene::TriangleMesh*) geom;
    const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(primID);
    const Vec3fa& p0 = mesh->vertex(tri.v[0]);
    const Vec3fa& p1 = mesh->vertex(tri.v[1]);
    const Vec3fa& p2 = mesh->vertex(tri.v[2]);
    new (dst) Triangle1(p0,p1,p2,mesh->id,primID,mesh->mask);
    prims++;
  }
  
  BBox3f TriangleMeshTriangle1::update(char* prim, size_t num, void* geom) const 
  {
    BBox3f bounds = empty;
    const TriangleMeshScene::TriangleMesh* mesh = (const TriangleMeshScene::TriangleMesh*) geom;
    
    for (size_t j=0; j<num; j++) 
    {
      Triangle1& dst = ((Triangle1*) prim)[j];
      const unsigned geomID = dst.geomID();
      const unsigned primID = dst.primID();
      const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(primID);
      const Vec3fa v0 = mesh->vertex(tri.v[0]);
      const Vec3fa v1 = mesh->vertex(tri.v[1]);
      const Vec3fa v2 = mesh->vertex(tri.v[2]);
      new (&dst) Triangle1(v0,v1,v2,geomID,primID,mesh->mask);
      bounds.extend(merge(BBox3f(v0),BBox3f(v1),BBox3f(v2)));
    }
    return bounds; 
  }
}
