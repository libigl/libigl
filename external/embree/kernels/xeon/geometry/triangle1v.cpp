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

#include "triangle1v.h"
#include "common/scene.h"

namespace embree
{
  SceneTriangle1v SceneTriangle1v::type;
  TriangleMeshTriangle1v TriangleMeshTriangle1v::type;

  SceneTriangle1vMB SceneTriangle1vMB::type;
  TriangleMeshTriangle1vMB TriangleMeshTriangle1vMB::type;

  Triangle1vType::Triangle1vType () 
      : PrimitiveType("triangle1v",sizeof(Triangle1v),1,false,1) {} 
  
  size_t Triangle1vType::blocks(size_t x) const {
    return x;
  }
    
  size_t Triangle1vType::size(const char* This) const {
    return 1;
  }

  void SceneTriangle1v::pack(char* dst, atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, void* geom) const 
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
    new (dst) Triangle1v(p0,p1,p2,mesh->id,primID,mesh->mask);
    prims++;
  }
  
  void SceneTriangle1v::pack(char* dst, const PrimRef* prims, size_t num, void* geom) const 
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
    new (dst) Triangle1v(p0,p1,p2,mesh->id,primID,mesh->mask);
    prims++;
  }
    
  BBox3f SceneTriangle1v::update(char* prim, size_t num, void* geom) const 
  {
    BBox3f bounds = empty;
    Scene* scene = (Scene*) geom;
    
    for (size_t j=0; j<num; j++) 
    {
      Triangle1v& dst = ((Triangle1v*) prim)[j];
      const unsigned geomID = dst.geomID();
      const unsigned primID = dst.primID();
      const TriangleMeshScene::TriangleMesh* mesh = scene->getTriangleMesh(geomID);
      const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(primID);
      const Vec3fa v0 = mesh->vertex(tri.v[0]);
      const Vec3fa v1 = mesh->vertex(tri.v[1]);
      const Vec3fa v2 = mesh->vertex(tri.v[2]);
      new (&dst) Triangle1v(v0,v1,v2,geomID,primID,mesh->mask);
      bounds.extend(merge(BBox3f(v0),BBox3f(v1),BBox3f(v2)));
    }
    return bounds; 
  }

  void TriangleMeshTriangle1v::pack(char* dst, atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, void* geom) const 
  {
    const PrimRef& prim = *prims;
    const unsigned geomID = prim.geomID();
    const unsigned primID = prim.primID();
    const TriangleMeshScene::TriangleMesh* mesh = (TriangleMeshScene::TriangleMesh*) geom;
    const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(primID);
    const Vec3fa& p0 = mesh->vertex(tri.v[0]);
    const Vec3fa& p1 = mesh->vertex(tri.v[1]);
    const Vec3fa& p2 = mesh->vertex(tri.v[2]);
    new (dst) Triangle1v(p0,p1,p2,mesh->id,primID,mesh->mask);
    prims++;
  }
  
  void TriangleMeshTriangle1v::pack(char* dst, const PrimRef* prims, size_t num, void* geom) const 
  {
    const PrimRef& prim = *prims;
    const unsigned geomID = prim.geomID();
    const unsigned primID = prim.primID();
    const TriangleMeshScene::TriangleMesh* mesh = (TriangleMeshScene::TriangleMesh*) geom;
    const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(primID);
    const Vec3fa& p0 = mesh->vertex(tri.v[0]);
    const Vec3fa& p1 = mesh->vertex(tri.v[1]);
    const Vec3fa& p2 = mesh->vertex(tri.v[2]);
    new (dst) Triangle1v(p0,p1,p2,mesh->id,primID,mesh->mask);
    prims++;
  }
  
  BBox3f TriangleMeshTriangle1v::update(char* prim, size_t num, void* geom) const 
  {
    BBox3f bounds = empty;
    const TriangleMeshScene::TriangleMesh* mesh = (const TriangleMeshScene::TriangleMesh*) geom;
    
    for (size_t j=0; j<num; j++) 
    {
      Triangle1v& dst = ((Triangle1v*) prim)[j];
      const unsigned geomID = dst.geomID();
      const unsigned primID = dst.primID();
      const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(primID);
      const Vec3fa v0 = mesh->vertex(tri.v[0]);
      const Vec3fa v1 = mesh->vertex(tri.v[1]);
      const Vec3fa v2 = mesh->vertex(tri.v[2]);
      new (&dst) Triangle1v(v0,v1,v2,geomID,primID,mesh->mask);
      bounds.extend(merge(BBox3f(v0),BBox3f(v1),BBox3f(v2)));
    }
    return bounds; 
  }

  Triangle1vMBType::Triangle1vMBType () 
  : PrimitiveType("triangle1v_mb",sizeof(Triangle1vMB),1,false,1) {} 
  
  size_t Triangle1vMBType::blocks(size_t x) const {
    return x;
  }

  size_t Triangle1vMBType::size(const char* This) const {
    return 1;
  }

  void SceneTriangle1vMB::pack(char* dst, atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, void* geom) const 
  {
    Scene* scene = (Scene*) geom;
    const PrimRef& prim = *prims;
    const unsigned geomID = prim.geomID();
    const unsigned primID = prim.primID();
    const TriangleMeshScene::TriangleMesh* mesh = scene->getTriangleMesh(geomID);
    const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(primID);
    const Vec3fa& a0 = mesh->vertex(tri.v[0],0);
    const Vec3fa& a1 = mesh->vertex(tri.v[0],1);
    const Vec3fa& b0 = mesh->vertex(tri.v[1],0);
    const Vec3fa& b1 = mesh->vertex(tri.v[1],1);
    const Vec3fa& c0 = mesh->vertex(tri.v[2],0);
    const Vec3fa& c1 = mesh->vertex(tri.v[2],1);
    new (dst) Triangle1vMB(a0,a1,b0,b1,c0,c1,mesh->id,primID,mesh->mask);
    prims++;
  }
    
  void SceneTriangle1vMB::pack(char* dst, const PrimRef* prims, size_t num, void* geom) const 
  {
    Scene* scene = (Scene*) geom;
    const PrimRef& prim = *prims;
    const unsigned geomID = prim.geomID();
    const unsigned primID = prim.primID();
    const TriangleMeshScene::TriangleMesh* mesh = scene->getTriangleMesh(geomID);
    const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(primID);
    const Vec3fa& a0 = mesh->vertex(tri.v[0],0);
    const Vec3fa& a1 = mesh->vertex(tri.v[0],1);
    const Vec3fa& b0 = mesh->vertex(tri.v[1],0);
    const Vec3fa& b1 = mesh->vertex(tri.v[1],1);
    const Vec3fa& c0 = mesh->vertex(tri.v[2],0);
    const Vec3fa& c1 = mesh->vertex(tri.v[2],1);
    new (dst) Triangle1vMB(a0,a1,b0,b1,c0,c1,mesh->id,primID,mesh->mask);
    prims++;
  }
  
  std::pair<BBox3f,BBox3f> SceneTriangle1vMB::update2(char* prim, size_t num, void* geom) const 
  {
    BBox3f bounds0 = empty, bounds1 = empty;
    
    for (size_t j=0; j<num; j++) 
    {
      const Triangle1vMB& tri = ((Triangle1vMB*) prim)[j];
      bounds0.extend(merge(BBox3f(tri.v0),BBox3f(tri.v1),BBox3f(tri.v2)));
      bounds1.extend(merge(BBox3f(tri.v0+tri.d0),BBox3f(tri.v1+tri.d1),BBox3f(tri.v2+tri.d2)));
    }
    return std::pair<BBox3f,BBox3f>(bounds0,bounds1);
  }

  void TriangleMeshTriangle1vMB::pack(char* dst, atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, void* geom) const 
  {
    const PrimRef& prim = *prims;
    const unsigned geomID = prim.geomID();
    const unsigned primID = prim.primID();
    const TriangleMeshScene::TriangleMesh* mesh = (TriangleMeshScene::TriangleMesh*) geom;
    const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(primID);
    const Vec3fa& a0 = mesh->vertex(tri.v[0],0);
    const Vec3fa& a1 = mesh->vertex(tri.v[0],1);
    const Vec3fa& b0 = mesh->vertex(tri.v[1],0);
    const Vec3fa& b1 = mesh->vertex(tri.v[1],1);
    const Vec3fa& c0 = mesh->vertex(tri.v[2],0);
    const Vec3fa& c1 = mesh->vertex(tri.v[2],1);
    new (dst) Triangle1vMB(a0,a1,b0,b1,c0,c1,mesh->id,primID,mesh->mask);
    prims++;
  }
  
  void TriangleMeshTriangle1vMB::pack(char* dst, const PrimRef* prims, size_t num, void* geom) const 
  {
    const PrimRef& prim = *prims;
    const unsigned geomID = prim.geomID();
    const unsigned primID = prim.primID();
    const TriangleMeshScene::TriangleMesh* mesh = (TriangleMeshScene::TriangleMesh*) geom;
    const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(primID);
    const Vec3fa& a0 = mesh->vertex(tri.v[0],0);
    const Vec3fa& a1 = mesh->vertex(tri.v[0],1);
    const Vec3fa& b0 = mesh->vertex(tri.v[1],0);
    const Vec3fa& b1 = mesh->vertex(tri.v[1],1);
    const Vec3fa& c0 = mesh->vertex(tri.v[2],0);
    const Vec3fa& c1 = mesh->vertex(tri.v[2],1);
    new (dst) Triangle1vMB(a0,a1,b0,b1,c0,c1,mesh->id,primID,mesh->mask);
    prims++;
  }

  std::pair<BBox3f,BBox3f> TriangleMeshTriangle1vMB::update2(char* prim, size_t num, void* geom) const 
  {
    BBox3f bounds0 = empty, bounds1 = empty;
    
    for (size_t j=0; j<num; j++) 
    {
      const Triangle1vMB& tri = ((Triangle1vMB*) prim)[j];
      bounds0.extend(merge(BBox3f(tri.v0),BBox3f(tri.v1),BBox3f(tri.v2)));
      bounds1.extend(merge(BBox3f(tri.v0+tri.d0),BBox3f(tri.v1+tri.d1),BBox3f(tri.v2+tri.d2)));
    }
    return std::pair<BBox3f,BBox3f>(bounds0,bounds1);
  }
}
