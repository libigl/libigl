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

#include "triangle4v.h"
#include "common/scene.h"

namespace embree
{
  SceneTriangle4v SceneTriangle4v::type;
  TriangleMeshTriangle4v TriangleMeshTriangle4v::type;

  Triangle4vType::Triangle4vType () 
  : PrimitiveType("triangle4v",sizeof(Triangle4v),4,false,1) {} 
  
  size_t Triangle4vType::blocks(size_t x) const {
    return (x+3)/4;
  }
  
  size_t Triangle4vType::size(const char* This) const {
    return ((Triangle4v*)This)->size();
  }

  void SceneTriangle4v::pack(char* This, atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, void* geom) const 
  {
    Scene* scene = (Scene*) geom;
    
    ssei geomID = -1, primID = -1, mask = -1;
    sse3f v0 = zero, v1 = zero, v2 = zero;
    
    for (size_t i=0; i<4 && prims; i++, prims++)
    {
      const PrimRef& prim = *prims;
      const TriangleMeshScene::TriangleMesh* mesh = scene->getTriangleMesh(prim.geomID());
      const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(prim.primID());
      const Vec3fa& p0 = mesh->vertex(tri.v[0]);
      const Vec3fa& p1 = mesh->vertex(tri.v[1]);
      const Vec3fa& p2 = mesh->vertex(tri.v[2]);
      geomID [i] = prim.geomID();
      primID [i] = prim.primID();
      mask   [i] = mesh->mask;
      v0.x[i] = p0.x; v0.y[i] = p0.y; v0.z[i] = p0.z;
      v1.x[i] = p1.x; v1.y[i] = p1.y; v1.z[i] = p1.z;
      v2.x[i] = p2.x; v2.y[i] = p2.y; v2.z[i] = p2.z;
    }
    new (This) Triangle4v(v0,v1,v2,geomID,primID,mask);
  }
  
  void SceneTriangle4v::pack(char* dst, const PrimRef* prims, size_t num, void* geom) const 
  {
    Scene* scene = (Scene*) geom;
    
    size_t p = 0;
    while (p < num)
    {
      ssei geomID = -1, primID = -1, mask = -1;
      sse3f v0 = zero, v1 = zero, v2 = zero;
      
      for (size_t i=0; i<4 && p < num; i++, p++)
      {
        const PrimRef& prim = prims[p];
        const TriangleMeshScene::TriangleMesh* mesh = scene->getTriangleMesh(prim.geomID());
        const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(prim.primID());
        const Vec3fa& p0 = mesh->vertex(tri.v[0]);
        const Vec3fa& p1 = mesh->vertex(tri.v[1]);
        const Vec3fa& p2 = mesh->vertex(tri.v[2]);
        geomID [i] = prim.geomID();
        primID [i] = prim.primID();
        mask   [i] = mesh->mask;
        v0.x[i] = p0.x; v0.y[i] = p0.y; v0.z[i] = p0.z;
        v1.x[i] = p1.x; v1.y[i] = p1.y; v1.z[i] = p1.z;
        v2.x[i] = p2.x; v2.y[i] = p2.y; v2.z[i] = p2.z;
      }
      new (dst) Triangle4v(v0,v1,v2,geomID,primID,mask);
      dst += sizeof(Triangle4v);
    }
  }
  
  BBox3f SceneTriangle4v::update(char* prim, size_t num, void* geom) const 
  {
    BBox3f bounds = empty;
    Scene* scene = (Scene*) geom;
    
    for (size_t j=0; j<num; j++) 
    {
      Triangle4v& dst = ((Triangle4v*) prim)[j];
      
      ssei vgeomID = -1, vprimID = -1, vmask = -1;
      sse3f v0 = zero, v1 = zero, v2 = zero;
      
      for (size_t i=0; i<4; i++)
      {
        if (dst.primID[i] == -1) break;
        const unsigned geomID = dst.geomID[i];
        const unsigned primID = dst.primID[i];
        const TriangleMeshScene::TriangleMesh* mesh = scene->getTriangleMesh(geomID);
        const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(primID);
        const Vec3fa p0 = mesh->vertex(tri.v[0]);
        const Vec3fa p1 = mesh->vertex(tri.v[1]);
        const Vec3fa p2 = mesh->vertex(tri.v[2]);
        bounds.extend(merge(BBox3f(p0),BBox3f(p1),BBox3f(p2)));
        vgeomID [i] = geomID;
        vprimID [i] = primID;
        vmask   [i] = mesh->mask;
        v0.x[i] = p0.x; v0.y[i] = p0.y; v0.z[i] = p0.z;
        v1.x[i] = p1.x; v1.y[i] = p1.y; v1.z[i] = p1.z;
        v2.x[i] = p2.x; v2.y[i] = p2.y; v2.z[i] = p2.z;
      }
      new (&dst) Triangle4v(v0,v1,v2,vgeomID,vprimID,vmask);
    }
    return bounds; 
  }

  void TriangleMeshTriangle4v::pack(char* This, atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, void* geom) const 
  {
    TriangleMeshScene::TriangleMesh* mesh = (TriangleMeshScene::TriangleMesh*) geom;
    
    ssei geomID = -1, primID = -1, mask = -1;
    sse3f v0 = zero, v1 = zero, v2 = zero;
    
    for (size_t i=0; i<4 && prims; i++, prims++)
    {
      const PrimRef& prim = *prims;
      const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(prim.primID());
      const Vec3fa& p0 = mesh->vertex(tri.v[0]);
      const Vec3fa& p1 = mesh->vertex(tri.v[1]);
      const Vec3fa& p2 = mesh->vertex(tri.v[2]);
      geomID [i] = mesh->id;
      primID [i] = prim.primID();
      mask   [i] = mesh->mask;
      v0.x[i] = p0.x; v0.y[i] = p0.y; v0.z[i] = p0.z;
      v1.x[i] = p1.x; v1.y[i] = p1.y; v1.z[i] = p1.z;
      v2.x[i] = p2.x; v2.y[i] = p2.y; v2.z[i] = p2.z;
    }
    new (This) Triangle4v(v0,v1,v2,geomID,primID,mask);
  }
  
  BBox3f TriangleMeshTriangle4v::update(char* prim, size_t num, void* geom) const 
  {
    BBox3f bounds = empty;
    TriangleMeshScene::TriangleMesh* mesh = (TriangleMeshScene::TriangleMesh*) geom;
    
    for (size_t j=0; j<num; j++) 
    {
      Triangle4v& dst = ((Triangle4v*) prim)[j];
      
      ssei vgeomID = -1, vprimID = -1, vmask = -1;
      sse3f v0 = zero, v1 = zero, v2 = zero;

      for (size_t i=0; i<4; i++)
      {
        if (dst.primID[i] == -1) break;
        const unsigned geomID = dst.geomID[i];
        const unsigned primID = dst.primID[i];
        const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(primID);
        const Vec3fa p0 = mesh->vertex(tri.v[0]);
        const Vec3fa p1 = mesh->vertex(tri.v[1]);
        const Vec3fa p2 = mesh->vertex(tri.v[2]);
        bounds.extend(merge(BBox3f(p0),BBox3f(p1),BBox3f(p2)));
        vgeomID [i] = geomID;
        vprimID [i] = primID;
        vmask   [i] = mesh->mask;
        v0.x[i] = p0.x; v0.y[i] = p0.y; v0.z[i] = p0.z;
        v1.x[i] = p1.x; v1.y[i] = p1.y; v1.z[i] = p1.z;
        v2.x[i] = p2.x; v2.y[i] = p2.y; v2.z[i] = p2.z;
        }
      new (&dst) Triangle4v(v0,v1,v2,vgeomID,vprimID,vmask);
    }
    return bounds; 
  }
}
