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

#include "triangle1v.h"
#include "common/scene.h"

namespace embree
{
  Triangle1vType Triangle1vType::type;
  TriangleMeshTriangle1v TriangleMeshTriangle1v::type;

  Triangle1vMBType Triangle1vMBType::type;
  TriangleMeshTriangle1vMB TriangleMeshTriangle1vMB::type;

  Triangle1vType::Triangle1vType () 
      : PrimitiveType("triangle1v",sizeof(Triangle1v),1,false,1) {} 
  
  size_t Triangle1vType::blocks(size_t x) const {
    return x;
  }
    
  size_t Triangle1vType::size(const char* This) const {
    return 1;
  }
  
  BBox3fa TriangleMeshTriangle1v::update(char* prim_i, size_t num, void* geom) const 
  {
    BBox3fa bounds = empty;
    const TriangleMesh* mesh = (const TriangleMesh*) geom;
    Triangle1v* prim = (Triangle1v*) prim_i;

    if (num == -1)
    {
      while (true)
      {
	const unsigned geomID = prim->geomID<1>();
	const unsigned primID = prim->primID<1>();
	const TriangleMesh::Triangle& tri = mesh->triangle(primID);
	const Vec3fa v0 = mesh->vertex(tri.v[0]);
	const Vec3fa v1 = mesh->vertex(tri.v[1]);
	const Vec3fa v2 = mesh->vertex(tri.v[2]);
	const bool last = prim->last();
	new (prim) Triangle1v(v0,v1,v2,geomID,primID,mesh->mask,last);
	bounds.extend(merge(BBox3fa(v0),BBox3fa(v1),BBox3fa(v2)));
	if (last) break;
	prim++;
      }
    }
    else
    {
      for (size_t i=0; i<num; i++, prim++)
      {
	const unsigned geomID = prim->geomID<0>();
	const unsigned primID = prim->primID<0>();
	const TriangleMesh::Triangle& tri = mesh->triangle(primID);
	const Vec3fa v0 = mesh->vertex(tri.v[0]);
	const Vec3fa v1 = mesh->vertex(tri.v[1]);
	const Vec3fa v2 = mesh->vertex(tri.v[2]);
	new (prim) Triangle1v(v0,v1,v2,geomID,primID,mesh->mask,false);
	bounds.extend(merge(BBox3fa(v0),BBox3fa(v1),BBox3fa(v2)));
      }
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

  std::pair<BBox3fa,BBox3fa> Triangle1vMBType::update2(char* prim, size_t num, void* geom) const 
  {
    BBox3fa bounds0 = empty, bounds1 = empty;
    
    for (size_t j=0; j<num; j++) 
    {
      const Triangle1vMB& tri = ((Triangle1vMB*) prim)[j];
      bounds0.extend(merge(BBox3fa(tri.v0),BBox3fa(tri.v1),BBox3fa(tri.v2)));
      bounds1.extend(merge(BBox3fa(tri.v0+tri.d0),BBox3fa(tri.v1+tri.d1),BBox3fa(tri.v2+tri.d2)));
    }
    return std::pair<BBox3fa,BBox3fa>(bounds0,bounds1);
  }

  std::pair<BBox3fa,BBox3fa> TriangleMeshTriangle1vMB::update2(char* prim_i, size_t num, void* geom) const 
  {
    BBox3fa bounds0 = empty, bounds1 = empty;
    Triangle1vMB* prim = (Triangle1vMB*) prim_i;

    if (num == -1)
    {
      while (true)
      {
	bounds0.extend(merge(BBox3fa(prim->v0),BBox3fa(prim->v1),BBox3fa(prim->v2)));
	bounds1.extend(merge(BBox3fa(prim->v0+prim->d0),BBox3fa(prim->v1+prim->d1),BBox3fa(prim->v2+prim->d2)));
	const bool last = prim->last();
	if (last) break;
	prim++;
      }
    }
    else
    {
      for (size_t i=0; i<num; i++, prim++)
      {
	bounds0.extend(merge(BBox3fa(prim->v0),BBox3fa(prim->v1),BBox3fa(prim->v2)));
	bounds1.extend(merge(BBox3fa(prim->v0+prim->d0),BBox3fa(prim->v1+prim->d1),BBox3fa(prim->v2+prim->d2)));
      }
    }
    return std::pair<BBox3fa,BBox3fa>(bounds0,bounds1);
  }
}
