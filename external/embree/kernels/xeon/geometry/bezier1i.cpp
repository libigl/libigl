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

#include "bezier1i.h"
#include "common/scene.h"

namespace embree
{
  SceneBezier1i SceneBezier1i::type;

  Bezier1iType::Bezier1iType () 
    : PrimitiveType("bezier1i",sizeof(Bezier1i),1,true,1) {} 
  
  size_t Bezier1iType::blocks(size_t x) const {
    return x;
  }
    
  size_t Bezier1iType::size(const char* This) const {
    return 1;
  }

  Bezier1iMBType Bezier1iMBType::type;

  Bezier1iMBType::Bezier1iMBType () 
    : PrimitiveType("bezier1imb",sizeof(Bezier1iMB),1,true,1) {} 
  
  size_t Bezier1iMBType::blocks(size_t x) const {
    return x;
  }
    
  size_t Bezier1iMBType::size(const char* This) const {
    return 1;
  }

  BBox3fa SceneBezier1i::update(char* prim_i, size_t num, void* geom) const 
  {
    BBox3fa bounds = empty;
    Scene* scene = (Scene*) geom;
    Bezier1i* prim = (Bezier1i*) prim_i;

    if (num == -1)
    {
      while (true)
      {
	const unsigned geomID = prim->geomID<1>();
	const unsigned primID = prim->primID<1>();
	const BezierCurves* curves = scene->getBezierCurves(geomID);
	const int vtx = curves->curve(primID);
	bounds.extend(curves->vertex(vtx+0));
	bounds.extend(curves->vertex(vtx+1));
	bounds.extend(curves->vertex(vtx+2));
	bounds.extend(curves->vertex(vtx+3));
	//prim->mask = curves->mask;
	const bool last = prim->last();
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
	const BezierCurves* curves = scene->getBezierCurves(geomID);
	const int vtx = curves->curve(primID);
	bounds.extend(curves->vertex(vtx+0));
	bounds.extend(curves->vertex(vtx+1));
	bounds.extend(curves->vertex(vtx+2));
	bounds.extend(curves->vertex(vtx+3));
	//prim->mask = curves->mask;
      }
    }
    return bounds; 
  }
}
