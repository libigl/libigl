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

#ifndef __EMBREE_VIRTUAL_ACCEL_H__
#define __EMBREE_VIRTUAL_ACCEL_H__

#include "embree2/rtcore.h"
#include "common/accel.h"
#include "common/accelset.h"
#include "common/buildsource.h"
#include "geometry/primitive.h"

namespace embree
{
  struct VirtualAccel : public Accel
  {
    struct VirtualBuildSource : public BuildSource
    {
      VirtualBuildSource (std::vector<AccelSet*>& accels) 
        : accels(accels) {}
      
      bool isEmpty () const { 
        return accels.size() == 0;
      }
      
      size_t groups () const { 
        return accels.size();
      }
      
      size_t prims (size_t group, size_t* pNumVertices) const {
        if (pNumVertices) *pNumVertices = 0;
        return accels[group]->size();
      }
      
      const BBox3f bounds(size_t group, size_t prim) const {
        return accels[group]->bounds(prim);
      }

      void bounds(size_t group, size_t begin, size_t end, BBox3f* bounds_o) const {
        assert(false); // FIXME
      }
      
      std::vector<AccelSet*>& accels;
    };
    
    struct VirtualAccelObjectType : public PrimitiveType
    {
      VirtualAccelObjectType () 
        : PrimitiveType("object",sizeof(AccelSetItem),1,false,1) {} 
      
      size_t blocks(size_t x) const {
        return x;
      }
      
      size_t size(const char* This) const {
        return 1;
      }
      
      void pack(char* This, atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, void* geom) const 
      {
        const PrimRef& prim = *prims;
        std::vector<AccelSet*>* accels = (std::vector<AccelSet*>*) geom;
        AccelSetItem* dst = (AccelSetItem*)This;
        dst->accel = (*accels)[prim.geomID()];
        dst->item = prim.primID();
        prims++;
      }
    };
    
  public:
    VirtualAccel (const std::string& ty, std::vector<AccelSet*>& accels);
    ~VirtualAccel ();
    
  public:
    void build (size_t threadIndex, size_t threadCount);
    
  public:
    Bounded* accel;
    Builder* builder;
    VirtualBuildSource source;
  };
}

#endif
