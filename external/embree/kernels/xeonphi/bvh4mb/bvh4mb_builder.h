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

#pragma once

#include "bvh4i/bvh4i_builder.h"
#include "bvh4mb.h"

namespace embree
{

  /*! derived binned-SAH builder supporting virtual geometry */  
  class BVH4mbBuilder : public BVH4iBuilder
  {
  public:

    /*! creates the builder */
    static Builder* create (void* accel, void* geometry, size_t mode = BVH4I_BUILDER_DEFAULT);

    BVH4mbBuilder (BVH4mb* bvh, void* geometry) : BVH4iBuilder((BVH4i*)bvh,geometry,sizeof(BVH4mb::Node)) 
    {
    }
    virtual void computePrimRefs  (const size_t threadIndex, const size_t threadCount);
    virtual void allocateData     (const size_t threadCount, const size_t newNumPrimitives);
    virtual void finalize         (const size_t threadIndex, const size_t threadCount);
    virtual void createAccel      (const size_t threadIndex, const size_t threadCount);
    virtual void printBuilderName();
    virtual size_t getNumPrimitives();
    virtual std::string getStatistics();


    /* parallel refit bvh4mb tree */
    void generate_subtrees(const BVH4i::NodeRef &ref,const size_t depth, size_t &subtrees);
    BBox3fa refit_toplevel(const BVH4i::NodeRef &ref,const size_t depth);

    /* scalar refit */
    BBox3fa refit(const BVH4i::NodeRef &ref);

    /* check bvh4mb tree */
    BBox3fa check_tree(const BVH4i::NodeRef &ref);

  protected:
    AtomicCounter atomicID;
    size_t subtrees;

    TASK_FUNCTION(BVH4mbBuilder,refitBVH4MB);    
    TASK_FUNCTION(BVH4mbBuilder,createTriangle01AccelMB);    
    TASK_FUNCTION(BVH4mbBuilder,computePrimRefsTrianglesMB);    

  };

}
