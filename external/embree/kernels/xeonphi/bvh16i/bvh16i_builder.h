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

#ifndef __EMBREE_BVH16I_BUILDER_MIC_H__
#define __EMBREE_BVH16I_BUILDER_MIC_H__

#include "bvh4i/bvh4i_builder.h"
#include "bvh16i.h"

namespace embree
{


  /*! derived binned-SAH builder supporting virtual geometry */  
  class BVH16iBuilder : public BVH4iBuilder
  {
  public:

    /*! creates the builder */
    static Builder* create (void* accel, BuildSource* source, void* geometry, size_t mode = BVH4I_BUILDER_DEFAULT);
    
  BVH16iBuilder(BVH16i* bvh, BuildSource* source, void* geometry) : BVH4iBuilder((BVH4i*)bvh,source,geometry) 
      {
	numNodesToAllocate = BVH4i::N; 
      }

    virtual void convertQBVHLayout(const size_t threadIndex, const size_t threadCount);
    virtual void printBuilderName();


  protected:
    TASK_FUNCTION(BVH16iBuilder,convertToSOALayoutMB);    

    void countLeaves(const size_t index);
    void getLeaves(unsigned int bvh4_ext_min, 
		   unsigned int node_index,
		   BVHNode *leaves, unsigned int &numLeaves);

    void convertBVH4iToBVH16i(const BVHNode *const bvh4,
			      const unsigned int bvh4_ext_min, 
			      const unsigned int bvh4_ext_max, 
			      BVH16i::Node *const bvh16,
			      size_t &index16,
			      unsigned int &parent_offset);

  };

}

#endif
