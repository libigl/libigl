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

#include "bvh4hair.h"
#include "bvh4hair_builder.h"
#include "common/accelinstance.h"
#include "geometry/triangle1.h"

namespace embree
{
#define DBG(x) x

  DECLARE_SYMBOL(Accel::Intersector1 ,BVH4HairIntersector1Bezier1i);
  DECLARE_SYMBOL(Accel::Intersector1 ,BVH4HairIntersector1Bezier1iNoFilter);
  DECLARE_SYMBOL(Accel::Intersector16,BVH4HairIntersector16Bezier1i);
  DECLARE_SYMBOL(Accel::Intersector16,BVH4HairIntersector16Bezier1iNoFilter);

  void BVH4HairRegister () 
  {
    int features = getCPUFeatures();

    SELECT_SYMBOL_KNC(features,BVH4HairIntersector1Bezier1i);
    SELECT_SYMBOL_KNC(features,BVH4HairIntersector1Bezier1iNoFilter);
    SELECT_SYMBOL_KNC(features,BVH4HairIntersector16Bezier1i);
    SELECT_SYMBOL_KNC(features,BVH4HairIntersector16Bezier1iNoFilter);



  }

  Accel::Intersectors BVH4HairIntersectors(BVH4Hair* bvh)
  {
    Accel::Intersectors intersectors;
    intersectors.ptr = bvh;
    intersectors.intersector1           = BVH4HairIntersector1Bezier1i;
    intersectors.intersector16          = BVH4HairIntersector16Bezier1i;
    intersectors.intersector16_filter   = BVH4HairIntersector16Bezier1i;
    intersectors.intersector16_nofilter = BVH4HairIntersector16Bezier1iNoFilter;

    return intersectors;
  }

  Accel* BVH4Hair::BVH4HairBinnedSAH(Scene* scene)
  {
    BVH4Hair* accel  = new BVH4Hair(SceneTriangle1::type,scene);    
    Builder* builder = new BVH4HairBuilder(accel,scene);   
    Accel::Intersectors intersectors = BVH4HairIntersectors(accel);
    return new AccelInstance(accel,builder,intersectors);    
  }

  // =================================================================================
  // =================================================================================
  // =================================================================================


  __aligned(64) float BVH4Hair::UnalignedNode::identityMatrix[16] = {
    1,0,0,0,
    0,1,0,0,
    0,0,1,0,
    0,0,0,1
  };

  __aligned(64) float BVH4Hair::UnalignedNode::invalidMatrix[16] = {
    NAN,0,0,0,
    0,NAN,0,0,
    0,0,NAN,0,
    NAN,NAN,NAN,NAN
  };




};
