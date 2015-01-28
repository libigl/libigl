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


#include "bvh16i.h"
#include "bvh16i_builder.h"

#include "common/accelinstance.h"
#include "geometry/triangle1.h"

namespace embree
{

  DECLARE_SYMBOL(Accel::Intersector16,BVH16iTriangle1Intersector16SingleMoeller);
  DECLARE_SYMBOL(Accel::Intersector1,BVH16iTriangle1Intersector1);

  void BVH16iRegister () 
  {
    int features = getCPUFeatures();
    
    /* default target */
    SELECT_SYMBOL_KNC(features,BVH16iTriangle1Intersector16SingleMoeller);
    SELECT_SYMBOL_KNC(features,BVH16iTriangle1Intersector1);
  }


  Accel::Intersectors BVH16iTriangle1Intersectors(BVH4i* bvh)
  {
    Accel::Intersectors intersectors;
    intersectors.ptr = bvh;
    intersectors.intersector1  = BVH16iTriangle1Intersector1;
    intersectors.intersector16 = BVH16iTriangle1Intersector16SingleMoeller;
    return intersectors;
  }

  Accel* BVH16i::BVH16iTriangle1ObjectSplitBinnedSAH(Scene* scene)
  { 
    BVH16i* accel = new BVH16i(SceneTriangle1::type);   
    Builder* builder = BVH16iBuilder::create(accel,&scene->flat_triangle_source_1,scene,BVH4iBuilder::BVH4I_BUILDER_DEFAULT);    
    Accel::Intersectors intersectors = BVH16iTriangle1Intersectors(accel);
    return new AccelInstance(accel,builder,intersectors);
  }


};
