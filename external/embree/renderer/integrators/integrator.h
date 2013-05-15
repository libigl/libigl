// ======================================================================== //
// Copyright 2009-2012 Intel Corporation                                    //
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

#ifndef __EMBREE_INTEGRATOR_H__
#define __EMBREE_INTEGRATOR_H__

#include "rtcore/common/ray.h"
#include "api/scene.h"
#include "samplers/sampler.h"

namespace embree
{
  /*! Interface to different integrators. The task of the integrator
   *  is to compute the incoming radiance along some ray. */
  class Integrator : public RefCount {
    ALIGNED_CLASS
  public:

    /*! Integrators need a virtual destructor. */
    virtual ~Integrator() {}

    /*! Allows the integrator to register required samples at the sampler. */
    virtual void requestSamples(Ref<SamplerFactory>& samplerFactory, const Ref<BackendScene>& scene) {}

    /*! Computes the radiance arriving at the origin of the ray from
     *  the ray direction. */
    virtual Col3f Li(const Ray&               ray,     /*!< Ray to compute the radiance along.                */
                     const Ref<BackendScene>& scene,   /*!< Scene geometry and lights.                        */
                     const PrecomputedSample& sample,  /*!< Sampler used to generate (pseudo) random numbers. */
                     size_t&                  numRays  /*!< Used to count the number of rays shot.            */) = 0;
  };
}

#endif
