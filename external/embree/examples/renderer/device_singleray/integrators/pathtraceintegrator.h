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

#ifndef __EMBREE_PATH_TRACE_INTEGRATOR_H__
#define __EMBREE_PATH_TRACE_INTEGRATOR_H__

#include "integrators/integrator.h"
#include "renderers/renderer.h"
#include "image/image.h"

namespace embree
{
  /*! Path tracer integrator. The implementation follows a single path
   *  from the camera into the scene and connect the path at each
   *  diffuse or glossy surface to all light sources. Except for this
   *  the path is never split, also not at glass surfaces. */

  class PathTraceIntegrator : public Integrator
  {
    /*! Tracks the state of the path. */
    class LightPath
    {
    public:

      /*! Constructs a path. */
      __forceinline LightPath (const Ray& ray, const Medium& medium = Medium::Vacuum(), const int depth = 0,
                               const Color& throughput = one, const bool ignoreVisibleLights = false, const bool unbend = true)
        : lastRay(ray), lastMedium(medium), depth(depth), throughput(throughput), ignoreVisibleLights(ignoreVisibleLights), unbend(unbend) {}

      /*! Extends a light path. */
      __forceinline LightPath extended(const Ray& nextRay, const Medium& nextMedium, const Color& weight, const bool ignoreVL) const {
        return LightPath(nextRay, nextMedium, depth+1, throughput*weight, ignoreVL, unbend && (nextRay.dir == lastRay.dir));
      }

    public:
      Ray lastRay;                 /*! Last ray in the path. */
      Medium lastMedium;           /*! Medium the last ray travels inside. */
      uint32 depth;                /*! Recursion depth of path. */
      Color throughput;            /*! Determines the fraction of radiance reaches the pixel along the path. */
      bool ignoreVisibleLights;    /*! If the previous shade point used shadow rays we have to ignore the emission
                                       of geometrical lights to not double count them. */
      bool unbend;                 /*! True of the ray path is a straight line. */
    };

  public:

    /*! Construction of integrator from parameters. */
    PathTraceIntegrator(const Parms& parms);

    /*! Registers samples we need tom the sampler. */
    void requestSamples(Ref<SamplerFactory>& samplerFactory, const Ref<BackendScene>& scene);

    /*! Function that is recursively called to compute the path. */
    Color Li(LightPath& lightPath, const Ref<BackendScene>& scene, IntegratorState& state);

    /*! Computes the radiance arriving at the origin of the ray from the ray direction. */
    Color Li(Ray& ray, const Ref<BackendScene>& scene, IntegratorState& state);

    /* Configuration. */
  private:
    size_t maxDepth;               //!< Maximal recursion depth (1=primary ray only)
    float minContribution;         //!< Minimal contribution of a path to the pixel.
    float epsilon;                 //!< Epsilon to avoid self intersections.
    Ref<Image> backplate;          //!< High resolution background.

    /*! Random variables. */
  private:
    int lightSampleID;            //!< 2D random variable to sample the light source.
    int firstScatterSampleID;     //!< 2D random variable to sample the BRDF.
    int firstScatterTypeSampleID; //!< 1D random variable to sample the BRDF type to choose.
    std::vector<int> precomputedLightSampleID;  //!< ID of precomputed light samples for lights that need precomputations.
  };
}

#endif
