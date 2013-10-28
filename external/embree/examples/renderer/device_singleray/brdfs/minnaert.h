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

#ifndef __EMBREE_MINNAERT_BRDF_H__
#define __EMBREE_MINNAERT_BRDF_H__

#include "../brdfs/brdf.h"
#include "../brdfs/optics.h"

namespace embree
{
  /*! Minnaert BRDF. For backscattering materials (e.g. dust, fabric)
   *  The BRDF has a reflectance parameter that determines the color of
   *  the surface and a parameter to determine the amount of backscattering. */
  class Minnaert : public BRDF
  {
  public:

    /*! Minnaert BRDF constructor. This is a diffuse reflection BRDF. */
    __forceinline Minnaert(const Color& R, const float b) : BRDF(DIFFUSE_REFLECTION), R(R), b(b) {}

    __forceinline Color eval(const Vector3f& wo, const DifferentialGeometry& dg, const Vector3f& wi) const {
      float cosThetaI = clamp(dot(wi,dg.Ns));
      float backScatter = powf(clamp(dot(wo,wi)), b);
      return R * backScatter * cosThetaI / float(pi);
    }

    Color sample(const Vector3f& wo, const DifferentialGeometry& dg, Sample3f& wi, const Vec2f& s) const {
      return eval(wo, dg, wi = cosineSampleHemisphere(s.x,s.y,dg.Ns));
    }

    float pdf(const Vector3f& wo, const DifferentialGeometry& dg, const Vector3f& wi) const {
      return cosineSampleHemispherePDF(wi,dg.Ns);
    }

  private:

    /*! The reflectance parameter. The vale 0 means no reflection,
     *  and 1 means full reflection. */
    Color R;

    /*! The amount of backscattering. A value of 0 means lambertian
     *  diffuse, and inf means maximum backscattering. */
    float b;
  };
}

#endif
