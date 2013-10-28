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

#ifndef __EMBREE_DIELECTRIC_LAYER_BRDF_H__
#define __EMBREE_DIELECTRIC_LAYER_BRDF_H__

#include "../brdfs/brdf.h"
#include "../brdfs/optics.h"

namespace embree
{
  /*! Dielectric layer BRDF model. Models a ground BRDF covered by a
   *  dielectric layer. The model works by refracting the incoming and
   *  outgoing light direction into the medium and then evaluating the
   *  ground BRDF. The model also supports a transmission coefficient
   *  of the dielectricum. */
  template<typename Ground>
    class DielectricLayer : public BRDF
  {
  public:

    /*! Dielectric layer BRDF constructor. \param T is the
     *  transmission coefficient of the dielectricum \param etai is
     *  the refraction index of the medium the incident ray travels in
     *  \param etat is the refraction index of the opposite medium \param ground is the ground BRDF to use */
    __forceinline DielectricLayer(const Color& T, float etai, float etat, const Ground& ground)
      : BRDF(ground.type), T(T), etait(etai*rcp(etat)), etati(etat*rcp(etai)), ground(ground) {}

    __forceinline Color eval(const Vector3f& wo, const DifferentialGeometry& dg, const Vector3f& wi) const
    {
      float cosThetaO = dot(wo,dg.Ns);
      float cosThetaI = dot(wi,dg.Ns);
      if (cosThetaI <= 0.0f || cosThetaO <= 0.0f) return zero;
      float cosThetaO1; Vector3f wo1 = refract(wo,dg.Ns,etait,cosThetaO,cosThetaO1);
      float cosThetaI1; Vector3f wi1 = refract(wi,dg.Ns,etait,cosThetaI,cosThetaI1);
      float Fi = 1.0f - fresnelDielectric(cosThetaI,cosThetaI1,etait);
      Color Fg = ground.eval(-wo1,dg,-wi1);
      float Fo = 1.0f - fresnelDielectric(cosThetaO,cosThetaO1,etait);
      return Fo * T * Fg * T * Fi;
    }

    Color sample(const Vector3f& wo, const DifferentialGeometry& dg, Sample3f& wi, const Vec2f& s) const
    {
      /*! refract ray into medium */
      float cosThetaO = dot(wo,dg.Ns);
      if (cosThetaO <= 0.0f) return zero;
      float cosThetaO1; Sample3f wo1 = refract(wo,dg.Ns,etait,cosThetaO,cosThetaO1);

      /*! sample ground BRDF */
      Sample3f wi1(zero); Color Fg = ground.sample(-wo1.value,dg,wi1,s);

      /*! refract ray out of medium */
      float cosThetaI1 = dot(wi1.value,dg.Ns);
      if (cosThetaI1 <= 0.0f) return zero;

      float cosThetaI; Sample3f wi0 = refract(-wi1.value,-dg.Ns,etati,cosThetaI1,cosThetaI);
      if (wi0.pdf == 0.0f) return zero;

      /*! accumulate contribution of path */
      wi = Sample3f(wi0,wi1.pdf);
      float Fi = 1.0f - fresnelDielectric(cosThetaI,cosThetaI1,etait);
      float Fo = 1.0f - fresnelDielectric(cosThetaO,cosThetaO1,etait);
      return Fo * T * Fg * T * Fi;
    }

    float pdf(const Vector3f& wo, const DifferentialGeometry& dg, const Vector3f& wi) const
    {
      float cosThetaO = dot(wo,dg.Ns);
      if (cosThetaO <= 0.0f) return 0.0f;
      Sample3f wo1 = refract(wo,dg.Ns,etait,cosThetaO);
      Sample3f wi1(zero); float p = ground.pdf(-wo1.value,dg,wi1);
      Sample3f wi0 = refract(-wi1.value,-dg.Ns,etati);
      if (wi0.pdf == 0.0f) return zero;
      return p;
    }

  private:
    Color T;             //!< Transmission coefficient of dielectricum
    float etait;         //!< Relative refraction index etai/etat of both media
    float etati;         //!< relative refraction index etat/etai of both media
    Ground ground;       //!< the BRDF of the ground layer
  };
}

#endif
