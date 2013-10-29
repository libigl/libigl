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

#ifndef __EMBREE_OPTICS_H__
#define __EMBREE_OPTICS_H__

#include "../default.h"

/*! \file optics.h This file collects all operations of light waves
 *  relevant for rendering, such as: computing reflection and
 *  refraction directions, or fresnel terms for different kind of
 *  media. */

namespace embree
{
  /*! Reflects a viewing vector V at a normal N. */
  __forceinline Sample3f reflect(const Vector3f& V, const Vector3f& N) {
    float cosi = dot(V,N);
    return Sample3f(2.0f*cosi*N-V, 1.0f);
  }

  /*! Reflects a viewing vector V at a normal N. Cosine between V
   *  and N is given as input. */
  __forceinline Sample3f reflect(const Vector3f& V, const Vector3f& N, float cosi) {
    return Sample3f(2.0f*cosi*N-V, 1.0f);
  }

  /*! Mirror a viewing vector V at a normal N. */
  __forceinline Vector3f mirror(const Vector3f& V, const Vector3f& N) {
    return V-2.0f*dot(V,N)*N;
  }

  /*! Refracts a viewing vector V at a normal N using the relative
   *  refraction index eta. Eta is refraction index of outside medium
   *  (where N points into) divided by refraction index of the inside
   *  medium. The vectors V and N have to point towards the same side
   *  of the surface. */
  __forceinline Sample3f refract(const Vector3f& V, const Vector3f& N, float eta)
  {
    float cosi = dot(V,N);
    float k = 1.0f-eta*eta*(1.0f-cosi*cosi);
    if (k < 0.0f) return Sample3f(Vector3f(zero),0.0f);
    float cost = sqrt(k);
    return Sample3f(eta*(cosi*N - V) - cost*N, sqr(eta));
  }

  /*! Refracts a viewing vector V at a normal N using the relative
   *  refraction index eta. Eta is refraction index of outside medium
   *  (where N points into) divided by refraction index of the inside
   *  medium. The vectors V and N have to point towards the same side
   *  of the surface. The cosine between V and N is given as input. */
  __forceinline Sample3f refract(const Vector3f& V, const Vector3f& N, float eta, float cosi)
  {
    float k = 1.0f-eta*eta*(1.0f-cosi*cosi);
    if (k < 0.0f) return Sample3f(Vector3f(zero),0.0f);
    float cost = sqrt(k);
    return Sample3f(eta*(cosi*N - V) - cost*N, sqr(eta));
  }

  /*!Refracts a viewing vector V at a normal N using the relative
   *  refraction index eta. Eta is refraction index of outside medium
   *  (where N points into) divided by refraction index of the inside
   *  medium. The vectors V and N have to point towards the same side
   *  of the surface. The cosine between V and N is given as input and
   *  the cosine of -N and transmission ray is computed as output. */
  __forceinline Sample3f refract(const Vector3f& V, const Vector3f& N, float eta, float cosi, float& cost)
  {
    float k = 1.0f-eta*eta*(1.0f-cosi*cosi);
    if (k < 0.0f) { cost = 0.0f; return Sample3f(Vector3f(zero),0.0f); }
    cost = sqrt(k);
    return Sample3f(eta*(cosi*N - V) - cost*N, sqr(eta));
  }

  /*! Computes fresnel coefficient for a media interface with ouside
   *  refraction index of etai and inside refraction index of
   *  etat. Both cosines have to be positive. */
  __forceinline float fresnelDielectric(float cosi, float cost, float etai, float etat)
  {
    float Rper = (etai*cosi - etat*cost) * rcp(etai*cosi + etat*cost);
    float Rpar = (etat*cosi - etai*cost) * rcp(etat*cosi + etai*cost);
    return 0.5f*(Rpar*Rpar + Rper*Rper);
  }

  /*! Computes fresnel coefficient for media interface with relative
   *  refraction index eta. Eta is the outside refraction index
   *  divided by the inside refraction index. Both cosines have to be
   *  positive. */
  __forceinline float fresnelDielectric(float cosi, float cost, float eta)
  {
    float Rper = (eta*cosi -     cost) * rcp(eta*cosi +     cost);
    float Rpar = (    cosi - eta*cost) * rcp(    cosi + eta*cost);
    return 0.5f*(Rpar*Rpar + Rper*Rper);
  }

  /*! Computes fresnel coefficient for media interface with relative
   *  refraction index eta. Eta is the outside refraction index
   *  divided by the inside refraction index. The cosine has to be
   *  positive. */
  __forceinline float fresnelDielectric(float cosi, float eta)
  {
    float k = 1.0f-eta*eta*(1.0f-cosi*cosi);
    if (k < 0.0f) return 1.0f;
    float cost = sqrt(k);
    return fresnelDielectric(cosi, cost, eta);
  }

  /*! Computes fresnel coefficient for conductor medium with complex
   *  refraction index (eta,k). The cosine has to be positive. */
  __forceinline Color fresnelConductor(float cosi, const Color& eta, const Color& k)
  {
    Color tmp = eta*eta + k*k;
    Color Rpar = (tmp * cosi*cosi - 2.0f*eta*cosi + Color(one)) *
      rcp(tmp * cosi*cosi + 2.0f*eta*cosi + Color(one));
    Color Rper = (tmp - 2.0f*eta*cosi + Color(cosi*cosi)) *
      rcp(tmp + 2.0f*eta*cosi + Color(cosi*cosi));
    return 0.5f * (Rpar + Rper);
  }
}

#endif
