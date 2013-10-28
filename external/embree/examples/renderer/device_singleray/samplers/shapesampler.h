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

#ifndef __EMBREE_SHAPESAMPLER_H__
#define __EMBREE_SHAPESAMPLER_H__

#include "../default.h"

/*! \file shapesampler.h Implements sampling functions for different
 *  geometric shapes. */

namespace embree
{
  ////////////////////////////////////////////////////////////////////////////////
  /// Sampling of Sphere
  ////////////////////////////////////////////////////////////////////////////////

  /*! Uniform sphere sampling. */
  __forceinline Sample3f uniformSampleSphere(const float& u, const float& v) {
    const float phi = float(two_pi) * u;
    const float cosTheta = 1.0f - 2.0f * v, sinTheta = 2.0f * sqrt(v * (1.0f - v));
    return Sample3f(Vector3f(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta), float(one_over_four_pi));
  }

  /*! Computes the probability density for the uniform shere sampling. */
  __forceinline float uniformSampleSpherePDF() {
    return float(one_over_four_pi);
  }

  /*! Cosine weighted sphere sampling. Up direction is the z direction. */
  __forceinline Sample3f cosineSampleSphere(const float& u, const float& v) {
    const float phi = float(two_pi) * u;
    const float vv = 2.0f*(v-0.5f);
    const float cosTheta = sign(vv)*sqrt(abs(vv)), sinTheta = cos2sin(cosTheta);
    return Sample3f(Vector3f(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta), 2.0f*cosTheta*float(one_over_pi));
  }

  /*! Computes the probability density for the cosine weighted sphere sampling. */
  __forceinline float cosineSampleSpherePDF(const Vector3f& s) {
    return 2.0f*abs(s.z)*float(one_over_pi);
  }

  /*! Cosine weighted sphere sampling. Up direction is provided as argument. */
  __forceinline Sample3f cosineSampleSphere(const float& u, const float& v, const Vector3f& N) {
    Sample3f s = cosineSampleSphere(u,v);
    return Sample3f(frame(N)*Vector3f(s),s.pdf);
  }
  
  /*! Computes the probability density for the cosine weighted sphere sampling. */
  __forceinline float cosineSampleSpherePDF(const Vector3f& s, const Vector3f& N) {
    return 2.0f*abs(dot(s,N))*float(one_over_pi);
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Sampling of Hemisphere
  ////////////////////////////////////////////////////////////////////////////////

  /*! Uniform hemisphere sampling. Up direction is the z direction. */
  __forceinline Sample3f uniformSampleHemisphere(const float& u, const float& v) {
    const float phi = float(two_pi) * u;
    const float cosTheta = v, sinTheta = cos2sin(v);
    return Sample3f(Vector3f(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta), float(one_over_two_pi));
  }

  /*! Computes the probability density for the uniform hemisphere sampling. */
  __forceinline float uniformSampleHemispherePDF(const Vector3f& s) {
    return select(s.z < 0.0f, 0.0f, float(one_over_two_pi));
  }

  /*! Uniform hemisphere sampling. Up direction is provided as argument. */
  __forceinline Sample3f uniformSampleHemisphere(const float& u, const float& v, const Vector3f& N) {
    Sample3f s = uniformSampleHemisphere(u,v);
    return Sample3f(frame(N)*Vector3f(s),s.pdf);
  }

  /*! Computes the probability density for the uniform hemisphere sampling. */
  __forceinline float uniformSampleHemispherePDF(const Vector3f& s, const Vector3f& N) {
    return select(dot(s,N) < 0.0f, 0.0f, float(one_over_two_pi));
  }

  /*! Cosine weighted hemisphere sampling. Up direction is the z direction. */
  __forceinline Sample3f cosineSampleHemisphere(const float& u, const float& v) {
    const float phi = float(two_pi) * u;
    const float cosTheta = sqrt(v), sinTheta = sqrt(1.0f - v);
    return Sample3f(Vector3f(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta), cosTheta*float(one_over_pi));
  }

  /*! Computes the probability density for the cosine weighted hemisphere sampling. */
  __forceinline float cosineSampleHemispherePDF(const Vector3f& s) {
    return select(s.z < 0.0f, 0.0f, s.z*float(one_over_pi));
  }

  /*! Cosine weighted hemisphere sampling. Up direction is provided as argument. */
  __forceinline Sample3f cosineSampleHemisphere(const float& u, const float& v, const Vector3f& N) {
    Sample3f s = cosineSampleHemisphere(u,v);
    return Sample3f(frame(N)*Vector3f(s),s.pdf);
  }

  /*! Computes the probability density for the cosine weighted hemisphere sampling. */
  __forceinline float cosineSampleHemispherePDF(const Vector3f& s, const Vector3f& N) {
    return select(dot(s,N) < 0.0f, 0.0f, dot(s,N)*float(one_over_pi));
  }

  /*! Samples hemisphere with power cosine distribution. Up direction
   *  is the z direction. */
  __forceinline Sample3f powerCosineSampleHemisphere(const float& u, const float& v, float exp) {
    const float phi = float(two_pi) * u;
    const float cosTheta = pow(v,rcp(exp+1));
    const float sinTheta = cos2sin(cosTheta);
    return Sample3f(Vector3f(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta), (exp+1.0f)*pow(cosTheta,exp)*float(one_over_two_pi));
  }

  /*! Computes the probability density for the power cosine sampling of the hemisphere. */
  __forceinline float powerCosineSampleHemispherePDF(const Vector3f& s, float exp) {
    return select(s.z < 0.0f, 0.0f, (exp+1.0f)*pow(s.z,exp)*float(one_over_two_pi));
  }

  /*! Samples hemisphere with power cosine distribution. Up direction
   *  is provided as argument. */
  __forceinline Sample3f powerCosineSampleHemisphere(const float& u, const float& v, const Vector3f& N, float exp) {
    Sample3f s = powerCosineSampleHemisphere(u,v,exp);
    return Sample3f(frame(N)*Vector3f(s),s.pdf);
  }

  /*! Computes the probability density for the power cosine sampling of the hemisphere. */
  __forceinline float powerCosineSampleHemispherePDF(const Vector3f& s, const Vector3f& N, float exp) {
    return select(dot(s,N) < 0.0f, 0.0f, (exp+1.0f)*pow(dot(s,N),exp)*float(one_over_two_pi));
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Sampling of Spherical Cone
  ////////////////////////////////////////////////////////////////////////////////

  /*! Uniform sampling of spherical cone. Cone direction is the z
   *  direction. */
  __forceinline Sample3f uniformSampleCone(const float& u, const float& v, const float& angle) {
    const float phi = float(two_pi) * u;
    const float cosTheta = 1.0f - v*(1.0f - cos(angle));
    const float sinTheta = cos2sin(cosTheta);
    return Sample3f(Vector3f(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta), rcp(float(four_pi)*sqr(sin(0.5f*angle))));
  }

  /*! Computes the probability density of uniform spherical cone sampling. */
  __forceinline float uniformSampleConePDF(const Vector3f& s, const float& angle) {
    return select(s.z < cos(angle), 0.0f, rcp(float(four_pi)*sqr(sin(0.5f*angle))));
  }

  /*! Uniform sampling of spherical cone. Cone direction is provided as argument. */
  __forceinline Sample3f uniformSampleCone(const float& u, const float& v, const float& angle, const Vector3f& N) {
    Sample3f s = uniformSampleCone(u,v,angle);
    return Sample3f(frame(N)*Vector3f(s),s.pdf);
  }

  /*! Computes the probability density of uniform spherical cone sampling. */
  __forceinline float uniformSampleConePDF(const Vector3f& s, const float& angle, const Vector3f& N) {
    return select(dot(s,N) < cos(angle), 0.0f, rcp(float(four_pi)*sqr(sin(0.5f*angle))));
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Sampling of Triangle
  ////////////////////////////////////////////////////////////////////////////////

  /*! Uniform sampling of triangle. */
  __forceinline Vector3f uniformSampleTriangle(const float& u, const float& v, const Vector3f& A, const Vector3f& B, const Vector3f& C) {
    float su = sqrt(u);
    return Vector3f(C+(1.0f-su)*(A-C)+(v*su)*(B-C));
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Sampling of Disk
  ////////////////////////////////////////////////////////////////////////////////

  /*! Uniform sampling of disk. */
  __forceinline Vec2f uniformSampleDisk(const Vec2f& sample, const float radius) {
    const float r = sqrtf(sample.x);
    const float theta = float(two_pi) * sample.y;
    return Vec2f(radius*r*cosf(theta), radius*r*sinf(theta));
  }
}

#endif
