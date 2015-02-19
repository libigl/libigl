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

#pragma once

/*! \file shapesampler.isph Implements sampling functions for different
 *  geometric shapes. */

//inline float cos2sin(const float f) { return sqrt(max(0.f,1.f-f*f)); }
//inline float sin2cos(const float f) { return sqrt(max(0.f,1.f-f*f)); }

/*! Cosine weighted hemisphere sampling. Up direction is the z direction. */
inline Sample3f cosineSampleHemisphere(const float u, const float v) {
  const float phi = 2.0f * (float(pi)) * u;
  const float cosTheta = sqrt(v);
  const float sinTheta = sqrt(1.0f - v);
  return Sample3f(Vec3fa(cos(phi) * sinTheta, 
                                  sin(phi) * sinTheta, 
                                  cosTheta), 
                       cosTheta*(1.f/(float(pi))));
}

/*! Cosine weighted hemisphere sampling. Up direction is provided as argument. */
inline Sample3f cosineSampleHemisphere(const float  u, const float  v, const Vec3fa& N) 
{
  Sample3f s = cosineSampleHemisphere(u,v);
  return Sample3f(frame(N)*s.v,s.pdf);
}

  /*! Samples hemisphere with power cosine distribution. Up direction
   *  is the z direction. */
inline Sample3f powerCosineSampleHemisphere(const float u, const float v, const float _exp) 
{
  const float phi = 2.0f * (float(pi)) * u;
  const float cosTheta = pow(v,1.0f/(_exp+1.0f));
  const float sinTheta = cos2sin(cosTheta);
  return Sample3f(Vec3fa(cos(phi) * sinTheta, 
				   sin(phi) * sinTheta, 
				   cosTheta), 
                       (_exp+1.0f)*pow(cosTheta,_exp)*0.5f/(float(pi)));
}

/*! Computes the probability density for the power cosine sampling of the hemisphere. */
inline float powerCosineSampleHemispherePDF(const Vec3fa& s, const float _exp) {
  if (s.z < 0.f) return 0.f;
  return (_exp+1.0f)*pow(s.z,_exp)*0.5f/float(pi);
}

/*! Samples hemisphere with power cosine distribution. Up direction
 *  is provided as argument. */
inline Sample3f powerCosineSampleHemisphere(const float u, const float v, const Vec3fa& N, const float _exp) {
  Sample3f s = powerCosineSampleHemisphere(u,v,_exp);
  return Sample3f(frame(N)*s.v,s.pdf);
}

////////////////////////////////////////////////////////////////////////////////
/// Sampling of Spherical Cone
////////////////////////////////////////////////////////////////////////////////


/*! Uniform sampling of spherical cone. Cone direction is the z
 *  direction. */
inline Sample3f UniformSampleCone(const float u, const float v, const float angle) {
  const float phi = (float)(2.0f * float(pi)) * u;
  const float cosTheta = 1.0f - v*(1.0f - cos(angle));
  const float sinTheta = cos2sin(cosTheta);
  return Sample3f(Vec3fa(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta), 1.0f/((float)(4.0f*float(pi))*sqr(sin(0.5f*angle))));
}

/*! Computes the probability density of spherical cone sampling. */
inline float UniformSampleConePDF(const Vec3fa &s, const float angle) {
  return select(s.z < cos(angle), 0.0f, 1.0f/((float)(4.0f*float(pi))*sqr(sin(0.5f*angle))));
}

/*! Uniform sampling of spherical cone. Cone direction is provided as argument. */
inline Sample3f UniformSampleCone(const float u, const float v, const float angle, const Vec3fa& N) {
  Sample3f s = UniformSampleCone(u,v,angle);
  return Sample3f(frame(N)*s.v,s.pdf);
}

/*! Computes the probability density of spherical cone sampling. */
inline float UniformSampleConePDF(const Vec3fa &s, const float angle, const Vec3fa &N) {
  // return make_select(dot(s,N) < cos(angle), 0.0f, 1.0f/((float)(4.0f*float(pi))*sqr(sin(0.5f*angle))));
  if (dot(s,N) < cos(angle))
    return 0.f;
  else
    return 1.0f/((float)(4.0f*float(pi))*sqr(sin(0.5f*angle)));
}

////////////////////////////////////////////////////////////////////////////////
/// Sampling of Triangle
////////////////////////////////////////////////////////////////////////////////

/*! Uniform sampling of triangle. */
inline Vec3fa UniformSampleTriangle(const float u, const float v, const Vec3fa& A, const Vec3fa& B, const Vec3fa& C) {
  const float su = sqrt(u);
  return C + (1.0f-su)*(A-C) + (v*su)*(B-C);
}

////////////////////////////////////////////////////////////////////////////////
/// Sampling of Disk
////////////////////////////////////////////////////////////////////////////////

/*! Uniform sampling of disk. */
inline Vec2f UniformSampleDisk(const Vec2f &sample, const float radius) 
{
  const float r = sqrt(sample.x);
  const float theta = (2.f*float(pi)) * sample.y;
  return Vec2f(radius*r*cos(theta), radius*r*sin(theta));
}
