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

#ifndef __RTCORE_RAY_H__
#define __RTCORE_RAY_H__

/*! \ingroup embree_kernel_api */
/*! \{ */

/*! \brief Ray structure for an individual ray */
struct RTCORE_ALIGN(16)  RTCRay
{
  /* ray data */
public:
  float org[3];      //!< Ray origin
  int   align0;
  
  float dir[3];      //!< Ray direction
  int   align1;
  
  float tnear;       //!< Start of ray segment
  float tfar;        //!< End of ray segment (set to hit distance)

  float time;        //!< Time of this ray for motion blur
  int   mask;        //!< Used to mask out objects during traversal
  
  /* hit data */
public:
  float Ng[3];       //!< Unnormalized geometry normal
  int   align2;
  
  float u;           //!< Barycentric u coordinate of hit
  float v;           //!< Barycentric v coordinate of hit

  int   geomID;        //!< geometry ID
  int   primID;        //!< primitive ID
  int   instID;        //!< instance ID
};

/*! Ray structure for packets of 4 rays. */
struct RTCORE_ALIGN(16) RTCRay4
{
  /* ray data */
public:
  float orgx[4];  //!< x coordinate of ray origin
  float orgy[4];  //!< y coordinate of ray origin
  float orgz[4];  //!< z coordinate of ray origin
  
  float dirx[4];  //!< x coordinate of ray direction
  float diry[4];  //!< y coordinate of ray direction
  float dirz[4];  //!< z coordinate of ray direction
  
  float tnear[4]; //!< Start of ray segment 
  float tfar[4];  //!< End of ray segment (set to hit distance)

  float time[4];  //!< Time of this ray for motion blur
  int   mask[4];  //!< Used to mask out objects during traversal
  
  /* hit data */
public:
  float Ngx[4];   //!< x coordinate of geometry normal
  float Ngy[4];   //!< y coordinate of geometry normal
  float Ngz[4];   //!< z coordinate of geometry normal
  
  float u[4];     //!< Barycentric u coordinate of hit
  float v[4];     //!< Barycentric v coordinate of hit
  
  int   geomID[4];  //!< geometry ID
  int   primID[4];  //!< primitive ID
  int   instID[4];  //!< instance ID
};

/*! Ray structure for packets of 8 rays. */
struct RTCORE_ALIGN(32) RTCRay8
{
  /* ray data */
public:
  float orgx[8];  //!< x coordinate of ray origin
  float orgy[8];  //!< y coordinate of ray origin
  float orgz[8];  //!< z coordinate of ray origin
  
  float dirx[8];  //!< x coordinate of ray direction
  float diry[8];  //!< y coordinate of ray direction
  float dirz[8];  //!< z coordinate of ray direction
  
  float tnear[8]; //!< Start of ray segment 
  float tfar[8];  //!< End of ray segment (set to hit distance)

  float time[8];  //!< Time of this ray for motion blur
  int   mask[8];  //!< Used to mask out objects during traversal
  
  /* hit data */
public:
  float Ngx[8];   //!< x coordinate of geometry normal
  float Ngy[8];   //!< y coordinate of geometry normal
  float Ngz[8];   //!< z coordinate of geometry normal
  
  float u[8];     //!< Barycentric u coordinate of hit
  float v[8];     //!< Barycentric v coordinate of hit
  
  int   geomID[8];  //!< geometry ID
  int   primID[8];  //!< primitive ID
  int   instID[8];  //!< instance ID
};

/*! \brief Ray structure for packets of 16 rays. */
struct RTCORE_ALIGN(64) RTCRay16
{
  /* ray data */
public:
  float orgx[16];  //!< x coordinate of ray origin
  float orgy[16];  //!< y coordinate of ray origin
  float orgz[16];  //!< z coordinate of ray origin
  
  float dirx[16];  //!< x coordinate of ray direction
  float diry[16];  //!< y coordinate of ray direction
  float dirz[16];  //!< z coordinate of ray direction
  
  float tnear[16]; //!< Start of ray segment 
  float tfar[16];  //!< End of ray segment (set to hit distance)

  float time[16];  //!< Time of this ray for motion blur
  int   mask[16];  //!< Used to mask out objects during traversal
  
  /* hit data */
public:
  float Ngx[16];   //!< x coordinate of geometry normal
  float Ngy[16];   //!< y coordinate of geometry normal
  float Ngz[16];   //!< z coordinate of geometry normal
  
  float u[16];     //!< Barycentric u coordinate of hit
  float v[16];     //!< Barycentric v coordinate of hit
  
  int   geomID[16];  //!< geometry ID
  int   primID[16];  //!< primitive ID
  int   instID[16];  //!< instance ID
};

/*! @} */

#endif
