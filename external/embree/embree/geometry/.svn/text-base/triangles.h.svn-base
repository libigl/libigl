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

#ifndef __EMBREE_ACCEL_TRIANGLES_H__
#define __EMBREE_ACCEL_TRIANGLES_H__

/*! include all triangles */
#include "virtual_object.h"

#include "triangle1i.h"
#include "triangle1v.h"
#include "triangle1.h"

#if defined (__SSE__)
#include "triangle4i.h"
#include "triangle4v.h"
#include "triangle4.h"
#endif

#if defined (__AVX__)
#include "triangle8.h"
#endif

/*! include single ray triangle intersectors */
#include "virtual_object_intersector1.h"

#if defined (__SSE__)
#include "intersector1_moeller.h"
#include "intersector1_pluecker.h"
#include "triangle1_intersector1_moeller.h"
#include "triangle4_intersector1_moeller.h"
#endif

#if defined (__AVX__)
#include "triangle8_intersector1_moeller.h"
#endif

#if defined(__MIC__)
#include "../geometry/triangle1_intersector1_moeller_mic.h"
#endif

/* include all 4 ray packet triangle intersectors */
#if defined(__SSE__)
#include "virtual_object_intersector4.h"
#include "intersector4_moeller.h"
#include "intersector4_pluecker.h"
#include "triangle1_intersector4_moeller.h"
#include "triangle4_intersector4_moeller.h"
#endif

/* include all 8 ray packet triangle intersectors */
#if defined(__AVX__)
#include "virtual_object_intersector8.h"
#include "intersector8_moeller.h"
#include "intersector8_pluecker.h"
#include "triangle1_intersector8_moeller.h"
#include "triangle4_intersector8_moeller.h"
#endif

/* include all 16 ray packet triangle intersectors */
#if defined(__MIC__)
#include "virtual_object_intersector16.h"
#include "../geometry/triangle1_intersector16_moeller_mic.h"
#endif

#endif
