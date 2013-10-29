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

#ifndef __EMBREE_BVH4AOS_MIC_GLOBALS_H__
#define __EMBREE_BVH4AOS_MIC_GLOBALS_H__

#include <math.h>
#include "simd/mic.h"
#include "simd/mic_f.h"
#include "simd/mic_i.h"
#include "simd/mic_m.h"

#include "math/vec3.h"

#define MAX_MIC_THREADS 256
#define MAX_MIC_CORES MAX_MIC_THREADS/4
#define CACHELINE_SIZE 64

#define COMPILER_MEMORY_FENCE __asm__ __volatile__("" ::: "memory")

#define __ALIGN(x) __declspec(align(x))
#define _INLINE __forceinline
const float float_infinity = INFINITY; 
const float float_minus_infinity = -INFINITY; 
const float float_flt_max = FLT_MAX; 
const float float_minus_flt_max = -FLT_MAX; 


#define FATAL(x) { std::cout << "FATAL error in " << __PRETTY_FUNCTION__ << " : " << x << std::endl << std::flush; exit(0); }
#define DBG_PRINT(x) std::cout << #x << " = " << (x) << std::endl << std::flush


#endif
