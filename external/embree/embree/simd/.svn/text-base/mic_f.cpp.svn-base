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

#include "mic.h"

namespace embree
{
  MIC_ALIGN float mic_f::float_idMod4[16] = {0.f,1.f,2.f,3.f,0.f,1.f,2.f,3.f,0.f,1.f,2.f,3.f,0.f,1.f,2.f,3.f};
  MIC_ALIGN float mic_f::float_idDiv4[16] = {0.f,0.f,0.f,0.f,1.f,1.f,1.f,1.f,2.f,2.f,2.f,2.f,3.f,3.f,3.f,3.f};
  MIC_ALIGN float mic_f::float_identity[16] = {0.f,1.f,2.f,3.f,4.f,5.f,6.f,7.f,8.f,9.f,10.f,11.f,12.f,13.f,14.f,15.f};
  MIC_ALIGN float mic_f::float_cancelLE[4] = { 1.0f,1.0f,1.0f,0.0f };
  
  MIC_ALIGN float mic_f::float_one_over[32] = {
    0.0f,
    1.0f /  1.0f,
    1.0f /  2.0f,
    1.0f /  3.0f,
    1.0f /  4.0f,
    1.0f /  5.0f,
    1.0f /  6.0f,
    1.0f /  8.0f,
    1.0f /  9.0f,
    1.0f / 10.0f,
    1.0f / 11.0f,
    1.0f / 12.0f,
    1.0f / 13.0f,
    1.0f / 14.0f,
    1.0f / 15.0f,
    1.0f / 16.0f,
    1.0f / 17.0f,
    1.0f / 18.0f,
    1.0f / 19.0f,
    1.0f / 21.0f,
    1.0f / 22.0f,
    1.0f / 23.0f,
    1.0f / 24.0f,
    1.0f / 25.0f,
    1.0f / 26.0f,
    1.0f / 27.0f,
    1.0f / 28.0f,
    1.0f / 29.0f,
    1.0f / 30.0f,
    1.0f / 31.0f
  };
}


