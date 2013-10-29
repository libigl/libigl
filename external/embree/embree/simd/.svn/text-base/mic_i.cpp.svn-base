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
  MIC_ALIGN int mic_i::int_identity[16] = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 };
  MIC_ALIGN int mic_i::int_idMod4[16] = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};
  MIC_ALIGN int mic_i::int_idDiv4[16] = {0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3};
  MIC_ALIGN int mic_i::int_pow4[16] = {1,4,16,64,256,1024,4096,16384,65536,262144,1048576,4194304,16777216,67108864,268435456,1073741824};
  MIC_ALIGN int mic_i::int_zlc4[4] = {0xffffffff,0xffffffff,0xffffffff,0};
  MIC_ALIGN int mic_i::int_addtriID4[16] = {0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3};
  MIC_ALIGN unsigned int mic_i::uint_shift1[32] = {
    ((unsigned int)1 << 0),
    ((unsigned int)1 << 1),
    ((unsigned int)1 << 2),
    ((unsigned int)1 << 3),
    ((unsigned int)1 << 4),
    ((unsigned int)1 << 5),
    ((unsigned int)1 << 6),
    ((unsigned int)1 << 7),
    ((unsigned int)1 << 8),
    ((unsigned int)1 << 9),
    ((unsigned int)1 << 10),
    ((unsigned int)1 << 11),
    ((unsigned int)1 << 12),
    ((unsigned int)1 << 13),
    ((unsigned int)1 << 14),
    ((unsigned int)1 << 15),
    ((unsigned int)1 << 16),
    ((unsigned int)1 << 17),
    ((unsigned int)1 << 18),
    ((unsigned int)1 << 19),
    ((unsigned int)1 << 20),
    ((unsigned int)1 << 21),
    ((unsigned int)1 << 22),
    ((unsigned int)1 << 23),
    ((unsigned int)1 << 24),
    ((unsigned int)1 << 25),
    ((unsigned int)1 << 26),
    ((unsigned int)1 << 27),
    ((unsigned int)1 << 28),
    ((unsigned int)1 << 29),
    ((unsigned int)1 << 30),
    ((unsigned int)1 << 31)
  };
  
  MIC_ALIGN int mic_i::int_aos2soa[16] = {
    0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15
  };
  
  MIC_ALIGN int mic_i::int_reverse_identity[16] = { 15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0 };
}
