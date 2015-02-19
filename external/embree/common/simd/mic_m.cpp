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

#include "mic.h"

namespace embree
{

  __aligned(64) unsigned int mic_m::shift1[32] = {
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

};
