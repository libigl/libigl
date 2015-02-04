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

#pragma once

#include "common/default.h"

namespace embree
{
  struct Mailbox 
  {
    mic_i geomIDs;      
    mic_i primIDs;
    unsigned int index;
    
    __forceinline Mailbox() {
      geomIDs = -1;
      primIDs = -1;
      index = 0;
    };
    
    __forceinline bool hit(const unsigned int geomID, const unsigned int primID) const
      {
	const mic_m m_geomMask = geomIDs == geomID;
	const mic_m m_primMask = primIDs == primID;
	return any(m_geomMask & m_primMask);
      }
    
    __forceinline void add(const unsigned int geomID, const unsigned int primID)
      {
	*(unsigned int*)&geomIDs[index] = geomID;
	*(unsigned int*)&primIDs[index] = primID;
	index = (index + 1 ) % 16;
      }
  };
}
