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
  namespace isa
  {
    struct Mailbox 
    {
#if defined (__AVX__)
      avxi geomIDs;      
      avxi primIDs;
#else
      ssei geomIDs;   
      ssei primIDs;
#endif
      unsigned int index;
      
      __forceinline Mailbox() {
        geomIDs = -1;
        primIDs = -1;
        index = 0;
      };
      
      __forceinline bool hit(const unsigned geomID, const unsigned primID) const
      {
#if defined (__AVX__)
	const avxb geomMask = geomIDs == geomID;
	const avxb primMask = primIDs == primID;
#else
	const sseb geomMask = geomIDs == geomID;
	const sseb primMask = primIDs == primID;
#endif
	return any(geomMask & primMask);
      }
      
      __forceinline void add(const unsigned geomID, const unsigned primID)
      {
	*(unsigned int*)&geomIDs[index] = geomID;
	*(unsigned int*)&primIDs[index] = primID;
#if defined (__AVX__)
	index = (index + 1 ) % 8;
#else
	index = (index + 1 ) % 4;
#endif
      }
    };
  }
}
