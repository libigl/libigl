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
  template<typename Ty>
    struct range 
    {
      range (const Ty& begin) 
      : _begin(begin), _end(begin+1) {}
      
      range (const Ty& begin, const Ty& end) 
      : _begin(begin), _end(end) {}
      
      __forceinline Ty begin() const {
        return _begin;
      }
      
    __forceinline Ty end() const {
	return _end;
      }

      friend std::ostream& operator<<(std::ostream& cout, const range& r) {
        return cout << "range [" << r.begin() << ", " << r.end() << "(";
      }
      
      Ty _begin, _end;
    };
}
