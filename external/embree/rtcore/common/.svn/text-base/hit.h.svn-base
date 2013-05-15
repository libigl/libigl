// ======================================================================== //
// Copyright 2009-2012 Intel Corporation                                    //
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

#ifndef __EMBREE_HIT_H__
#define __EMBREE_HIT_H__

#include "../common/default.h"

namespace embree
{
  /*! Hit information. Stores all information about a hit. */
  struct Hit
  {
    /*! Default constructor creates invalid hit. */
    Hit () : id0(-1), id1(-1), t(inf) {};

    /*! Tests if we hit something. */
    __forceinline operator bool() const { return id0 != -1; }

  public:
    int id0;           //!< 1st primitive ID
    int id1;           //!< 2nd primitive ID
    float u;           //!< Barycentric u coordinate of hit
    float v;           //!< Barycentric v coordinate of hit
    float t;           //!< Distance of hit
  };

  /*! Outputs hit to to stream. */
  inline std::ostream& operator<<(std::ostream& cout, const Hit& hit) {
    return cout << "{ id0 = " << hit.id0 << ", id1 = " << hit.id1 <<  ", "
                << "u = " << hit.u <<  ", v = " << hit.v <<  ", t = " << hit.t << " }";
  }
}

#endif
