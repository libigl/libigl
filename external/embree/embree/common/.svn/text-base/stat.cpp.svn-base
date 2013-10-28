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

#include "stat.h"

namespace embree
{
  Stat Stat::instance; 
  
  Stat::Stat () {
  }

  Stat::~Stat () 
  {
#ifdef __USE_STAT_COUNTERS__
    Stat::print(std::cout);
#endif
  }

  void Stat::print(std::ostream& cout)
  {
    Counters& cntrs = instance.cntrs;

    /* print absolute numbers */
    cout << "--------- ABSOLUTE ---------" << std::endl;
    cout << "  #normal_travs   = " << float(cntrs.code.normal.travs            )*1E-6 << "M" << std::endl;
    cout << "    #nodes        = " << float(cntrs.code.normal.trav_nodes       )*1E-6 << "M" << std::endl;
    cout << "    #leaves       = " << float(cntrs.code.normal.trav_leaves      )*1E-6 << "M" << std::endl;
    cout << "    #tris         = " << float(cntrs.code.normal.trav_tris        )*1E-6 << "M" << std::endl;
    if (cntrs.code.shadow.travs) {
      cout << "  #shadow_travs = " << float(cntrs.code.shadow.travs      )*1E-6 << "M" << std::endl;
      cout << "    #nodes      = " << float(cntrs.code.shadow.trav_nodes )*1E-6 << "M" << std::endl;
      cout << "    #leaves     = " << float(cntrs.code.shadow.trav_leaves)*1E-6 << "M" << std::endl;
      cout << "    #tris       = " << float(cntrs.code.shadow.trav_tris  )*1E-6 << "M" << std::endl;
    }
    cout << std::endl;

    /* print per traversal numbers */
    cout << "--------- PER TRAVERSAL ---------" << std::endl;
    float active_normal_travs       = float(cntrs.active.normal.travs      )/float(cntrs.all.normal.travs      );
    float active_normal_trav_nodes  = float(cntrs.active.normal.trav_nodes )/float(cntrs.all.normal.trav_nodes );
    float active_normal_trav_leaves = float(cntrs.active.normal.trav_leaves)/float(cntrs.all.normal.trav_leaves);
    float active_normal_trav_tris   = float(cntrs.active.normal.trav_tris  )/float(cntrs.all.normal.trav_tris  );
    cout << "  #normal_travs   = " << float(cntrs.all.normal.travs      )/float(cntrs.all.normal.travs) << ", " << 100.0f*active_normal_travs       << "% active" << std::endl;
    cout << "    #nodes        = " << float(cntrs.all.normal.trav_nodes )/float(cntrs.all.normal.travs) << ", " << 100.0f*active_normal_trav_nodes  << "% active" << std::endl;
    cout << "    #leaves       = " << float(cntrs.all.normal.trav_leaves)/float(cntrs.all.normal.travs) << ", " << 100.0f*active_normal_trav_leaves << "% active" << std::endl;
    cout << "    #tris         = " << float(cntrs.all.normal.trav_tris  )/float(cntrs.all.normal.travs) << ", " << 100.0f*active_normal_trav_tris   << "% active" << std::endl;
    if (cntrs.all.shadow.travs) {
      float active_shadow_travs       = float(cntrs.active.shadow.travs      )/float(cntrs.all.shadow.travs      );
      float active_shadow_trav_nodes  = float(cntrs.active.shadow.trav_nodes )/float(cntrs.all.shadow.trav_nodes );
      float active_shadow_trav_leaves = float(cntrs.active.shadow.trav_leaves)/float(cntrs.all.shadow.trav_leaves);
      float active_shadow_trav_tris   = float(cntrs.active.shadow.trav_tris  )/float(cntrs.all.shadow.trav_tris  );
      cout << "  #shadow_travs = " << float(cntrs.all.shadow.travs      )/float(cntrs.all.shadow.travs) << ", " << 100.0f*active_shadow_travs       << "% active" << std::endl;
      cout << "    #nodes      = " << float(cntrs.all.shadow.trav_nodes )/float(cntrs.all.shadow.travs) << ", " << 100.0f*active_shadow_trav_nodes  << "% active" << std::endl;
      cout << "    #leaves     = " << float(cntrs.all.shadow.trav_leaves)/float(cntrs.all.shadow.travs) << ", " << 100.0f*active_shadow_trav_leaves << "% active" << std::endl;
      cout << "    #tris       = " << float(cntrs.all.shadow.trav_tris  )/float(cntrs.all.shadow.travs) << ", " << 100.0f*active_shadow_trav_tris   << "% active" << std::endl;
    }
    cout << std::endl;
  }
}
