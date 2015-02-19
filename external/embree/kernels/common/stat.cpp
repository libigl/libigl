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

#include "stat.h"

namespace embree
{
  Stat Stat::instance; 
  
  Stat::Stat () {
  }

  Stat::~Stat () 
  {
#ifdef RTCORE_STAT_COUNTERS
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
    cout << "    #prims        = " << float(cntrs.code.normal.trav_prims       )*1E-6 << "M" << std::endl;
    cout << "    #prim_hits    = " << float(cntrs.code.normal.trav_prim_hits   )*1E-6 << "M" << std::endl;

#if defined(__MIC__)
    cout << "    #stack nodes  = " << float(cntrs.code.normal.trav_stack_nodes )*1E-6 << "M" << std::endl;

    size_t normal_box_hits = 0;
    for (size_t i=0;i<=16;i++) normal_box_hits += cntrs.code.normal.trav_hit_boxes[i];
    cout << "    #hit_boxes    = ";
    for (size_t i=0;i<=16;i++) cout << "[" << i << "] " << 100.0f * cntrs.code.normal.trav_hit_boxes[i] / normal_box_hits << " ";
    cout << std::endl;
#endif
    if (cntrs.code.shadow.travs) {
      cout << "  #shadow_travs = " << float(cntrs.code.shadow.travs         )*1E-6 << "M" << std::endl;
      cout << "    #nodes      = " << float(cntrs.code.shadow.trav_nodes    )*1E-6 << "M" << std::endl;
      cout << "    #leaves     = " << float(cntrs.code.shadow.trav_leaves   )*1E-6 << "M" << std::endl;
      cout << "    #prims      = " << float(cntrs.code.shadow.trav_prims    )*1E-6 << "M" << std::endl;
      cout << "    #prim_hits  = " << float(cntrs.code.shadow.trav_prim_hits)*1E-6 << "M" << std::endl;

#if defined(__MIC__)
      cout << "    #stack nodes = " << float(cntrs.code.shadow.trav_stack_nodes )*1E-6 << "M" << std::endl;

      size_t shadow_box_hits = 0;
      for (size_t i=0;i<=16;i++) shadow_box_hits += cntrs.code.shadow.trav_hit_boxes[i];
      cout << "    #hit_boxes    = ";
      for (size_t i=0;i<=16;i++) cout << "[" << i << "] " << 100.0f * cntrs.code.shadow.trav_hit_boxes[i] / shadow_box_hits << " ";
      cout << std::endl;
#endif

    }
    cout << std::endl;

    /* print per traversal numbers */
    cout << "--------- PER TRAVERSAL ---------" << std::endl;
    float active_normal_travs       = float(cntrs.active.normal.travs      )/float(cntrs.all.normal.travs      );
    float active_normal_trav_nodes  = float(cntrs.active.normal.trav_nodes )/float(cntrs.all.normal.trav_nodes );
    float active_normal_trav_leaves = float(cntrs.active.normal.trav_leaves)/float(cntrs.all.normal.trav_leaves);
    float active_normal_trav_prims   = float(cntrs.active.normal.trav_prims  )/float(cntrs.all.normal.trav_prims  );
    float active_normal_trav_prim_hits = float(cntrs.active.normal.trav_prim_hits  )/float(cntrs.all.normal.trav_prim_hits  );

    cout << "  #normal_travs   = " << float(cntrs.all.normal.travs      )/float(cntrs.all.normal.travs) << ", " << 100.0f*active_normal_travs       << "% active" << std::endl;
    cout << "    #nodes        = " << float(cntrs.all.normal.trav_nodes )/float(cntrs.all.normal.travs) << ", " << 100.0f*active_normal_trav_nodes  << "% active" << std::endl;
    cout << "    #leaves       = " << float(cntrs.all.normal.trav_leaves)/float(cntrs.all.normal.travs) << ", " << 100.0f*active_normal_trav_leaves << "% active" << std::endl;
    cout << "    #prims        = " << float(cntrs.all.normal.trav_prims  )/float(cntrs.all.normal.travs) << ", " << 100.0f*active_normal_trav_prims   << "% active" << std::endl;
    cout << "    #prim_hits    = " << float(cntrs.all.normal.trav_prim_hits  )/float(cntrs.all.normal.travs) << ", " << 100.0f*active_normal_trav_prim_hits   << "% active" << std::endl;

    if (cntrs.all.shadow.travs) {
      float active_shadow_travs       = float(cntrs.active.shadow.travs      )/float(cntrs.all.shadow.travs      );
      float active_shadow_trav_nodes  = float(cntrs.active.shadow.trav_nodes )/float(cntrs.all.shadow.trav_nodes );
      float active_shadow_trav_leaves = float(cntrs.active.shadow.trav_leaves)/float(cntrs.all.shadow.trav_leaves);
      float active_shadow_trav_prims   = float(cntrs.active.shadow.trav_prims  )/float(cntrs.all.shadow.trav_prims  );
      float active_shadow_trav_prim_hits = float(cntrs.active.shadow.trav_prim_hits  )/float(cntrs.all.shadow.trav_prim_hits  );

      cout << "  #shadow_travs = " << float(cntrs.all.shadow.travs      )/float(cntrs.all.shadow.travs) << ", " << 100.0f*active_shadow_travs       << "% active" << std::endl;
      cout << "    #nodes      = " << float(cntrs.all.shadow.trav_nodes )/float(cntrs.all.shadow.travs) << ", " << 100.0f*active_shadow_trav_nodes  << "% active" << std::endl;
      cout << "    #leaves     = " << float(cntrs.all.shadow.trav_leaves)/float(cntrs.all.shadow.travs) << ", " << 100.0f*active_shadow_trav_leaves << "% active" << std::endl;
      cout << "    #prims      = " << float(cntrs.all.shadow.trav_prims  )/float(cntrs.all.shadow.travs) << ", " << 100.0f*active_shadow_trav_prims   << "% active" << std::endl;
      cout << "    #prim_hits  = " << float(cntrs.all.shadow.trav_prim_hits  )/float(cntrs.all.shadow.travs) << ", " << 100.0f*active_shadow_trav_prim_hits   << "% active" << std::endl;

    }
    cout << std::endl;
  }
}
