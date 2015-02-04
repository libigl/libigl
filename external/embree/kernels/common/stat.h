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

#include "default.h"

/* Makros to gather statistics */
#ifdef RTCORE_STAT_COUNTERS
#define STAT(x) x
#define STAT3(s,x,y,z) \
  STAT(Stat::get().code  .s+=x);               \
  STAT(Stat::get().active.s+=y);               \
  STAT(Stat::get().all   .s+=z);
#else
#define STAT(x)
#define STAT3(s,x,y,z)
#endif

namespace embree
{
  /*! Gathers ray tracing statistics. */
  class Stat
  { 
  public:
    
    Stat ();
    ~Stat ();

    class Counters 
    {
    public:
      Counters () { 
        clear(); 
      }
      
      void clear() { 
        memset(this,0,sizeof(Counters)); 
      }

    public:

	/* per packet and per ray stastics */
	struct {
	  /* normal and shadow ray statistics */
	  struct {
	    AtomicCounter travs;
	    AtomicCounter trav_nodes;
	    AtomicCounter trav_leaves;
	    AtomicCounter trav_prims;
	    AtomicCounter trav_prim_hits;
#if defined(__MIC__)
	    AtomicCounter trav_hit_boxes[16+1];
	    AtomicCounter trav_stack_nodes;

#endif

	  } normal, shadow;
	} all, active, code;

    };

  public:

    static __forceinline Counters& get() {
      return instance.cntrs;
    }
    
    static void clear() {
      instance.cntrs.clear();
    }
    
    static void print(std::ostream& cout);

  private: 
    Counters cntrs;
  private:
    static Stat instance;
  };
}
