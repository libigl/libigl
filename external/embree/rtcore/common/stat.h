/* ************************************************************************* *\
 INTEL CORPORATION PROPRIETARY INFORMATION
 This software is supplied under the terms of a license agreement or 
 nondisclosure agreement with Intel Corporation and may not be copied 
 or disclosed except in accordance with the terms of that agreement. 
 Copyright (C) 2009 Intel Corporation. All Rights Reserved.
\* ************************************************************************* */

#ifndef __EMBREE_STAT_H__
#define __EMBREE_STAT_H__

#include "default.h"

/* Makros to gather statistics */
#ifdef __USE_STAT_COUNTERS__
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

      union {
        
        struct {
          
          /* per packet and per ray stastics */
          struct {
            /* normal and shadow ray statistics */
            struct {
              size_t travs;
              size_t trav_nodes;
              size_t trav_leaves;
              size_t trav_tris;
            } normal, shadow;
          } all, active, code;
        };

        size_t data[32];
      };
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

#endif
