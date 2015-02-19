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

#include "tessellation_cache.h"

namespace embree
{

  void TessellationCache::printStats()
  {
    CACHE_STATS(
                assert(cache_hits + cache_misses == cache_accesses);
                DBG_PRINT(cache_accesses);
                DBG_PRINT(cache_misses);
                DBG_PRINT(cache_hits);
                DBG_PRINT(cache_evictions);
                DBG_PRINT(100.0f * cache_hits / cache_accesses);
                DBG_PRINT(cache_clears);
                );
  }

  void TessellationCache::clearStats()
  {
  CACHE_STATS(
              TessellationCache::cache_accesses  = 0;
              TessellationCache::cache_hits      = 0;
              TessellationCache::cache_misses    = 0;
              TessellationCache::cache_clears    = 0;
              TessellationCache::cache_evictions = 0;          
              );
  }

  CACHE_STATS(
              AtomicCounter TessellationCache::cache_accesses  = 0;
              AtomicCounter TessellationCache::cache_hits      = 0;
              AtomicCounter TessellationCache::cache_misses    = 0;
              AtomicCounter TessellationCache::cache_clears    = 0;
              AtomicCounter TessellationCache::cache_evictions = 0;                
              );           

};

extern "C" void printTessCacheStats()
{
  embree::TessellationCache::printStats();
  embree::TessellationCache::clearStats();
}
