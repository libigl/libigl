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

#include "acceln.h"
#include "embree2/rtcore_ray.h"

namespace embree
{
  AccelN::AccelN () : N(0), M(0) { 
    memset(accels,0,sizeof(accels));
  }

  AccelN::~AccelN() 
  {
    for (size_t i=0; i<N; i++)
      delete accels[i];
  }

  void AccelN::add(Accel* accel) 
  {
    assert(N<16);
    assert(accel);
    if (N<16) accels[N++] = accel;
  }
  
  void AccelN::intersect (void* ptr, RTCRay& ray) 
  {
    AccelN* This = (AccelN*)ptr;
    for (size_t i=0; i<This->M; i++)
      This->validAccels[i]->intersect(ray);
  }

  void AccelN::intersect4 (const void* valid, void* ptr, RTCRay4& ray) 
  {
    AccelN* This = (AccelN*)ptr;
    for (size_t i=0; i<This->M; i++)
      This->validAccels[i]->intersect4(valid,ray);
  }

  void AccelN::intersect8 (const void* valid, void* ptr, RTCRay8& ray) 
  {
    AccelN* This = (AccelN*)ptr;
    for (size_t i=0; i<This->M; i++)
      This->validAccels[i]->intersect8(valid,ray);
  }

  void AccelN::intersect16 (const void* valid, void* ptr, RTCRay16& ray) 
  {
    AccelN* This = (AccelN*)ptr;
    for (size_t i=0; i<This->M; i++)
      This->validAccels[i]->intersect16(valid,ray);
  }

  void AccelN::occluded (void* ptr, RTCRay& ray) 
  {
    AccelN* This = (AccelN*)ptr;
    for (size_t i=0; i<This->M; i++) {
      This->validAccels[i]->occluded(ray); 
      if (ray.geomID == 0) break;
    }
  }

  void AccelN::occluded4 (const void* valid, void* ptr, RTCRay4& ray) 
  {
    AccelN* This = (AccelN*)ptr;
    for (size_t i=0; i<This->M; i++) {
      This->validAccels[i]->occluded4(valid,ray);
#if defined(__SSE2__)
      sseb valid0 = ((sseb*)valid)[0];
      sseb hit0   = ((ssei*)ray.geomID)[0] == ssei(0);
      if (all(valid0,hit0)) break;
#endif
    }
  }

  void AccelN::occluded8 (const void* valid, void* ptr, RTCRay8& ray) 
  {
    AccelN* This = (AccelN*)ptr;
    for (size_t i=0; i<This->M; i++) {
      This->validAccels[i]->occluded8(valid,ray);
#if defined(__SSE2__)
      sseb valid0 = ((sseb*)valid)[0];
      sseb hit0   = ((ssei*)ray.geomID)[0] == ssei(0);
      sseb valid1 = ((sseb*)valid)[1];
      sseb hit1   = ((ssei*)ray.geomID)[1] == ssei(0);
      if (all(valid0,hit0) && all(valid1,hit1)) break;
#endif
    }
  }

  void AccelN::occluded16 (const void* valid, void* ptr, RTCRay16& ray) 
  {
    AccelN* This = (AccelN*)ptr;
    for (size_t i=0; i<This->M; i++) {
      This->validAccels[i]->occluded16(valid,ray);
#if defined(__MIC__)
      mic_m valid0 = ((mic_m*)valid)[0];
      mic_m hit0   = ((mic_i*)ray.geomID)[0] == mic_i(0);
      if (all(valid0,hit0)) break;
#endif
    }
  }

  void AccelN::print(size_t ident)
  {
    for (size_t i=0; i<M; i++)
    {
      for (size_t j=0; j<ident; j++) std::cout << " "; 
      std::cout << "accels[" << i << "]" << std::endl;
      validAccels[i]->intersectors.print(ident+2);
    }
  }

  void AccelN::immutable()
  {
    for (size_t i=0; i<N; i++)
      accels[i]->immutable();
  }

  void AccelN::build (size_t threadIndex, size_t threadCount) 
  {
    /* build all acceleration structures */
    M = 0;

    for (size_t i=0; i<N; i++) 
    {
      accels[i]->build(threadIndex,threadCount);

      if (accels[i]->bounds.empty()) continue;
      validAccels[M++] = accels[i];
    }

    if (M == 1) {
      intersectors = validAccels[0]->intersectors;
    }
    else 
    {
      intersectors.ptr = this;
      intersectors.intersector1 = Intersector1(&intersect,&occluded,"AccelN::intersector1");
      intersectors.intersector4 = Intersector4(&intersect4,&occluded4,"AccelN::intersector4");
      intersectors.intersector8 = Intersector8(&intersect8,&occluded8,"AccelN::intersector8");
      intersectors.intersector16= Intersector16(&intersect16,&occluded16,"AccelN::intersector16");
    }
    
    /*! calculate bounds */
    bounds = empty;
    for (size_t i=0; i<M; i++) 
      bounds.extend(validAccels[i]->bounds);
  }

  void AccelN::select(bool filter4, bool filter8, bool filter16)
  {
    for (size_t i=0; i<N; i++) 
      accels[i]->intersectors.select(filter4,filter8,filter16);
  }
}
