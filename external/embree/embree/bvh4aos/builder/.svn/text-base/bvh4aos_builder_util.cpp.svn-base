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


#include "bvh4aos_globals.h"
#include "bvh4aos_builder_common.h"
#include "bvh4aos_builder_util.h"

#include "simd/mic_i.h"

namespace embree
{

  MIC_ALIGN float Centroid_Scene_AABB::initCentroidScene[16] = { 
    float_flt_max       ,float_flt_max      ,float_flt_max       ,0.0f,
    float_minus_flt_max ,float_minus_flt_max,float_minus_flt_max ,0.0f,
    float_flt_max       ,float_flt_max      ,float_flt_max       ,0.0f,
    float_minus_flt_max ,float_minus_flt_max,float_minus_flt_max ,0.0f
  };

  void QuadTreeBarrier::sync(const unsigned int threadID, const unsigned int MAX_THREADS_SYNC)
  {
    if (unlikely(MAX_THREADS_SYNC == 1)) return;

    const unsigned int MAX_CORES_SYNC = MAX_THREADS_SYNC >> 2;
    const unsigned int coreID = threadID >> 2;
    const unsigned int MODE = data[coreID].mode;

    // Drain store buffer for NGO stores
    atomic_add((atomic_t*)&data[coreID].data[0],0);

    data[coreID].prefetchL1Ex(); 

    if (threadID == 0)
      {		

	data[0].setThreadStateToDone(MODE,threadID);

	// == wait for core 0 ==
	data[0].waitForAllThreadsOnCore(MODE);

	// == wait for (possible) two children cores
	const unsigned int nextCoreID0 = 1;
	const unsigned int nextCoreID1 = 2;

	data[nextCoreID0].prefetchL1(); 
	data[nextCoreID1].prefetchL1(); 

	if (nextCoreID0 < MAX_CORES_SYNC) data[nextCoreID0].waitForAllThreadsOnCore(MODE);
	if (nextCoreID1 < MAX_CORES_SYNC) data[nextCoreID1].waitForAllThreadsOnCore(MODE);

	data[nextCoreID0].prefetchL1Ex(); 
	data[nextCoreID1].prefetchL1Ex(); 

	// == run signal to core 0 ==
	data[0].switchModeAndSendRunSignal(MODE);

	// == propagate run signal to core 1,2 ==
	if (nextCoreID0 < MAX_CORES_SYNC) data[nextCoreID0].switchModeAndSendRunSignal(MODE);
	if (nextCoreID1 < MAX_CORES_SYNC) data[nextCoreID1].switchModeAndSendRunSignal(MODE);
      }			
    else
      {				

	const unsigned int nextCoreID0 = 2*coreID + 1;
	const unsigned int nextCoreID1 = 2*coreID + 2;
	data[nextCoreID0].prefetchL1(); 
	data[nextCoreID1].prefetchL1(); 

	if (threadID % 4 == 0)
	  {	

	    if (nextCoreID0 < MAX_CORES_SYNC) data[nextCoreID0].waitForAllThreadsOnCore(MODE);
	    if (nextCoreID1 < MAX_CORES_SYNC) data[nextCoreID1].waitForAllThreadsOnCore(MODE);
	  }

	data[coreID].setThreadStateToDone(MODE,threadID % 4);
	data[coreID].waitForThreadReceivesRunSignal(MODE,threadID % 4);	

	// == propagte run signal to the two children ==

	if (threadID % 4 == 0)
	  {
	    data[nextCoreID0].prefetchL1Ex(); 
	    data[nextCoreID1].prefetchL1Ex(); 

	    if (nextCoreID0 < MAX_CORES_SYNC) data[nextCoreID0].switchModeAndSendRunSignal(MODE);
	    if (nextCoreID1 < MAX_CORES_SYNC) data[nextCoreID1].switchModeAndSendRunSignal(MODE);
	  }
      }

  }


  void QuadTreeBarrier::syncWithReduction(const unsigned int threadID, 
					  const unsigned int MAX_THREADS_SYNC,
					  void (* reductionFct)(const unsigned int currentThreadID,
								const unsigned int childThreadID))
  {
    if (unlikely(MAX_THREADS_SYNC == 1)) return;

    const unsigned int MAX_CORES_SYNC = MAX_THREADS_SYNC >> 2;
    const unsigned int coreID = threadID >> 2;
    const unsigned int MODE = data[coreID].mode;

    // Drain store buffer for NGO stores
    atomic_add((atomic_t*)&data[coreID].data[0],0);

    data[coreID].prefetchL1Ex(); 

    if (threadID == 0)
      {		

	data[0].setThreadStateToDone(MODE,threadID);

	// == wait for core 0 ==
	data[0].waitForAllThreadsOnCore(MODE);

	(*reductionFct)(threadID,threadID+1);
	(*reductionFct)(threadID,threadID+2);
	(*reductionFct)(threadID,threadID+3);

	// == wait for (possible) two children cores
	const unsigned int nextCoreID0 = 1;
	const unsigned int nextCoreID1 = 2;
	const unsigned int nextThreadID0 = nextCoreID0 * 4;
	const unsigned int nextThreadID1 = nextCoreID1 * 4;

	data[nextCoreID0].prefetchL1(); 
	data[nextCoreID1].prefetchL1(); 

	if (nextCoreID0 < MAX_CORES_SYNC) 
	  {
	    data[nextCoreID0].waitForAllThreadsOnCore(MODE);
	    (*reductionFct)(threadID,nextThreadID0);
	  }

	if (nextCoreID1 < MAX_CORES_SYNC) 
	  {
	    data[nextCoreID1].waitForAllThreadsOnCore(MODE);
	    (*reductionFct)(threadID,nextThreadID1);
	  }

	data[nextCoreID0].prefetchL1Ex(); 
	data[nextCoreID1].prefetchL1Ex(); 

	// == run signal to core 0 ==
	data[0].switchModeAndSendRunSignal(MODE);

	// == propagate run signal to core 1,2 ==
	if (nextCoreID0 < MAX_CORES_SYNC) data[nextCoreID0].switchModeAndSendRunSignal(MODE);
	if (nextCoreID1 < MAX_CORES_SYNC) data[nextCoreID1].switchModeAndSendRunSignal(MODE);
      }			
    else
      {				

	const unsigned int nextCoreID0 = 2*coreID + 1;
	const unsigned int nextCoreID1 = 2*coreID + 2;
	const unsigned int nextThreadID0 = nextCoreID0 * 4;
	const unsigned int nextThreadID1 = nextCoreID1 * 4;

	data[nextCoreID0].prefetchL1(); 
	data[nextCoreID1].prefetchL1(); 

	if (threadID % 4 == 0)
	  {	
	    data[coreID].waitForAllOtherThreadsOnCore(MODE,threadID);
	    (*reductionFct)(threadID,threadID+1);
	    (*reductionFct)(threadID,threadID+2);
	    (*reductionFct)(threadID,threadID+3);

	    if (nextCoreID0 < MAX_CORES_SYNC) 
	      {
		data[nextCoreID0].waitForAllThreadsOnCore(MODE);
		(*reductionFct)(threadID,nextThreadID0);
	      }
	  

	    if (nextCoreID1 < MAX_CORES_SYNC) 
	      {
		data[nextCoreID1].waitForAllThreadsOnCore(MODE);
		(*reductionFct)(threadID,nextThreadID1);
	      }
	  }

	data[coreID].setThreadStateToDone(MODE,threadID % 4);

	data[coreID].waitForThreadReceivesRunSignal(MODE,threadID % 4);	

	// == propagte run signal to the two children ==

	if (threadID % 4 == 0)
	  {
	    data[nextCoreID0].prefetchL1Ex(); 
	    data[nextCoreID1].prefetchL1Ex(); 

	    if (nextCoreID0 < MAX_CORES_SYNC) data[nextCoreID0].switchModeAndSendRunSignal(MODE);
	    if (nextCoreID1 < MAX_CORES_SYNC) data[nextCoreID1].switchModeAndSendRunSignal(MODE);
	  }
      }  
  }


  void SimpleBarrier::sync(const unsigned int ID, const unsigned int MAX)
  {
    if (mode == 0)
      {			
	if (ID == 0)
	  {	
	    store16i((char*)&count1[  0],mic_i::zero());
	    store16i((char*)&count1[ 64],mic_i::zero());
	    store16i((char*)&count1[128],mic_i::zero());
	    store16i((char*)&count1[192],mic_i::zero());
	    for (unsigned int i=1; i<MAX; i++)
	      {
		unsigned int delay = MIN_DELAY_CYCLES;
		while(count0[i] == 0) 
		  {
		    delayThreadExpFalloff(delay);
		  }
	      }

	    mode  = 1;
	    flag1 = 0;
	    flag0 = 1;
	  }			
	else
	  {					
	    unsigned int delay = MIN_DELAY_CYCLES;
	    count0[ID] = 1;
	    while(unlikely(flag0 == 0)) delayThreadExpFalloff(delay);
	  }		
      }					
    else						
      {
	if (ID == 0)
	  {	
	    store16i((char*)&count0[  0],mic_i::zero());
	    store16i((char*)&count0[ 64],mic_i::zero());
	    store16i((char*)&count0[128],mic_i::zero());
	    store16i((char*)&count0[192],mic_i::zero());
	    for (unsigned int i=1; i<MAX; i++)
	      {
		unsigned int delay = MIN_DELAY_CYCLES;
		while(count1[i] == 0) delayThreadExpFalloff(delay);
	      }

	    mode  = 0;
	    flag0 = 0;
	    flag1 = 1;
	  }			
	else
	  {					
	    unsigned int delay = MIN_DELAY_CYCLES;
	    count1[ID] = 1;
	    while(unlikely(flag1 == 0)) 
	      {
		delayThreadExpFalloff(delay);
	      }
	  }		
      }					
  }


  static double getFrequencyInMHz()
  {
    struct timeval tvstart, tvstop;
    unsigned long long int cycles[2];
    
    gettimeofday(&tvstart, NULL);
    cycles[0] = rdtsc();
    gettimeofday(&tvstart, NULL);
    usleep(250000);
    gettimeofday(&tvstop, NULL);
    cycles[1] = rdtsc();
    gettimeofday(&tvstop, NULL);
  
    const unsigned long microseconds = ((tvstop.tv_sec-tvstart.tv_sec)*1000000) + (tvstop.tv_usec-tvstart.tv_usec);
    unsigned long mhz = (unsigned long) (cycles[1]-cycles[0]) / microseconds;

    //PRINT(mhz);
    return (double)mhz;
  }

  float Timer::freqInMHz = getFrequencyInMHz();

};
