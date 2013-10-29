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

#ifndef __EMBREE_BVH4AOS_BUILDER_UTIL_H__
#define __EMBREE_BVH4AOS_BUILDER_UTIL_H__

#include "bvh4aos_box.h"

#include <unistd.h>
#include <sys/time.h>

namespace embree
{
  typedef volatile int atomic_t;

  static __inline int atomic_add(atomic_t *p, int v)
  {
    __asm __volatile("	lock;	"
		     "	xaddl	%0, %1 ;	"
		     : "+r" (v),			/* 0 (result) */
		     "=m" (*p)			/* 1 */
		     : "m" (*p));			/* 2 */
    return (v);
  }

  static __inline long atomic_add64(volatile long *p, long v)
  {
    __asm __volatile("	lock;	"
		     "	xadd	%0, %1 ;	"
		     : "+r" (v),			/* 0 (result) */
		     "=m" (*p)			/* 1 */
		     : "m" (*p));			/* 2 */
    return (v);
  }


  /* atomic compare and exchange */

  _INLINE int atomic_cmpxchg(volatile int *dst, int exp, int src)
  {
    unsigned int res;


    __asm __volatile(
		     "	lock; 		        "
		     "	cmpxchgl %2,%1 ;	"
		     : "=a" (res),			/* 0 */
		     "=m" (*dst)			/* 1 */
		     : "r" (src),			/* 2 */
		     "a" (exp),			/* 3 */
		     "m" (*dst)			/* 4 */
		     : "memory");

    return (res);
  }
 
  /* atomic test and set bit */
  _INLINE unsigned int test_and_set_bit(int nr, volatile unsigned int* addr)
  {
    unsigned int oldbit;
    asm volatile( "lock; "
		  "btsl %2,%1\n\tsbbl %0,%0"
		  :"=r" (oldbit), "+m" (*addr)
		  :"dIr" (nr) : "memory");
    return oldbit;
  }


  /* atomic clear of a bit */
  _INLINE void clear_bit(int nr, volatile unsigned int * addr)
  {
    asm volatile( "lock; "
		  "btrl %1,%0"
		  :"+m" (*addr)
		  :"dIr" (nr));
  }

  _INLINE void nop() 
  {
    asm ("nop");
  }


  _INLINE void atomic_min(volatile float *__restrict__ ptr, const float b)
  {
    const int int_b = *(int*)&b;
    while (1)
      {
	float a = *ptr;
	if (a <= b) { break; }
	const int int_a = *(int*)&a;
	const int result = atomic_cmpxchg((int*)ptr,int_a,int_b);
	if (result == int_a) break;
      }
  }

  _INLINE void atomic_max(volatile float *__restrict__ ptr, const float b)
  {
    const int int_b = *(int*)&b;
    while (1)
      {
	float a = *ptr;
	if (a >= b) { break; }
	const int int_a = *(int*)&a;
	const int result = atomic_cmpxchg((int*)ptr,int_a,int_b);
	if (result == int_a) break;
      }
  }

  _INLINE void atomic_min_i32(volatile int *__restrict__ ptr, const int b)
  {
    while (1)
      {
	int a = *ptr;
	if (a <= b) { break; }
	const int int_a = *(int*)&a;
	const int result = atomic_cmpxchg((int*)ptr,int_a,b);
	if (result == int_a) break;
      }
  }

  _INLINE void atomic_max_i32(volatile int *__restrict__ ptr, const int b)
  {
    while (1)
      {
	int a = *ptr;
	if (a >= b) { break; }
	const int int_a = *(int*)&a;
	const int result = atomic_cmpxchg((int*)ptr,int_a,b);
	if (result == int_a) break;
      }
  }

  _INLINE void atomic_min_ui32(volatile unsigned int *__restrict__ ptr, const unsigned int b)
  {
    while (1)
      {
	unsigned int a = *ptr;
	if (a <= b) { break; }
	const unsigned int int_a = *(unsigned int*)&a;
	const unsigned int result = atomic_cmpxchg((int*)ptr,int_a,b);
	if (result == int_a) break;
      }
  }

  _INLINE void atomic_max_ui32(volatile unsigned int *__restrict__ ptr, const unsigned int b)
  {
    while (1)
      {
	unsigned int a = *ptr;
	if (a >= b) { break; }
	const unsigned int int_a = *(unsigned int*)&a;
	const unsigned int result = atomic_cmpxchg((int*)ptr,int_a,b);
	if (result == int_a) break;
      }
  }



  /* class should have the right alignment to prevent cache trashing */
  __ALIGN(64)
    class AtomicCounter
    {
    private:
      atomic_t m_counter;

      char dummy[64-sizeof(atomic_t)]; // one counter per cache line
   
    public:

      AtomicCounter() {
	reset();
      }

      AtomicCounter(const unsigned int v) {
	m_counter = v;
      }

      _INLINE void prefetchEx() {
	prefetch<PFHINT_L1EX>((void*)&m_counter);
      }

      _INLINE void reset(unsigned int i = 0) {
	*(volatile unsigned int*)&m_counter = i;

      }

      _INLINE unsigned int inc() {
	return atomic_add(&m_counter,1);
      }

      _INLINE unsigned int dec() {
	return atomic_add(&m_counter,-1);
      }

      _INLINE unsigned int add(const int i) {
	return atomic_add(&m_counter, i);
      }

      _INLINE unsigned int val() {
	return m_counter;
      };

      _INLINE AtomicCounter operator++(int)
      {
	return inc();
      }

    };

  __ALIGN(64)
    class AtomicCounter64
    {
    private:
      volatile long m_counter;

      char dummy[64-sizeof(long)]; // one counter per cache line
   
    public:

      AtomicCounter64() {
	reset();
      }

      AtomicCounter64(const long v) {
	m_counter = v;
      }

      _INLINE operator long () {
	return m_counter;
      };

      _INLINE void reset(long i = 0) {
	m_counter = i;
      }

      _INLINE long inc() {
	return atomic_add64(&m_counter,1);
      }

      _INLINE long dec() {
	return atomic_add64(&m_counter,-1);
      }

      _INLINE long add(const long i) {
	return atomic_add64(&m_counter, i);
      }

      _INLINE long val() {
	return m_counter;
      };

      _INLINE AtomicCounter64 operator++(int)
      {
	return inc();
      }


      _INLINE AtomicCounter64 operator+=(long n)
      {
	return add(n);
      }

    };


#define MIN_DELAY_CYCLES  8
#define MAX_DELAY_CYCLES 64
#define MAX_NUM_CORES    64

  _INLINE void delayThread(unsigned int cycles)
  {
    _mm_delay_32(cycles);
  }


  _INLINE void delayThreadExpFalloff(unsigned int &cycles, const unsigned int max_delay_cycles = MAX_DELAY_CYCLES)
  {
    _mm_delay_32(cycles);
    cycles += cycles;
    cycles = min(cycles,max_delay_cycles);
  }


  // =============================================================================
  // === fast memory barrier with tree-based reduction and signal propagation, ===
  // === assumes four threads per core which is given by the tasking system   ====
  // =============================================================================

  class QuadTreeBarrier
  {
  public:

    __ALIGN(64)
      class CoreSyncData {
    public:
      volatile unsigned char threadState[2][4];
      volatile unsigned int mode;
      volatile unsigned int data[16-3];

      _INLINE void init() 
      {
	*(volatile unsigned int*)&threadState[0][0] = 0; 
	*(volatile unsigned int*)&threadState[1][0] = 0; 
	mode = 0;
	data[0] = 0;
      }

      _INLINE void prefetchL1Ex()
      { 
	prefetch<PFHINT_L1EX>((char*)threadState); 
      }

      _INLINE void prefetchL1()
      { 
	prefetch<PFHINT_L1>((char*)threadState); 
      }
    
      _INLINE void switchModeAndSendRunSignal(const unsigned int m)
      { 
	mode = 1 - mode;
	*(volatile unsigned int*)&threadState[m][0] = 0; 
      }


      _INLINE void setThreadStateToDone(const unsigned int m, const unsigned int threadID)
      {
	threadState[m][threadID % 4] = 1;
      }

      _INLINE bool allFourThreadsDone(const unsigned int m, const unsigned int orMask = 0) 
      { return (*(volatile unsigned int*)&threadState[m][0] | orMask)== 0x01010101; }

      _INLINE bool threadDone(const unsigned int m, const unsigned int threadID) 
      { return threadState[m][threadID % 4] == 1; }
 

      _INLINE void waitForAllThreadsOnCore(const unsigned int m)
      {
	unsigned int count = MIN_DELAY_CYCLES;
	while(likely(allFourThreadsDone(m) == false)) 
	  delayThread(count);
      }

      _INLINE void waitForAllOtherThreadsOnCore(const unsigned int m, const unsigned int threadID)
      {
	unsigned int count = MIN_DELAY_CYCLES;
	const unsigned int orMask = (unsigned int)1 << ((threadID % 4) * 8);
	while(likely(allFourThreadsDone(m,orMask) == false)) 
	  delayThread(count);
      }



      _INLINE void waitForThreadReceivesRunSignal(const unsigned int m, const unsigned int threadID)
      {
	unsigned int count = MIN_DELAY_CYCLES;
	while(likely(threadDone(m,threadID) == true)) 
	  delayThread(count);
      }

    };

    void sync(const unsigned int threadID, const unsigned int MAX_THREADS_SYNC);

    void syncWithReduction(const unsigned int threadID, 
			   const unsigned int MAX_THREADS_SYNC,
			   void (* reductionFct)(const unsigned int currentThreadID,
						 const unsigned int childThreadID));

  
    MIC_ALIGN CoreSyncData data[MAX_NUM_CORES]; // == one cacheline for each core ==
  
    QuadTreeBarrier()
      {    
	assert(sizeof(CoreSyncData) == 64);
	for (unsigned int i=0;i<MAX_NUM_CORES;i++) data[i].init();
      }

  };


  // =================================================
  // === fast memory barrier with linear reduction ===
  // =================================================

  __ALIGN(64)
    class SimpleBarrier
    {
    public:
      volatile unsigned int mode;
      volatile unsigned int flag0;
      volatile unsigned int flag1;
      volatile unsigned int fill[16-3];
      volatile unsigned char count0[MAX_MIC_THREADS]; 
      volatile unsigned char count1[MAX_MIC_THREADS]; 

      SimpleBarrier()
	{    
	  mode      = 0;
	  flag0     = 0;
	  flag1     = 0;
	  for (size_t i=0;i<MAX_MIC_THREADS;i++)
	    {
	      count0[i] = 0;
	      count1[i] = 0;
	    }
	}

      void sync(const unsigned int ID, const unsigned int MAX);

    };




  __ALIGN(64)
    class AtomicMutex
    {
    public:

      // ----------------
      // -- mutex flag --
      // ----------------

      volatile int flag;

      // -------------------------
      // -- additional elements --
      // -------------------------
 
      atomic_t m_counter;
      volatile unsigned int index;
      volatile char threadState[4];

      // -------------------------
      // -------------------------
      // -------------------------
  
      volatile int data[12]; 

      _INLINE void init()
      {
	flag = 0;
	m_counter = 0;
	index = 0;
	threadState[0] = 0;
	threadState[1] = 0;
	threadState[2] = 0;
	threadState[3] = 0;
      }

      AtomicMutex()
	{
	  init();
	}

      _INLINE void lock()
      {
	while(1) {
	  unsigned int numDelayCycles = 1024;
	  while( flag == 1) delayThread(numDelayCycles); // read without atomic op first
	  if ( atomic_cmpxchg(&flag,0,1) == 0) break;
	}
      }

      _INLINE void unlock()
      {
	flag = 0;
      }

      _INLINE void reset(int i = 0) {
	flag = i;
      }


      _INLINE void resetCounter(unsigned int i = 0) {
	*(volatile unsigned int*)&m_counter = i;
      }


      _INLINE unsigned int inc() {
	return atomic_add(&m_counter,1);
      }

      _INLINE unsigned int val() {
	return m_counter;
      };


    };

  template<class T, unsigned int SIZE>
    class WorkStack 
  {
  public:
    MIC_ALIGN AtomicMutex mutex;
    MIC_ALIGN T t[SIZE];

    _INLINE void init() {
      mutex.reset();
      mutex.index = 0;
      mutex.threadState[0] = 0;
      mutex.threadState[1] = 0;
      mutex.threadState[2] = 0;
      mutex.threadState[3] = 0;
    }


    _INLINE void prefetchL1Ex()
    {
      prefetch<PFHINT_L1EX>((int*)&mutex);
    }

    _INLINE void setThreadState(const unsigned int threadID, 
				const unsigned char val)
    {
      mutex.threadState[threadID % 4] = val;
    }

    _INLINE unsigned int getThreadState(const unsigned int threadID)
    {
      return mutex.threadState[threadID % 4];
    }


    _INLINE bool allThreadsDone()
    {
      return *(volatile unsigned int*)mutex.threadState == 0; 
    }

    _INLINE void waitUntilAllThreadsAreActive()
    {
      while (*(volatile unsigned int*)mutex.threadState != 0x01010101 ) delayThread(MAX_DELAY_CYCLES); 
    }


    _INLINE unsigned int getAllThreadStates()
    {
      return *(volatile unsigned int*)mutex.threadState; 
    }

    _INLINE WorkStack() {
      init();
    }

    _INLINE bool isFull(const unsigned int plus = 0)
    {
      return mutex.index + plus >= SIZE;
    }

    _INLINE bool isEmpty()
    {
      return mutex.index == 0;
    }
  
    _INLINE unsigned int usedSlots()
    {
      return mutex.index;
    }

    _INLINE void lock()
    {
      mutex.lock();
    }

    _INLINE void unlock()
    {
      mutex.unlock();
    }

    _INLINE bool request_lock(const unsigned int num)
    {
      lock();
      if (isFull(num)) { unlock(); return false; }    
      return true;
    }

    _INLINE void push_nolock(const T &v)
    {
      if (unlikely(isFull())) FATAL("stack is full");
      t[mutex.index++] = v;
    }

    _INLINE T *getNextElementPtr()
    {
      return &t[mutex.index];
    }

    _INLINE void incUsedSlots(const unsigned int num)
    {
      mutex.index += num;
    }

    _INLINE bool push(const T &v)
    {
      lock();
      if (isFull()) { unlock(); return false; }
      t[mutex.index++] = v;
      unlock();
      return true;
    }

    _INLINE bool pop_nolock(T &v)
    {
      if (isEmpty()) return false;
      v = t[--mutex.index];
      return true;
    }

    _INLINE bool pop_nolock_largest(T &v, const unsigned int threadID)
    {
      mutex.threadState[threadID % 4] = 0;
      if (isEmpty()) return false;
      unsigned int largest = 0;
      for (unsigned int i=1;i<mutex.index;i++)
	if (t[i] > t[largest])
	  largest = i;
      v = t[largest];
      t[largest] = t[mutex.index-1];
      mutex.index--;
      mutex.threadState[threadID % 4] = 1;

      return true;
    }



    _INLINE bool pop(T &v)
    {
      if (isEmpty()) { return false; }
      lock();
      if (isEmpty()) 
	{	
	  unlock(); return false; 
	}
      v = t[--mutex.index];
      unlock();
      return true;
    }


    _INLINE bool pop(T &v, const unsigned int threadID)
    {
      mutex.threadState[threadID % 4] = 0;
      if (isEmpty()) { return false; }
      lock();
      if (isEmpty()) 
	{	
	  unlock(); return false; 
	}
      v = t[--mutex.index];
      mutex.threadState[threadID % 4] = 1;
      unlock();
      return true;
    }

    _INLINE bool pop_largest(T &v, const unsigned int threadID)
    {
      mutex.threadState[threadID % 4] = 0;
      if (isEmpty()) { return false; }
      lock();
      if (isEmpty()) 
	{	
	  unlock(); return false; 
	}

      unsigned int largest = 0;
      for (unsigned int i=1;i<mutex.index;i++)
	if (t[i] > t[largest])
	  largest = i;
      v = t[largest];
      t[largest] = t[mutex.index-1];
      mutex.index--;

      mutex.threadState[threadID % 4] = 1;
      unlock();
      return true;
    }


    _INLINE bool pop_smallest(T &v, const unsigned int threadID)
    {
      mutex.threadState[threadID % 4] = 0;
      if (isEmpty()) { return false; }
      lock();
      if (isEmpty()) 
	{	
	  unlock(); return false; 
	}

      unsigned int smallest = 0;
      for (unsigned int i=1;i<mutex.index;i++)
	if (t[i] < t[smallest])
	  smallest = i;
      v = t[smallest];
      t[smallest] = t[mutex.index-1];
      mutex.index--;

      mutex.threadState[threadID % 4] = 1;
      unlock();
      return true;
    }


    _INLINE bool pop_first(T &v)
    {

      lock();
      if (isEmpty()) { unlock(); return false; }
      v = t[0];
      for (unsigned int i=0;i<mutex.index-1;i++)
	t[i] = t[i+1];
      mutex.index--;
      unlock();
      return true;
    }

    _INLINE T *get_addr(const unsigned int i)
    {
      return &t[i];
    }

    _INLINE const T &get(const unsigned int i) const
    {
      return t[i];
    }

  };


  _INLINE unsigned long long rdtsc()
  {
    return _rdtsc();
  }

  _INLINE void __cpuid(int *res, int in)
  {
    unsigned long ax, bx, cx, dx;

    ax = in;   
    __asm__ __volatile__ ("cpuid"
			  : "=a" (ax),
			    "=b" (bx),
			    "=c" (cx),
			    "=d" (dx)
			  : "a" (ax));
    res[0] = ax;
    res[1] = bx;
    res[2] = cx;
    res[3] = dx;
  }

  class Timer
  {
  public:
    unsigned long long tstart, toverhead;
    static float freqInMHz;

    Timer() {
      toverhead = 0;
      start();
      toverhead = stop();
    }

    _INLINE void start() {
      tstart = rdtsc();
    }

    _INLINE unsigned long long stop() const {
      unsigned long long tend = rdtsc();
      return tend - tstart - toverhead;
    }
  };

  template<class T> 
    void quicksort_ascending(T *__restrict__ t, 
			     const int& begin, 
			     const int& end)
    {
      if (likely(begin < end)) 
	{

	  const T pivotvalue = t[begin];
	  int left  = begin - 1;
	  int right = end   + 1;

	  while(1) 
	    {
	      while (t[--right] > pivotvalue);
	      while (t[++left] < pivotvalue);

	      if (left >= right) break;

	      const T temp = t[right];
	      t[right] = t[left];
	      t[left] = temp;
	    }

	  const int pivot = right;
	  quicksort_ascending(t, begin, pivot);
	  quicksort_ascending(t, pivot + 1, end);
	}
    }

  template<class T, bool less_fct(const T &a, const T &b)> 
    void quicksort_ascending_ext(T *__restrict__ t, 
			     const int& begin, 
			     const int& end)
    {
      if (likely(begin < end)) 
	{

	  const T pivotvalue = t[begin];
	  int left  = begin - 1;
	  int right = end   + 1;

	  while(1) 
	    {
	      while (less_fct(pivotvalue,t[--right]));
	      while (less_fct(t[++left],pivotvalue));

	      if (left >= right) break;

	      const T temp = t[right];
	      t[right] = t[left];
	      t[left] = temp;
	    }

	  const int pivot = right;
	  quicksort_ascending_ext<T,less_fct>(t, begin, pivot);
	  quicksort_ascending_ext<T,less_fct>(t, pivot + 1, end);
	}
    }


  template<class T> 
    void quicksort_decending(T *__restrict__ t, 
			     const int& begin, 
			     const int& end)
    {
      if (likely(begin < end)) 
	{

	  const T pivotvalue = t[begin];
	  int left  = begin - 1;
	  int right = end   + 1;

	  while(1) 
	    {
	      while (t[--right] < pivotvalue);
	      while (t[++left] > pivotvalue);

	      if (left >= right) break;

	      const T temp = t[right];
	      t[right] = t[left];
	      t[left] = temp;
	    }

	  const int pivot = right;
	  quicksort_decending(t, begin, pivot);
	  quicksort_decending(t, pivot + 1, end);
	}
    }


  template<class T>
    _INLINE T bitInterleave(T x, 
			    T y, 
			    T z)
    {
      x = (x | (x << 16)) & 0x030000FF; 
      x = (x | (x <<  8)) & 0x0300F00F; 
      x = (x | (x <<  4)) & 0x030C30C3; 
      x = (x | (x <<  2)) & 0x09249249; 


      y = (y | (y << 16)) & 0x030000FF; 
      y = (y | (y <<  8)) & 0x0300F00F; 
      y = (y | (y <<  4)) & 0x030C30C3; 
      y = (y | (y <<  2)) & 0x09249249; 


      z = (z | (z << 16)) & 0x030000FF; 
      z = (z | (z <<  8)) & 0x0300F00F; 
      z = (z | (z <<  4)) & 0x030C30C3; 
      z = (z | (z <<  2)) & 0x09249249; 


      return x | (y << 1) | (z << 2);
    }

  template<class T>
    _INLINE T bitInterleave(T u, T v)
    {
      u = (u | (u << 8)) & 0x00FF00FF;
      u = (u | (u << 4)) & 0x0F0F0F0F;
      u = (u | (u << 2)) & 0x33333333;
      u = (u | (u << 1)) & 0x55555555;
      v = (v | (v << 8)) & 0x00FF00FF;
      v = (v | (v << 4)) & 0x0F0F0F0F;
      v = (v | (v << 2)) & 0x33333333;
      v = (v | (v << 1)) & 0x55555555;
      return (u | (v << 1));
    }

  template<class T>
    _INLINE T bitReverse(T v)
    {
      v = ((v >> 1) & 0x55555555) | ((v & 0x55555555) << 1);
      v = ((v >> 2) & 0x33333333) | ((v & 0x33333333) << 2);
      v = ((v >> 4) & 0x0F0F0F0F) | ((v & 0x0F0F0F0F) << 4);
      v = ((v >> 8) & 0x00FF00FF) | ((v & 0x00FF00FF) << 8);
      v = ( v >> 16             ) | ( v               << 16);
      return v;
    }

  _INLINE unsigned int clz(const unsigned int x)
  {
    return _lzcnt_u32(x);
  }


  _INLINE void do_cpuid(unsigned int eax, unsigned int *p)
  {
    __asm __volatile("cpuid"
		     : "=a" (p[0]), "=b" (p[1]), "=c" (p[2]), "=d" (p[3])
		     :  "0" (eax));
  }

};
#endif
