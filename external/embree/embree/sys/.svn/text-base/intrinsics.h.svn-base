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

#ifndef __EMBREE_INTRINSICS_H__
#define __EMBREE_INTRINSICS_H__

#include "platform.h"
#include <xmmintrin.h>

#if defined __MIC__
#include <immintrin.h>
#endif

////////////////////////////////////////////////////////////////////////////////
/// Windows Platform
////////////////////////////////////////////////////////////////////////////////

#if defined(__WIN32__)

#include <intrin.h>

__forceinline uint64 __rdpmc(int i) {
  return __readpmc(i);
}

__forceinline int __popcnt(int in) {
  return _mm_popcnt_u32(in);
}

#if !defined(_MSC_VER)
__forceinline unsigned int __popcnt(unsigned int in) {
  return _mm_popcnt_u32(in);
}
#endif

#if defined(__X86_64__)
__forceinline long long __popcnt(long long in) {
  return _mm_popcnt_u64(in);
}
__forceinline size_t __popcnt(size_t in) {
  return _mm_popcnt_u64(in);
}
#endif

__forceinline int __bsf(int v) {
  unsigned long r = 0; _BitScanForward(&r,v); return r;
}

__forceinline int __bsr(int v) {
  unsigned long r = 0; _BitScanReverse(&r,v); return r;
}

__forceinline int __btc(int v, int i) {
  long r = v; _bittestandcomplement(&r,i); return r;
}

__forceinline int __bts(int v, int i) {
  long r = v; _bittestandset(&r,i); return r;
}

__forceinline int __btr(int v, int i) {
  long r = v; _bittestandreset(&r,i); return r;
}

#if defined(__X86_64__)

__forceinline size_t __bsf(size_t v) {
  size_t r = 0; _BitScanForward64((unsigned long*)&r,v); return r;
}

__forceinline size_t __bsr(size_t v) {
  size_t r = 0; _BitScanReverse64((unsigned long*)&r,v); return r;
}

__forceinline size_t __btc(size_t v, size_t i) {
  //size_t r = v; _bittestandcomplement64((__int64*)&r,i); return r;
  return v ^ (size_t(1) << i); // faster than using intrinsics, as intrinsic goes through memory
}

__forceinline size_t __bts(size_t v, size_t i) {
  __int64 r = v; _bittestandset64(&r,i); return r;
}

__forceinline size_t __btr(size_t v, size_t i) {
  __int64 r = v; _bittestandreset64(&r,i); return r;
}

#endif

#if defined(__X86_64__)

typedef int64 atomic_t;

__forceinline int64 atomic_add(volatile int64* m, const int64 v) {
  return _InterlockedExchangeAdd64(m,v);
}

__forceinline int64 atomic_xchg(volatile int64 *p, int64 v) {
  return _InterlockedExchange64((volatile long long *)p, v);
}

__forceinline int64 atomic_cmpxchg(volatile int64* m, const int64 v, const int64 c) {
  return _InterlockedCompareExchange64(m,v,c);
}

#else

typedef int32 atomic_t;

__forceinline int32 atomic_add(volatile int32* p, const int32 v) {
  return _InterlockedExchangeAdd((volatile long*)p,v);
}

__forceinline int32 atomic_xchg(volatile int32 *p, int32 v) {
  return _InterlockedExchange((volatile long*)p, v);
}

__forceinline int32 atomic_cmpxchg(volatile int32* p, const int32 v, const int32 c) {
  return _InterlockedCompareExchange((volatile long*)p,v,c);
}

#endif

////////////////////////////////////////////////////////////////////////////////
/// Unix Platform
////////////////////////////////////////////////////////////////////////////////

#else

__forceinline void __cpuid(int out[4], int op) {
  asm volatile ("cpuid" : "=a"(out[0]), "=b"(out[1]), "=c"(out[2]), "=d"(out[3]) : "a"(op)); 
}

__forceinline uint64 __rdtsc()  {
  uint32 high,low;
  asm volatile ("rdtsc" : "=d"(high), "=a"(low));
  return (((uint64)high) << 32) + (uint64)low;
}

__forceinline uint64 __rdpmc(int i) {
  uint32 high,low;
  asm volatile ("rdpmc" : "=d"(high), "=a"(low) : "c"(i));
  return (((uint64)high) << 32) + (uint64)low;
}

__forceinline unsigned int __popcnt(unsigned int in) {
  int r = 0; asm ("popcnt %1,%0" : "=r"(r) : "r"(in)); return r;
}

__forceinline int __bsf(int v) {
  int r = 0; asm ("bsf %1,%0" : "=r"(r) : "r"(v)); return r;
}

__forceinline int __bsr(int v) {
  int r = 0; asm ("bsr %1,%0" : "=r"(r) : "r"(v)); return r;
}

__forceinline int __btc(int v, int i) {
  int r = 0; asm ("btc %1,%0" : "=r"(r) : "r"(i), "0"(v) : "flags" ); return r;
}

__forceinline int __bts(int v, int i) {
  int r = 0; asm ("bts %1,%0" : "=r"(r) : "r"(i), "0"(v) : "flags"); return r;
}

__forceinline int __btr(int v, int i) {
  int r = 0; asm ("btr %1,%0" : "=r"(r) : "r"(i), "0"(v) : "flags"); return r;
}

__forceinline size_t __bsf(size_t v) {
  size_t r = 0; asm ("bsf %1,%0" : "=r"(r) : "r"(v)); return r;
}

__forceinline size_t __bsr(size_t v) {
  size_t r = 0; asm ("bsr %1,%0" : "=r"(r) : "r"(v)); return r;
}

__forceinline size_t __btc(size_t v, size_t i) {
  size_t r = 0; asm ("btc %1,%0" : "=r"(r) : "r"(i), "0"(v) : "flags" ); return r;
}

__forceinline size_t __bts(size_t v, size_t i) {
  size_t r = 0; asm ("bts %1,%0" : "=r"(r) : "r"(i), "0"(v) : "flags"); return r;
}

__forceinline size_t __btr(size_t v, size_t i) {
  size_t r = 0; asm ("btr %1,%0" : "=r"(r) : "r"(i), "0"(v) : "flags"); return r;
}

#if defined(__X86_64__)

typedef int64 atomic_t;

__forceinline int64 atomic_add( int64 volatile* value, int64 input ) {
  return __sync_fetch_and_add(value, input);
}

__forceinline int64 atomic_xchg( int64 volatile* value, int64 input ) {
  return __sync_lock_test_and_set(value, input);
}

__forceinline int64 atomic_cmpxchg( int64 volatile* value, const int64 input, int64 comparand ) {
  return __sync_val_compare_and_swap(value, comparand, input);
}

#else

typedef int32 atomic_t;

__forceinline int32 atomic_add( int32 volatile* value, int32 input ) {
  return __sync_fetch_and_add(value, input);
}

__forceinline int32 atomic_xchg( int32 volatile* value, int32 input ) {
  return __sync_lock_test_and_set(value, input);
}

__forceinline int32 atomic_cmpxchg( int32 volatile* value, const int32 input, int32 comparand ) {
  return __sync_val_compare_and_swap(value, comparand, input);
}

#endif

#endif

////////////////////////////////////////////////////////////////////////////////
/// All Platforms
////////////////////////////////////////////////////////////////////////////////

#if defined(__X86_64__)

template<typename T>
__forceinline T* atomic_xchg_ptr( T* volatile* value, const T* input)
{  return (T*)atomic_xchg((int64*)value,(int64)input); }

template<typename T>
__forceinline T* atomic_cmpxchg_ptr( T* volatile* value, const T* input, T* comparand )
{  return (T*)atomic_cmpxchg((int64*)value,(int64)input,(int64)comparand); }

#else

template<typename T>
__forceinline T* atomic_xchg_ptr( T* volatile* value, const T* input)
{  return (T*)atomic_xchg((int32*)value,(int32)input); }

template<typename T>
__forceinline T* atomic_cmpxchg_ptr( T* volatile* value, const T* input, T* comparand )
{  return (T*)atomic_cmpxchg((int32*)value,(int32)input,(int32)comparand); }

#endif

__forceinline uint64 rdtsc()
{
  int dummy[4]; 
  __cpuid(dummy,0); 
  uint64 clock = __rdtsc(); 
  __cpuid(dummy,0); 
  return clock;
}

__forceinline int cast_f2i(float f) {
  union { float f; int i; } v; v.f = f; return v.i;
}

__forceinline float cast_i2f(int i) {
  union { float f; int i; } v; v.i = i; return v.f;
}

#if defined(__MIC__)
__forceinline void __pause (const int cycles = 16) { _mm_delay_32(cycles); }
#else
__forceinline void __pause (const int cycles = 0) {  _mm_pause(); }
#endif

/* prefetches */
__forceinline void prefetchL1 (const void* ptr) { _mm_prefetch((const char*)ptr,_MM_HINT_T0); }
__forceinline void prefetchL2 (const void* ptr) { _mm_prefetch((const char*)ptr,_MM_HINT_T1); }
__forceinline void prefetchL3 (const void* ptr) { _mm_prefetch((const char*)ptr,_MM_HINT_T2); }
__forceinline void prefetchNTA(const void* ptr) { _mm_prefetch((const char*)ptr,_MM_HINT_NTA); }

#endif
