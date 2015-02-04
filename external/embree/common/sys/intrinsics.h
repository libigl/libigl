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

#include "platform.h"

#if defined(__SSE__)
#include <xmmintrin.h>
#endif

#if defined(__SSE2__)
#include <emmintrin.h>
#endif

#if defined(__SSE3__)
#include <pmmintrin.h>
#endif

#if defined(__SSSE3__)
#include <tmmintrin.h>
#endif

#if defined (__SSE4_1__)
#include <smmintrin.h>
#endif

#if defined (__SSE4_2__)
#include <nmmintrin.h>
#endif

#if defined(__AVX__) || defined(__MIC__)
#include <immintrin.h>
#endif

#if defined(__BMI__) && defined(__GNUC__)
#define _tzcnt_u32 __tzcnt_u32
#define _tzcnt_u64 __tzcnt_u64
#endif

#if defined(__LZCNT__)
#define _lzcnt_u32 __lzcnt32
#define _lzcnt_u64 __lzcnt64
#endif

#if defined(__MIC__)
#include <immintrin.h>
#endif

////////////////////////////////////////////////////////////////////////////////
/// Windows Platform
////////////////////////////////////////////////////////////////////////////////

#if defined(__WIN32__)

#include <intrin.h>

//FIXME
__forceinline size_t read_tsc()  {
	return 0;
}

#if defined(__SSE4_2__)

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

#endif

__forceinline int __bsf(int v) {
#if defined(__AVX2__) 
  return _tzcnt_u32(v);
#else
  unsigned long r = 0; _BitScanForward(&r,v); return r;
#endif
}

__forceinline unsigned int __bsf(unsigned int v) {
#if defined(__AVX2__) 
  return _tzcnt_u32(v);
#else
  unsigned long r = 0; _BitScanForward(&r,v); return r;
#endif
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

__forceinline int bitscan(int v) {
#if defined(__AVX2__) 
  return _tzcnt_u32(v);
#else
  return __bsf(v);
#endif
}

__forceinline int clz(const int x)
{
#if defined(__AVX2__)
  return _lzcnt_u32(x);
#else
  if (unlikely(x == 0)) return 32;
  return 31 - __bsr(x);    
#endif
}

__forceinline int __bscf(int& v) 
{
  int i = __bsf(v);
  v &= v-1;
  return i;
}

__forceinline unsigned int __bscf(unsigned int& v) 
{
  unsigned int i = __bsf(v);
  v &= v-1;
  return i;
}

#if defined(__X86_64__)

__forceinline size_t __bsf(size_t v) {
#if defined(__AVX2__) 
  return _tzcnt_u64(v);
#else
  unsigned long r = 0; _BitScanForward64(&r,v); return r;
#endif
}

__forceinline size_t __bsr(size_t v) {
  unsigned long r = 0; _BitScanReverse64(&r,v); return r;
}

__forceinline size_t __btc(size_t v, size_t i) {
  size_t r = v; _bittestandcomplement64((__int64*)&r,i); return r;
}

__forceinline size_t __bts(size_t v, size_t i) {
  __int64 r = v; _bittestandset64(&r,i); return r;
}

__forceinline size_t __btr(size_t v, size_t i) {
  __int64 r = v; _bittestandreset64(&r,i); return r;
}

__forceinline size_t bitscan(size_t v) {
#if defined(__AVX2__)
#if defined(__X86_64__)
  return _tzcnt_u64(v);
#else
  return _tzcnt_u32(v);
#endif
#else
  return __bsf(v);
#endif
}

__forceinline size_t __bscf(size_t& v) 
{
  size_t i = __bsf(v);
  v &= v-1;
  return i;
}

#endif

typedef int32 atomic32_t;

__forceinline int32 atomic_add(volatile int32* p, const int32 v) {
  return _InterlockedExchangeAdd((volatile long*)p,v);
}

__forceinline int32 atomic_sub(volatile int32* p, const int32 v) {
  return _InterlockedExchangeAdd((volatile long*)p,-v);
}

__forceinline int32 atomic_xchg(volatile int32 *p, int32 v) {
  return _InterlockedExchange((volatile long*)p, v);
}

__forceinline int32 atomic_cmpxchg(volatile int32* p, const int32 c, const int32 v) {
  return _InterlockedCompareExchange((volatile long*)p,v,c);
}

#if defined(__X86_64__)

typedef int64 atomic64_t;

__forceinline int64 atomic_add(volatile int64* m, const int64 v) {
  return _InterlockedExchangeAdd64(m,v);
}

__forceinline int64 atomic_sub(volatile int64* m, const int64 v) {
  return _InterlockedExchangeAdd64(m,-v);
}

__forceinline int64 atomic_xchg(volatile int64 *p, int64 v) {
  return _InterlockedExchange64((volatile long long *)p, v);
}

__forceinline int64 atomic_cmpxchg(volatile int64* m, const int64 c, const int64 v) {
  return _InterlockedCompareExchange64(m,v,c);
}

#endif

////////////////////////////////////////////////////////////////////////////////
/// Unix Platform
////////////////////////////////////////////////////////////////////////////////

#else

#if defined(__i386__) && defined(__PIC__)

__forceinline void __cpuid(int out[4], int op) 
{
  asm volatile ("xchg{l}\t{%%}ebx, %1\n\t"
                "cpuid\n\t"
                "xchg{l}\t{%%}ebx, %1\n\t"
                : "=a"(out[0]), "=r"(out[1]), "=c"(out[2]), "=d"(out[3]) 
                : "0"(op)); 
}

__forceinline void __cpuid_count(int out[4], int op1, int op2) 
{
  asm volatile ("xchg{l}\t{%%}ebx, %1\n\t"
                "cpuid\n\t"
                "xchg{l}\t{%%}ebx, %1\n\t"
                : "=a" (out[0]), "=r" (out[1]), "=c" (out[2]), "=d" (out[3])
                : "0" (op1), "2" (op2)); 
}

#else

__forceinline void __cpuid(int out[4], int op) {
  asm volatile ("cpuid" : "=a"(out[0]), "=b"(out[1]), "=c"(out[2]), "=d"(out[3]) : "a"(op)); 
}

__forceinline void __cpuid_count(int out[4], int op1, int op2) {
  asm volatile ("cpuid" : "=a"(out[0]), "=b"(out[1]), "=c"(out[2]), "=d"(out[3]) : "a"(op1), "c"(op2)); 
}

#endif

__forceinline uint64 read_tsc()  {
  uint32 high,low;
  asm volatile ("rdtsc" : "=d"(high), "=a"(low));
  return (((uint64)high) << 32) + (uint64)low;
}

#if !defined(__MIC__)

#if defined(__SSE4_2__)
__forceinline unsigned int __popcnt(unsigned int in) {
  int r = 0; asm ("popcnt %1,%0" : "=r"(r) : "r"(in)); return r;
}
#endif

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

__forceinline unsigned int __bsf(unsigned int v) {
  unsigned int r = 0; asm ("bsf %1,%0" : "=r"(r) : "r"(v)); return r;
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

__forceinline int bitscan(int v) {
#if defined(__AVX2__) 
  return _tzcnt_u32(v);
#else
  return __bsf(v);
#endif
}

__forceinline unsigned int bitscan(unsigned int v) {
#if defined(__AVX2__) 
  return _tzcnt_u32(v);
#else
  return __bsf(v);
#endif
}

__forceinline size_t bitscan(size_t v) {
#if defined(__AVX2__)
#if defined(__X86_64__)
  return _tzcnt_u64(v);
#else
  return _tzcnt_u32(v);
#endif
#else
  return __bsf(v);
#endif
}

__forceinline int clz(const int x)
{
#if defined(__AVX2__)
  return _lzcnt_u32(x);
#else
  if (unlikely(x == 0)) return 32;
  return 31 - __bsr(x);    
#endif
}

__forceinline int __bscf(int& v) 
{
  int i = bitscan(v);
#if defined(__AVX2__)
  v &= v-1;
#else
  v = __btc(v,i);
#endif
  return i;
}

__forceinline unsigned int __bscf(unsigned int& v) 
{
  unsigned int i = bitscan(v);
  v &= v-1;
  return i;
}

__forceinline size_t __bscf(size_t& v) 
{
  size_t i = bitscan(v);
#if defined(__AVX2__)
  v &= v-1;
#else
  v = __btc(v,i);
#endif
  return i;
}

#else

__forceinline unsigned int clz(const unsigned int x) {
  return _lzcnt_u32(x); 
}

__forceinline size_t clz(const size_t x) {
  return _lzcnt_u64(x); 
}

__forceinline unsigned int bitscan(unsigned int v) {
  return _mm_tzcnt_32(v); 
}

__forceinline size_t bitscan64(size_t v) {
  return _mm_tzcnt_64(v); 
}

__forceinline unsigned int bitscan(const int index, const unsigned int v) { 
  return _mm_tzcnti_32(index,v); 
};

__forceinline size_t bitscan64(const ssize_t index, const size_t v) { 
  return _mm_tzcnti_64(index,v); 
};

__forceinline int __popcnt(int v) {
  return _mm_countbits_32(v); 
}

__forceinline unsigned int __popcnt(unsigned int v) {
  return _mm_countbits_32(v); 
}

__forceinline unsigned int countbits(unsigned int v) {
  return _mm_countbits_32(v); 
};

__forceinline size_t __popcnt(size_t v) {
  return _mm_countbits_64(v); 
}

__forceinline size_t countbits64(size_t v) { 
  return _mm_countbits_64(v); 
};

__forceinline int __bsf(int v) {
  return bitscan(v); 
}

__forceinline unsigned int __bsf(unsigned int v) {
  return bitscan(v); 
}

__forceinline size_t __bsf(size_t v) {
  return bitscan(v); 
}

__forceinline size_t __btc(size_t v, size_t i) {
  return v ^ (size_t(1) << i); 
}

__forceinline unsigned int __bsr(unsigned int v) {
  return 31 - _lzcnt_u32(v); 
}

__forceinline size_t __bsr(size_t v) {
  return 63 - _lzcnt_u64(v); 
}

#endif

#if defined(__X86_64__)

typedef int64 atomic64_t;

__forceinline int64 atomic_add( int64 volatile* value, int64 input ) {
  return __sync_fetch_and_add(value, input);
}

__forceinline int64 atomic_sub( int64 volatile* value, int64 input ) {
  return __sync_fetch_and_add(value, -input);
}

__forceinline int64 atomic_xchg( int64 volatile* value, int64 input ) {
  return __sync_lock_test_and_set(value, input);
}

__forceinline int64 atomic_cmpxchg( int64 volatile* value, int64 comparand, const int64 input) {
  return __sync_val_compare_and_swap(value, comparand, input);
}

#endif

typedef int32 atomic32_t;

__forceinline int32 atomic_add( int32 volatile* value, int32 input ) {
  return __sync_fetch_and_add(value, input);
}

__forceinline int32 atomic_sub( int32 volatile* value, int32 input ) {
  return __sync_fetch_and_add(value, -input);
}

__forceinline int32 atomic_xchg( int32 volatile* value, int32 input ) {
  return __sync_lock_test_and_set(value, input);
}

__forceinline int32 atomic_cmpxchg( int32 volatile* value, int32 comparand, const int32 input ) {
  return __sync_val_compare_and_swap(value, comparand, input);
}

typedef int8 atomic8_t;

__forceinline int8 atomic_add( int8 volatile* value, int8 input ) {
  return __sync_fetch_and_add(value, input);
}

__forceinline int8 atomic_sub( int8 volatile* value, int8 input ) {
  return __sync_fetch_and_add(value, -input);
}

__forceinline int8 atomic_xchg( int8 volatile* value, int8 input ) {
  return __sync_lock_test_and_set(value, input);
}

__forceinline int8 atomic_cmpxchg( int8 volatile* value, int8 comparand, const int8 input ) {
  return __sync_val_compare_and_swap(value, comparand, input);
}

#endif

////////////////////////////////////////////////////////////////////////////////
/// All Platforms
////////////////////////////////////////////////////////////////////////////////

#if defined(__X86_64__)
  typedef atomic64_t atomic_t;
#else
  typedef atomic32_t atomic_t;
#endif

#if defined(__X86_64__)

template<typename T>
__forceinline T* atomic_xchg_ptr( T* volatile* value, const T* input)
{  return (T*)atomic_xchg((int64*)value,(int64)input); }

template<typename T>
__forceinline T* atomic_cmpxchg_ptr( T* volatile* value, T* comparand, const T* input )
{  return (T*)atomic_cmpxchg((int64*)value,(int64)comparand,(int64)input); }

#else

template<typename T>
__forceinline T* atomic_xchg_ptr( T* volatile* value, const T* input)
{  return (T*)atomic_xchg((int32*)value,(int32)input); }

template<typename T>
__forceinline T* atomic_cmpxchg_ptr( T* volatile* value,  T* comparand, const T* input )
{  return (T*)atomic_cmpxchg((int32*)value,(int32)comparand,(int32)input); }

#endif

__forceinline void atomic_min_f32(volatile float *__restrict__ ptr, const float b)
{
  const int int_b = *(int*)&b;
  while (1)
  {
    float a = *ptr;
    if (a <= b) break;
    const int int_a = *(int*)&a;
    const int result = atomic_cmpxchg((int*)ptr,int_a,int_b);
    if (result == int_a) break;
  }
}

__forceinline void atomic_max_f32(volatile float *__restrict__ ptr, const float b)
{
  const int int_b = *(int*)&b;
  while (1)
  {
    float a = *ptr;
    if (a >= b) break;
    const int int_a = *(int*)&a;
    const int result = atomic_cmpxchg((int*)ptr,int_a,int_b);
    if (result == int_a) break;
  }
}

__forceinline void atomic_min_i32(volatile int *__restrict__ ptr, const int b)
{
  while (1)
  {
    int a = *ptr;
    if (a <= b) break;
    const int int_a = *(int*)&a;
    const int result = atomic_cmpxchg((int*)ptr,int_a,b);
    if (result == int_a) break;
  }
}

__forceinline void atomic_max_i32(volatile int *__restrict__ ptr, const int b)
{
  while (1)
  {
    int a = *ptr;
    if (a >= b) break;
    const int int_a = *(int*)&a;
    const int result = atomic_cmpxchg((int*)ptr,int_a,b);
    if (result == int_a) break;
  }
}

__forceinline void atomic_min_ui32(volatile unsigned int *__restrict__ ptr, const unsigned int b)
{
  while (1)
  {
    unsigned int a = *ptr;
    if (a <= b) break;
    const unsigned int int_a = *(unsigned int*)&a;
    const unsigned int result = atomic_cmpxchg((int*)ptr,int_a,b);
    if (result == int_a) break;
  }
}

__forceinline void atomic_max_ui32(volatile unsigned int *__restrict__ ptr, const unsigned int b)
{
  while (1)
  {
    unsigned int a = *ptr;
    if (a >= b) break;
    const unsigned int int_a = *(unsigned int*)&a;
    const unsigned int result = atomic_cmpxchg((int*)ptr,int_a,b);
    if (result == int_a) break;
  }
}

static const unsigned int BITSCAN_NO_BIT_SET_32 = 32;
static const size_t       BITSCAN_NO_BIT_SET_64 = 64;


__forceinline uint64 rdtsc()
{
#if !defined(__MIC__)
  int dummy[4]; 
  __cpuid(dummy,0); 
  uint64 clock = read_tsc(); 
  __cpuid(dummy,0); 
  return clock;
#else
  return read_tsc(); 
#endif
}

__forceinline int cast_f2i(float f) {
  union { float f; int i; } v; v.f = f; return v.i;
}

__forceinline float cast_i2f(int i) {
  union { float f; int i; } v; v.i = i; return v.f;
}

#if defined(__MIC__)
__forceinline void __pause_cpu (const unsigned int cycles = 256) { 
  _mm_delay_32(cycles); 
}
#else
__forceinline void __pause_cpu (const int cycles = 0) {
  for (size_t i=0; i<8; i++)
    _mm_pause();    
}
#endif

__forceinline void __pause_cpu_expfalloff(unsigned int &cycles, const unsigned int max_cycles) 
{ 
  __pause_cpu(cycles);
  cycles += cycles;
  if (cycles > max_cycles) 
    cycles = max_cycles;
}

/* prefetches */
__forceinline void prefetchL1 (const void* ptr) { _mm_prefetch((const char*)ptr,_MM_HINT_T0); }
__forceinline void prefetchL2 (const void* ptr) { _mm_prefetch((const char*)ptr,_MM_HINT_T1); }
__forceinline void prefetchL3 (const void* ptr) { _mm_prefetch((const char*)ptr,_MM_HINT_T2); }
__forceinline void prefetchNTA(const void* ptr) { _mm_prefetch((const char*)ptr,_MM_HINT_NTA); }

__forceinline void prefetchL1EX(const void* ptr) { 
#if defined(__MIC__)
  _mm_prefetch((const char*)ptr,_MM_HINT_ET0); 
#else
  _mm_prefetch((const char*)ptr,_MM_HINT_T0); 
#endif
}

__forceinline void prefetchL2EX(const void* ptr) { 
#if defined(__MIC__)
  _mm_prefetch((const char*)ptr,_MM_HINT_ET2); 
#else
  _mm_prefetch((const char*)ptr,_MM_HINT_T1); 
#endif
}
