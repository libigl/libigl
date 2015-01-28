/*
 Copyright (c) 2010, Intel Corporation. All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of Intel Corporation nor the names of its contributors may
   be used to endorse or promote products derived from this software without
   specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
   THE POSSIBILITY OF SUCH DAMAGE.
*/

/***

 Provide feedback to: maxim.locktyukhin intel com, phil.j.kerly intel com

 Version 1.0 - Initial release.

    This AVX intrinsics emulation header file designed to work with Intel C/C++
    as well as GCC compilers.

    Known Issues and limitations:

    - does not support immediate values higher than 0x7 for _mm[256]_cmp_[ps|pd]
      intrinsics, UD2 instruction will be generated instead

    - -O0 optimization level may _sometimes_ result with compile time errors due
      to failed forced inline and compiler not being able to generate instruction
      with constant immediate operand becasue of it, compiling with -O1 and/or
      -finline-functions should help.

***/


#ifndef __EMU_M256_AVXIMMINTRIN_EMU_H__
#define __EMU_M256_AVXIMMINTRIN_EMU_H__

#ifndef __GNUC__
#pragma message (" --- Intel remark: AVX intrinsics are emulated with SSE ---")
#endif

/*
 * Intel(R) AVX compiler intrinsics.
 */

#ifdef __cplusplus
extern "C" {
#endif

/*
 * This is an emulation of Intel AVX
 */

#if defined( _MSC_VER ) || defined( __INTEL_COMPILER )
    #define __EMU_M256_ALIGN( a ) __declspec(align(a))
    #define __emu_inline          __forceinline
    #define __emu_int64_t         __int64
#elif defined( __GNUC__ )
    #define __EMU_M256_ALIGN( a ) __attribute__((__aligned__(a)))
    #define __emu_inline          __inline __attribute__((__always_inline__))
    #define __emu_int64_t         long long
#else
    #error "unsupported platform"
#endif

typedef union __EMU_M256_ALIGN(16) __emu__m256
{
    float       __emu_arr[8];
    __m128      __emu_m128[2];
} __emu__m256;

typedef union __EMU_M256_ALIGN(16) __emu__m256d
{
    double      __emu_arr[4];
    __m128d     __emu_m128[2];
} __emu__m256d;

typedef union __EMU_M256_ALIGN(16) __emu__m256i
{
    int         __emu_arr[8];
    __m128i     __emu_m128[2];
} __emu__m256i;

static __emu_inline __emu__m256  __emu_set_m128( const __m128 arr[] ) { __emu__m256 ret;  ret.__emu_m128[0] = arr[0]; ret.__emu_m128[1] = arr[1]; return (ret); }
static __emu_inline __emu__m256d __emu_set_m128d( const __m128d arr[] ) { __emu__m256d ret; ret.__emu_m128[0] = arr[0]; ret.__emu_m128[1] = arr[1]; return (ret); }
static __emu_inline __emu__m256i __emu_set_m128i( const __m128i arr[] ) { __emu__m256i ret; ret.__emu_m128[0] = arr[0]; ret.__emu_m128[1] = arr[1]; return (ret); }


#define __EMU_M256_IMPL_M1( type, func ) \
static __emu_inline __emu##type __emu_mm256_##func( __emu##type m256_param1 ) \
{   __emu##type res; \
    res.__emu_m128[0] = _mm_##func( m256_param1.__emu_m128[0] ); \
    res.__emu_m128[1] = _mm_##func( m256_param1.__emu_m128[1] ); \
    return ( res ); \
}

#define __EMU_M256_IMPL_M1_RET( ret_type, type, func ) \
static __emu_inline __emu##ret_type __emu_mm256_##func( __emu##type m256_param1 ) \
{   __emu##ret_type res; \
    res.__emu_m128[0] = _mm_##func( m256_param1.__emu_m128[0] ); \
    res.__emu_m128[1] = _mm_##func( m256_param1.__emu_m128[1] ); \
    return ( res ); \
}

#define __EMU_M256_IMPL_M1_RET_NAME( ret_type, type, func, name ) \
    static __emu_inline __emu##ret_type __emu_mm256_##name( __emu##type m256_param1 ) \
{   __emu##ret_type res; \
    res.__emu_m128[0] = _mm_##func( m256_param1.__emu_m128[0] ); \
    res.__emu_m128[1] = _mm_##func( m256_param1.__emu_m128[1] ); \
    return ( res ); \
}

#define __EMU_M256_IMPL_M1_LH( type, type_128, func ) \
static __emu_inline __emu##type __emu_mm256_##func( type_128 m128_param ) \
{   __emu##type res; \
    res.__emu_m128[0] = _mm_##func( m128_param ); \
    __m128 m128_param_high = _mm_movehl_ps( *(__m128*)&m128_param, *(__m128*)&m128_param ); \
    res.__emu_m128[1] = _mm_##func( *(type_128*)&m128_param_high ); \
    return ( res ); \
}

#define __EMU_M256_IMPL_M1_HL( type_128, type, func ) \
static __emu_inline type_128 __emu_mm256_##func( __emu##type m256_param1 ) \
{   type_128 res, tmp; \
    res = _mm_##func( m256_param1.__emu_m128[0] ); \
    tmp = _mm_##func( m256_param1.__emu_m128[1] ); \
    *(((__emu_int64_t*)&res)+1) = *(__emu_int64_t*)&tmp; \
    return ( res ); \
}

#define __EMU_M256_IMPL_M1P_DUP( type, type_param, func ) \
static __emu_inline __emu##type __emu_mm256_##func( type_param param ) \
{   __emu##type res; \
    res.__emu_m128[0] = _mm_##func( param ); \
    res.__emu_m128[1] = _mm_##func( param ); \
    return ( res ); \
}

#define __EMU_M256_IMPL_M1I_DUP( type, func ) \
    static __emu_inline __emu##type __emu_mm256_##func( __emu##type m256_param1, const int param2 ) \
{   __emu##type res; \
    res.__emu_m128[0] = _mm_##func( m256_param1.__emu_m128[0], param2 ); \
    res.__emu_m128[1] = _mm_##func( m256_param1.__emu_m128[1], param2 ); \
    return ( res ); \
}

#define __EMU_M256_IMPL2_M1I_DUP( type, func ) \
static __emu_inline __emu##type __emu_mm256_##func( __emu##type m256_param1, const int param2 ) \
{   __emu##type res; \
    res.__emu_m128[0] = __emu_mm_##func( m256_param1.__emu_m128[0], param2 ); \
    res.__emu_m128[1] = __emu_mm_##func( m256_param1.__emu_m128[1], param2 ); \
    return ( res ); \
}

#define __EMU_M256_IMPL2_M1I_SHIFT( type, func, shift_for_hi ) \
static __emu_inline __emu##type __emu_mm256_##func( __emu##type m256_param1, const int param2 ) \
{   __emu##type res; \
    res.__emu_m128[0] = __emu_mm_##func( m256_param1.__emu_m128[0], param2 & ((1<<shift_for_hi)-1) ); \
    res.__emu_m128[1] = __emu_mm_##func( m256_param1.__emu_m128[1], param2 >> shift_for_hi); \
    return ( res ); \
}

#define __EMU_M256_IMPL_M2( type, func ) \
static __emu_inline __emu##type __emu_mm256_##func( __emu##type m256_param1, __emu##type m256_param2 ) \
{   __emu##type res; \
    res.__emu_m128[0] = _mm_##func( m256_param1.__emu_m128[0], m256_param2.__emu_m128[0] ); \
    res.__emu_m128[1] = _mm_##func( m256_param1.__emu_m128[1], m256_param2.__emu_m128[1] ); \
    return ( res ); \
}

#define __EMU_M256_IMPL2_M2T( type, type_2, func ) \
static __emu_inline __emu##type __emu_mm256_##func( __emu##type m256_param1, __emu##type_2 m256_param2 ) \
{   __emu##type res; \
    res.__emu_m128[0] = __emu_mm_##func( m256_param1.__emu_m128[0], m256_param2.__emu_m128[0] ); \
    res.__emu_m128[1] = __emu_mm_##func( m256_param1.__emu_m128[1], m256_param2.__emu_m128[1] ); \
    return ( res ); \
}

#define __EMU_M256_IMPL_M2I_DUP( type, func ) \
static __emu_inline __emu##type __emu_mm256_##func( __emu##type m256_param1, __emu##type m256_param2, const int param3 ) \
{   __emu##type res; \
    res.__emu_m128[0] = _mm_##func( m256_param1.__emu_m128[0], m256_param2.__emu_m128[0], param3 ); \
    res.__emu_m128[1] = _mm_##func( m256_param1.__emu_m128[1], m256_param2.__emu_m128[1], param3 ); \
    return ( res ); \
}

#define __EMU_M256_IMPL2_M2I_DUP( type, func ) \
static __emu_inline __emu##type __emu_mm256_##func( __emu##type m256_param1, __emu##type m256_param2, const int param3 ) \
{   __emu##type res; \
    res.__emu_m128[0] = __emu_mm_##func( m256_param1.__emu_m128[0], m256_param2.__emu_m128[0], param3 ); \
    res.__emu_m128[1] = __emu_mm_##func( m256_param1.__emu_m128[1], m256_param2.__emu_m128[1], param3 ); \
    return ( res ); \
}

#define __EMU_M256_IMPL_M2I_SHIFT( type, func, shift_for_hi ) \
static __emu_inline  __emu##type __emu_mm256_##func( __emu##type m256_param1, __emu##type m256_param2, const int param3 ) \
{   __emu##type res; \
    res.__emu_m128[0] = _mm_##func( m256_param1.__emu_m128[0], m256_param2.__emu_m128[0], param3 & ((1<<shift_for_hi)-1) ); \
    res.__emu_m128[1] = _mm_##func( m256_param1.__emu_m128[1], m256_param2.__emu_m128[1], param3 >> shift_for_hi ); \
    return ( res ); \
}

#define __EMU_M256_IMPL_M3( type, func ) \
static __emu_inline __emu##type __emu_mm256_##func( __emu##type m256_param1, __emu##type m256_param2, __emu##type m256_param3 ) \
{   __emu##type res; \
    res.__emu_m128[0] = _mm_##func( m256_param1.__emu_m128[0], m256_param2.__emu_m128[0], m256_param3.__emu_m128[0] ); \
    res.__emu_m128[1] = _mm_##func( m256_param1.__emu_m128[1], m256_param2.__emu_m128[1], m256_param3.__emu_m128[1] ); \
    return ( res ); \
}


/*
 * Compare predicates for scalar and packed compare intrinsics
 */
#define _CMP_EQ_OQ     0x00  /* Equal (ordered, nonsignaling)                       */
#define _CMP_LT_OS     0x01  /* Less-than (ordered, signaling)                      */
#define _CMP_LE_OS     0x02  /* Less-than-or-equal (ordered, signaling)             */
#define _CMP_UNORD_Q   0x03  /* Unordered (nonsignaling)                            */
#define _CMP_NEQ_UQ    0x04  /* Not-equal (unordered, nonsignaling)                 */
#define _CMP_NLT_US    0x05  /* Not-less-than (unordered, signaling)                */
#define _CMP_NLE_US    0x06  /* Not-less-than-or-equal (unordered, signaling)       */
#define _CMP_ORD_Q     0x07  /* Ordered (nonsignaling)                              */

#define _CMP_EQ_UQ     _CMP_EQ_OQ  /* Equal (unordered, non-signaling)                    */
#define _CMP_NGE_US    _CMP_LT_OS  /* Not-greater-than-or-equal (unordered, signaling)    */
#define _CMP_NGT_US    _CMP_LE_OS  /* Not-greater-than (unordered, signaling)             */
#define _CMP_FALSE_OQ  0x0B  /* False (ordered, nonsignaling)                       */
#define _CMP_NEQ_OQ    _CMP_NEQ_UQ  /* Not-equal (ordered, non-signaling)                  */
#define _CMP_GE_OS     _CMP_NLT_US  /* Greater-than-or-equal (ordered, signaling)          */
#define _CMP_GT_OS     _CMP_NLE_US  /* Greater-than (ordered, signaling)                   */
#define _CMP_TRUE_UQ   0x0F  /* True (unordered, non-signaling)                     */

#define _CMP_EQ_OS     _CMP_EQ_OQ  /* Equal (ordered, signaling)                          */
#define _CMP_LT_OQ     _CMP_LT_OS  /* Less-than (ordered, nonsignaling)                   */
#define _CMP_LE_OQ     _CMP_LE_OS  /* Less-than-or-equal (ordered, nonsignaling)          */
#define _CMP_UNORD_S   _CMP_UNORD_Q  /* Unordered (signaling)                               */
#define _CMP_NEQ_US    _CMP_NEQ_UQ  /* Not-equal (unordered, signaling)                    */
#define _CMP_NLT_UQ    _CMP_NLT_US  /* Not-less-than (unordered, nonsignaling)             */
#define _CMP_NLE_UQ    _CMP_NLE_US  /* Not-less-than-or-equal (unordered, nonsignaling)    */
#define _CMP_ORD_S     _CMP_ORD_Q  /* Ordered (signaling)                                 */

#define _CMP_EQ_US     _CMP_EQ_OQ  /* Equal (unordered, signaling)                        */
#define _CMP_NGE_UQ    _CMP_LT_OS  /* Not-greater-than-or-equal (unordered, nonsignaling) */
#define _CMP_NGT_UQ    _CMP_LE_OS  /* Not-greater-than (unordered, nonsignaling)          */
#define _CMP_FALSE_OS  0x1B  /* False (ordered, signaling)                          */
#define _CMP_NEQ_OS    _CMP_NEQ_UQ  /* Not-equal (ordered, signaling)                      */
#define _CMP_GE_OQ     _CMP_NLT_US  /* Greater-than-or-equal (ordered, nonsignaling)       */
#define _CMP_GT_OQ     _CMP_NLE_US  /* Greater-than (ordered, nonsignaling)                */
#define _CMP_TRUE_US   0x1F  /* True (unordered, signaling)                         */

__EMU_M256_IMPL_M2( __m256d, add_pd )
__EMU_M256_IMPL_M2( __m256, add_ps )

__EMU_M256_IMPL_M2( __m256d, addsub_pd )
__EMU_M256_IMPL_M2( __m256, addsub_ps  )

__EMU_M256_IMPL_M2( __m256d, and_pd  )
__EMU_M256_IMPL_M2( __m256, and_ps )

__EMU_M256_IMPL_M2( __m256d, andnot_pd )
__EMU_M256_IMPL_M2( __m256, andnot_ps )

__EMU_M256_IMPL_M2( __m256d, div_pd )
__EMU_M256_IMPL_M2( __m256, div_ps )

__EMU_M256_IMPL_M2( __m256d, hadd_pd )
__EMU_M256_IMPL_M2( __m256, hadd_ps )

__EMU_M256_IMPL_M2( __m256d, hsub_pd )
__EMU_M256_IMPL_M2( __m256, hsub_ps )

__EMU_M256_IMPL_M2( __m256d, max_pd )
__EMU_M256_IMPL_M2( __m256, max_ps )

__EMU_M256_IMPL_M2( __m256d, min_pd )
__EMU_M256_IMPL_M2( __m256, min_ps )

__EMU_M256_IMPL_M2( __m256d, mul_pd )
__EMU_M256_IMPL_M2( __m256, mul_ps )

__EMU_M256_IMPL_M2( __m256d, or_pd )
__EMU_M256_IMPL_M2( __m256, or_ps )

__EMU_M256_IMPL_M2I_SHIFT( __m256d, shuffle_pd, 2 )
__EMU_M256_IMPL_M2I_DUP( __m256, shuffle_ps )

__EMU_M256_IMPL_M2( __m256d, sub_pd )
__EMU_M256_IMPL_M2( __m256, sub_ps )

__EMU_M256_IMPL_M2( __m256d, xor_pd )
__EMU_M256_IMPL_M2( __m256, xor_ps )

#if defined (__SSE4_2__) || defined (__SSE4_1__) || defined(__EMBREE_SMMINTRIN_EMU_H__)

__EMU_M256_IMPL_M2I_SHIFT( __m256d, blend_pd, 2 )
__EMU_M256_IMPL_M2I_SHIFT( __m256, blend_ps, 4 )

__EMU_M256_IMPL_M3( __m256d, blendv_pd )
__EMU_M256_IMPL_M3( __m256, blendv_ps )

__EMU_M256_IMPL_M2I_DUP( __m256, dp_ps )

__EMU_M256_IMPL_M1I_DUP( __m256d, round_pd )
#define _mm256_ceil_pd(val)   _mm256_round_pd((val), 0x0A);
#define _mm256_floor_pd(val)  _mm256_round_pd((val), 0x09);

__EMU_M256_IMPL_M1I_DUP( __m256, round_ps )
#define _mm256_ceil_ps(val)   _mm256_round_ps((val), 0x0A);
#define _mm256_floor_ps(val)  _mm256_round_ps((val), 0x09);

#define __emu_mm_test_impl_ps( op, sfx, vec_type ) \
static __emu_inline int     __emu_mm_test##op##_##sfx(vec_type s1, vec_type s2) {               \
    __m128  sign_bits_ps = _mm_castsi128_ps( _mm_set1_epi32( 1 << 31 ) );                       \
    s1 = _mm_and_##sfx( s1,  sign_bits_##sfx );                                                 \
    s2 = _mm_and_##sfx( s2,  sign_bits_##sfx );                                                 \
    return _mm_test##op##_si128( _mm_cast##sfx##_si128( s1 ), _mm_cast##sfx##_si128( s2 ) );        \
}

#define __emu_mm_test_impl_pd( op, sfx, vec_type ) \
static __emu_inline int     __emu_mm_test##op##_##sfx(vec_type s1, vec_type s2) {               \
    __m128d sign_bits_pd = _mm_castsi128_pd( _mm_set_epi32( 1 << 31, 0, 1 << 31, 0 ) );         \
    s1 = _mm_and_##sfx( s1,  sign_bits_##sfx );                                                 \
    s2 = _mm_and_##sfx( s2,  sign_bits_##sfx );                                                 \
    return _mm_test##op##_si128( _mm_cast##sfx##_si128( s1 ), _mm_cast##sfx##_si128( s2 ) );        \
}

__emu_mm_test_impl_pd( z, pd, __m128d )
__emu_mm_test_impl_pd( c, pd, __m128d )
__emu_mm_test_impl_pd( nzc, pd, __m128d )

__emu_mm_test_impl_ps( z, ps, __m128 )
__emu_mm_test_impl_ps( c, ps, __m128 )
__emu_mm_test_impl_ps( nzc, ps, __m128 )



#define __emu_mm256_test_impl( prfx, op, sfx, sfx_impl, vec_type ) \
static __emu_inline int     __emu_mm256_test##op##_##sfx(vec_type s1, vec_type s2) { \
    int ret1 = prfx##_test##op##_##sfx_impl( s1.__emu_m128[0], s2.__emu_m128[0] );           \
    int ret2 = prfx##_test##op##_##sfx_impl( s1.__emu_m128[1], s2.__emu_m128[1] );           \
    return ( ret1 && ret2 );                                                         \
}

__emu_mm256_test_impl( _mm, z,   si256, si128, __emu__m256i )
__emu_mm256_test_impl( _mm, c,   si256, si128, __emu__m256i )
__emu_mm256_test_impl( _mm, nzc, si256, si128, __emu__m256i )

__emu_mm256_test_impl( __emu_mm, z,   pd, pd, __emu__m256d )
__emu_mm256_test_impl( __emu_mm, c,   pd, pd, __emu__m256d )
__emu_mm256_test_impl( __emu_mm, nzc, pd, pd, __emu__m256d )

__emu_mm256_test_impl( __emu_mm, z,   ps, ps, __emu__m256 )
__emu_mm256_test_impl( __emu_mm, c,   ps, ps, __emu__m256 )
__emu_mm256_test_impl( __emu_mm, nzc, ps, ps, __emu__m256 )

#endif

#if 0 //defined( __GNUC__ ) && ( __GNUC__ == 4 ) && (__GNUC_MINOR__ < 4 )
/* use macro implementation instead of inline functions to allow -O0 for GCC pre 4.4 */

#pragma message ("Using macro for GCC <4.4" )

#define __emu_mm_cmp_ps(m1, m2, predicate) \
({ \
    __m128 res_ = (m1), m2_ = (m2); \
    if ( 7 < (unsigned)predicate ) __asm__ __volatile__ ( "ud2" : : : "memory" ); \
    __asm__ ( "cmpps %[pred_], %[m2_], %[res_]" : [res_] "+x" (res_) : [m2_] "xm" (m2_), [pred_] "i" (predicate) ); \
    res_; })

#define __emu_mm256_cmp_ps(m1, m2, predicate) \
({                                     \
    __emu__m256 res_ = (m1), m2_ = (m2); \
    if ( 7 < (unsigned)predicate ) __asm__ __volatile__ ( "ud2" : : : "memory" ); /* not supported yet */ \
    __asm__ ( "cmpps %[pred_], %[m2_], %[res_]" : [res_] "+x" (res_.__emu_m128[0]) : [m2_] "xm" (m2_.__emu_m128[0]), [pred_] "i" (predicate) ); \
    __asm__ ( "cmpps %[pred_], %[m2_], %[res_]" : [res_] "+x" (res_.__emu_m128[1]) : [m2_] "xm" (m2_.__emu_m128[1]), [pred_] "i" (predicate) ); \
    res_; })

#define __emu_mm_cmp_pd(m1, m2, predicate) \
({ \
    __m128 res_ = (m1), m2_ = (m2); \
    if ( 7 < (unsigned)predicate ) __asm__ __volatile__ ( "ud2" : : : "memory" ); /* not supported yet */ \
    __asm__ ( "cmppd %[pred_], %[m2_], %[res_]" : [res_] "+x" (res_) : [m2_] "xm" (m2_), [pred_] "i" (predicate) ); \
    res_; })

#define __emu_mm256_cmp_pd(m1, m2, predicate) \
({ \
    __emu__m256 res_ = (m1), m2_ = (m2); \
    if ( 7 < (unsigned)predicate ) __asm__ __volatile__ ( "ud2" : : : "memory" ); /* not supported yet */ \
    __asm__ ( "cmppd %[pred_], %[m2_], %[res_]" : [res_] "+x" (res_.__emu_m128[0]) : [m2_] "xm" (m2_.__emu_m128[0]), [pred_] "i" (predicate) ); \
    __asm__ ( "cmppd %[pred_], %[m2_], %[res_]" : [res_] "+x" (res_.__emu_m128[1]) : [m2_] "xm" (m2_.__emu_m128[1]), [pred_] "i" (predicate) ); \
    res_; })


#define __emu_mm_cmp_ss(m1, m2, predicate) \
({ \
    __m128 res_ = (m1), m2_ = (m2); \
    if ( 7 < (unsigned)predicate ) __asm__ __volatile__ ( "ud2" : : : "memory" ); /* not supported yet */ \
    __asm__ ( "cmpss %[pred_], %[m2_], %[res_]" : [res_] "+x" (res_) : [m2_] "xm" (m2_), [pred_] "i" (predicate) ); \
    res_; })

#define __emu_mm_cmp_sd(m1, m2, predicate) \
({ \
    __m128 res_ = (m1), m2_ = (m2); \
    if ( 7 < (unsigned)predicate ) __asm__ __volatile__ ( "ud2" : : : "memory" ); /* not supported yet */ \
    __asm__ ( "cmpsd %[pred_], %[m2_], %[res_]" : [res_] "+x" (res_) : [m2_] "xm" (m2_), [pred_] "i" (predicate) ); \
    res_; })



#else /* __GNUC__==4 && __GNUC_MINOR__ <4 */


static __emu_inline __m128 __emu_mm_cmp_ps(__m128 m1, __m128 m2, const int predicate)
{
    __m128 res;

    if ( predicate >= 0 && predicate <= 7 ) {
        res = m1;
        __asm__ ( "cmpps %[pred_], %[m2_], %[res_]" : [res_] "+x" (res) : [m2_] "xm" (m2), [pred_] "i" (predicate) );
    } else {
        __asm__ __volatile__ ( "ud2" : : : "memory" ); /* not supported yet */
    }

    return ( res );
}
__EMU_M256_IMPL2_M2I_DUP( __m256, cmp_ps )

static __emu_inline __m128d __emu_mm_cmp_pd(__m128d m1, __m128d m2, const int predicate)
{
    __m128d res;

    if ( predicate >= 0 && predicate <= 7 ) {
        res = m1;
        __asm__ ( "cmppd %[pred_], %[m2_], %[res_]" : [res_] "+x" (res) : [m2_] "xm" (m2), [pred_] "i" (predicate) );
    } else {
        __asm__ __volatile__ ( "ud2" : : : "memory" ); /* not supported yet */
    }

    return ( res );
}
__EMU_M256_IMPL2_M2I_DUP( __m256d, cmp_pd )


static __emu_inline __m128d __emu_mm_cmp_sd(__m128d m1, __m128d m2, const int predicate)
{
    __m128d res;

    if ( predicate >= 0 && predicate <= 7 ) {
        res = m1;
        __asm__ ( "cmpsd %[pred_], %[m2_], %[res_]" : [res_] "+x" (res) : [m2_] "xm" (m2), [pred_] "i" (predicate) );
    } else {
        __asm__ __volatile__ ( "ud2" : : : "memory" ); /* not supported yet */
    }

    return ( res );
}

static __emu_inline __m128 __emu_mm_cmp_ss(__m128 m1, __m128 m2, const int predicate)
{
    __m128 res;

    if ( predicate >= 0 && predicate <= 7 ) {
        res = m1;
        __asm__ ( "cmpss %[pred_], %[m2_], %[res_]" : [res_] "+x" (res) : [m2_] "xm" (m2), [pred_] "i" (predicate) );
    } else {
        __asm__ __volatile__ ( "ud2" : : : "memory" ); /* not supported yet */
    }

    return ( res );
}

#endif


__EMU_M256_IMPL_M1_LH( __m256d, __m128i, cvtepi32_pd )
__EMU_M256_IMPL_M1_RET( __m256, __m256i, cvtepi32_ps )
//__EMU_M256_IMPL_M1_HL( __m128, __m256d, cvtpd_ps )
__EMU_M256_IMPL_M1_RET( __m256i, __m256, cvtps_epi32 )
__EMU_M256_IMPL_M1_LH( __m256d, __m128, cvtps_pd )
//__EMU_M256_IMPL_M1_HL( __m128i, __m256d, cvttpd_epi32 )
//__EMU_M256_IMPL_M1_HL( __m128i, __m256d, cvtpd_epi32 )
__EMU_M256_IMPL_M1_RET( __m256i, __m256, cvttps_epi32 )

static __emu_inline __m128  __emu_mm256_extractf128_ps(__emu__m256 m1, const int offset) { return m1.__emu_m128[ offset ]; }
static __emu_inline __m128d __emu_mm256_extractf128_pd(__emu__m256d m1, const int offset) { return m1.__emu_m128[ offset ]; }
static __emu_inline __m128i __emu_mm256_extractf128_si256(__emu__m256i m1, const int offset) { return m1.__emu_m128[ offset ]; }

static __emu_inline void __emu_mm256_zeroall(void) {}
static __emu_inline void __emu_mm256_zeroupper(void) {}

static __emu_inline __m128  __emu_mm_permutevar_ps(__m128 a, __m128i control)
{
    int const* sel = (int const*)&control;
    float const* src = (float const*)&a;
    __EMU_M256_ALIGN(16) float dest[4];
    int i=0;

    for (; i<4; ++i)
        dest[i] = src[ 3 & sel[i] ];

    return ( *(__m128*)dest );
}
__EMU_M256_IMPL2_M2T( __m256, __m256i, permutevar_ps )

static __emu_inline __m128  __emu_mm_permute_ps(__m128 a, int control) { return _mm_castsi128_ps( _mm_shuffle_epi32( *(__m128i*)&a, control ) ); }
__EMU_M256_IMPL2_M1I_DUP( __m256, permute_ps )


static __emu_inline __m128d __emu_mm_permutevar_pd(__m128d a, __m128i control)
{
    __emu_int64_t const* sel = (__emu_int64_t const*)&control;
    double const* src = (double const*)&a;
    __EMU_M256_ALIGN(16) double dest[2];
    int i=0;

    for (; i<2; ++i)
        dest[i] = src[ (2 & sel[i]) >> 1 ];

    return ( *(__m128d*)dest );
}
__EMU_M256_IMPL2_M2T( __m256d, __m256i, permutevar_pd )

static __emu_inline __m128d __emu_mm_permute_pd(__m128d a, int control)
{
    double const* src = (double const*)&a;
    __EMU_M256_ALIGN(16) double dest[2];
    int i=0;

    for (; i<2; ++i)
        dest[i] = src[ 1 & (control >> i) ];

    return ( *(__m128d*)dest );
}
__EMU_M256_IMPL2_M1I_SHIFT( __m256d, permute_pd, 2 )


#define __emu_mm256_permute2f128_impl( name, m128_type, m256_type ) \
static __emu_inline m256_type name( m256_type m1, m256_type m2, int control) { \
    m256_type res; \
    __m128 zero = _mm_setzero_ps(); \
    const m128_type param[4] = { m1.__emu_m128[0], m1.__emu_m128[1], m2.__emu_m128[0], m2.__emu_m128[1] }; \
    res.__emu_m128[0] = (control & 8) ? *(m128_type*)&zero : param[ control & 0x3 ]; control >>= 4; \
    res.__emu_m128[1] = (control & 8) ? *(m128_type*)&zero : param[ control & 0x3 ]; \
    return ( res ); \
}

__emu_mm256_permute2f128_impl( __emu_mm256_permute2f128_ps, __m128, __emu__m256 )
__emu_mm256_permute2f128_impl( __emu_mm256_permute2f128_pd, __m128d, __emu__m256d )
__emu_mm256_permute2f128_impl( __emu_mm256_permute2f128_si256, __m128i, __emu__m256i )


#define __emu_mm_broadcast_impl( name, res_type, type )     \
static __emu_inline res_type  name(type const *a) {         \
    const size_t size = sizeof( res_type ) / sizeof( type );\
    __EMU_M256_ALIGN(16) type res[ size ];                  \
    size_t i = 0;                                           \
    for ( ; i < size; ++i )                                 \
        res[ i ] = *a;                                      \
    return (*(res_type*)&res);                              \
}

__emu_mm_broadcast_impl( __emu_mm_broadcast_ss, __m128, float )
__emu_mm_broadcast_impl( __emu_mm256_broadcast_ss, __emu__m256, float )

__emu_mm_broadcast_impl( __emu_mm_broadcast_sd, __m128, double )
__emu_mm_broadcast_impl( __emu_mm256_broadcast_sd, __emu__m256d, double )

__emu_mm_broadcast_impl( __emu_mm256_broadcast_ps, __emu__m256, __m128 )
__emu_mm_broadcast_impl( __emu_mm256_broadcast_pd, __emu__m256d, __m128d )


static __emu_inline __emu__m256  __emu_mm256_insertf128_ps(__emu__m256 a, __m128 b, int offset)  { a.__emu_m128[ offset ] = b; return a; }
static __emu_inline __emu__m256d __emu_mm256_insertf128_pd(__emu__m256d a, __m128d b, int offset)  { a.__emu_m128[ offset ] = b; return a; }
static __emu_inline __emu__m256i __emu_mm256_insertf128_si256(__emu__m256i a, __m128i b, int offset)  { a.__emu_m128[ offset ] = b; return a; }


#define __emu_mm_load_impl( name, sfx, m256_sfx, m256_type, type_128, type )           \
static __emu_inline __emu##m256_type  __emu_mm256_##name##_##m256_sfx(const type* a) { \
    __emu##m256_type res;                                                              \
    res.__emu_m128[0] = _mm_##name##_##sfx( (const type_128 *)a );                     \
    res.__emu_m128[1] = _mm_##name##_##sfx( (const type_128 *)(1+(const __m128 *)a) ); \
    return (res);                                                                      \
}

#define __emu_mm_store_impl( name, sfx, m256_sfx, m256_type, type_128, type ) \
static __emu_inline void __emu_mm256_##name##_##m256_sfx(type *a, __emu##m256_type b) {  \
    _mm_##name##_##sfx( (type_128*)a, b.__emu_m128[0] );                                 \
    _mm_##name##_##sfx( (type_128*)(1+(__m128*)a), b.__emu_m128[1] );   \
}

__emu_mm_load_impl( load, pd, pd, __m256d, double, double )
__emu_mm_store_impl( store, pd, pd, __m256d, double, double )

__emu_mm_load_impl( load, ps, ps, __m256, float, float )
__emu_mm_store_impl( store, ps, ps, __m256, float, float )

__emu_mm_load_impl( loadu, pd, pd, __m256d, double, double )
__emu_mm_store_impl( storeu, pd, pd, __m256d, double, double )

__emu_mm_load_impl( loadu, ps, ps, __m256, float, float )
__emu_mm_store_impl( storeu, ps, ps, __m256, float, float )

__emu_mm_load_impl( load, si128, si256, __m256i, __m128i, __emu__m256i )
__emu_mm_store_impl( store, si128, si256, __m256i, __m128i, __emu__m256i )

__emu_mm_load_impl( loadu, si128, si256, __m256i, __m128i, __emu__m256i )
__emu_mm_store_impl( storeu, si128, si256, __m256i, __m128i, __emu__m256i )


#define __emu_maskload_impl( name, vec_type, mask_vec_type, type, mask_type ) \
static __emu_inline vec_type  name(type const *a, mask_vec_type mask) {   \
    const size_t size_type = sizeof( type );                          \
    const size_t size = sizeof( vec_type ) / size_type;               \
    __EMU_M256_ALIGN(16) type res[ size ];                            \
    const mask_type* p_mask = (const mask_type*)&mask;                \
    size_t i = 0;                                                     \
    mask_type sign_bit = 1;                                           \
    sign_bit <<= (8*size_type - 1);                                   \
    for ( ; i < size; ++i )                                           \
        res[ i ] = (sign_bit & *(p_mask + i)) ? *(a+i) : 0;           \
    return (*(vec_type*)&res);                                        \
}

#define __emu_maskstore_impl( name, vec_type, mask_vec_type, type, mask_type ) \
static __emu_inline void  name(type *a, mask_vec_type mask, vec_type data) { \
    const size_t size_type = sizeof( type );                          \
    const size_t size = sizeof( vec_type ) / sizeof( type );          \
    type* p_data = (type*)&data;                                      \
    const mask_type* p_mask = (const mask_type*)&mask;                \
    size_t i = 0;                                                     \
    mask_type sign_bit = 1;                                           \
    sign_bit <<= (8*size_type - 1);                                   \
    for ( ; i < size; ++i )                                           \
        if ( *(p_mask + i ) & sign_bit)                               \
            *(a + i) = *(p_data + i);                                 \
}

__emu_maskload_impl( __emu_mm256_maskload_pd, __emu__m256d, __emu__m256i, double, __emu_int64_t )
__emu_maskstore_impl( __emu_mm256_maskstore_pd, __emu__m256d, __emu__m256i, double, __emu_int64_t )

__emu_maskload_impl( __emu_mm_maskload_pd, __m128d, __m128i, double, __emu_int64_t )
__emu_maskstore_impl( __emu_mm_maskstore_pd, __m128d, __m128i, double, __emu_int64_t )

__emu_maskload_impl( __emu_mm256_maskload_ps, __emu__m256, __emu__m256i, float, int )
__emu_maskstore_impl( __emu_mm256_maskstore_ps, __emu__m256, __emu__m256i, float, int )

__emu_maskload_impl( __emu_mm_maskload_ps, __m128, __m128i, float, int )
__emu_maskstore_impl( __emu_mm_maskstore_ps, __m128, __m128i, float, int )


__EMU_M256_IMPL_M1( __m256, movehdup_ps )
__EMU_M256_IMPL_M1( __m256, moveldup_ps )
__EMU_M256_IMPL_M1( __m256d, movedup_pd )

__emu_mm_load_impl( lddqu, si128, si256, __m256i, __m128i, __emu__m256i )

__emu_mm_store_impl( stream, si128, si256, __m256i, __m128i, __emu__m256i )
__emu_mm_store_impl( stream, pd, pd, __m256d, double, double )
__emu_mm_store_impl( stream, ps, ps, __m256, float, float )


__EMU_M256_IMPL_M1( __m256, rcp_ps )
__EMU_M256_IMPL_M1( __m256, rsqrt_ps )

__EMU_M256_IMPL_M1( __m256d, sqrt_pd )
__EMU_M256_IMPL_M1( __m256, sqrt_ps )

__EMU_M256_IMPL_M2( __m256d, unpackhi_pd )
__EMU_M256_IMPL_M2( __m256, unpackhi_ps )
__EMU_M256_IMPL_M2( __m256d, unpacklo_pd )
__EMU_M256_IMPL_M2( __m256, unpacklo_ps )


static __emu_inline int     __emu_mm256_movemask_pd(__emu__m256d a)
{
    return
        (_mm_movemask_pd( a.__emu_m128[1] ) << 2) |
        _mm_movemask_pd( a.__emu_m128[0] );
}

static __emu_inline int     __emu_mm256_movemask_ps(__emu__m256 a)
{
    return
        (_mm_movemask_ps( a.__emu_m128[1] ) << 4) |
         _mm_movemask_ps( a.__emu_m128[0] );
}

static __emu_inline __emu__m256d __emu_mm256_setzero_pd(void) { __m128d ret[2] = { _mm_setzero_pd(), _mm_setzero_pd() }; return __emu_set_m128d( ret ); }
static __emu_inline __emu__m256  __emu_mm256_setzero_ps(void) { __m128  ret[2] = { _mm_setzero_ps(), _mm_setzero_ps() }; return __emu_set_m128( ret ); }
static __emu_inline __emu__m256i __emu_mm256_setzero_si256(void) { __m128i ret[2] = { _mm_setzero_si128(), _mm_setzero_si128() }; return __emu_set_m128i( ret ); }

static __emu_inline __emu__m256d __emu_mm256_set_pd(double a1, double a2, double a3, double a4)
{ __m128d ret[2] = { _mm_set_pd( a3, a4 ), _mm_set_pd( a1, a2 ) }; return __emu_set_m128d( ret ); }

static __emu_inline __emu__m256  __emu_mm256_set_ps(float a1, float a2, float a3, float a4, float a5, float a6, float a7, float a8)
{ __m128 ret[2] = { _mm_set_ps( a5, a6, a7, a8 ), _mm_set_ps( a1, a2, a3, a4 ) }; return __emu_set_m128( ret ); }

static __emu_inline __emu__m256i __emu_mm256_set_epi8(char a1, char a2, char a3, char a4, char a5, char a6, char a7, char a8,
                                       char a9, char a10, char a11, char a12, char a13, char a14, char a15, char a16,
                                       char a17, char a18, char a19, char a20, char a21, char a22, char a23, char a24,
                                       char a25, char a26, char a27, char a28, char a29, char a30, char a31, char a32)
{   __m128i ret[2] = { _mm_set_epi8( a17, a18, a19, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29, a30, a31, a32 ),
                       _mm_set_epi8( a1,   a2,  a3,  a4,  a5,  a6,  a7,  a8,  a9, a10, a11, a12, a13, a14, a15, a16 ) };
    return __emu_set_m128i( ret );
}

static __emu_inline __emu__m256i __emu_mm256_set_epi16(short a1, short a2, short a3, short a4, short a5, short a6, short a7, short a8,
                                                       short a9, short a10, short a11, short a12, short a13, short a14, short a15, short a16)
{   __m128i ret[2] = { _mm_set_epi16( a9, a10, a11, a12, a13, a14, a15, a16 ),
                       _mm_set_epi16( a1,  a2,  a3,  a4,  a5,  a6,  a7,  a8 ) };
    return __emu_set_m128i( ret );
}

static __emu_inline __emu__m256i __emu_mm256_set_epi32(int a1, int a2, int a3, int a4, int a5, int a6, int a7, int a8)
{ __m128i ret[2] = { _mm_set_epi32( a5, a6, a7, a8 ), _mm_set_epi32( a1, a2, a3, a4 ) }; return __emu_set_m128i( ret ); }

static __emu_inline __m128i __emu_mm_set_epi64x( __emu_int64_t a, __emu_int64_t b ) { return _mm_set_epi64( *(__m64*)&a, *(__m64*)&b ); }

static __emu_inline __emu__m256i __emu_mm256_set_epi64x(__emu_int64_t a1, __emu_int64_t a2, __emu_int64_t a3, __emu_int64_t a4)
{ __m128i ret[2] = { __emu_mm_set_epi64x( a3, a4 ), __emu_mm_set_epi64x( a1, a2 ) }; return __emu_set_m128i( ret ); }


static __emu_inline __emu__m256d __emu_mm256_setr_pd(double a1, double a2, double a3, double a4)
{ __m128d ret[2] = { _mm_setr_pd( a1, a2 ), _mm_setr_pd( a3, a4 ) }; return __emu_set_m128d( ret ); }

static __emu_inline __emu__m256  __emu_mm256_setr_ps(float a1, float a2, float a3, float a4, float a5, float a6, float a7, float a8)
{ __m128 ret[2] = { _mm_setr_ps( a1, a2, a3, a4 ), _mm_setr_ps( a5, a6, a7, a8 ) }; return __emu_set_m128( ret ); }

static __emu_inline __emu__m256i __emu_mm256_setr_epi8(char a1, char a2, char a3, char a4, char a5, char a6, char a7, char a8,
                                                      char a9, char a10, char a11, char a12, char a13, char a14, char a15, char a16,
                                                      char a17, char a18, char a19, char a20, char a21, char a22, char a23, char a24,
                                                      char a25, char a26, char a27, char a28, char a29, char a30, char a31, char a32)
{   __m128i ret[2] = { _mm_setr_epi8( a1,   a2,  a3,  a4,  a5,  a6,  a7,  a8,  a9, a10, a11, a12, a13, a14, a15, a16 ),
                       _mm_setr_epi8( a17, a18, a19, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29, a30, a31, a32 ) };
    return __emu_set_m128i( ret );
}

static __emu_inline __emu__m256i __emu_mm256_setr_epi16(short a1, short a2, short a3, short a4, short a5, short a6, short a7, short a8,
                                                       short a9, short a10, short a11, short a12, short a13, short a14, short a15, short a16)
{   __m128i ret[2] = { _mm_setr_epi16( a1,  a2,  a3,  a4,  a5,  a6,  a7,  a8 ),
                       _mm_setr_epi16( a9, a10, a11, a12, a13, a14, a15, a16 ) }; return __emu_set_m128i( ret );
}

static __emu_inline __emu__m256i __emu_mm256_setr_epi32(int a1, int a2, int a3, int a4, int a5, int a6, int a7, int a8)
{   __m128i ret[2] = { _mm_setr_epi32( a1, a2, a3, a4 ), _mm_setr_epi32( a5, a6, a7, a8 ),  }; return __emu_set_m128i( ret ); }

static __emu_inline __emu__m256i __emu_mm256_setr_epi64x(__emu_int64_t a1, __emu_int64_t a2, __emu_int64_t a3, __emu_int64_t a4)
{   __m128i ret[2] = { __emu_mm_set_epi64x( a2, a1 ), __emu_mm_set_epi64x( a4, a3 ) }; return __emu_set_m128i( ret ); }



__EMU_M256_IMPL_M1P_DUP( __m256d, double, set1_pd )
__EMU_M256_IMPL_M1P_DUP( __m256, float, set1_ps )
__EMU_M256_IMPL_M1P_DUP( __m256i, char, set1_epi8 )
__EMU_M256_IMPL_M1P_DUP( __m256i, short, set1_epi16 )
__EMU_M256_IMPL_M1P_DUP( __m256i, int, set1_epi32 )

inline __emu__m256i __emu_mm256_set1_epi64x(__emu_int64_t a)
{
    __emu_int64_t res[4] = { a, a, a, a };
    return *((__emu__m256i*)res);
}

/*
 * Support intrinsics to do vector type casts. These intrinsics do not introduce
 * extra moves to generated code. When cast is done from a 128 to 256-bit type
 * the low 128 bits of the 256-bit result contain source parameter value; the
 * upper 128 bits of the result are undefined
 */
__EMU_M256_IMPL_M1_RET( __m256, __m256d, castpd_ps )
__EMU_M256_IMPL_M1_RET( __m256d, __m256, castps_pd )

__EMU_M256_IMPL_M1_RET_NAME( __m256i, __m256, castps_si128, castps_si256 )
__EMU_M256_IMPL_M1_RET_NAME( __m256i, __m256d, castpd_si128, castpd_si256 )

__EMU_M256_IMPL_M1_RET_NAME( __m256, __m256i, castsi128_ps, castsi256_ps )
__EMU_M256_IMPL_M1_RET_NAME( __m256d, __m256i, castsi128_pd, castsi256_pd )

static __emu_inline __m128  __emu_mm256_castps256_ps128(__emu__m256 a) { return ( a.__emu_m128[0] ); }
static __emu_inline __m128d __emu_mm256_castpd256_pd128(__emu__m256d a) { return ( a.__emu_m128[0] ); }
static __emu_inline __m128i __emu_mm256_castsi256_si128(__emu__m256i a) { return ( a.__emu_m128[0] ); }

static __emu_inline __emu__m256  __emu_mm256_castps128_ps256(__m128 a) { __m128 ret[2] = { a, _mm_setzero_ps() }; return __emu_set_m128( ret ); }
static __emu_inline __emu__m256d __emu_mm256_castpd128_pd256(__m128d a) { __m128d ret[2] = { a, _mm_setzero_pd() }; return __emu_set_m128d( ret ); }
static __emu_inline __emu__m256i __emu_mm256_castsi128_si256(__m128i a) { __m128i ret[2] = { a, _mm_setzero_si128() }; return __emu_set_m128i( ret ); }

#if defined __cplusplus
} /* End "C" */
#endif /* __cplusplus */






#ifndef __EMU_M256_NOMAP

#define __m256  __emu__m256
#define __m256i __emu__m256i
#define __m256d __emu__m256d

#define _mm256_add_pd    __emu_mm256_add_pd
#define _mm256_add_ps    __emu_mm256_add_ps

#define _mm256_addsub_pd __emu_mm256_addsub_pd
#define _mm256_addsub_ps __emu_mm256_addsub_ps

#define _mm256_and_pd __emu_mm256_and_pd
#define _mm256_and_ps __emu_mm256_and_ps

#define _mm256_andnot_pd __emu_mm256_andnot_pd
#define _mm256_andnot_ps __emu_mm256_andnot_ps

#define _mm256_blend_pd __emu_mm256_blend_pd
#define _mm256_blend_ps __emu_mm256_blend_ps

#define _mm256_blendv_pd __emu_mm256_blendv_pd
#define _mm256_blendv_ps __emu_mm256_blendv_ps

#define _mm256_div_pd __emu_mm256_div_pd
#define _mm256_div_ps __emu_mm256_div_ps

#define _mm256_dp_ps __emu_mm256_dp_ps

#define _mm256_hadd_pd __emu_mm256_hadd_pd
#define _mm256_hadd_ps __emu_mm256_hadd_ps

#define _mm256_hsub_pd __emu_mm256_hsub_pd
#define _mm256_hsub_ps __emu_mm256_hsub_ps

#define _mm256_max_pd __emu_mm256_max_pd
#define _mm256_max_ps __emu_mm256_max_ps

#define _mm256_min_pd __emu_mm256_min_pd
#define _mm256_min_ps __emu_mm256_min_ps

#define _mm256_mul_pd __emu_mm256_mul_pd
#define _mm256_mul_ps __emu_mm256_mul_ps

#define _mm256_or_pd __emu_mm256_or_pd
#define _mm256_or_ps __emu_mm256_or_ps

#define _mm256_shuffle_pd __emu_mm256_shuffle_pd
#define _mm256_shuffle_ps __emu_mm256_shuffle_ps

#define _mm256_sub_pd __emu_mm256_sub_pd
#define _mm256_sub_ps __emu_mm256_sub_ps

#define _mm256_xor_pd __emu_mm256_xor_pd
#define _mm256_xor_ps __emu_mm256_xor_ps


#define _mm_cmp_pd __emu_mm_cmp_pd
#define _mm256_cmp_pd __emu_mm256_cmp_pd

#define _mm_cmp_ps __emu_mm_cmp_ps
#define _mm256_cmp_ps __emu_mm256_cmp_ps

#define _mm_cmp_sd __emu_mm_cmp_sd
#define _mm_cmp_ss __emu_mm_cmp_ss

#define _mm256_cvtepi32_pd __emu_mm256_cvtepi32_pd
#define _mm256_cvtepi32_ps __emu_mm256_cvtepi32_ps

#define _mm256_cvtpd_ps __emu_mm256_cvtpd_ps
#define _mm256_cvtps_epi32 __emu_mm256_cvtps_epi32
#define _mm256_cvtps_pd __emu_mm256_cvtps_pd

#define _mm256_cvttpd_epi32 __emu_mm256_cvttpd_epi32
#define _mm256_cvtpd_epi32 __emu_mm256_cvtpd_epi32
#define _mm256_cvttps_epi32 __emu_mm256_cvttps_epi32

#define _mm256_extractf128_ps __emu_mm256_extractf128_ps
#define _mm256_extractf128_pd __emu_mm256_extractf128_pd
#define _mm256_extractf128_si256 __emu_mm256_extractf128_si256

#define _mm256_zeroall __emu_mm256_zeroall
#define _mm256_zeroupper __emu_mm256_zeroupper

#define _mm256_permutevar_ps __emu_mm256_permutevar_ps
#define _mm_permutevar_ps __emu_mm_permutevar_ps

#define _mm256_permute_ps __emu_mm256_permute_ps
#define _mm_permute_ps __emu_mm_permute_ps

#define _mm256_permutevar_pd __emu_mm256_permutevar_pd
#define _mm_permutevar_pd __emu_mm_permutevar_pd

#define _mm256_permute_pd __emu_mm256_permute_pd
#define _mm_permute_pd __emu_mm_permute_pd

#define _mm256_permute2f128_ps __emu_mm256_permute2f128_ps
#define _mm256_permute2f128_pd __emu_mm256_permute2f128_pd
#define _mm256_permute2f128_si256 __emu_mm256_permute2f128_si256

#define _mm256_broadcast_ss __emu_mm256_broadcast_ss
#define _mm_broadcast_ss __emu_mm_broadcast_ss

#define _mm256_broadcast_sd __emu_mm256_broadcast_sd

#define _mm256_broadcast_ps __emu_mm256_broadcast_ps
#define _mm256_broadcast_pd __emu_mm256_broadcast_pd

#define _mm256_insertf128_ps __emu_mm256_insertf128_ps
#define _mm256_insertf128_pd __emu_mm256_insertf128_pd
#define _mm256_insertf128_si256 __emu_mm256_insertf128_si256

#define _mm256_load_pd __emu_mm256_load_pd
#define _mm256_store_pd __emu_mm256_store_pd
#define _mm256_load_ps __emu_mm256_load_ps
#define _mm256_store_ps __emu_mm256_store_ps

#define _mm256_loadu_pd __emu_mm256_loadu_pd
#define _mm256_storeu_pd __emu_mm256_storeu_pd
#define _mm256_loadu_ps __emu_mm256_loadu_ps
#define _mm256_storeu_ps __emu_mm256_storeu_ps

#define _mm256_load_si256 __emu_mm256_load_si256
#define _mm256_store_si256 __emu_mm256_store_si256
#define _mm256_loadu_si256 __emu_mm256_loadu_si256
#define _mm256_storeu_si256 __emu_mm256_storeu_si256

#define _mm256_maskload_pd __emu_mm256_maskload_pd
#define _mm256_maskstore_pd __emu_mm256_maskstore_pd
#define _mm_maskload_pd __emu_mm_maskload_pd
#define _mm_maskstore_pd __emu_mm_maskstore_pd

#define _mm256_maskload_ps __emu_mm256_maskload_ps
#define _mm256_maskstore_ps __emu_mm256_maskstore_ps
#define _mm_maskload_ps __emu_mm_maskload_ps
#define _mm_maskstore_ps __emu_mm_maskstore_ps

#define _mm256_movehdup_ps __emu_mm256_movehdup_ps
#define _mm256_moveldup_ps __emu_mm256_moveldup_ps

#define _mm256_movedup_pd __emu_mm256_movedup_pd
#define _mm256_lddqu_si256 __emu_mm256_lddqu_si256

#define _mm256_stream_si256 __emu_mm256_stream_si256
#define _mm256_stream_pd __emu_mm256_stream_pd
#define _mm256_stream_ps __emu_mm256_stream_ps

#define _mm256_rcp_ps __emu_mm256_rcp_ps
#define _mm256_rsqrt_ps __emu_mm256_rsqrt_ps

#define _mm256_sqrt_pd __emu_mm256_sqrt_pd
#define _mm256_sqrt_ps __emu_mm256_sqrt_ps

#define _mm256_round_pd __emu_mm256_round_pd

#define _mm256_round_ps __emu_mm256_round_ps

#define _mm256_unpackhi_pd __emu_mm256_unpackhi_pd
#define _mm256_unpackhi_ps __emu_mm256_unpackhi_ps

#define _mm256_unpacklo_pd __emu_mm256_unpacklo_pd
#define _mm256_unpacklo_ps __emu_mm256_unpacklo_ps

#define _mm256_testz_si256 __emu_mm256_testz_si256
#define _mm256_testc_si256 __emu_mm256_testc_si256
#define _mm256_testnzc_si256 __emu_mm256_testnzc_si256

#define _mm256_testz_pd __emu_mm256_testz_pd
#define _mm256_testc_pd __emu_mm256_testc_pd
#define _mm256_testnzc_pd __emu_mm256_testnzc_pd
#define _mm_testz_pd __emu_mm_testz_pd
#define _mm_testc_pd __emu_mm_testc_pd
#define _mm_testnzc_pd __emu_mm_testnzc_pd

#define _mm256_testz_ps __emu_mm256_testz_ps
#define _mm256_testc_ps __emu_mm256_testc_ps
#define _mm256_testnzc_ps __emu_mm256_testnzc_ps
#define _mm_testz_ps __emu_mm_testz_ps
#define _mm_testc_ps __emu_mm_testc_ps
#define _mm_testnzc_ps __emu_mm_testnzc_ps

#define _mm256_movemask_pd __emu_mm256_movemask_pd
#define _mm256_movemask_ps __emu_mm256_movemask_ps

#define _mm256_setzero_pd __emu_mm256_setzero_pd
#define _mm256_setzero_ps __emu_mm256_setzero_ps
#define _mm256_setzero_si256 __emu_mm256_setzero_si256

#define _mm256_set_pd __emu_mm256_set_pd
#define _mm256_set_ps __emu_mm256_set_ps
#define _mm256_set_epi8 __emu_mm256_set_epi8
#define _mm256_set_epi16 __emu_mm256_set_epi16
#define _mm256_set_epi32 __emu_mm256_set_epi32
#define _mm256_set_epi64x __emu_mm256_set_epi64x

#define _mm256_setr_pd __emu_mm256_setr_pd
#define _mm256_setr_ps __emu_mm256_setr_ps
#define _mm256_setr_epi8 __emu_mm256_setr_epi8
#define _mm256_setr_epi16 __emu_mm256_setr_epi16
#define _mm256_setr_epi32 __emu_mm256_setr_epi32
#define _mm256_setr_epi64x __emu_mm256_setr_epi64x

#define _mm256_set1_pd __emu_mm256_set1_pd
#define _mm256_set1_ps __emu_mm256_set1_ps
#define _mm256_set1_epi8 __emu_mm256_set1_epi8
#define _mm256_set1_epi16 __emu_mm256_set1_epi16
#define _mm256_set1_epi32 __emu_mm256_set1_epi32
#define _mm256_set1_epi64x __emu_mm256_set1_epi64x

#define _mm256_castpd_ps __emu_mm256_castpd_ps
#define _mm256_castps_pd __emu_mm256_castps_pd
#define _mm256_castps_si256 __emu_mm256_castps_si256
#define _mm256_castpd_si256 __emu_mm256_castpd_si256
#define _mm256_castsi256_ps __emu_mm256_castsi256_ps
#define _mm256_castsi256_pd __emu_mm256_castsi256_pd
#define _mm256_castps256_ps128 __emu_mm256_castps256_ps128
#define _mm256_castpd256_pd128 __emu_mm256_castpd256_pd128
#define _mm256_castsi256_si128 __emu_mm256_castsi256_si128
#define _mm256_castps128_ps256 __emu_mm256_castps128_ps256
#define _mm256_castpd128_pd256 __emu_mm256_castpd128_pd256
#define _mm256_castsi128_si256 __emu_mm256_castsi128_si256

#endif /* __EMU_M256_NOMAP */



#endif /* __EMU_M256_AVXIMMINTRIN_EMU_H__ */
