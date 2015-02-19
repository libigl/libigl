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

namespace embree
{
  /*! swap that is forced to inline */
  template<typename T> __forceinline void swap(T& a, T& b) { T t = b; b = a; a = t; }

  /*! An item on the stack holds the node ID and distance of that node. */
  template<typename T>
  struct StackItemT
  {
    /*! Sort 2 stack items. */
    __forceinline friend void sort(StackItemT& s1, StackItemT& s2) {
      if (s2.dist < s1.dist) swap(s2,s1);
    }
    
    /*! Sort 3 stack items. */
    __forceinline friend void sort(StackItemT& s1, StackItemT& s2, StackItemT& s3)
    {
      if (s2.dist < s1.dist) swap(s2,s1);
      if (s3.dist < s2.dist) swap(s3,s2);
      if (s2.dist < s1.dist) swap(s2,s1);
    }
    
    /*! Sort 4 stack items. */
    __forceinline friend void sort(StackItemT& s1, StackItemT& s2, StackItemT& s3, StackItemT& s4)
    {
      if (s2.dist < s1.dist) swap(s2,s1);
      if (s4.dist < s3.dist) swap(s4,s3);
      if (s3.dist < s1.dist) swap(s3,s1);
      if (s4.dist < s2.dist) swap(s4,s2);
      if (s3.dist < s2.dist) swap(s3,s2);
    }
    
  public:
    T ptr; 
    float dist;
  };

  /*! An item on the stack holds the node ID and distance of that node. */
  template<typename T>
    struct __aligned(16) StackItemInt32
  {
    __forceinline static void swap2(StackItemInt32<T>& a, StackItemInt32<T>& b) 
    { 
      const ssef sse_a = load4f(&a); 
      const ssef sse_b = load4f(&b);
      store4f(&a,sse_b);
      store4f(&b,sse_a);
    }

    /*! Sort 2 stack items. */
    __forceinline friend void sort(StackItemInt32& s1, StackItemInt32& s2) {
      if (s2.dist < s1.dist) swap2(s2,s1);
    }
    
    /*! Sort 3 stack items. */
    __forceinline friend void sort(StackItemInt32& s1, StackItemInt32& s2, StackItemInt32& s3)
    {
      if (s2.dist < s1.dist) swap2(s2,s1);
      if (s3.dist < s2.dist) swap2(s3,s2);
      if (s2.dist < s1.dist) swap2(s2,s1);
    }
    
    /*! Sort 4 stack items. */
    __forceinline friend void sort(StackItemInt32& s1, StackItemInt32& s2, StackItemInt32& s3, StackItemInt32& s4)
    {
      if (s2.dist < s1.dist) swap2(s2,s1);
      if (s4.dist < s3.dist) swap2(s4,s3);
      if (s3.dist < s1.dist) swap2(s3,s1);
      if (s4.dist < s2.dist) swap2(s4,s2);
      if (s3.dist < s2.dist) swap2(s3,s2);
    }
    
  public:
    T ptr; 
    unsigned int dist;
  };


  /*! An item on the stack holds the node ID and distance of that node. */
  template<>
    struct StackItemT<unsigned int>
  {
    /*! Sort 2 stack items. */
    __forceinline friend void sort(StackItemT& a, StackItemT& b)
    {
      int64 s1 = a.all;
      int64 s2 = b.all;
      if (s2 < s1) swap(s2,s1);
      a.all = s1;
      b.all = s2;
    }

    /*! Sort 3 stack items. */
    __forceinline friend void sort(StackItemT& a, StackItemT& b, StackItemT& c)
    {
      int64 s1 = a.all;
      int64 s2 = b.all;
      int64 s3 = c.all;
      if (s2 < s1) swap(s2,s1);
      if (s3 < s2) swap(s3,s2);
      if (s2 < s1) swap(s2,s1);
      a.all = s1;
      b.all = s2;
      c.all = s3;
    }
    
    /*! Sort 4 stack items. */
    __forceinline friend void sort(StackItemT& a, StackItemT& b, StackItemT& c, StackItemT& d)
    {
      int64 s1 = a.all;
      int64 s2 = b.all;
      int64 s3 = c.all;
      int64 s4 = d.all;
      if (s2 < s1) swap(s2,s1);
      if (s4 < s3) swap(s4,s3);
      if (s3 < s1) swap(s3,s1);
      if (s4 < s2) swap(s4,s2);
      if (s3 < s2) swap(s3,s2);
      a.all = s1;
      b.all = s2;
      c.all = s3;
      d.all = s4;
    }

  public:
    union {
      struct { unsigned ptr; float dist; };
      int64 all;
    };
  };


  struct __aligned(8) StackItemInt64
  {

#if 1
    /*! Sort 2 stack items. */
    __forceinline friend void sort(StackItemInt64& s1, StackItemInt64& s2) {
      if (s2.dist < s1.dist) swap(s2,s1);
    }
    
    /*! Sort 3 stack items. */
    __forceinline friend void sort(StackItemInt64& s1, StackItemInt64& s2, StackItemInt64& s3)
    {
      if (s2.dist < s1.dist) swap(s2,s1);
      if (s3.dist < s2.dist) swap(s3,s2);
      if (s2.dist < s1.dist) swap(s2,s1);
    }
    
    /*! Sort 4 stack items. */
    __forceinline friend void sort(StackItemInt64& s1, StackItemInt64& s2, StackItemInt64& s3, StackItemInt64& s4)
    {
      if (s2.dist < s1.dist) swap(s2,s1);
      if (s4.dist < s3.dist) swap(s4,s3);
      if (s3.dist < s1.dist) swap(s3,s1);
      if (s4.dist < s2.dist) swap(s4,s2);
      if (s3.dist < s2.dist) swap(s3,s2);
    }
#else

    /*! Sort 2 stack items. */
    __forceinline friend void sort(StackItemInt64& s1, StackItemInt64& s2) {
      if (s2.i64 < s1.i64) swap(s2.i64,s1.i64);
    }
    
    /*! Sort 3 stack items. */
    __forceinline friend void sort(StackItemInt64& s1, StackItemInt64& s2, StackItemInt64& s3)
    {
      if (s2.i64 < s1.i64) swap(s2.i64,s1.i64);
      if (s3.i64 < s2.i64) swap(s3.i64,s2.i64);
      if (s2.i64 < s1.i64) swap(s2.i64,s1.i64);
   }
    
    /*! Sort 4 stack items. */
    __forceinline friend void sort(StackItemInt64& s1, StackItemInt64& s2, StackItemInt64& s3, StackItemInt64& s4)
    {
      if (s2.i64 < s1.i64) swap(s2.i64,s1.i64);
      if (s4.i64 < s3.i64) swap(s4.i64,s3.i64);
      if (s3.i64 < s1.i64) swap(s3.i64,s1.i64);
      if (s4.i64 < s2.i64) swap(s4.i64,s2.i64);
      if (s3.i64 < s2.i64) swap(s3.i64,s2.i64);
    }

#endif
    
  public:
    union {
      struct {
	unsigned int ptr; 
	unsigned int dist;
      };
      uint64 i64;
    };

  };

 
#if defined(__MIC__)

  /*! An item on the stack holds the node and distance. */
  struct __aligned(16) StackItem 
  {
    /*! pointer to the node */
    void* ptr;

    /*! distance of stack node */
    float dist;

    /*! Sort 2 stack items. */
    __forceinline friend void sort(StackItem& s1, StackItem& s2) 
    {
      if (s2.dist < s1.dist) swap(s2,s1);
    }

    /*! Sort 3 stack items. */
    __forceinline friend void sort(StackItem& s1, StackItem& s2, StackItem& s3)
    {
      if (s2.dist < s1.dist) swap(s2,s1);
      if (s3.dist < s2.dist) swap(s3,s2);
      if (s2.dist < s1.dist) swap(s2,s1);
    }
    
    /*! Sort 4 stack items. */
    __forceinline friend void sort(StackItem& s1, StackItem& s2, StackItem& s3, StackItem& s4)
    {
      if (s2.dist < s1.dist) swap(s2,s1);
      if (s4.dist < s3.dist) swap(s4,s3);
      if (s3.dist < s1.dist) swap(s3,s1);
      if (s4.dist < s2.dist) swap(s4,s2);
      if (s3.dist < s2.dist) swap(s3,s2);
    }
  };

#else

  struct __aligned(16) StackItem  
  {  
    /*! Copy operator */  
    StackItem& operator=(const StackItem& other) { all = other.all; return *this; }  
    
    union {  
      struct { float dist; float pad; void* ptr; };  
      __m128 all;  
    };  
  };  
  
  /*! Sort 2 stack items. */  
  __forceinline void sort(StackItem& a, StackItem& b)  
  {  
    __m128 s1 = a.all;  
    __m128 s2 = b.all;  
    if (_mm_comilt_ss(s2, s1)) swap(s2, s1);  
    a.all = s1;  
    b.all = s2;  
  }  

  /*! Sort 3 stack items. */  
  __forceinline void sort(StackItem& a, StackItem& b, StackItem& c)  
  {  
    __m128 s1 = a.all;  
    __m128 s2 = b.all;  
    __m128 s3 = c.all;  
    if (_mm_comilt_ss(s2, s1)) swap(s2, s1);  
    if (_mm_comilt_ss(s3, s2)) swap(s3, s2);  
    if (_mm_comilt_ss(s2, s1)) swap(s2, s1);  
    a.all = s1;  
    b.all = s2;  
    c.all = s3;  
  }  
  
  /*! Sort 4 stack items. */  
  __forceinline void sort(StackItem& a, StackItem& b, StackItem& c, StackItem& d)  
  {  
    __m128 s1 = a.all;  
    __m128 s2 = b.all;  
    __m128 s3 = c.all;  
    __m128 s4 = d.all;  
    if (_mm_comilt_ss(s2, s1)) swap(s2,s1);  
    if (_mm_comilt_ss(s4, s3)) swap(s4,s3);  
    if (_mm_comilt_ss(s3, s1)) swap(s3,s1);  
    if (_mm_comilt_ss(s4, s2)) swap(s4,s2);  
    if (_mm_comilt_ss(s3, s2)) swap(s3,s2);  
    a.all = s1;  
    b.all = s2;  
    c.all = s3;  
    d.all = s4;  
  }  
#endif

  struct __aligned(16) StackItemNearFar 
  {
  public:
    __forceinline static void swap2(StackItemNearFar& a, StackItemNearFar& b) 
    { 
#if defined(__AVX__)
      ssef sse_a = load4f(&a);
      ssef sse_b = load4f(&b);
      store4f(&a,sse_b);
      store4f(&b,sse_a);
#else
      StackItemNearFar t = b; b = a; a = t;
#endif
    }
    
    __forceinline friend bool operator<(const StackItemNearFar& s1, const StackItemNearFar& s2) {
      return s1.tNear > s2.tNear;
    }
    
    /*! Sort 2 stack items. */
    __forceinline friend void sort(StackItemNearFar& s1, StackItemNearFar& s2) {
      if (s2.tNear < s1.tNear) swap2(s2,s1);
    }
    
    /*! Sort 3 stack items. */
    __forceinline friend void sort(StackItemNearFar& s1, StackItemNearFar& s2, StackItemNearFar& s3)
    {
      if (s2.tNear < s1.tNear) swap2(s2,s1);
      if (s3.tNear < s2.tNear) swap2(s3,s2);
      if (s2.tNear < s1.tNear) swap2(s2,s1);
    }
    
    /*! Sort 4 stack items. */
    __forceinline friend void sort(StackItemNearFar& s1, StackItemNearFar& s2, StackItemNearFar& s3, StackItemNearFar& s4)
    {
      if (s2.tNear < s1.tNear) swap2(s2,s1);
      if (s4.tNear < s3.tNear) swap2(s4,s3);
      if (s3.tNear < s1.tNear) swap2(s3,s1);
      if (s4.tNear < s2.tNear) swap2(s4,s2);
      if (s3.tNear < s2.tNear) swap2(s3,s2);
    }
    
  public:
    size_t ref;
    float tNear,tFar;
  };
}
