// ======================================================================== //
// Copyright 2009-2012 Intel Corporation                                    //
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

#ifndef __EMBREE_STACK_ITEM_H__
#define __EMBREE_STACK_ITEM_H__

namespace embree
{
  /*! swap that is forced to inline */
  template<typename T> __forceinline void swap(T& a, T& b) { T t = b; b = a; a = t; }

#if 0
  /*! An item on the stack holds the node and distance. */
  struct __align(16) StackItem 
  {
    /*! pointer to the node */
    void* ptr;

    /*! distance of stack node */
    float dist;

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
  struct __align(16) StackItem  
  {  
    /*! Copy operator */  
    StackItem& operator=(const StackItem& other) { all = other.all; return *this; }  
    
    union {  
      struct { float dist; float pad; void* ptr; };  
      __m128 all;  
    };  
  };  
  
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
}

#endif
