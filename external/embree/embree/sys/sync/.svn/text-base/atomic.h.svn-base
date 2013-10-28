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

#ifndef __EMBREE_ATOMIC_H__
#define __EMBREE_ATOMIC_H__

#include "../intrinsics.h"

namespace embree
{
  struct Atomic {
  protected:
    Atomic( const Atomic& ); // don't implement
    Atomic& operator =( const Atomic& ); // don't implement

  public:
    __forceinline Atomic( void ) : data(0) {}
    __forceinline Atomic( const atomic_t data ) : data(data) {}
    __forceinline Atomic& operator =( const atomic_t input ) { data = input; return *this; }
    __forceinline operator atomic_t() const { return data; }

  public:
    __forceinline atomic_t sub( const atomic_t input ) { return atomic_add(&data,-input) - input; }
    __forceinline atomic_t add( const atomic_t input ) { return atomic_add(&data, input) + input; }
    __forceinline friend atomic_t operator +=( Atomic& value, const atomic_t input ) { return atomic_add(&value.data, input) + input; }
    __forceinline friend atomic_t operator ++( Atomic& value ) { return atomic_add(&value.data,  1) + 1; }
    __forceinline friend atomic_t operator --( Atomic& value ) { return atomic_add(&value.data, -1) - 1; }
    __forceinline friend atomic_t operator ++( Atomic& value, int ) { return atomic_add(&value.data,  1); }
    __forceinline friend atomic_t operator --( Atomic& value, int ) { return atomic_add(&value.data, -1); }
    __forceinline friend atomic_t cmpxchg    ( Atomic& value, const atomic_t v, const atomic_t c) { return atomic_cmpxchg(&value.data,v,c); }

  private:
    volatile atomic_t data;
  };
 }

#endif
