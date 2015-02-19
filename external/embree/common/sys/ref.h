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

#include "constants.h"
#include "sync/atomic.h"

namespace embree
{
  class RefCount
  {
  public:
    RefCount(int val = 0) : refCounter(val) {}
    virtual ~RefCount() {};
  
    virtual void refInc() { atomic_add(&refCounter,1); }
    virtual void refDec() { if (atomic_add(&refCounter,-1) == 1) delete this; }
  private:
    atomic_t refCounter;
  };
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Reference to single object
  ////////////////////////////////////////////////////////////////////////////////

  template<typename Type>
  class Ref {
  public:
    Type* const ptr;

    ////////////////////////////////////////////////////////////////////////////////
    /// Constructors, Assignment & Cast Operators
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Ref( void ) : ptr(NULL) {}
    __forceinline Ref(NullTy) : ptr(NULL) {}
    __forceinline Ref( const Ref& input ) : ptr(input.ptr) { if ( ptr ) ptr->refInc(); }

    __forceinline Ref( Type* const input ) : ptr(input) {
      if ( ptr )
        ptr->refInc();
    }

    __forceinline ~Ref( void ) {
      if (ptr) ptr->refDec();
    }

    __forceinline Ref& operator =( const Ref& input )
    {
      if ( input.ptr ) input.ptr->refInc();
      if (ptr) ptr->refDec();
      *(Type**)&ptr = input.ptr;
      return *this;
    }

    __forceinline Ref& operator =( NullTy ) {
      if (ptr) ptr->refDec();
      *(Type**)&ptr = NULL;
      return *this;
    }

    __forceinline operator bool( void ) const { return ptr != NULL; }

    __forceinline const Type& operator  *( void ) const { return *ptr; }
    __forceinline       Type& operator  *( void )       { return *ptr; }
    __forceinline const Type* operator ->( void ) const { return  ptr; }
    __forceinline       Type* operator ->( void )       { return  ptr; }

    template<typename TypeOut>
    __forceinline       Ref<TypeOut> cast()       { return Ref<TypeOut>(static_cast<TypeOut*>(ptr)); }
    template<typename TypeOut>
    __forceinline const Ref<TypeOut> cast() const { return Ref<TypeOut>(static_cast<TypeOut*>(ptr)); }

    template<typename TypeOut>
    __forceinline       Ref<TypeOut> dynamicCast()       { return Ref<TypeOut>(dynamic_cast<TypeOut*>(ptr)); }
    template<typename TypeOut>
    __forceinline const Ref<TypeOut> dynamicCast() const { return Ref<TypeOut>(dynamic_cast<TypeOut*>(ptr)); }
  };

  template<typename Type> __forceinline  bool operator < ( const Ref<Type>& a, const Ref<Type>& b ) { return a.ptr <  b.ptr ; }

  template<typename Type> __forceinline  bool operator ==( const Ref<Type>& a, NullTy             ) { return a.ptr == NULL  ; }
  template<typename Type> __forceinline  bool operator ==( NullTy            , const Ref<Type>& b ) { return NULL  == b.ptr ; }
  template<typename Type> __forceinline  bool operator ==( const Ref<Type>& a, const Ref<Type>& b ) { return a.ptr == b.ptr ; }

  template<typename Type> __forceinline  bool operator !=( const Ref<Type>& a, NullTy             ) { return a.ptr != NULL  ; }
  template<typename Type> __forceinline  bool operator !=( NullTy            , const Ref<Type>& b ) { return NULL  != b.ptr ; }
  template<typename Type> __forceinline  bool operator !=( const Ref<Type>& a, const Ref<Type>& b ) { return a.ptr != b.ptr ; }
}
