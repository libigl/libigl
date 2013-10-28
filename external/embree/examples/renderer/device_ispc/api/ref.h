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

#ifndef __EMBREE_ISPC_REF_H__
#define __EMBREE_ISPC_REF_H__

#include "default.h"
#include "ref_ispc.h"

namespace embree
{
  ////////////////////////////////////////////////////////////////////////////////
  /// Reference to ISPC object
  ////////////////////////////////////////////////////////////////////////////////

  class ISPCRef {
  public:
    void* ptr;

  public:
    __forceinline ISPCRef( void ) : ptr(NULL) {}
    __forceinline ISPCRef(NullTy) : ptr(NULL) {}

    __forceinline ISPCRef( const ISPCRef& input ) : ptr(input.ptr) { 
      ispc::RefCount__incref(ptr); 
    }

    __forceinline ISPCRef( void* input ) : ptr(input) {
      ispc::RefCount__incref(ptr); 
    }

    __forceinline ~ISPCRef( void ) {
      ispc::RefCount__decref(ptr); 
    }

    __forceinline ISPCRef& operator =( const ISPCRef& input ) {
      ispc::RefCount__incref(input.ptr); 
      ispc::RefCount__decref(ptr); 
      ptr = input.ptr;
      return *this;
    }

    __forceinline ISPCRef& operator =( NullTy ) {
      ispc::RefCount__decref(ptr); 
      ptr = NULL;
      return *this;
    }

    __forceinline operator bool( void ) const { return ptr != NULL; }
  };
}

#endif

