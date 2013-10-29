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

#ifndef __EMBREE_HANDLE_H__
#define __EMBREE_HANDLE_H__

#include "sys/platform.h"
#include "sys/constants.h"
#include "device.h"

namespace embree
{
  extern Device* g_device;
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Automatic reference counting for Embree Handles
  ////////////////////////////////////////////////////////////////////////////////

  template<typename Type>
  class Handle 
  {
  public:
    
    __forceinline Handle ( void ) : handle(NULL) {}
    __forceinline Handle (NullTy) : handle(NULL) {}
    __forceinline Handle( Type const input ) : handle(input) {}

    __forceinline Handle( const Handle& input ) : handle(input.handle) { 
      if (handle) g_device->rtIncRef(handle); 
    }

    __forceinline ~Handle( void ) {
      if (handle) g_device->rtDecRef(handle);
    }

    __forceinline Handle& operator =( const Handle& input )
    {
      if (input.handle) g_device->rtIncRef(input.handle);
      if (handle) g_device->rtDecRef(handle);
      handle = input.handle;
      return *this;
    }

    __forceinline Handle& operator =( NullTy ) {
      if (handle) g_device->rtDecRef(handle);
      handle = NULL;
      return *this;
    }

    __forceinline operator bool() const { return handle != NULL; }
    __forceinline operator Type() const { return handle; }

  private:
    Type handle;
  };
}

#endif

