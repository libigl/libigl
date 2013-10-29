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

#ifndef __EMBREE_MUTEX_H__
#define __EMBREE_MUTEX_H__

#include "../platform.h"
#include "atomic.h"

namespace embree
{
  /*! system mutex */
  class MutexSys {
    friend class ConditionSys;
  public:

    MutexSys( void );
    ~MutexSys( void );

    void lock( void );
    void unlock( void );

  protected:
    void* mutex;

    MutexSys( const MutexSys& );             // don't implement
    MutexSys& operator =( const MutexSys& ); // don't implement
  };

  /*! spinning mutex */
  class MutexActive {
  public:
    __forceinline MutexActive( void ) : flag(0) {}
    void lock  ( void );
    void unlock( void );
  protected:
    Atomic flag;

    MutexActive( const MutexActive& );             // don't implement
    MutexActive& operator =( const MutexActive& ); // don't implement
  };

  /*! safe mutex lock and unlock helper */
  template<typename Mutex> class Lock {
  public:
    Lock (Mutex& mutex) : mutex(mutex) { mutex.lock(); }
    ~Lock() { mutex.unlock(); }
  protected:
    Mutex& mutex;
    Lock( const Lock& );             // don't implement
    Lock& operator =( const Lock& ); // don't implement
  };
}

#endif
