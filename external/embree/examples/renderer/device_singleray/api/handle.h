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

#include "device/device.h"
#include "parms.h"

namespace embree
{
  /* generic handle interface */
  struct _RTHandle : public RefCount {
  public:

    /*! Start counting with 1, as each created handle is owned by the
        application. */
    _RTHandle () : RefCount(1) {}

    /*! Recreated the object the handle refers to. */
    virtual void create() = 0;

    /*! Clear handle. */
    virtual void clear() { }

    /*! Sets a parameter of the handle. */
    virtual void set(const std::string& property, const embree::Variant& data) = 0;
  };

  /*******************************************************************
                 generic handle implementations
  *******************************************************************/

  /*! Constant handle. A special handle type that does not allow
   *  setting parameters. */
  template<typename T> class ConstHandle : public _RTHandle {
  public:

    /*! Creates a constant handle from the object to reference. */
    ConstHandle(const Ref<T>& ptr) : instance(ptr) {}

    /*! Recreating the underlying object is not allowed. */
    void create() { throw std::runtime_error("cannot modify constant handle"); }

    /*! Setting parameters is not allowed. */
    void set(const std::string& property, const Variant& data) { throw std::runtime_error("cannot modify constant handle"); }

    Ref<T> getInstance() { return instance; }

  public:
    Ref<T> instance;  //!< Referenced object.
  };

  /*! Base Handle */
  template<typename B> class InstanceHandle : public _RTHandle {
  public:
    InstanceHandle () {}

    Ref<B> getInstance() { return instance; }

    /*! checks if object is newer than other object */
    template<typename A>
    bool is_newer_than(const InstanceHandle<A>* other) const {
      return time > other->time;
    }

  public:
    Ref<B> instance; //!< Referenced object.
  };

  /*! Constructor handle. A normal handle type that buffers parameters and
   *  can create the underlying object form constructor. */
  template<typename T, typename B> class ConstructorHandle : public InstanceHandle<B> {
  public:

    ConstructorHandle () 
      : modified(true) {}

    /*! Creates a new object. */
    void create() { 
      if (!modified) return;
      this->instance = new T(this->parms); 
      modified = false; 
    }

    /*! Sets a new parameter. */
    void set(const std::string& property, const Variant& data) { 
      this->parms.add(property,data);
      modified = true; 
    }
    
  private:
    bool modified;   //!< tells if the object got modified
    Parms parms;     //!< Parameter container for handle.
  };

  /*! Create handle. A normal handle type that buffers parameters and
   *  can create the underlying object from a create function. */
  template<typename T, typename B> class CreateHandle : public InstanceHandle<B> {
  public:

    CreateHandle () 
      : modified(true) {}

    /*! Creates a new object. */
    void create() { 
      if (!modified) return;
      this->instance = T::create(this->parms); 
      modified = false; 
    }

    /*! Clear handle. */
    void clear() {
      this->parms.clear();
      modified = true; 
    }

    /*! Sets a new parameter. */
    void set(const std::string& property, const Variant& data) { 
      this->parms.add(property,data);
      modified = true; 
    }
    
  private:
    bool modified;   //!< tells if the object got modified
    Parms parms;     //!< Parameter container for handle.
  };
  
  /*! Safe handle cast. Casts safely to a specified output type.  */
  template<typename T> static Ref<T> castHandle(Device::RTHandle handle_i, const char* name) {
    Ref<T> handle = dynamic_cast<T*>((_RTHandle*)handle_i);
    if (!handle          ) throw std::runtime_error("invalid "+std::string(name)+" handle");
    //if (!handle->instance) throw std::runtime_error("invalid "+std::string(name)+" value");
    return handle;
  }
}

#endif
