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

#include "sys/intrinsics.h"

namespace embree
{
  /*! An atomic set. Insert calls are atomic and take calls are atomic
      but only when not called intermixed. Thus one can only
      atomically remove or insert items. */
  template<typename T>
    class __aligned(64) atomic_set 
  {
    ALIGNED_CLASS;
  public:

    /*! set item */
    class item : public T
    {
    public:
      
      /*! default constructor */
      item () : next(NULL) {}

    public:
      item* next;
    };

    /*! Iterator for the set. */
    class iterator
    {
    public:

      /*! default constructor */
     __forceinline  iterator () 
        : root(NULL) {}

      /*! initialize the iterator from a set */
      __forceinline iterator (atomic_set& other) 
        : root(other.root) {}


      /*! return next element */
      __forceinline item* next()
      {
        item* ptr; 
        while (!try_take(ptr));
        return ptr;
      }

    private:
      __forceinline bool try_take(item*& ptr) 
      {
        ptr = root;
        if (ptr == NULL) return true;
        return atomic_cmpxchg_ptr(&root,ptr,ptr->next) == ptr;
      }
  
    private:
      item* root;
    };

    /*! Not thread safe iterator for iterating over elements of a list of blocks. */
    class block_iterator_unsafe
    {
      typedef typename T::T Type;
    public:
      __forceinline block_iterator_unsafe (atomic_set& other) 
        : root(other.root), pos(0) 
      {
	next();
      }
      
      __forceinline void operator++(int)
      {
        pos++;
        next();
      }

      size_t size()
      { 
        size_t s = 0;
        while (root) {
          s += root->size();
          root = root->next;
        }
        return s;
      }

      __forceinline operator bool( ) const { return root; }
      __forceinline const Type& operator*( ) const { return (*root)[pos]; }
      __forceinline       Type& operator*( )       { return (*root)[pos]; }
      __forceinline const Type* operator->( ) const { return &(*root)[pos]; }
      __forceinline       Type* operator->( )       { return &(*root)[pos]; }
  
    private:

      __forceinline void next()
      {
	while (root && pos >= root->size()) {
          root = root->next;
          pos = 0;
        }
      }

    private:
      item* root;
      size_t pos;
    };

  public:

     /*! default constructor */
    __forceinline atomic_set () : root(NULL) {}

    /*! copy constructor */
    __forceinline atomic_set (const atomic_set& other) {
      this->root = other.root; other.root = NULL;
    }

    /*! assignment operator */
    __forceinline atomic_set& operator=(const atomic_set& other) {
      this->root = other.root; other.root = NULL;
      return *this;
    }
 
    /*! add element to front of list */
    __forceinline item* insert(item* ptr) {
      while (!try_insert(ptr)) {
      } 
      return ptr;
    }

    /*! remove element from front of list */
    __forceinline item* take()
    {
      item* ptr; 
      while (!try_take(ptr)) {
      }
      return ptr;
    }

    /*! add element to front of list */
    __forceinline item* insert_unsafe(item* ptr) {
      ptr->next = root;
      root = ptr;
      return ptr;
    }

    /*! remove element from front of list */
    __forceinline item* take_unsafe()
    {
      if (root == NULL) return NULL;
      item* cur = root;
      root = cur->next;
      return cur;
    }

    size_t length() const
    { 
      size_t s = 0;
      item* i = root;
      while (i) {
        s++;
        i = i->next;
      }
      return s;
    }

    __forceinline item* head() {
      return root;
    }

  private:

    __forceinline bool try_insert(item* ptr) 
    {
      if (ptr == NULL) return true;
      item* cur = root;
      ptr->next = cur;
      return atomic_cmpxchg_ptr(&root,cur,ptr) == cur;
    }

    __forceinline bool try_take(item*& ptr) 
    {
      ptr = root;
      if (ptr == NULL) return true;
      return atomic_cmpxchg_ptr(&root,ptr,ptr->next) == ptr;
    }

  private:
    mutable item* root;
  };
}
