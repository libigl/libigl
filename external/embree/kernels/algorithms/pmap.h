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

#include "common/default.h"
#include "common/buffer.h"
#include "algorithms/sort.h"
#include "algorithms/parallel_for.h"

namespace embree
{
  /*! implementation of a map key/value map with parallel construction */
  template<typename Key, typename Val>
  class pmap
  {
    struct KeyValue
    {
      __forceinline KeyValue () {}

      __forceinline KeyValue (const Key key, const Val val)
	: key(key), val(val) {}

      __forceinline operator Key() const {
	return key;
      }

    public:
      Key key;
      Val val;
    };

  public:
    
    /*! parallel map constructors */
    pmap () {}

    template<typename SourceKey>
      pmap (const std::vector<SourceKey>& keys, const std::vector<Val>& values) { init(keys,values); }

    template<typename SourceKey>
      pmap (const BufferT<SourceKey>& keys, const BufferT<Val>& values) { init(keys,values); }

    /*! initialized the parallel map from a vector with keys and values */
    template<typename SourceKey>
      void init(const std::vector<SourceKey>& keys, const std::vector<Val>& values) 
    {
      /* reserve sufficient space for all data */
      assert(keys.size() == values.size());
      vec.resize(keys.size());
      temp.resize(keys.size());

      /* generate key/value pairs */
      parallel_for( size_t(0), keys.size(), size_t(4*4096), [&](const range<size_t>& r) {
	for (size_t i=r.begin(); i<r.end(); i++)
	  vec[i] = KeyValue((Key)keys[i],values[i]);
      });

      /* perform parallel radix sort of the key/value pairs */
      radix_sort<KeyValue,Key>(vec.data(),temp.data(),keys.size());
    }

    /*! initialized the parallel map from user buffers with keys and values */
    template<typename SourceKey>
      void init(const BufferT<SourceKey>& keys, const BufferT<Val>& values) 
    {
      /* reserve sufficient space for all data */
      assert(keys.size() == values.size());
      vec.resize(keys.size());
      temp.resize(keys.size());
      
      /* generate key/value pairs */
      parallel_for( size_t(0), keys.size(), size_t(4*4096), [&](const range<size_t>& r) {
	for (size_t i=r.begin(); i<r.end(); i++)
	  vec[i] = KeyValue((Key)keys[i],values[i]);
      });

      /* perform parallel radix sort of the key/value pairs */
      radix_sort<KeyValue,Key>(vec.data(),temp.data(),keys.size());
    }

    /*! Returns a pointer to the value associated with the specified key. The pointer will be NULL of the key is not contained in the map. */
    __forceinline const Val* lookup(const Key& key) const 
    {
      typename std::vector<KeyValue>::const_iterator i = std::lower_bound(vec.begin(), vec.end(), key);
      if (i == vec.end()) return NULL;
      if (i->key != key) return NULL;
      return &i->val;
    }

    /*! If the key is in the map, the function returns the value associated with the key, otherwise it returns the default value. */
    __forceinline Val lookup(const Key& key, const Val& def) const 
    {
      typename std::vector<KeyValue>::const_iterator i = std::lower_bound(vec.begin(), vec.end(), key);
      if (i == vec.end()) return def;
      if (i->key != key) return def;
      return i->val;
    }

    /*! cleans temporary state required for re-construction */
    void cleanup() {
      temp.clear();
    }

    /*! clears all state */
    void clear() {
      vec.clear();
      temp.clear();
    }

  private:
    std::vector<KeyValue> vec;    //!< vector containing sorted elements
    std::vector<KeyValue> temp;   //!< temporary vector required during construction only
  };
}
