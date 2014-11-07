// +-------------------------------------------------------------------------
// | shortVec.h
// | 
// | Author: Gilbert Bernstein
// +-------------------------------------------------------------------------
// | COPYRIGHT:
// |    Copyright Gilbert Bernstein 2013
// |    See the included COPYRIGHT file for further details.
// |    
// |    This file is part of the Cork library.
// |
// |    Cork is free software: you can redistribute it and/or modify
// |    it under the terms of the GNU Lesser General Public License as
// |    published by the Free Software Foundation, either version 3 of
// |    the License, or (at your option) any later version.
// |
// |    Cork is distributed in the hope that it will be useful,
// |    but WITHOUT ANY WARRANTY; without even the implied warranty of
// |    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// |    GNU Lesser General Public License for more details.
// |
// |    You should have received a copy 
// |    of the GNU Lesser General Public License
// |    along with Cork.  If not, see <http://www.gnu.org/licenses/>.
// +-------------------------------------------------------------------------
#pragma once

#include "prelude.h"
#include "memPool.h"

#include <algorithm>

// Don't know why, but can't get this to compile
// when I make the datablock def. a member of ShortVec<T,LEN>
template<class T, uint LEN>
struct ShortVecBlock_Private {
    byte data[sizeof(T)*LEN];
};
// We allocate these blocks instead of typed arrays
// in order to ensure that we get control of allocation/deallocation
// rather than the allocator attempting to do so.

template<class T, uint LEN>
class ShortVec
{
public: // constructor/destructor
    ShortVec(uint size = 0);
    ShortVec(uint size, const T &fill_val);
    ShortVec(const ShortVec<T,LEN> &cp);
    //ShortVec(ShortVec<T,LEN> &&cp); // move constructor
    ~ShortVec();
    
    ShortVec<T,LEN>& operator=(const ShortVec<T,LEN> &vec);
    
public: // index accessors
    inline       T& operator[](uint i)       { return data[i]; }
    inline const T& operator[](uint i) const { return data[i]; }
    
public: // iterators
    typedef       T* iterator;
    typedef const T* const_iterator;
    iterator        begin()       { return data; }
    const_iterator  begin() const { return data; }
    iterator        end()         { return data + user_size; }
    const_iterator  end()   const { return data + user_size; }
    
public: // inspectors
    inline uint size() const { return user_size; }
    
public: // modifiers
    void resize(uint newsize);
    void push_back(const T &datum);
    void erase(const T &val); // erase if it can be found
    
private: // helper functions
    T*   allocData(uint space, uint &allocated);
    void deallocData(T* data_ptr, uint allocated);
    
    void constructRange(T* array, int begin, int end);
    void copyConstructRange(T* src, T* dest, int begin, int end);
    void destructRange(T* array, int begin, int end);
    
    // resize, manage allocation/deallocation,
    // but not construction/destruction
    void resizeHelper(uint newsize);
    
public: // shared data structures and data.
    static MemPool< ShortVecBlock_Private<T,LEN> > pool;
    
private: // instance data
    uint user_size;     // actual number of entries from client perspective
    uint internal_size; // number of entries allocated, if greater than LEN;
                        // if allocated from the memory pool, this is 0
    T* data;
};

template<class T, uint LEN>
MemPool< ShortVecBlock_Private<T,LEN> > ShortVec<T,LEN>::pool;

template<class T, uint LEN> inline
T* ShortVec<T,LEN>::allocData(uint space, uint &allocated)
{
    T* result;
    if(space <= LEN) {
        allocated = LEN;
        result =  reinterpret_cast<T*>(pool.alloc());
    } else {
        allocated = space;
        result = reinterpret_cast<T*>(new byte[sizeof(T)*space]);
    }
    //if(LEN == 2) std::cout << "        Allocing:   " << result << std::endl;
    return result;
}
template<class T, uint LEN> inline
void ShortVec<T,LEN>::deallocData(T* data_ptr, uint allocated)
{
    //if(LEN == 2) std::cout << "        Deallocing: " << data_ptr << std::endl;
    if(allocated <= LEN)
        pool.free(reinterpret_cast< ShortVecBlock_Private<T,LEN>* >(data_ptr));
    else
        delete[] reinterpret_cast<byte*>(data_ptr);
}

template<class T, uint LEN> inline
void ShortVec<T,LEN>::constructRange(T* array, int begin, int end)
{
    for(int i=begin; i<end; i++)
        new (&(array[i])) T();
}
template<class T, uint LEN> inline
void ShortVec<T,LEN>::copyConstructRange(T* src, T* dest, int begin, int end)
{
    // copy actual data over
    for(int i=begin; i<end; i++)
        new (&(dest[i])) T(src[i]);
}
template<class T, uint LEN> inline
void ShortVec<T,LEN>::destructRange(T* array, int begin, int end)
{
    for(int i=begin; i<end; i++)
        (&(array[i]))->~T();
}

// we use a strictly increasing allocation size policy
// with array length doubling to ensure that the cost of
// copying array entries on a reallocation has
// constant amortized cost.  This is important when the
// vector is used to accumulate a list of values.
template<class T, uint LEN> inline
void ShortVec<T,LEN>::resizeHelper(uint newsize) {
    if(newsize > internal_size) { // we need more space!
        // setup the new data block with at least twice as much space
        uint new_space;
        T *new_data = allocData(std::max(newsize, internal_size*2), new_space);
        // copy data and destroy old copies
        copyConstructRange(data, new_data, 0, user_size);
        destructRange(data, 0, user_size);
        // free old data
        deallocData(data, internal_size);
        data = new_data;
        internal_size = new_space;
    }
    user_size = newsize;
}




template<class T, uint LEN> inline
ShortVec<T,LEN>::ShortVec(uint size) : user_size(size)
{
    data = allocData(user_size, internal_size);
    constructRange(data, 0, user_size);
}
template<class T, uint LEN> inline
ShortVec<T,LEN>::ShortVec(uint size, const T &fill_val) : user_size(size)
{
    data = allocData(user_size, internal_size);
    // then fill
    for(uint i=0; i<user_size; i++)
        new (&data[i]) T(fill_val);
}
template<class T, uint LEN> inline
ShortVec<T,LEN>::ShortVec(const ShortVec<T,LEN> &cp) : user_size(cp.user_size)
{
    data = allocData(user_size, internal_size);
    // copy actual data over
    copyConstructRange(cp.data, data, 0, user_size);
}
/*template<class T, uint LEN> inline
ShortVec<T,LEN>::ShortVec(ShortVec<T,LEN> &&cp)
{
    user_size = cp.user_size;
    internal_size = cp.internal_size;
    data = cp.data;
    cp.user_size = 0; // ensure that no destructors are called
    cp.internal_size = 0; // ensure that pool free is called
    cp.data = NULL; // on a null pointer, which will do nothing
}*/
template<class T, uint LEN> inline
ShortVec<T,LEN>::~ShortVec()
{
    destructRange(data, 0, user_size);
    deallocData(data, internal_size);
}

template<class T, uint LEN> inline
ShortVec<T,LEN>& ShortVec<T,LEN>::operator=(const ShortVec<T,LEN> &vec)
{
    uint old_size = user_size;
    
    // ensure there is enough space allocated at the destination
    resizeHelper(vec.user_size);
    
    // copy assignment for all data in range overlap
    for(uint i=0; i<std::min(vec.user_size, old_size); i++)
        data[i] = vec.data[i];
    
    // if the new range is larger, copy construct the portion
    // outside of the old range
    if(vec.user_size > old_size) {
        copyConstructRange(vec.data, data, old_size, vec.user_size);
    }
    
    // if the new range is smaller, destruct old unused entries
    if(vec.user_size < old_size)
        destructRange(data, vec.user_size, old_size);
    
    return *this;
}

template<class T, uint LEN> inline
void ShortVec<T,LEN>::resize(uint newsize) {
    uint oldsize = user_size;
    
    resizeHelper(newsize);
    
    if(oldsize < newsize)
        constructRange(data, oldsize, newsize);
    
    if(newsize < oldsize)
        destructRange(data, newsize, oldsize);
}

template<class T, uint LEN> inline
void ShortVec<T,LEN>::push_back(const T &datum)
{
    uint i = user_size;
    resizeHelper(user_size+1); // make room
    new (&(data[i])) T(datum);
}

template<class T, uint LEN> inline
void ShortVec<T,LEN>::erase(const T &val)
{
    for(uint i=0; i<user_size; i++) {
        if(data[i] == val) {
            std::swap(data[i], data[user_size-1]);
            resize(user_size-1);
            break;
        }
    }
}



