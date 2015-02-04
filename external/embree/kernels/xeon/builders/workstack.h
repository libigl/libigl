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

namespace embree
{
  template<class T, unsigned int SIZE>
    class WorkStack 
  {
    ALIGNED_CLASS;
  public:
    AlignedAtomicMutex __aligned(64) mutex;
    __aligned(64) T t[SIZE];
    
    __forceinline void init() {
      mutex.reset();
      mutex.index = 0;
    }

    __forceinline void reset() {
      mutex.index = 0;
    }

    __forceinline WorkStack() {
      init();
    }

    __forceinline bool isFull(const size_t plus = 0) {
      return mutex.index + plus >= SIZE;
    }

    __forceinline bool isEmpty() {
      return mutex.index == 0;
    }
  
    __forceinline size_t size() {
      return mutex.index;
    }

    __forceinline bool islocked() {
      return mutex.islocked();
    }

    __forceinline void lock() {
      mutex.lock();
    }

    __forceinline void unlock() {
      mutex.unlock();
    }

    __forceinline bool request_lock(const size_t num)
    {
      lock();
      if (isFull(num)) { unlock(); return false; }    
      return true;
    }

    __forceinline void push_nolock(const T &v)
    {
      if (unlikely(isFull())) FATAL("stack is full");
      t[mutex.index++] = v;
    }

    __forceinline bool push(const T &v)
    {
      lock();
      if (isFull()) { unlock(); return false; }
      t[mutex.index++] = v;
      unlock();
      return true;
    }

    __forceinline bool pop_nolock(T &v)
    {
      if (isEmpty()) return false;
      v = t[--mutex.index];
      return true;
    }

    __forceinline bool pop_nolock_largest(T &v)
    {
      if (isEmpty()) return false;
      size_t largest = 0;
      for (size_t i=1;i<mutex.index;i++)
	if (t[i] > t[largest])
	  largest = i;
      v = t[largest];
      t[largest] = t[mutex.index-1];
      mutex.index--;

      return true;
    }

    __forceinline bool try_pop(T &v)
    {
      if (isEmpty() || islocked()) { return false; }
      lock();
      if (isEmpty()) 
	{	
	  unlock(); return false; 
	}
      mutex.index--;
      v = t[mutex.index];
      unlock();
      return true;
    }

    __forceinline bool pop(T &v)
    {
      if (isEmpty()) { return false; }
      lock();
      if (isEmpty()) 
	{	
	  unlock(); return false; 
	}
      mutex.index--;
      v = t[mutex.index];
      unlock();
      return true;
    }

    __forceinline bool pop_largest(T &v)
    {
      if (isEmpty()) { return false; }
      lock();
      if (isEmpty()) 
	{	
	  unlock(); return false; 
	}

      size_t largest = 0;
      for (size_t i=1;i<mutex.index;i++)
	if (t[i] > t[largest])
	  largest = i;
      v = t[largest];
      t[largest] = t[mutex.index-1];
      mutex.index--;
      unlock();
      return true;
    }

    __forceinline bool pop_smallest(T &v)
    {
      if (isEmpty()) { return false; }
      lock();
      if (isEmpty()) 
	{	
	  unlock(); return false; 
	}
      
      size_t smallest = 0;
      for (size_t i=1;i<mutex.index;i++)
	if (t[i] < t[smallest])
	  smallest = i;
      v = t[smallest];
      t[smallest] = t[mutex.index-1];
      mutex.index--;
      unlock();
      return true;
    }

    __forceinline bool pop_first(T &v)
    {

      lock();
      if (isEmpty()) { unlock(); return false; }
      v = t[0];
      for (size_t i=0;i<mutex.index-1;i++)
	t[i] = t[i+1];
      mutex.index--;
      unlock();
      return true;
    }

    __forceinline T* begin()
    {
      return t;
    }

    __forceinline T* end()
    {
      return &t[mutex.index];
    }

    __forceinline T& get(const size_t index)
    {
      return t[index];
    }

    __forceinline void increaseSize(const size_t plus)
    {
      mutex.index += plus;
    }
  };

  template<class T>
    class WorkHeap 
  {
    ALIGNED_CLASS;
  public:

    WorkHeap() {
      mutex.reset();
    }

    void reset() {
      heap.clear();
    }

    size_t size() const { 
      return heap.size();
    }

    T& front() { return heap[0]; }
    T* begin() { return &heap[0]; }
    T* end  () { return &heap[0]+heap.size(); }

    void push(T& br)
    {
      heap.push_back(br);
      std::push_heap(heap.begin(),heap.end());
    }

    bool pop(T& br)
    {
      mutex.lock();
      if  (heap.size() == 0) {
	mutex.unlock();
	return false;
      }
      br = heap.front();
      std::pop_heap(heap.begin(),heap.end());
      heap.pop_back();
      mutex.unlock();
      return true;
    }
    
  private:
    AlignedAtomicMutex __aligned(64) mutex;
    vector_t<T> heap;
  };
}
