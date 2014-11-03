// +-------------------------------------------------------------------------
// | iterPool.h
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
// +-------------------------------------------------------------------------
// | WHAT IS THIS?
// | 
// | An IterPool is
// |    *)  a memory pool -- new items can be requested, old items released
// |    *)  iterable -- all of the allocated items can be iterated over
// | It was built primarily to support dynamic connectivity structures
// | for triangle meshes, such as those used during re-meshing.
// |
// | These processes often require the frequent allocation and
// | deallocation of mesh elements and become easier to write when
// | elements do not move around in memory, since valid pointers can be
// | maintained.  std::vector does not provide guaranteed static locations.
// |
// | Using a memory pool--threaded with a linked list--provides fast
// | allocation/deallocation and relatively speedy enumeration of all
// | currently allocated elements.  While not as memory access friendly
// | as an array or std::vector, the cost of linked list enumeration
// | comes out in the wash compared to the pointer chasing performed
// | for each element.
// +-------------------------------------------------------------------------

#include "memPool.h"
#include <utility>

template<class T>
class IterPool
{
private:
    uint numAlloced;
public:
    IterPool(int minInitBlocks=10) :
        numAlloced(0),
        block_list(nullptr),
        pool(minInitBlocks)
    {}
    IterPool(IterPool &&src)
        : numAlloced(src.numAlloced),
          block_list(src.block_list),
          pool(std::move(src.pool))
    {
        src.block_list = nullptr;
    }
    ~IterPool() {
        // run through and destruct all remaining elements
        for_each([](T* obj) {
            obj->~T();
        });
    }
    
    void clear() {
        for_each([](T* obj) {
            obj->~T();
        });
        numAlloced = 0;
        block_list = nullptr;
        pool.clear();
    }
    
    void operator=(IterPool &&src)
    {
        for_each([](T* obj) {
            obj->~T();
        });
        block_list = src.block_list;
        src.block_list = nullptr;
        pool = std::move(src.pool);
    }
    
    
private:
    struct Block {
        T       datum;
        Block   *next;
        Block   *prev;
    };
    Block *block_list;
    
public: // allocation/deallocation support
    T* alloc() {
        Block *new_block = pool.alloc();
        if(block_list) block_list->prev = new_block;
        new_block->next = block_list;
        new_block->prev = NULL;
        block_list = new_block;
        
        T* obj = (T*)new_block;
        new (obj) T(); // invoke default constructor when allocating
        
        numAlloced++;
        
        return obj;
    }
    void free(T* item) {
        if(item == NULL)   return;
        item->~T(); // invoke destructor before releasing
        
        numAlloced--;
        
        Block *ptr = (Block*)(item);
        if(ptr->next)   ptr->next->prev = ptr->prev;
        if(ptr->prev)   ptr->prev->next = ptr->next;
        if(ptr == block_list)   block_list = ptr->next;
        pool.free(ptr);
    }
    
public:
    inline void for_each(std::function<void(T*)> func) const {
        for(Block *block = block_list;
          block != NULL;
          block = block->next) {
            func((T*)(block));
        }
    }
    inline bool contains(T* tptr) const {
        for(Block *block = block_list;
          block != NULL;
          block = block->next) {
            if(tptr == (T*)(block))
                return true;
        }
        return false;
    }
    inline uint size() const {
        return numAlloced;
    }
public: // iteration support
    class iterator {
    public:
        iterator() : ptr(NULL) {}
        iterator(Block *init) : ptr(init) {}
        iterator(const iterator &cp) : ptr(cp.ptr) {}
        
        iterator& operator++() { // prefix version
            ptr = ptr->next;
            return *this;
        }
        iterator operator++(int) {
            iterator it(ptr);
            ptr = ptr->next;
            return it;
        } // postfix version
        T& operator*() {
            return ptr->datum;
        }
        T* operator->() {
            return (T*)(ptr);
        }
        bool operator==(const iterator &rhs) {
            return ptr == rhs.ptr;
        }
        bool operator!=(const iterator &rhs) {
            return ptr != rhs.ptr;
        }
    private:
        Block *ptr;
    };
    
    iterator begin() {
        return iterator(block_list);
    }
    iterator end() {
        return iterator(NULL);
    }
    
private:
    MemPool<Block> pool;
};



