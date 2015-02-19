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

#include "builders/builder_util.h"
#include "builders/binning.h"

#include "bvh4i.h"
#include "bvh4i_statistics.h"

namespace embree
{


  class __aligned(16) SmallBuildRecord 
  {
  public:
    unsigned int begin;
    unsigned int end;
    unsigned int depth;
    unsigned int parentNodeID;
    unsigned int parentLocalID;
      
    __forceinline unsigned int size() const {
      return end - begin;
    }

    __forceinline unsigned int items() const {
      return end - begin;
    }
      
    __forceinline void init(const unsigned int _begin, const unsigned int _end)			 
    {
      begin         = _begin;
      end           = _end;
      depth         = 1;
      parentNodeID  = 0;
      parentLocalID = 0;

      assert(begin < end);
    }
      
    __forceinline bool operator<(const SmallBuildRecord& br) const { return size() < br.size(); } 
    __forceinline bool operator>(const SmallBuildRecord& br) const { return size() > br.size(); } 

    __forceinline friend std::ostream &operator<<(std::ostream &o, const SmallBuildRecord &br)
      {
	o << "begin " << br.begin << " end " << br.end << " size " << br.size() << " depth " << br.depth << " parentNodeID " << br.parentNodeID << " parentLocalID " << br.parentLocalID;
	return o;
      }

  };

  struct __aligned(8) MortonID32Bit
  {
    unsigned int code;
    unsigned int index;
            
    __forceinline unsigned int get(const unsigned int shift, const unsigned and_mask) const {
      return (code >> shift) & and_mask;
    }

    __forceinline unsigned int getByte(const size_t b) const {
      assert(b < 4);
      const unsigned char *__restrict const ptr = (const unsigned char*)&code;
      return ptr[b];
    }
      
    __forceinline void operator=(const MortonID32Bit& v) {
      *(size_t*)this = *(size_t*)&v;
    };
      
    __forceinline friend std::ostream &operator<<(std::ostream &o, const MortonID32Bit& mc) {
      o << "index " << mc.index << " code = " << mc.code;
      return o;
    }

    __forceinline bool operator<(const MortonID32Bit &m) const { return code < m.code; } 
    __forceinline bool operator>(const MortonID32Bit &m) const { return code > m.code; } 
  };

  class BVH4iBuilderMorton : public Builder
  {
    ALIGNED_CLASS;

    enum { RECURSE = 1, CREATE_TOP_LEVEL = 2 };

    static const size_t GLOBAL_WORK_STACK_ENTRIES = 512;
    static const size_t MAX_TOP_LEVEL_BINS = 1024;
    static const size_t ALLOCATOR_NODE_BLOCK_SIZE = 64;
    static const size_t MORTON_LEAF_THRESHOLD = 4;
    static const size_t LATTICE_BITS_PER_DIM = 10;
    static const size_t LATTICE_SIZE_PER_DIM = size_t(1) << LATTICE_BITS_PER_DIM;
    typedef AtomicIDBlock<ALLOCATOR_NODE_BLOCK_SIZE> NodeAllocator;

  public:    

    /*! Constructor. */
    BVH4iBuilderMorton (BVH4i* bvh, void* geometry, const bool tree_rotations = false);

    /*! Destructor. */
    ~BVH4iBuilderMorton();

    /*! creates the builder */
    static Builder* create (void* accel, void* geometry, bool tree_rotations) { 
      return new BVH4iBuilderMorton((BVH4i*)accel,geometry,tree_rotations);
    }

    /* build function */
    void build(size_t threadIndex, size_t threadCount);

    /*! initialized morton encoding */
    void initEncodingAllocateData();
    
    /*! allocate data arrays */
    void allocateData(size_t threadCount);

    /*! main build task */
    void build_main(const size_t threadID, const size_t numThreads);

    /*! precalculate some per thread data */
    TASK_FUNCTION(BVH4iBuilderMorton,initThreadState);

    /*! task that calculates the bounding box of the scene */
    TASK_FUNCTION(BVH4iBuilderMorton,computeBounds);

    /*! task that calculates the morton codes for each primitive in the scene */
    TASK_FUNCTION(BVH4iBuilderMorton,computeMortonCodes);
    
    /*! parallel sort of the morton codes */
    TASK_FUNCTION(BVH4iBuilderMorton,radixsort);

    /*! builds top of the tree in parallel */
    TASK_FUNCTION(BVH4iBuilderMorton,createTopLevelTree);

    /*! task that builds a list of sub-trees */
    TASK_FUNCTION(BVH4iBuilderMorton,recurseSubMortonTrees);

  public:


    /*! creates a leaf node */
    BBox3fa createSmallLeaf(SmallBuildRecord& current);


    BBox3fa createLeaf(SmallBuildRecord& current, NodeAllocator& alloc);

    /*! fallback split mode */
    void split_fallback(SmallBuildRecord& current, SmallBuildRecord& leftChild, SmallBuildRecord& rightChild);

    /*! split a build record into two */
    bool split(SmallBuildRecord& current, SmallBuildRecord& left, SmallBuildRecord& right);

    /*! create the top-levels of the tree */
    size_t createSingleBVH4iNode(SmallBuildRecord& current, SmallBuildRecord *__restrict__ const children);

    /*! main recursive build function */
    BBox3fa recurse(SmallBuildRecord& current, 
		    NodeAllocator& alloc);
    
    /*! refit the toplevel part of the BVH */
    BBox3fa refit_toplevel(const BVH4i::NodeRef &ref);

    /*! refit the sub-BVHs */
    BBox3fa refit(const BVH4i::NodeRef &ref);
    
    /*! recreates morton codes when reaching a region where all codes are identical */
    void recreateMortonCodes(SmallBuildRecord& current) const;

  public:
    BVH4i      * bvh;         //!< Output BVH
    Scene      * scene;

    unsigned int encodeShift;
    unsigned int encodeMask;

    __aligned(64) SmallBuildRecord buildRecords[MAX_TOP_LEVEL_BINS];    
    __aligned(64) unsigned int thread_startGroup[MAX_MIC_THREADS];      
    __aligned(64) unsigned int thread_startGroupOffset[MAX_MIC_THREADS];


    /*! state for radix sort */
  public:
    static const size_t RADIX_BITS = 8;
    static const size_t RADIX_BUCKETS = (1 << RADIX_BITS);
    static const size_t RADIX_BUCKETS_MASK = (RADIX_BUCKETS-1);
    __aligned(64) unsigned int radixCount[MAX_MIC_THREADS][RADIX_BUCKETS];

  protected:
    MortonID32Bit * __restrict__ morton;
    BVH4i::Node   * __restrict__ node;
    Triangle1     * __restrict__ accel;

    size_t numPrimitives;
    size_t numGroups;
    size_t numNodes;
    size_t numAllocatedNodes;
    size_t size_morton;
    size_t size_node;
    size_t size_accel;
    size_t topLevelItemThreshold;
    size_t numBuildRecords;
    size_t numPrimitivesOld;
    bool enableTreeRotations;

    __aligned(64) Centroid_Scene_AABB global_bounds;

    /*! node allocator */
    __aligned(64) AlignedAtomicCounter32  atomicID;
    __aligned(64) AlignedAtomicCounter32  numBuildRecordCounter;

    __forceinline unsigned int allocGlobalNode(int size)
    {
      const unsigned int currentIndex = atomicID.add(size);
      if (unlikely(currentIndex >= numAllocatedNodes)) {
        FATAL("not enough nodes allocated");
      }
      return currentIndex;
    }


  };



  struct __aligned(16) MortonID64Bit
  {
    size_t code;
    unsigned int groupID;
    unsigned int primID;
            
    __forceinline unsigned int get(const unsigned int shift, const unsigned and_mask) const {
      return (code >> shift) & and_mask;
    }

    __forceinline unsigned int getByte(const size_t b) const {
      assert(b < 8);
      const unsigned char *const ptr = (const unsigned char*)&code;
      return ptr[b];
    }
      
    __forceinline void operator=(const MortonID64Bit& v) {
      ((size_t*)this)[0] = ((size_t*)&v)[0];
      ((size_t*)this)[1] = ((size_t*)&v)[1];
    };
      
    __forceinline friend std::ostream &operator<<(std::ostream &o, const MortonID64Bit& mc) {
      o << "primID " << mc.primID << " groupID " << mc.groupID << " code = " << mc.code;
      return o;
    }

    __forceinline bool operator<(const MortonID64Bit &m) const { return code < m.code; } 
    __forceinline bool operator>(const MortonID64Bit &m) const { return code > m.code; } 
  };


  __forceinline void xchg(MortonID64Bit &a, MortonID64Bit &b)
  {
    const mic_f ai = broadcast4to16f((float*)&a);
    const mic_f bi = broadcast4to16f((float*)&b);
    store4f((float*)&a,bi);
    store4f((float*)&b,ai);
  }


  class BVH4iBuilderMorton64Bit : public Builder
  {
    ALIGNED_CLASS;

    enum { RECURSE = 1, CREATE_TOP_LEVEL = 2 };

    static const size_t GLOBAL_WORK_STACK_ENTRIES = 512;
    static const size_t MAX_TOP_LEVEL_BINS = 1024;
    static const size_t ALLOCATOR_NODE_BLOCK_SIZE = 64;
    static const size_t MORTON_LEAF_THRESHOLD = 4;
    static const size_t LATTICE_BITS_PER_DIM = 21;
    static const size_t LATTICE_SIZE_PER_DIM = size_t(1) << LATTICE_BITS_PER_DIM;
    typedef AtomicIDBlock<ALLOCATOR_NODE_BLOCK_SIZE> NodeAllocator;

  public:

    /*! Constructor. */
    BVH4iBuilderMorton64Bit(BVH4i* bvh, void* geometry);

    /*! Destructor. */
    ~BVH4iBuilderMorton64Bit();

    /*! creates the builder */
    static Builder* create (void* accel, void* geometry) { 
      return new BVH4iBuilderMorton64Bit((BVH4i*)accel,geometry);
    }

    /* build function */
    void build(size_t threadIndex, size_t threadCount);
    
    /*! allocate data arrays */
    void allocateData(size_t threadCount);

    /*! main build function */
    void build_main(const size_t threadID, const size_t numThreads);

    /*! precalculate some per thread data */
    TASK_FUNCTION(BVH4iBuilderMorton64Bit,initThreadState);
    
    /*! task that calculates the bounding box of the scene */
    TASK_FUNCTION(BVH4iBuilderMorton64Bit,computeBounds);

    /*! task that calculates the morton codes for each primitive in the scene */
    TASK_FUNCTION(BVH4iBuilderMorton64Bit,computeMortonCodes);
    
    /*! parallel sort of the morton codes */
    TASK_FUNCTION(BVH4iBuilderMorton64Bit,radixsort);

    /*! builds top of the tree in parallel */
    TASK_FUNCTION(BVH4iBuilderMorton64Bit,createTopLevelTree);

    /*! task that builds a list of sub-trees */
    TASK_FUNCTION(BVH4iBuilderMorton64Bit,recurseSubMortonTrees);

    /*! task that does tree roations on a list of sub-trees */
    TASK_FUNCTION(BVH4iBuilderMorton64Bit,doTreeRotationsOnSubTrees);


    /*! creates a leaf node */
    BBox3fa createSmallLeaf(SmallBuildRecord& current);

    /*! creates a leaf node */
    BBox3fa createLeaf(SmallBuildRecord& current, NodeAllocator& alloc);

    /*! fallback split mode */
    void split_fallback(SmallBuildRecord& current, SmallBuildRecord& leftChild, SmallBuildRecord& rightChild);

    /*! split a build record into two */
    bool split(SmallBuildRecord& current, SmallBuildRecord& left, SmallBuildRecord& right);

    /*! create the top-levels of the tree */
    size_t createSingleBVH4iNode(SmallBuildRecord& current, SmallBuildRecord *__restrict__ const children);

    /*! main recursive build function */
    BBox3fa recurse(SmallBuildRecord& current, 
		    NodeAllocator& alloc);
    
    /*! refit the toplevel part of the BVH */
    BBox3fa refit_toplevel(const BVH4i::NodeRef &ref);

    /*! refit the sub-BVHs */
    BBox3fa refit(const BVH4i::NodeRef &ref);


  public:



  public:
    BVH4i      * bvh;         //!< Output BVH
    Scene      * scene;

    __aligned(64) SmallBuildRecord buildRecords[MAX_TOP_LEVEL_BINS];    
    __aligned(64) unsigned int thread_startGroup[MAX_MIC_THREADS];      
    __aligned(64) unsigned int thread_startGroupOffset[MAX_MIC_THREADS];


    /*! state for radix sort */
  public:
    static const size_t RADIX_BITS = 8;
    static const size_t RADIX_BUCKETS = (1 << RADIX_BITS);
    static const size_t RADIX_BUCKETS_MASK = (RADIX_BUCKETS-1);
    __aligned(64) unsigned int radixCount[MAX_MIC_THREADS][RADIX_BUCKETS];

  protected:
    MortonID64Bit * __restrict__ morton;
    BVH4i::Node   * __restrict__ node;
    Triangle1     * __restrict__ accel;

    size_t numPrimitives;
    size_t numGroups;
    size_t numNodes;
    size_t numAllocatedNodes;
    size_t size_morton;
    size_t size_node;
    size_t size_accel;
    size_t topLevelItemThreshold;
    size_t numBuildRecords;
    size_t numPrimitivesOld;


    __aligned(64) Centroid_Scene_AABB global_bounds;

    /*! node allocator */
    __aligned(64) AlignedAtomicCounter32  atomicID;
    __aligned(64) AlignedAtomicCounter32  numBuildRecordCounter;

    __forceinline unsigned int allocGlobalNode(int size)
    {
      const unsigned int currentIndex = atomicID.add(size);
      if (unlikely(currentIndex >= numAllocatedNodes)) {
        FATAL("not enough nodes allocated");
      }
      return currentIndex;
    }

  };

}
