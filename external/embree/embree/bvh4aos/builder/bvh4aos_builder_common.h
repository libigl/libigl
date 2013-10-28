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

#ifndef __EMBREE_BVH4AOS_BUILDER_COMMON_H__
#define __EMBREE_BVH4AOS_BUILDER_COMMON_H__

#include "bvh4aos_globals.h"
#include "bvh4aos_box.h"
#include "bvh4aos_builder_util.h"
#include "bvh4aos_task_scheduler.h"


#define LOCK   ParallelQBVHBuilder::threadMutex.lock()
#define UNLOCK ParallelQBVHBuilder::threadMutex.unlock()

#define QBVH_BUILDER_LEAF_ITEM_THRESHOLD 4
#define NUM_ATOMIC_BUILD_RECORDS 512
#define LOCAL_NODE_IDS 64
#define NUM_LOCAL_ATOMIC_BUILD_RECORDS 16

#define BVH_INDEX_SHIFT 6
#define BVH_ITEMS_MASK   (((unsigned int)1 << BVH_INDEX_SHIFT)-1)
#define BVH_LEAF_MASK    ((unsigned int)1 << 31)
#define BVH_OFFSET_MASK  (~(BVH_ITEMS_MASK | BVH_LEAF_MASK))

#define QBVH_MAX_STACK_DEPTH 64
#define QBVH_INDEX_SHIFT      7
#define QBVH_LEAF_BIT_SHIFT   5
#define QBVH_ITEMS_MASK      (((unsigned int)1 << QBVH_LEAF_BIT_SHIFT)-1)
#define QBVH_LEAF_MASK       ((unsigned int)1 << QBVH_LEAF_BIT_SHIFT)
#define QBVH_OFFSET_MASK     (~(QBVH_ITEMS_MASK | QBVH_LEAF_MASK))
#define QBVH_TERMINAL_TOKEN  QBVH_LEAF_MASK


namespace embree
{

  /* ---------------------- */
  /* --- Triangle Accel --- */
  /* ---------------------- */

  __ALIGN(64)
    class TriangleAccel // === if layout changes => update computeTriangleAccel,getLeafBoundsAndComputeAccelerationData ===
    {
    public:
      float a[3];
      unsigned int id0;     
      float b[3]; // edgeAB
      unsigned int id1;
      float c[3]; // edgeAC
      unsigned int data;
      float normal[3]; // edgeAB x edgeAC
      unsigned int reserved;
    };


  _INLINE std::ostream &operator<<(std::ostream &o, const TriangleAccel &t)
    {
      o << "v0 " << t.a[0] << " " << t.a[1] << " " << t.a[2] << std::endl;
      o << "v1 " << t.b[0] << " " << t.b[1] << " " << t.b[2] << std::endl;
      o << "v2 " << t.c[0] << " " << t.c[1] << " " << t.c[2] << std::endl;
      o << "id0 " << t.id0 << std::endl;
      o << "id1 " << t.id1 << std::endl;
      return o;
    }

  class BuilderTriangle
  {
  public:
    unsigned int v0,v1,v2;   
    unsigned int id0,id1;
  };

  _INLINE std::ostream &operator<<(std::ostream &o, const BuilderTriangle &t)
    {
      o << "v0  " << t.v0 << std::endl;
      o << "v1  " << t.v1 << std::endl;
      o << "v2  " << t.v2 << std::endl;
      o << "id0 " << t.id0 << std::endl;
      o << "id1 " << t.id1 << std::endl;
      return o;
    }


  /* ----------- */
  /* --- BVH --- */
  /* ----------- */

  _INLINE unsigned int bvhItemOffset(const unsigned int children)
  {
    return (children & ~BVH_LEAF_MASK) >> BVH_INDEX_SHIFT;
  }

  _INLINE unsigned int bvhItems(const unsigned int children)
  {
    return children & BVH_ITEMS_MASK;
  }

  _INLINE unsigned int bvhChildren(const unsigned int children)
  {
    return children & BVH_ITEMS_MASK;
  }

  _INLINE unsigned int bvhChildID(const unsigned int children)
  {
    return (children & BVH_OFFSET_MASK) >> BVH_INDEX_SHIFT;
  };

  _INLINE unsigned int bvhLeaf(const unsigned int children) {
    return (children & BVH_LEAF_MASK);
  };

  /* ---------------- */
  /* --- QUAD BVH --- */
  /* ---------------- */


  class QBVHNode
  {
  public:
    struct
    {
      float x,y,z;
      unsigned int d;
    } m_min[4];

    struct
    {
      float x,y,z;
      unsigned int d;
    } m_max[4];
  };

  _INLINE unsigned int qbvhItemOffset(const unsigned int children)
  {
    return children & BVH_OFFSET_MASK; // 6 bits instead of 7
  }

  _INLINE unsigned int qbvhItemOffsetToID(const unsigned int children)
  {
    return children >> BVH_INDEX_SHIFT; // 6 bits instead of 7
  }

  _INLINE unsigned int qbvhItems(const unsigned int children)
  {
    return children & QBVH_ITEMS_MASK; // 6 bits instead of 7
  }

  _INLINE unsigned int qbvhChildID(const unsigned int node)
  {
    return (node & QBVH_OFFSET_MASK) >> QBVH_INDEX_SHIFT;
  };

  _INLINE QBVHNode *qbvhChildPtr(const QBVHNode * __restrict__ const ptr, const unsigned int node)
  {
    const unsigned int offset = node & QBVH_OFFSET_MASK;
    return (QBVHNode*)((char*)ptr + offset);
  };

  _INLINE QBVHNode *qbvhChildPtrNoMask(const QBVHNode * __restrict__ const ptr, const unsigned int node)
  {
    return (QBVHNode*)((char*)ptr + (unsigned long)node);
  };

  _INLINE unsigned int qbvhLeaf(const unsigned int node) {
    return (node & QBVH_LEAF_MASK);
  };

  _INLINE unsigned int qbvhLeaf(const unsigned int node, const unsigned int mask) {
    return (node & mask);
  };

  _INLINE unsigned int qbvhChildren(const unsigned int node) {
    return (node & QBVH_ITEMS_MASK);
  };

  _INLINE unsigned int qbvhCreateNode(const unsigned int nodeID,
				      const unsigned int children) {
    return (nodeID << QBVH_INDEX_SHIFT) | children;
  };


  _INLINE std::ostream &operator<<(std::ostream &o, const QBVHNode &v)
    {
      o << std::endl;
      for (int i=0;i<4;i++)
	{
	  o << "[" << i << "]" << std::endl;
	  o << "min [" << v.m_min[i].x << "," << v.m_min[i].y << "," << v.m_min[i].z << "," << v.m_min[i].d << "]" << std::endl;
	  o << "max [" << v.m_max[i].x << "," << v.m_max[i].y << "," << v.m_max[i].z << "," << v.m_max[i].d << "]" << std::endl;

	}
      return o;
    }

  /* ------------------ */
  /* --- Binary BVH --- */
  /* ------------------ */

  enum
  {
    MIN_X = 0,
    MIN_Y = 1,
    MIN_Z = 2,
    MAX_X = 4,
    MAX_Y = 5,
    MAX_Z = 6
  };

  class BVHNode : public AABB
  {
  public:
    _INLINE unsigned int isLeaf() const {
      return bvhLeaf(ext_min.children);
    };

    _INLINE int firstChildID() const {
      return bvhChildID(ext_min.children);
    };
    _INLINE int items() const {
      return bvhItems(ext_min.children);
    }
    _INLINE unsigned int itemListOfs() const {
      return bvhItemOffset(ext_min.children);
    }

    _INLINE unsigned int getData() const {
      return ext_max.t;
    }

    _INLINE void createLeaf(const unsigned int offset,
			    const unsigned int entries,
			    const unsigned int data = 0) {
      ext_min.children = (offset << BVH_INDEX_SHIFT) | BVH_LEAF_MASK | entries;
      ext_max.t    = data;
    }

    _INLINE void createNode(const unsigned int index,			  
			    const unsigned short children = 0,
			    const unsigned int items_subtree = 0) {
      assert((index %2) == 0);
      ext_min.children = (index << BVH_INDEX_SHIFT) | children;
      ext_max.t = items_subtree;
    }

  };

  _INLINE std::ostream &operator<<(std::ostream &o, const BVHNode &v)
    {
      if (v.isLeaf())
	{
	  o << "LEAF" << " ";
	  o << "offset " << v.itemListOfs() << " ";
	  o << "items  " << v.items() << " ";
	}
      else
	{
	  o << "NODE" << " ";
	  o << "firstChildID " << v.firstChildID() << " ";
	}  
      o << "min [" << v.m_min[0] << ", " << v.m_min[1] << ", " << v.m_min[2] << ", " << v.ext_min.t <<"] ";
      o << "max [" << v.m_max[0] << ", " << v.m_max[1] << ", " << v.m_max[2] << ", " << v.ext_max.t <<"] ";

      return o;
    } 



  /* ----------------------------------------------------------------------- */
  /* ----------------------------------------------------------------------- */
  /* ----------------------------------------------------------------------- */

  __ALIGN(8)
    struct MortonID32Bit
    {
      unsigned int code;
      unsigned int index;

      _INLINE bool operator<(const MortonID32Bit &m) const 
      { 
	return code < m.code; 
      } 

      _INLINE bool operator>(const MortonID32Bit &m) const 
      { 
	return code > m.code; 
      } 

      _INLINE unsigned char getByte(const unsigned long i) const
      {
	const unsigned char *__restrict const ptr = (const unsigned char*__restrict )&code;
	return ptr[i];
      }

      _INLINE mic_i getByteAndBroadCast(const unsigned long i) const
      {
	return upconv1i_uint8( ((unsigned char*__restrict )&code) + i );
      }

      _INLINE void operator=(const MortonID32Bit& v) { 
	*(long*)this = *(long*)&v;
      }

    };

  _INLINE std::ostream &operator<<(std::ostream &o, const MortonID32Bit &mc)
    {
      o << "index " << mc.index << " code = " << mc.code;
      return o;
    }

  /* ----------------------------------------------------------------------- */
  /* ------------------------- Centroid_Scene_AABB ------------------------- */
  /* ----------------------------------------------------------------------- */

  __ALIGN(64)
    class Centroid_Scene_AABB
    {
    public:
      static float initCentroidScene[16];

      mic_f aabb;

      _INLINE void reset() {
	aabb = upconv16f(initCentroidScene);
      }  

      _INLINE void extend(const mic_f &b_min, const mic_f &b_max, const mic_f &centroid)
      {
	aabb = mask_min(0x700,aabb,aabb,b_min);
	aabb = mask_max(0x7000,aabb,aabb,b_max);
	aabb = mask_min(0x7,aabb,aabb,centroid);
	aabb = mask_max(0x70,aabb,aabb,centroid);
      }


      _INLINE void extend_centroid(const mic_f &centroid)
      {
	aabb = mask_min(0x7,aabb,aabb,centroid);
	aabb = mask_max(0x70,aabb,aabb,centroid);
      }

      _INLINE void extend_scene(const mic_f &b_min, const mic_f &b_max)
      {
	aabb = mask_min(0x700,aabb,aabb,b_min);
	aabb = mask_max(0x7000,aabb,aabb,b_max);
      }

      _INLINE void extend(const AABB &b)
      {
	const mic_f b_min = upconv4f((float*)&b.m_min);
	const mic_f b_max = upconv4f((float*)&b.m_max);
	const mic_f centroid = (b_min+b_max) * 0.5f;
	extend(b_min,b_max,centroid);
      }


      _INLINE void extend(const Centroid_Scene_AABB &v)
      {
	aabb = mask_min(0x700,aabb,aabb,v.aabb);
	aabb = mask_max(0x7000,aabb,aabb,v.aabb);
	aabb = mask_min(0x7,aabb,aabb,v.aabb);
	aabb = mask_max(0x70,aabb,aabb,v.aabb);    
      }


      _INLINE void extend_atomic(const Centroid_Scene_AABB &v)
      {
	atomic_min(&aabb[0] ,v.aabb[0]);
	atomic_min(&aabb[1] ,v.aabb[1]);
	atomic_min(&aabb[2] ,v.aabb[2]);

	atomic_max(&aabb[4] ,v.aabb[4]);
	atomic_max(&aabb[5] ,v.aabb[5]);
	atomic_max(&aabb[6] ,v.aabb[6]);

	atomic_min(&aabb[8] ,v.aabb[8]);
	atomic_min(&aabb[9] ,v.aabb[9]);
	atomic_min(&aabb[10],v.aabb[10]);

	atomic_max(&aabb[12],v.aabb[12]);
	atomic_max(&aabb[13],v.aabb[13]);
	atomic_max(&aabb[14],v.aabb[14]);
      }


      _INLINE void extend_atomic_centroid_bounds(const Centroid_Scene_AABB &v)
      {
	atomic_min(&aabb[0] ,v.aabb[0]);
	atomic_min(&aabb[1] ,v.aabb[1]);
	atomic_min(&aabb[2] ,v.aabb[2]);

	atomic_max(&aabb[4] ,v.aabb[4]);
	atomic_max(&aabb[5] ,v.aabb[5]);
	atomic_max(&aabb[6] ,v.aabb[6]);
      }


      _INLINE mic_f minCentroidAABB() {
	return upconv4f(&aabb[0]);
      }

      _INLINE mic_f maxCentroidAABB() {
	return upconv4f(&aabb[4]);
      }

      _INLINE mic_f minSceneAABB() {
	return upconv4f(&aabb[8]);
      }

      _INLINE mic_f maxSceneAABB() {
	return upconv4f(&aabb[12]);
      }


      _INLINE bool encloseCentroid(const mic_f& c)
      {
	const mic_m m_min = ge(0x7777,c,minCentroidAABB());
	const mic_m m_max = le(0x7777,c,maxCentroidAABB());
	return (m_max & m_min) == 0x7777;
      }

      _INLINE void setMinCentroidAABB(const mic_f &v)
      {
	store4f(&aabb[0],v);
      }

      _INLINE void setMaxCentroidAABB(const mic_f &v)
      {
	store4f(&aabb[4],v);
      }

      _INLINE void setMinSceneAABB(const mic_f &v)
      {
	store4f(&aabb[8],v);
      }

      _INLINE void setMaxSceneAABB(const mic_f &v)
      {
	store4f(&aabb[12],v);
      }


      _INLINE mic_f centroidDiagonal() {
	return maxCentroidAABB() - minCentroidAABB();
      }

      _INLINE mic_f sceneDiagonal() {
	return maxSceneAABB() - minSceneAABB();
      }

      _INLINE Vec3fa &sceneAABB_min() {
	return *(Vec3fa*)&aabb[8];
      };

      _INLINE Vec3fa &sceneAABB_max() {
	return *(Vec3fa*)&aabb[12];
      };

      _INLINE Vec3fa &centroidAABB_min() {
	return *(Vec3fa*)&aabb[0];
      };

      _INLINE Vec3fa &centroidAABB_max() {
	return *(Vec3fa*)&aabb[4];
      };

      _INLINE unsigned int getMaxDimCentroidDiagonal() 
      {
	const mic_f diag = sel(0x7777,centroidDiagonal(),mic_f::zero());
	const mic_f diag_max = set_max4(diag);
	return bsf64(eq(diag_max,diag));
      }

      _INLINE float sceneArea() {
	const Vec3fa d = sceneAABB_max() - sceneAABB_min();
	return 2.0f * (d.x*d.y+d.x*d.z+d.y*d.z);
      };

      _INLINE float centroidArea() {
	const Vec3fa d = centroidAABB_max() - centroidAABB_min();
	return 2.0f * (d.x*d.y+d.x*d.z+d.y*d.z);
      };

      _INLINE void storeSceneAABB(void *s) const
      {
	compactustore16f(0xff00,(float*)s,aabb);
      }

      _INLINE void storeCentroidAABB(void *s) const
      {
	compactustore16f(0x00ff,(float*)s,aabb);
      }

      _INLINE unsigned int add_atomic(const unsigned int index, const unsigned int add)
      {
	return atomic_add(((atomic_t*)&aabb) + index, add);
      }

      _INLINE unsigned int extract_uint(const unsigned int index)
      {
	return ((unsigned int*)&aabb)[index];
      }

    };

  _INLINE std::ostream &operator<<(std::ostream &o, const Centroid_Scene_AABB &v)
    {
      o << "centroid: ";
      o << "min [" << v.aabb[0] << ", " << v.aabb[1] << ", " << v.aabb[2] << "] ";
      o << "max [" << v.aabb[4] << ", " << v.aabb[5] << ", " << v.aabb[6] << "] ";
      o << "scene: ";
      o << "min [" << v.aabb[8]  << ", " << v.aabb[9]  << ", " << v.aabb[10] << "] ";
      o << "max [" << v.aabb[12] << ", " << v.aabb[13] << ", " << v.aabb[14] << "] ";
      return o;
    } 

  _INLINE bool operator==(const Centroid_Scene_AABB &a, const Centroid_Scene_AABB &b) 
  {
    return eq(0x7777,a.aabb,b.aabb) == 0x7777;
  }

  /* ----------------------------------------------------------------------- */
  /* ----------------------------------------------------------------------- */
  /* ----------------------------------------------------------------------- */


  struct BestSplitData {
    int index;
    int dim;
    int left;
    float cost;

    _INLINE void init(const float c)
    {
      index = -1;
      dim   = -1;
      left  = -1;
      cost  =  c;
    }
  };

  _INLINE std::ostream &operator<<(std::ostream &o, const BestSplitData &b)
    {
      o << "index " << b.index << std::endl;
      o << "dim   " << b.dim   << std::endl;
      o << "left  " << b.left  << std::endl;
      o << "cost  " << b.cost  << std::endl;
      return o;
    }

  /* --------------------------------------------------------------- */
  /* ------------------------- BuildRecord ------------------------- */
  /* --------------------------------------------------------------- */


  enum {
    BUILD_RECORD_INIT  = 0,
    BUILD_RECORD_NODE  = 1,
    BUILD_RECORD_LEAF  = 2
  };

  __ALIGN(64)
    class BuildRecord {
  public:
    Centroid_Scene_AABB cs_AABB;
    unsigned int begin;
    unsigned int end;
    unsigned int childrenID;
    unsigned int nodeID;

    volatile unsigned int flags;
    float sArea;
    int splitIndex;
    int splitNumLeft;
    
    unsigned int sourceID;
    unsigned int data;
    char dummy[8+16];

    BuildRecord() { assert(sizeof(BuildRecord) == 128); }


    _INLINE void prefetchL1Ex() const {
      mic_f *p = (mic_f*)this;
      prefetch<PFHINT_L1EX>(p+0);
      prefetch<PFHINT_L1EX>(p+1);
    }

    _INLINE void prefetchL1() const {
      mic_f *p = (mic_f*)this;
      prefetch<PFHINT_L1>(p+0);
      prefetch<PFHINT_L1>(p+1);
    }

    _INLINE void prefetchL2Ex() const {
      mic_f *p = (mic_f*)this;
      prefetch<PFHINT_L2EX>(p+0);
      prefetch<PFHINT_L2EX>(p+1);
    }

    _INLINE void prefetchL2() const {
      mic_f *p = (mic_f*)this;
      prefetch<PFHINT_L2>(p+0);
      prefetch<PFHINT_L2>(p+1);
    }

    _INLINE void initFLags() {
      flags = BUILD_RECORD_INIT;
    }

    _INLINE void createNode() {
      flags = BUILD_RECORD_NODE;
    }


    _INLINE void createLeaf() {
      flags = BUILD_RECORD_LEAF;
    }

    _INLINE bool isLeaf() {
      return flags == BUILD_RECORD_LEAF;
    }

    _INLINE float sceneArea() {
      return sArea;
    }

    _INLINE float centroidArea() {
      return cs_AABB.centroidArea();
    }

    _INLINE unsigned int items() const {
      return end - begin;
    }

    _INLINE void init(const Centroid_Scene_AABB & _cs_AABB, 
		      const unsigned int _begin, 
		      const unsigned int _end)
    {
      cs_AABB        = _cs_AABB;
      begin          = _begin;
      end            = _end;
      childrenID     = (unsigned int)-1;
      nodeID         = (unsigned int)-1;
      sArea          = cs_AABB.sceneArea();
      flags          = BUILD_RECORD_NODE;
    }

    _INLINE void init(const Centroid_Scene_AABB & _cs_AABB, 
		      const unsigned int _begin, 
		      const unsigned int _end,
		      const unsigned int _sourceID)			 
    {
      cs_AABB        = _cs_AABB;
      begin          = _begin;
      end            = _end;
      childrenID     = (unsigned int)-1;
      nodeID         = (unsigned int)-1;
      flags          = BUILD_RECORD_NODE;
      sourceID       = _sourceID;
    }


    _INLINE void init(const unsigned int _begin, 
		      const unsigned int _end,
		      const unsigned int _sourceID)			 
    {
      begin          = _begin;
      end            = _end;
      childrenID     = (unsigned int)-1;
      nodeID         = (unsigned int)-1;
      flags          = BUILD_RECORD_NODE;
      sourceID       = _sourceID;
    }


    _INLINE void combine(BuildRecord &right)
    {
      assert(right.begin == end);
      end = right.end;
      cs_AABB.extend(right.cs_AABB);
    }

    _INLINE bool operator<(const BuildRecord &br) const 
    { return items() < br.items(); } 
   
    _INLINE bool operator>(const BuildRecord &br) const 
    { return items() > br.items(); } 

    _INLINE void operator=(const BuildRecord& v) { 
      const mic_f b0 = upconv16f((float*)&v);
      const mic_f b1 = upconv16f((float*)&v + 16);
      store16f((float*)this +  0, b0);
      store16f((float*)this + 16, b1);
    };

  };

  _INLINE bool less_begin(const BuildRecord &a, const BuildRecord &b) 
  { return a.begin < b.begin; } 


  _INLINE std::ostream &operator<<(std::ostream &o, const BuildRecord &br)
    {
      o << br.cs_AABB << std::endl;
      o << "begin      " << br.begin << std::endl;
      o << "end        " << br.end << std::endl;
      o << "items      " << br.end-br.begin << std::endl;
      o << "childrenID " << br.childrenID << std::endl;
      o << "nodeID     " << br.nodeID << std::endl;
      o << "flags      " << br.flags << std::endl;
      o << "sArea      " << br.sArea << std::endl;
      return o;
    };


  /* ----------------------------------------------------------------------- */
  /* ------------------------- Parallel QVHBuilder ------------------------- */
  /* ----------------------------------------------------------------------- */

  static _INLINE mic_m lt_split(const AABB *__restrict__ const aabb,
				const unsigned int dim,
				const mic_f &c,
				const mic_f &s,
				const mic_f &bestSplit_f)
  {
    const mic_f b_min = mic_f(aabb->m_min[dim]);
    const mic_f b_max = mic_f(aabb->m_max[dim]);
    prefetch<PFHINT_NT>(aabb + 2);
    const mic_f centroid_2 = b_min + b_max;
    const mic_f binID = (centroid_2 - c)*s;
    return lt(binID,bestSplit_f);    
  }


  static _INLINE mic_m ge_split(const AABB *__restrict__ const aabb,
				const unsigned int dim,
				const mic_f &c,
				const mic_f &s,
				const mic_f &bestSplit_f)
  {
    const mic_f b_min = mic_f(aabb->m_min[dim]);
    const mic_f b_max = mic_f(aabb->m_max[dim]);
    prefetch<PFHINT_NT>(aabb-2);
    const mic_f centroid_2 = b_min + b_max;
    const mic_f binID = (centroid_2 - c)*s;
    return ge(binID,bestSplit_f);    
  }

  /* ----------------------------------------------------------------------- */
  /* ----------------------------------------------------------------------- */
  /* ----------------------------------------------------------------------- */

  __ALIGN(64)
    class Bin16
    {
    public:
      mic_f min_x[3];
      mic_f min_y[3];
      mic_f min_z[3];
      mic_f max_x[3];
      mic_f max_y[3];
      mic_f max_z[3];
      mic_i count[3];
      mic_i thread_count[3];

      //mic_i dummy[9+2];

      Bin16() {}

      _INLINE void prefetchL1()
      {
#pragma unroll
	for (unsigned int i=0;i<sizeof(Bin16);i+=64)
	  prefetch<PFHINT_L1>((char*)this + i);
      }

      _INLINE void prefetchL1EX()
      {
	prefetch<PFHINT_L1EX>(&min_x[0]);
	prefetch<PFHINT_L1EX>(&min_x[1]);
	prefetch<PFHINT_L1EX>(&min_x[2]);

	prefetch<PFHINT_L1EX>(&min_y[0]);
	prefetch<PFHINT_L1EX>(&min_y[1]);
	prefetch<PFHINT_L1EX>(&min_y[2]);

	prefetch<PFHINT_L1EX>(&min_z[0]);
	prefetch<PFHINT_L1EX>(&min_z[1]);
	prefetch<PFHINT_L1EX>(&min_z[2]);

	prefetch<PFHINT_L1EX>(&max_x[0]);
	prefetch<PFHINT_L1EX>(&max_x[1]);
	prefetch<PFHINT_L1EX>(&max_x[2]);

	prefetch<PFHINT_L1EX>(&max_y[0]);
	prefetch<PFHINT_L1EX>(&max_y[1]);
	prefetch<PFHINT_L1EX>(&max_y[2]);

	prefetch<PFHINT_L1EX>(&max_z[0]);
	prefetch<PFHINT_L1EX>(&max_z[1]);
	prefetch<PFHINT_L1EX>(&max_z[2]);

	prefetch<PFHINT_L1EX>(&count[0]);
	prefetch<PFHINT_L1EX>(&count[1]);
	prefetch<PFHINT_L1EX>(&count[2]);
      }


      _INLINE void reset()
      {
	const mic_f init_min = mic_f::inf();
	const mic_f init_max = mic_f::minus_inf();
	const mic_i zero     = mic_i::zero();

	min_x[0] = init_min;
	min_x[1] = init_min;
	min_x[2] = init_min;

	min_y[0] = init_min;
	min_y[1] = init_min;
	min_y[2] = init_min;

	min_z[0] = init_min;
	min_z[1] = init_min;
	min_z[2] = init_min;

	max_x[0] = init_max;
	max_x[1] = init_max;
	max_x[2] = init_max;

	max_y[0] = init_max;
	max_y[1] = init_max;
	max_y[2] = init_max;

	max_z[0] = init_max;
	max_z[1] = init_max;
	max_z[2] = init_max;

	count[0] = zero;
	count[1] = zero;
	count[2] = zero;
      }


      _INLINE void merge(const Bin16& b)
      {
#pragma unroll(3)
	for (unsigned int i=0;i<3;i++)
	  {
	    min_x[i] = _min(min_x[i],b.min_x[i]);
	    min_y[i] = _min(min_y[i],b.min_y[i]);
	    min_z[i] = _min(min_z[i],b.min_z[i]);

	    max_x[i] = _max(max_x[i],b.max_x[i]);
	    max_y[i] = _max(max_y[i],b.max_y[i]);
	    max_z[i] = _max(max_z[i],b.max_z[i]);

	    count[i] += b.count[i];
	  }

      } 

      _INLINE void merge_with_prefetch(const Bin16& b)
      {

	for (unsigned int i=0;i<3;i++)
	  {
	    prefetch<PFHINT_L2>(&b.min_x[i]);  
	    prefetch<PFHINT_L2>(&b.min_y[i]);  
	    prefetch<PFHINT_L2>(&b.min_z[i]);  

	    prefetch<PFHINT_L2>(&b.max_x[i]);  
	    prefetch<PFHINT_L2>(&b.max_y[i]);  
	    prefetch<PFHINT_L2>(&b.max_z[i]);  
	    prefetch<PFHINT_L2>(&b.count[i]);  
	  }

#pragma unroll(3)
	for (unsigned int i=0;i<3;i++)
	  {
	    prefetch<PFHINT_NT>(&b.min_x[i]);  
	    prefetch<PFHINT_NT>(&b.min_y[i]);  
	    prefetch<PFHINT_NT>(&b.min_z[i]);  

	    prefetch<PFHINT_NT>(&b.max_x[i]);  
	    prefetch<PFHINT_NT>(&b.max_y[i]);  
	    prefetch<PFHINT_NT>(&b.max_z[i]);  
	    prefetch<PFHINT_NT>(&b.count[i]);  

	    min_x[i] = _min(min_x[i],b.min_x[i]);
	    min_y[i] = _min(min_y[i],b.min_y[i]);
	    min_z[i] = _min(min_z[i],b.min_z[i]);

	    max_x[i] = _max(max_x[i],b.max_x[i]);
	    max_y[i] = _max(max_y[i],b.max_y[i]);
	    max_z[i] = _max(max_z[i],b.max_z[i]);
	    count[i] += b.count[i];
	  }

      } 

      _INLINE mic_f prefix_area_rl(const unsigned int i)
      {
	prefetch<PFHINT_NT>(&min_x[i]);  
	prefetch<PFHINT_NT>(&min_y[i]);  
	prefetch<PFHINT_NT>(&min_z[i]);  

	prefetch<PFHINT_NT>(&max_x[i]);  
	prefetch<PFHINT_NT>(&max_y[i]);  
	prefetch<PFHINT_NT>(&max_z[i]);  

	const mic_f r_min_x = prefix_min(reverse(min_x[i]));
	const mic_f r_min_y = prefix_min(reverse(min_y[i]));
	const mic_f r_min_z = prefix_min(reverse(min_z[i]));
	const mic_f r_max_x = prefix_max(reverse(max_x[i]));
	const mic_f r_max_y = prefix_max(reverse(max_y[i]));
	const mic_f r_max_z = prefix_max(reverse(max_z[i]));

	const mic_f dx = r_max_x - r_min_x;
	const mic_f dy = r_max_y - r_min_y;
	const mic_f dz = r_max_z - r_min_z;
   
	const mic_f area_rl = (dx*dy+dx*dz+dy*dz) * mic_f::const2();
	return reverse(shl1_zero_extend(area_rl));
	//return reverse(area_rl);
      }

      static _INLINE mic_f prefix_area_rl(const mic_f min_x,
					  const mic_f min_y,
					  const mic_f min_z,
					  const mic_f max_x,
					  const mic_f max_y,
					  const mic_f max_z)
      {
	const mic_f r_min_x = prefix_min(reverse(min_x));
	const mic_f r_min_y = prefix_min(reverse(min_y));
	const mic_f r_min_z = prefix_min(reverse(min_z));
	const mic_f r_max_x = prefix_max(reverse(max_x));
	const mic_f r_max_y = prefix_max(reverse(max_y));
	const mic_f r_max_z = prefix_max(reverse(max_z));

	const mic_f dx = r_max_x - r_min_x;
	const mic_f dy = r_max_y - r_min_y;
	const mic_f dz = r_max_z - r_min_z;
   
	const mic_f area_rl = (dx*dy+dx*dz+dy*dz) * mic_f::const2();
	return reverse(shl1_zero_extend(area_rl));
      }


      _INLINE mic_f prefix_area_lr(const unsigned int i)
      {
	const mic_f r_min_x = prefix_min(min_x[i]);
	const mic_f r_min_y = prefix_min(min_y[i]);
	const mic_f r_min_z = prefix_min(min_z[i]);
	const mic_f r_max_x = prefix_max(max_x[i]);
	const mic_f r_max_y = prefix_max(max_y[i]);
	const mic_f r_max_z = prefix_max(max_z[i]);

	const mic_f dx = r_max_x - r_min_x;
	const mic_f dy = r_max_y - r_min_y;
	const mic_f dz = r_max_z - r_min_z;
  
	const mic_f area_lr = (dx*dy+dx*dz+dy*dz) * mic_f::const2();
	return area_lr;
      }

      static _INLINE mic_f prefix_area_lr(const mic_f min_x,
					  const mic_f min_y,
					  const mic_f min_z,
					  const mic_f max_x,
					  const mic_f max_y,
					  const mic_f max_z)
      {
	const mic_f r_min_x = prefix_min(min_x);
	const mic_f r_min_y = prefix_min(min_y);
	const mic_f r_min_z = prefix_min(min_z);
	const mic_f r_max_x = prefix_max(max_x);
	const mic_f r_max_y = prefix_max(max_y);
	const mic_f r_max_z = prefix_max(max_z);

	const mic_f dx = r_max_x - r_min_x;
	const mic_f dy = r_max_y - r_min_y;
	const mic_f dz = r_max_z - r_min_z;
  
	const mic_f area_lr = (dx*dy+dx*dz+dy*dz) * mic_f::const2();
	return area_lr;
      }


      _INLINE mic_i prefix_count(const unsigned int i)
      {
	return prefix_sum(count[i]);
      }

      static _INLINE mic_i prefix_count(const mic_i c)
      {
	return prefix_sum(c);
      }

    };

  _INLINE bool operator==(const Bin16 &a, const Bin16 &b) { 
    mic_m mask = MIC_M_ALL;
#pragma unroll(3)
    for (unsigned int i=0;i<3;i++)
      {
	mask &= eq(a.min_x[i],b.min_x[i]);
	mask &= eq(a.min_y[i],b.min_y[i]);
	mask &= eq(a.min_z[i],b.min_z[i]);

	mask &= eq(a.max_x[i],b.max_x[i]);
	mask &= eq(a.max_y[i],b.max_y[i]);
	mask &= eq(a.max_z[i],b.max_z[i]);

	mask &= eq(a.count[i],b.count[i]);
      }
    return mask == MIC_M_ALL;
  };

  _INLINE bool operator!=(const Bin16 &a, const Bin16 &b) { 
    return !(a==b);
  }

  // =====================================================================================
  // =====================================================================================
  // =====================================================================================

  _INLINE void fastbin(const AABB * __restrict__ const aabb,
		       const unsigned int start,
		       const unsigned int end,
		       const mic_f &centroidBoundsMin_2,
		       const mic_f &scale,
		       mic_f lArea[3],
		       mic_f rArea[3],
		       mic_i lNum[3])
  {
    const AABB * __restrict__ aptr = aabb + start;

    prefetch<PFHINT_NT>(aptr);
    prefetch<PFHINT_L2>(aptr+2);
    prefetch<PFHINT_L2>(aptr+4);
    prefetch<PFHINT_L2>(aptr+6);
    prefetch<PFHINT_L2>(aptr+8);
    prefetch<PFHINT_L2>(aptr+10);

    const mic_f init_min = mic_f::inf();
    const mic_f init_max = mic_f::minus_inf();
    const mic_i zero     = mic_i::zero();

    mic_f min_x0,min_x1,min_x2;
    mic_f min_y0,min_y1,min_y2;
    mic_f min_z0,min_z1,min_z2;
    mic_f max_x0,max_x1,max_x2;
    mic_f max_y0,max_y1,max_y2;
    mic_f max_z0,max_z1,max_z2;
    mic_i count0,count1,count2;

    min_x0 = init_min;
    min_x1 = init_min;
    min_x2 = init_min;
    min_y0 = init_min;
    min_y1 = init_min;
    min_y2 = init_min;
    min_z0 = init_min;
    min_z1 = init_min;
    min_z2 = init_min;

    max_x0 = init_max;
    max_x1 = init_max;
    max_x2 = init_max;
    max_y0 = init_max;
    max_y1 = init_max;
    max_y2 = init_max;
    max_z0 = init_max;
    max_z1 = init_max;
    max_z2 = init_max;

    count0 = zero;
    count1 = zero;
    count2 = zero;

    for (unsigned int j = start;j < end;j++,aptr++)
      {
	prefetch<PFHINT_NT>(aptr+2);
	prefetch<PFHINT_L2>(aptr+8);

	const mic_f b_min = upconv4f((float*)aptr->m_min);
	const mic_f b_max = upconv4f((float*)aptr->m_max);    

	const mic_f centroid_2 = b_min + b_max;
	const mic_i binID = mic_i((centroid_2 - centroidBoundsMin_2)*scale);

	assert(0 <= binID[0] && binID[0] < 16);
	assert(0 <= binID[1] && binID[1] < 16);
	assert(0 <= binID[2] && binID[2] < 16);

	const mic_i id = mic_i::identity();
	const mic_m m_update_x = eq(id,swAAAA_i(binID));
	const mic_m m_update_y = eq(id,swBBBB_i(binID));
	const mic_m m_update_z = eq(id,swCCCC_i(binID));

	min_x0 = mask_min(m_update_x,min_x0,min_x0,swAAAA(b_min));
	min_y0 = mask_min(m_update_x,min_y0,min_y0,swBBBB(b_min));
	min_z0 = mask_min(m_update_x,min_z0,min_z0,swCCCC(b_min));
	// ------------------------------------------------------------------------      
	max_x0 = mask_max(m_update_x,max_x0,max_x0,swAAAA(b_max));
	max_y0 = mask_max(m_update_x,max_y0,max_y0,swBBBB(b_max));
	max_z0 = mask_max(m_update_x,max_z0,max_z0,swCCCC(b_max));
	// ------------------------------------------------------------------------
	min_x1 = mask_min(m_update_y,min_x1,min_x1,swAAAA(b_min));
	min_y1 = mask_min(m_update_y,min_y1,min_y1,swBBBB(b_min));
	min_z1 = mask_min(m_update_y,min_z1,min_z1,swCCCC(b_min));      
	// ------------------------------------------------------------------------      
	max_x1 = mask_max(m_update_y,max_x1,max_x1,swAAAA(b_max));
	max_y1 = mask_max(m_update_y,max_y1,max_y1,swBBBB(b_max));
	max_z1 = mask_max(m_update_y,max_z1,max_z1,swCCCC(b_max));
	// ------------------------------------------------------------------------
	min_x2 = mask_min(m_update_z,min_x2,min_x2,swAAAA(b_min));
	min_y2 = mask_min(m_update_z,min_y2,min_y2,swBBBB(b_min));
	min_z2 = mask_min(m_update_z,min_z2,min_z2,swCCCC(b_min));
	// ------------------------------------------------------------------------      
	max_x2 = mask_max(m_update_z,max_x2,max_x2,swAAAA(b_max));
	max_y2 = mask_max(m_update_z,max_y2,max_y2,swBBBB(b_max));
	max_z2 = mask_max(m_update_z,max_z2,max_z2,swCCCC(b_max));
	// ------------------------------------------------------------------------
	count0 = mask_add(m_update_x,count0,count0,mic_i::one());
	count1 = mask_add(m_update_y,count1,count1,mic_i::one());
	count2 = mask_add(m_update_z,count2,count2,mic_i::one());      
	//evictL1(aptr-2);
      }
    prefetch<PFHINT_L1EX>(&rArea[0]);
    prefetch<PFHINT_L1EX>(&lArea[0]);
    prefetch<PFHINT_L1EX>(&lNum[0]);
    rArea[0] = Bin16::prefix_area_rl(min_x0,min_y0,min_z0,max_x0,max_y0,max_z0);
    lArea[0] = Bin16::prefix_area_lr(min_x0,min_y0,min_z0,max_x0,max_y0,max_z0);
    lNum[0]  = Bin16::prefix_count(count0);

    prefetch<PFHINT_L1EX>(&rArea[1]);
    prefetch<PFHINT_L1EX>(&lArea[1]);
    prefetch<PFHINT_L1EX>(&lNum[1]);
    rArea[1] = Bin16::prefix_area_rl(min_x1,min_y1,min_z1,max_x1,max_y1,max_z1);
    lArea[1] = Bin16::prefix_area_lr(min_x1,min_y1,min_z1,max_x1,max_y1,max_z1);
    lNum[1]  = Bin16::prefix_count(count1);

    prefetch<PFHINT_L1EX>(&rArea[2]);
    prefetch<PFHINT_L1EX>(&lArea[2]);
    prefetch<PFHINT_L1EX>(&lNum[2]);
    rArea[2] = Bin16::prefix_area_rl(min_x2,min_y2,min_z2,max_x2,max_y2,max_z2);
    lArea[2] = Bin16::prefix_area_lr(min_x2,min_y2,min_z2,max_x2,max_y2,max_z2);
    lNum[2]  = Bin16::prefix_count(count2);

  }


  _INLINE std::ostream &operator<<(std::ostream &o, const Bin16 &v)
    {
#pragma unroll(3)
      for (unsigned int i=0;i<3;i++)
	{
	  DBG_PRINT(v.min_x[i]);
	  DBG_PRINT(v.min_y[i]);
	  DBG_PRINT(v.min_z[i]);

	  DBG_PRINT(v.max_x[i]);
	  DBG_PRINT(v.max_y[i]);
	  DBG_PRINT(v.max_z[i]);

	  DBG_PRINT(v.count[i]);
	}

      return o;
    }

  // =====================================================================================
  // =====================================================================================
  // =====================================================================================


  _INLINE mic_f computeTriangleAccel(const unsigned int primID,				  
				     const BuilderTriangle * __restrict__ const tptr,
				     const Vec3fa          * __restrict__ const vptr)
  {
    const unsigned int i0 = tptr[primID].v0;
    const unsigned int i1 = tptr[primID].v1;
    const unsigned int i2 = tptr[primID].v2;

    prefetch<PFHINT_L1>(&vptr[i0]);
    prefetch<PFHINT_L1>(&vptr[i1]);
    prefetch<PFHINT_L1>(&vptr[i2]);

    const mic_f vtxA = upconv4f((float*)&vptr[i0]);
    const mic_f vtxB = upconv4f((float*)&vptr[i1]);
    const mic_f vtxC = upconv4f((float*)&vptr[i2]);
    const mic_f e1 = vtxA - vtxB;
    const mic_f e2 = vtxC - vtxA;	     
    const mic_f normal = lcross_xyz(e1,e2);
    const mic_f v0 = sel(0x8888,cast_to_mic_f(mic_i(primID)),vtxA);
    const mic_f v1 = sel(0x8888,cast_to_mic_f(mic_i(tptr[primID].id0)),e1);
    const mic_f v2 = sel(0x8888,cast_to_mic_f(mic_i(tptr[primID].id1)),e2);
    const mic_f v3 = sel(0x8888,cast_to_mic_f(0x0),normal);
    return lane_shuffle_gather<0>(v0,v1,v2,v3);
  }



  class ParallelQBVHBuilder : public QBVHTaskScheduler
  {
  public:

    enum {
      RECURSE     = 1,
      PUSH_LOCAL  = 2,
      PUSH_GLOBAL = 3
    };


    _INLINE static unsigned int getAtomicID(const unsigned int mode, 
					    unsigned int &localNodeID,
					    unsigned int &localNodeIDs,
					    const unsigned int nodes = 4)
    {
      assert(nodes % 4 == 0);
      if (unlikely(mode == PUSH_GLOBAL))
	{
	  const unsigned int currentIndex = ParallelQBVHBuilder::atomicID.add(nodes);
	  if (unlikely(currentIndex >= max_nodes))
	    {
	      LOCK;
	      std::cout << "not enough BVH nodes allocated: " << currentIndex << " <-> " << max_nodes << std::endl << std::flush;
	      FATAL("BVH memory pre-allocation error");
	      UNLOCK;
	    }
	  return currentIndex;
	}
      else
	{
	  const unsigned int currentIndex = localNodeID + localNodeIDs;

	  if (unlikely(currentIndex >= max_nodes))
	    {
	      LOCK;
	      std::cout << "not enough BVH nodes allocated: " << currentIndex << " <-> " << max_nodes << std::endl << std::flush;
	      FATAL("BVH memory pre-allocation error");
	      UNLOCK;
	    }

	  localNodeIDs += nodes;
	  if (unlikely(localNodeIDs == LOCAL_NODE_IDS))
	    {
	      localNodeID = ParallelQBVHBuilder::atomicID.add(LOCAL_NODE_IDS);
	      localNodeIDs = 0;	
	    }
	  return currentIndex;
	}   
    }

    static void fastbinning(const unsigned int start, const unsigned int end, Bin16 &bin16,const unsigned int threadID);

    static void thread_createTriangleAccel(const unsigned int threadID);
    static void thread_convertToOptimizedQBVH(const unsigned int threadID);
    static void thread_createAABBs(const unsigned int threadID);

    static void thread_build_local(const unsigned int threadID);
    static void thread_binning(const unsigned int threadID);
    static void thread_partition_parallel_tmp_copy(const unsigned int threadID);
  
    static void computeSAHfromReduceBins();
  
    static unsigned int (* recursePtr)(BuildRecord &current, 
				       const unsigned int mode, 
				       const unsigned int threadID,
				       unsigned int &localNodeID,
				       unsigned int &localNodeIDs);


    static unsigned int recurseSAH(BuildRecord &current, 
				   const unsigned int mode, 
				   const unsigned int threadID,
				   unsigned int &localNodeID,
				   unsigned int &localNodeIDs);


    static bool split(BuildRecord &current,
		      BuildRecord &leftChild,
		      BuildRecord &rightChild);

    static bool split_parallel(BuildRecord &current,
			       BuildRecord &leftChild,
			       BuildRecord &rightChild);

    static unsigned int partition_parallel(const BuildRecord &current,
					   const unsigned int bestSplit,
					   const unsigned int bestSplitDim,
					   const mic_f &centroidBoundsMin_2,
					   const mic_f &scale,
					   Centroid_Scene_AABB &left,
					   Centroid_Scene_AABB &right);

    static void reduceBin16(const unsigned int currentThreadID,
			    const unsigned int childThreadID);

    struct SplitData {
      BuildRecord rec;
      Centroid_Scene_AABB left;
      Centroid_Scene_AABB right;
      mic_f centroidBoundsMin_2;
      mic_f scale;

      unsigned int worker_threads;
      unsigned int work_items_per_thread;
      unsigned int bestSplit;
      unsigned int bestSplitDim;

      _INLINE void prefetchL2()
      {
#pragma unroll
	for (unsigned int i=0;i<sizeof(SplitData);i+=64)
	  prefetch<PFHINT_L2>(((char*)this) + 64 * i);
      }

    };
  
    static Bin16 threadBin16[MAX_MIC_THREADS];

    static AABB          * __restrict__ aabb;
    static AABB          * __restrict__ tmp_aabb;

    static BVHNode         * __restrict__ globalNodePtr;
    static Vec3fa          * __restrict__ globalVertexPtr;  
    static BuilderTriangle * __restrict__ globalBuildTrianglePtr;
    static TriangleAccel   * __restrict__ globalAccelPtr;

    static unsigned int numActiveLocalWorkQueues;
  
    static size_t triangles;
    static size_t vertices;
    static size_t qbvh_nodes;
    static size_t max_nodes;
    static bool enableTaskStealing;

    static AABB initQBVHNode[2];
    static Centroid_Scene_AABB global_cs_AABB;
    static AtomicMutex threadMutex;
    static AtomicCounter atomicID;
    static WorkStack<BuildRecord,NUM_ATOMIC_BUILD_RECORDS> atomicBuildRecordStack;
    static SplitData splitData;
    static WorkStack<BuildRecord,NUM_LOCAL_ATOMIC_BUILD_RECORDS> localAtomicBuildRecordStack[MAX_MIC_CORES];
 


    static _INLINE unsigned int getBinSlot(const AABB *__restrict__ const aabb,
					   const unsigned int dim,
					   const mic_f &c,
					   const mic_f &s)
    {
      const mic_f b_min = upconv4f((float*)&aabb->m_min);
      const mic_f b_max = upconv4f((float*)&aabb->m_max);
      const mic_f centroid_2 = b_min + b_max;
      const mic_i binID = mic_i((centroid_2 - c)*s);
      return binID[dim];    
    }  


    static _INLINE void pushLargestElementToTop()
    {
      int recordID = -1;
      int items = -1;
      for (unsigned int i=0;i<atomicBuildRecordStack.usedSlots();i++)
	if ((int)atomicBuildRecordStack.t[i].items() > items)
	  {
	    mic_f *ptr = (mic_f*)&atomicBuildRecordStack.t[i+1];
	    prefetch<PFHINT_L1EX>(ptr + 0);
	    prefetch<PFHINT_L1EX>(ptr + 1);
	    recordID = i;
	    items = atomicBuildRecordStack.t[i].items();	 
	  }
      if (recordID != -1)
	{
	  MIC_ALIGN BuildRecord tmp;
	  tmp = atomicBuildRecordStack.t[recordID];
	  atomicBuildRecordStack.t[recordID] = atomicBuildRecordStack.t[atomicBuildRecordStack.usedSlots()-1];
	  atomicBuildRecordStack.t[atomicBuildRecordStack.usedSlots()-1] = tmp;
	}
    }



    static unsigned int handleNonSplitableLeaves(BuildRecord &current, BuildRecord record[4]);

    static void transferBuildRecordsFromGlobalToLocalStack();


    static bool checkRange(const unsigned int begin,
			   const unsigned int end,
			   const unsigned int side);

    static void checkBuildRecord(const BuildRecord &current);

    // ==================================
    // ==================================
    // ==================================

  public:

    static void build();
    static void convertToOptimizedQBVH();

    _INLINE static TriangleAccel *getTriangleAccelPtr() { return globalAccelPtr; }
    _INLINE static Bin16 *getThreadBin16Ptr(const unsigned int index) { return &threadBin16[index]; }
    _INLINE static size_t getQBVHNodes() { return qbvh_nodes; }
    _INLINE static size_t getMaxAllocatedNodes() { return max_nodes; }
    static void check_tree(const bool checkPrimBounds = false);


  };


  __ALIGN(16)
    struct RadixNode {
      unsigned int start;
      unsigned int end;
      unsigned int leftIsLeaf  : 1;
      unsigned int leftChild   : 31;
      unsigned int rightIsLeaf : 1;
      unsigned int rightChild  : 31;

      _INLINE unsigned int items()
      {
	return end-start+1;
      }

      _INLINE bool contains(const unsigned int ID)
      {
	return ID >= start && ID <= end;
      }

      _INLINE bool isLeaf()
      {
	return leftChild == 0;
      }
    };


  _INLINE std::ostream& operator<<(std::ostream& o, const RadixNode& n) {
    o << " start " << n.start 
      << " end " << n.end 
      << " leftIsLeaf " << n.leftIsLeaf
      << " leftChild " << n.leftChild
      << " rightIsLeaf " << n.rightIsLeaf
      << " rightChild " << n.rightChild;
    return o;
  }

  class ParallelQBVHBuilderMortonCode : public ParallelQBVHBuilder
  {
  public:

    enum {
      EXTRACT_SUBTREES,
      REFIT_FROM_SUBTREES
    };


    static void thread_refit_qbvh(const unsigned int threadID); 
    static void refit_qbvh(const unsigned int index);


    static void refit_subtrees(const unsigned int index, 
			       const unsigned int threshold, 
			       const unsigned int mode);

    static MortonID32Bit *mortonID[2];
    static RadixNode *radixTree;
    static unsigned int numMortonIDBlocks;

    static void thread_radixsort_blocked(const unsigned int threadID);

    static void thread_getSceneBounds(const unsigned int threadID);
    static void thread_computeMortonIDs(const unsigned int threadID);
    static void thread_getTaskIntervals(const unsigned int threadID);
    static void thread_createBuildRecordsFromTaskIntervals(const unsigned int threadID);

    static void createBuildRecordsFromTaskIntervals();

    static void createTopLevel(const unsigned long start,
			       const unsigned long end,
			       const unsigned long fatherID,
			       const MortonID32Bit *__restrict__ const m,
			       unsigned long &i_index);

    static unsigned int recurseMortonID(BuildRecord &current, 
					const unsigned int mode, 
					const unsigned int threadID,
					unsigned int &localNodeID,
					unsigned int &localNodeIDs);

    static void handleLargeLeaves(const unsigned int nodeID,
				  const unsigned int start,
				  const unsigned int items);

    static _INLINE void getLeafBoundsAndComputeAccelerationData(const unsigned int start,
								const unsigned int items,
								const unsigned int nodeID)
    {
      assert(items<=4);

      prefetch<PFHINT_L1EX>((float*)&globalNodePtr[nodeID]);		

      MIC_ALIGN Centroid_Scene_AABB bounds;
      bounds.reset();

      const Vec3fa          * __restrict__ const vptr = globalVertexPtr;
      const BuilderTriangle * __restrict__ const tptr = globalBuildTrianglePtr;
      const MortonID32Bit   *__restrict__  const    m = mortonID[0];

      assert(tptr != NULL);

      for (long i=0;i<items;i++)
	{	
	  const unsigned int primID = m[start+i].index;

	  const unsigned int i0 = tptr[primID].v0;
	  const unsigned int i1 = tptr[primID].v1;
	  const unsigned int i2 = tptr[primID].v2;

	  prefetch<PFHINT_L1>(&vptr[i0]);
	  prefetch<PFHINT_L1>(&vptr[i1]);
	  prefetch<PFHINT_L1>(&vptr[i2]);

	  const mic_f vtxA = upconv4f((float*)&vptr[i0]);
	  const mic_f vtxB = upconv4f((float*)&vptr[i1]);
	  const mic_f vtxC = upconv4f((float*)&vptr[i2]);
	  const mic_f e1 = vtxA - vtxB;
	  const mic_f e2 = vtxC - vtxA;	     
	  const mic_f normal = lcross_xyz(e1,e2);
	  const mic_f v0 = sel(0x8888,cast_to_mic_f(mic_i(primID)),vtxA);
	  const mic_f v1 = sel(0x8888,cast_to_mic_f(mic_i(tptr[primID].id0)),e1);
	  const mic_f v2 = sel(0x8888,cast_to_mic_f(mic_i(tptr[primID].id1)),e2);
	  const mic_f v3 = sel(0x8888,cast_to_mic_f(0x0),normal);
	  const mic_f acc_f = lane_shuffle_gather<0>(v0,v1,v2,v3);
	  const mic_f b_min = _min(vtxA,_min(vtxB,vtxC));
	  const mic_f b_max = _max(vtxA,_max(vtxB,vtxC));
	  bounds.extend_scene(b_min,b_max);
	  store16f_ngo(&globalAccelPtr[start + i],acc_f);
	}
      bounds.storeSceneAABB((float*)&globalNodePtr[nodeID]);
    }

  
  public:

    static void build();

  
  };


};

#endif
