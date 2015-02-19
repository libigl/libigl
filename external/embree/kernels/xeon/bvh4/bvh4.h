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

#include "embree2/rtcore.h"
#include "common/alloc.h"
#include "common/accel.h"
#include "common/scene.h"
#include "geometry/primitive.h"
#include "common/ray.h"

namespace embree
{
  /*! Multi BVH with 4 children. Each node stores the bounding box of
   * it's 4 children as well as 4 child pointers. */
  class BVH4 : public AccelData
  {
    ALIGNED_CLASS;
  public:
    
    /*! forward declaration of node type */
    struct BaseNode;
    struct Node;
    struct NodeMB;
    struct UnalignedNode;
    struct NodeSingleSpaceMB;
    struct NodeDualSpaceMB;
    struct NodeConeMB;
#define BVH4HAIR_MB_VERSION 0

#if BVH4HAIR_MB_VERSION == 0
    typedef NodeSingleSpaceMB UnalignedNodeMB;
#elif BVH4HAIR_MB_VERSION == 1
    typedef NodeDualSpaceMB UnalignedNodeMB;
#elif BVH4HAIR_MB_VERSION == 2
    typedef NodeConeMB UnalignedNodeMB;
#endif

    /*! branching width of the tree */
    static const size_t N = 4;

    /*! Number of address bits the Node and primitives are aligned
        to. Maximally 2^alignment-1 many primitive blocks per leaf are
        supported. */
    static const size_t alignment = 4;

    /*! highest address bit is used as barrier for some algorithms */
    static const size_t barrier_mask = (1LL << (8*sizeof(size_t)-1));

    /*! Masks the bits that store the number of items per leaf. */
    static const size_t align_mask = (1 << alignment)-1;  
    static const size_t items_mask = (1 << alignment)-1;  

    /*! different supported node types */
    static const size_t tyNode = 0;
    static const size_t tyNodeMB = 1;
    static const size_t tyUnalignedNode = 2;
    static const size_t tyUnalignedNodeMB = 3;
    static const size_t tyLeaf = 8;

    /*! Empty node */
    static const size_t emptyNode = tyLeaf;

    /*! Invalid node, used as marker in traversal */
    static const size_t invalidNode = (((size_t)-1) & (~items_mask)) | tyLeaf;
      
    /*! Maximal depth of the BVH. */
    static const size_t maxBuildDepth = 32;
    static const size_t maxBuildDepthLeaf = maxBuildDepth+16;
    static const size_t maxDepth = maxBuildDepthLeaf+maxBuildDepthLeaf+maxBuildDepth;
    
    /*! Maximal number of primitive blocks in a leaf. */
    static const size_t maxLeafBlocks = items_mask-tyLeaf;

    /*! Cost of one traversal step. */
    static const int travCost = 1;
    static const int travCostAligned = 2;
    static const int travCostUnaligned = 3; // FIXME: find best cost
    static const int intCost = 6;

    /*! Pointer that points to a node or a list of primitives */
    struct NodeRef
    {
      /*! Default constructor */
      __forceinline NodeRef () {}

      /*! Construction from integer */
      __forceinline NodeRef (size_t ptr) : ptr(ptr) {}

      /*! Cast to size_t */
      __forceinline operator size_t() const { return ptr; }

      /*! Prefetches the node this reference points to */
      __forceinline void prefetch(int types) const {
	prefetchL1(((char*)ptr)+0*64);
	prefetchL1(((char*)ptr)+1*64);
	if (types > 0x1) {
	  prefetchL1(((char*)ptr)+2*64);
	  prefetchL1(((char*)ptr)+3*64);
	  /*prefetchL1(((char*)ptr)+4*64);
	  prefetchL1(((char*)ptr)+5*64);
	  prefetchL1(((char*)ptr)+6*64);
	  prefetchL1(((char*)ptr)+7*64);*/
	}
      }

      /*! Sets the barrier bit. */
      __forceinline void setBarrier() { ptr |= barrier_mask; }
      
      /*! Clears the barrier bit. */
      __forceinline void clearBarrier() { ptr &= ~barrier_mask; }

      /*! Checks if this is an barrier. A barrier tells the top level tree rotations how deep to enter the tree. */
      __forceinline bool isBarrier() const { return (ptr & barrier_mask) != 0; }

      /*! checks if this is a leaf */
      __forceinline size_t isLeaf() const { return ptr & tyLeaf; }

      /*! checks if this is a leaf */
      __forceinline int isLeaf(int types) const { 
	if      (types == 0x0001) return !isNode();
	/*else if (types == 0x0010) return !isNodeMB();
	else if (types == 0x0100) return !isUnalignedNode();
	else if (types == 0x1000) return !isUnalignedNodeMB();*/
	else return isLeaf();
      }
      
      /*! checks if this is a node */
      __forceinline int isNode() const { return (ptr & (size_t)align_mask) == tyNode; }
      __forceinline int isNode(int types) const { return (types == 0x1) || ((types & 0x1) && isNode()); }

      /*! checks if this is a motion blur node */
      __forceinline int isNodeMB() const { return (ptr & (size_t)align_mask) == tyNodeMB; }
      __forceinline int isNodeMB(int types) const { return (types == 0x10) || ((types & 0x10) && isNodeMB()); }

      /*! checks if this is a node with unaligned bounding boxes */
      __forceinline int isUnalignedNode() const { return (ptr & (size_t)align_mask) == tyUnalignedNode; }
      __forceinline int isUnalignedNode(int types) const { return (types == 0x100) || ((types & 0x100) && isUnalignedNode()); }

      /*! checks if this is a motion blur node with unaligned bounding boxes */
      __forceinline int isUnalignedNodeMB() const { return (ptr & (size_t)align_mask) == tyUnalignedNodeMB; }
      __forceinline int isUnalignedNodeMB(int types) const { return (types == 0x1000) || ((types & 0x1000) && isUnalignedNodeMB()); }

      /*! returns base node pointer */
      __forceinline BaseNode* baseNode(int types) { 
	assert(!isLeaf()); 
	if (types == 0x1) return (BaseNode*)ptr; 
	else              return (BaseNode*)(ptr & ~(size_t)align_mask); 
      }
      __forceinline const BaseNode* baseNode(int types) const { 
	assert(!isLeaf()); 
	if (types == 0x1) return (const BaseNode*)ptr; 
	else              return (const BaseNode*)(ptr & ~(size_t)align_mask); 
      }

      /*! returns node pointer */
      __forceinline       Node* node()       { assert(isNode()); return (      Node*)ptr; }
      __forceinline const Node* node() const { assert(isNode()); return (const Node*)ptr; }

      /*! returns motion blur node pointer */
      __forceinline       NodeMB* nodeMB()       { assert(isNodeMB()); return (      NodeMB*)(ptr & ~(size_t)align_mask); }
      __forceinline const NodeMB* nodeMB() const { assert(isNodeMB()); return (const NodeMB*)(ptr & ~(size_t)align_mask); }

      /*! returns unaligned node pointer */
      __forceinline       UnalignedNode* unalignedNode()       { assert(isUnalignedNode()); return (      UnalignedNode*)(ptr & ~(size_t)align_mask); }
      __forceinline const UnalignedNode* unalignedNode() const { assert(isUnalignedNode()); return (const UnalignedNode*)(ptr & ~(size_t)align_mask); }

      /*! returns unaligned motion blur node pointer */
      __forceinline       UnalignedNodeMB* unalignedNodeMB()       { assert(isUnalignedNodeMB()); return (      UnalignedNodeMB*)(ptr & ~(size_t)align_mask); }
      __forceinline const UnalignedNodeMB* unalignedNodeMB() const { assert(isUnalignedNodeMB()); return (const UnalignedNodeMB*)(ptr & ~(size_t)align_mask); }
            
      /*! returns leaf pointer */
      __forceinline char* leaf(size_t& num) const {
        assert(isLeaf());
        num = (ptr & (size_t)items_mask)-tyLeaf;
        return (char*)(ptr & ~(size_t)align_mask);
      }

      /*! clear all bit flags */
      __forceinline void clearFlags() {
        ptr &= ~(size_t)align_mask;
      }

    private:
      size_t ptr;
    };
    
    /*! BVH4 Base Node */
    struct BaseNode
    {
      /*! Clears the node. */
      __forceinline void clear() {
	for (size_t i=0; i<N; i++) children[i] = emptyNode;
      }

        /*! Returns reference to specified child */
      __forceinline       NodeRef& child(size_t i)       { assert(i<N); return children[i]; }
      __forceinline const NodeRef& child(size_t i) const { assert(i<N); return children[i]; }

      /*! verifies the node */
      __forceinline bool verify() const  // FIXME: call in statistics
      {
	for (size_t i=0; i<BVH4::N; i++) {
	  if (child(i) == BVH4::emptyNode) {
	    for (; i<BVH4::N; i++) {
	      if (child(i) != BVH4::emptyNode)
		return false;
	    }
	    break;
	  }
	}
	return true;
      }

      NodeRef children[N];    //!< Pointer to the 4 children (can be a node or leaf)
    };
    
    /*! BVH4 Node */
    struct Node : public BaseNode
    {
      /*! Clears the node. */
      __forceinline void clear() {
        lower_x = lower_y = lower_z = pos_inf; 
        upper_x = upper_y = upper_z = neg_inf;
	BaseNode::clear();
      }

      /*! Sets bounding box and ID of child. */
      __forceinline void set(size_t i, const NodeRef& childID) {
	assert(i < N);
        children[i] = childID;
      }

      /*! Sets bounding box of child. */
      __forceinline void set(size_t i, const BBox3fa& bounds) 
      {
        assert(i < N);
        lower_x[i] = bounds.lower.x; lower_y[i] = bounds.lower.y; lower_z[i] = bounds.lower.z;
        upper_x[i] = bounds.upper.x; upper_y[i] = bounds.upper.y; upper_z[i] = bounds.upper.z;
      }

      /*! Sets bounding box and ID of child. */
      __forceinline void set(size_t i, const BBox3fa& bounds, const NodeRef& childID) {
        set(i,bounds);
        children[i] = childID;
      }

      /*! Returns bounds of node. */
      __forceinline BBox3fa bounds() const {
        const Vec3fa lower(reduce_min(lower_x),reduce_min(lower_y),reduce_min(lower_z));
        const Vec3fa upper(reduce_max(upper_x),reduce_max(upper_y),reduce_max(upper_z));
        return BBox3fa(lower,upper);
      }

      /*! Returns bounds of specified child. */
      __forceinline BBox3fa bounds(size_t i) const 
      {
        assert(i < N);
        const Vec3fa lower(lower_x[i],lower_y[i],lower_z[i]);
        const Vec3fa upper(upper_x[i],upper_y[i],upper_z[i]);
        return BBox3fa(lower,upper);
      }

      /*! Returns extent of bounds of specified child. */
      __forceinline BBox3fa extend(size_t i) const {
	return bounds(i).size();
      }

      /*! Returns bounds of all children */
      __forceinline void bounds(BBox<ssef>& bounds0, BBox<ssef>& bounds1, BBox<ssef>& bounds2, BBox<ssef>& bounds3) const {
        transpose(lower_x,lower_y,lower_z,ssef(zero),bounds0.lower,bounds1.lower,bounds2.lower,bounds3.lower);
        transpose(upper_x,upper_y,upper_z,ssef(zero),bounds0.upper,bounds1.upper,bounds2.upper,bounds3.upper);
      }

      /*! swap two children of the node */
      __forceinline void swap(size_t i, size_t j)
      {
	assert(i<N && j<N);
	std::swap(children[i],children[j]);
	std::swap(lower_x[i],lower_x[j]);
	std::swap(lower_y[i],lower_y[j]);
	std::swap(lower_z[i],lower_z[j]);
	std::swap(upper_x[i],upper_x[j]);
	std::swap(upper_y[i],upper_y[j]);
	std::swap(upper_z[i],upper_z[j]);
      }

      /*! Returns reference to specified child */
      __forceinline       NodeRef& child(size_t i)       { assert(i<N); return children[i]; }
      __forceinline const NodeRef& child(size_t i) const { assert(i<N); return children[i]; }

      /*! intersection with single rays */
      template<bool robust>
      __forceinline size_t intersect(size_t nearX, size_t nearY, size_t nearZ,
				     const sse3f& org, const sse3f& rdir, const sse3f& org_rdir, const ssef& tnear, const ssef& tfar, 
				     ssef& dist) const
      {
	const size_t farX  = nearX ^ sizeof(ssef), farY  = nearY ^ sizeof(ssef), farZ  = nearZ ^ sizeof(ssef);
#if defined (__AVX2__)
	const ssef tNearX = msub(load4f((const char*)&lower_x+nearX), rdir.x, org_rdir.x);
	const ssef tNearY = msub(load4f((const char*)&lower_x+nearY), rdir.y, org_rdir.y);
	const ssef tNearZ = msub(load4f((const char*)&lower_x+nearZ), rdir.z, org_rdir.z);
	const ssef tFarX  = msub(load4f((const char*)&lower_x+farX ), rdir.x, org_rdir.x);
	const ssef tFarY  = msub(load4f((const char*)&lower_x+farY ), rdir.y, org_rdir.y);
	const ssef tFarZ  = msub(load4f((const char*)&lower_x+farZ ), rdir.z, org_rdir.z);
#else
	const ssef tNearX = (load4f((const char*)&lower_x+nearX) - org.x) * rdir.x;
	const ssef tNearY = (load4f((const char*)&lower_x+nearY) - org.y) * rdir.y;
	const ssef tNearZ = (load4f((const char*)&lower_x+nearZ) - org.z) * rdir.z;
	const ssef tFarX  = (load4f((const char*)&lower_x+farX ) - org.x) * rdir.x;
	const ssef tFarY  = (load4f((const char*)&lower_x+farY ) - org.y) * rdir.y;
	const ssef tFarZ  = (load4f((const char*)&lower_x+farZ ) - org.z) * rdir.z;
#endif

        if (robust) {
          const float round_down = 1.0f-2.0f*float(ulp);
          const float round_up   = 1.0f+2.0f*float(ulp);
          const ssef tNear = max(tNearX,tNearY,tNearZ,tnear);
          const ssef tFar  = min(tFarX ,tFarY ,tFarZ ,tfar);
          const sseb vmask = round_down*tNear <= round_up*tFar;
          const size_t mask = movemask(vmask);
          dist = tNear;
          return mask;
        }

#if defined(__SSE4_1__)
	const ssef tNear = maxi(maxi(tNearX,tNearY),maxi(tNearZ,tnear));
	const ssef tFar  = mini(mini(tFarX ,tFarY ),mini(tFarZ ,tfar ));
	const sseb vmask = cast(tNear) > cast(tFar);
	const size_t mask = movemask(vmask)^0xf;
#else
	const ssef tNear = max(tNearX,tNearY,tNearZ,tnear);
	const ssef tFar  = min(tFarX ,tFarY ,tFarZ ,tfar);
	const sseb vmask = tNear <= tFar;
	const size_t mask = movemask(vmask);
#endif
	dist = tNear;
	return mask;
      }

      /*! intersection with ray packet of size 4 */
      template<bool robust>
      __forceinline sseb intersect(size_t i, const sse3f& org, const sse3f& rdir, const sse3f& org_rdir, const ssef& tnear, const ssef& tfar, ssef& dist) const
      {
#if defined(__AVX2__)
	const ssef lclipMinX = msub(lower_x[i],rdir.x,org_rdir.x);
	const ssef lclipMinY = msub(lower_y[i],rdir.y,org_rdir.y);
	const ssef lclipMinZ = msub(lower_z[i],rdir.z,org_rdir.z);
	const ssef lclipMaxX = msub(upper_x[i],rdir.x,org_rdir.x);
	const ssef lclipMaxY = msub(upper_y[i],rdir.y,org_rdir.y);
	const ssef lclipMaxZ = msub(upper_z[i],rdir.z,org_rdir.z);
#else
	const ssef lclipMinX = (lower_x[i] - org.x) * rdir.x;
	const ssef lclipMinY = (lower_y[i] - org.y) * rdir.y;
	const ssef lclipMinZ = (lower_z[i] - org.z) * rdir.z;
	const ssef lclipMaxX = (upper_x[i] - org.x) * rdir.x;
	const ssef lclipMaxY = (upper_y[i] - org.y) * rdir.y;
	const ssef lclipMaxZ = (upper_z[i] - org.z) * rdir.z;
#endif

        if (robust) {
          const float round_down = 1.0f-2.0f*float(ulp);
          const float round_up   = 1.0f+2.0f*float(ulp);
          const ssef lnearP = max(max(min(lclipMinX, lclipMaxX), min(lclipMinY, lclipMaxY)), min(lclipMinZ, lclipMaxZ));
          const ssef lfarP  = min(min(max(lclipMinX, lclipMaxX), max(lclipMinY, lclipMaxY)), max(lclipMinZ, lclipMaxZ));
          const sseb lhit   = round_down*max(lnearP,tnear) <= round_up*min(lfarP,tfar);      
          dist = lnearP;
          return lhit;
        }

#if defined(__SSE4_1__)
	const ssef lnearP = maxi(maxi(mini(lclipMinX, lclipMaxX), mini(lclipMinY, lclipMaxY)), mini(lclipMinZ, lclipMaxZ));
	const ssef lfarP  = mini(mini(maxi(lclipMinX, lclipMaxX), maxi(lclipMinY, lclipMaxY)), maxi(lclipMinZ, lclipMaxZ));
	const sseb lhit   = maxi(lnearP,tnear) <= mini(lfarP,tfar);      
#else
	
	const ssef lnearP = max(max(min(lclipMinX, lclipMaxX), min(lclipMinY, lclipMaxY)), min(lclipMinZ, lclipMaxZ));
	const ssef lfarP  = min(min(max(lclipMinX, lclipMaxX), max(lclipMinY, lclipMaxY)), max(lclipMinZ, lclipMaxZ));
	const sseb lhit   = max(lnearP,tnear) <= min(lfarP,tfar);      
#endif
	dist = lnearP;
	return lhit;
      }
      
      /*! intersection with ray packet of size 8 */
#if defined(__AVX__)
      template<bool robust>
      __forceinline avxb intersect8(size_t i, const avx3f& org, const avx3f& rdir, const avx3f& org_rdir, const avxf& tnear, const avxf& tfar, avxf& dist) const
      {
#if defined(__AVX2__)
	const avxf lclipMinX = msub(lower_x[i],rdir.x,org_rdir.x);
	const avxf lclipMinY = msub(lower_y[i],rdir.y,org_rdir.y);
	const avxf lclipMinZ = msub(lower_z[i],rdir.z,org_rdir.z);
	const avxf lclipMaxX = msub(upper_x[i],rdir.x,org_rdir.x);
	const avxf lclipMaxY = msub(upper_y[i],rdir.y,org_rdir.y);
	const avxf lclipMaxZ = msub(upper_z[i],rdir.z,org_rdir.z);
#else
	const avxf lclipMinX = (lower_x[i] - org.x) * rdir.x;
	const avxf lclipMinY = (lower_y[i] - org.y) * rdir.y;
	const avxf lclipMinZ = (lower_z[i] - org.z) * rdir.z;
	const avxf lclipMaxX = (upper_x[i] - org.x) * rdir.x;
	const avxf lclipMaxY = (upper_y[i] - org.y) * rdir.y;
	const avxf lclipMaxZ = (upper_z[i] - org.z) * rdir.z;
#endif

        if (robust) {
          const float round_down = 1.0f-2.0f*float(ulp);
          const float round_up   = 1.0f+2.0f*float(ulp);
          const avxf lnearP = max(max(min(lclipMinX, lclipMaxX), min(lclipMinY, lclipMaxY)), min(lclipMinZ, lclipMaxZ));
          const avxf lfarP  = min(min(max(lclipMinX, lclipMaxX), max(lclipMinY, lclipMaxY)), max(lclipMinZ, lclipMaxZ));
          const avxb lhit   = round_down*max(lnearP,tnear) <= round_up*min(lfarP,tfar);      
          dist = lnearP;
          return lhit;
        }

#if defined(__AVX2__)
	const avxf lnearP = maxi(maxi(mini(lclipMinX, lclipMaxX), mini(lclipMinY, lclipMaxY)), mini(lclipMinZ, lclipMaxZ));
	const avxf lfarP  = mini(mini(maxi(lclipMinX, lclipMaxX), maxi(lclipMinY, lclipMaxY)), maxi(lclipMinZ, lclipMaxZ));
	const avxb lhit   = maxi(lnearP,tnear) <= mini(lfarP,tfar);      
#else
	const avxf lnearP = max(max(min(lclipMinX, lclipMaxX), min(lclipMinY, lclipMaxY)), min(lclipMinZ, lclipMaxZ));
	const avxf lfarP  = min(min(max(lclipMinX, lclipMaxX), max(lclipMinY, lclipMaxY)), max(lclipMinZ, lclipMaxZ));
	const avxb lhit   = max(lnearP,tnear) <= min(lfarP,tfar);      
#endif
	dist = lnearP;
	return lhit;
      }
#endif
      
    public:
      ssef lower_x;           //!< X dimension of lower bounds of all 4 children.
      ssef upper_x;           //!< X dimension of upper bounds of all 4 children.
      ssef lower_y;           //!< Y dimension of lower bounds of all 4 children.
      ssef upper_y;           //!< Y dimension of upper bounds of all 4 children.
      ssef lower_z;           //!< Z dimension of lower bounds of all 4 children.
      ssef upper_z;           //!< Z dimension of upper bounds of all 4 children.
    };

    /*! Motion Blur Node */
    struct NodeMB : public BaseNode
    {
      /*! Clears the node. */
      __forceinline void clear()  {
        lower_x = lower_y = lower_z = ssef(nan);
        upper_x = upper_y = upper_z = ssef(nan);
        lower_dx = lower_dy = lower_dz = ssef(nan); // initialize with NAN and update during refit
        upper_dx = upper_dy = upper_dz = ssef(nan);
	BaseNode::clear();
      }

      /*! Sets bounding box and ID of child. */
      __forceinline void set(size_t i, const BBox3fa& bounds, NodeRef childID) {
        lower_x[i] = bounds.lower.x; lower_y[i] = bounds.lower.y; lower_z[i] = bounds.lower.z;
        upper_x[i] = bounds.upper.x; upper_y[i] = bounds.upper.y; upper_z[i] = bounds.upper.z;
        children[i] = childID;
      }

      /*! Sets bounding box and ID of child. */
      __forceinline void set(size_t i, NodeRef childID) {
	children[i] = childID;
      }

      /*! Sets bounding box and ID of child. */
      __forceinline void set(size_t i, const BBox3fa& bounds) {
        lower_x[i] = bounds.lower.x; lower_y[i] = bounds.lower.y; lower_z[i] = bounds.lower.z;
        upper_x[i] = bounds.upper.x; upper_y[i] = bounds.upper.y; upper_z[i] = bounds.upper.z;
      }

      /*! Sets bounding box and ID of child. */
      __forceinline void set(size_t i, const BBox3fa& bounds0, const BBox3fa& bounds1) 
      {
        lower_x[i] = bounds0.lower.x; lower_y[i] = bounds0.lower.y; lower_z[i] = bounds0.lower.z;
        upper_x[i] = bounds0.upper.x; upper_y[i] = bounds0.upper.y; upper_z[i] = bounds0.upper.z;

        /*! for empty bounds we have to avoid inf-inf=nan */
        if (unlikely(bounds0.empty())) { 
          lower_dx[i] = lower_dy[i] = lower_dz[i] = zero;
          upper_dx[i] = upper_dy[i] = upper_dz[i] = zero;
        } 
        /*! standard case */
        else {
          const Vec3fa dlower = bounds1.lower-bounds0.lower;
          const Vec3fa dupper = bounds1.upper-bounds0.upper;
          lower_dx[i] = dlower.x; lower_dy[i] = dlower.y; lower_dz[i] = dlower.z;
          upper_dx[i] = dupper.x; upper_dy[i] = dupper.y; upper_dz[i] = dupper.z;
        }
      }

      /*! tests if the node has valid bounds */
      __forceinline bool hasBounds() const {
        return lower_dx.i[0] != cast_f2i(float(nan));
      }

      /*! Return bounding box for time 0 */
      __forceinline BBox3fa bounds0(size_t i) const {
        return BBox3fa(Vec3fa(lower_x[i],lower_y[i],lower_z[i]),
                      Vec3fa(upper_x[i],upper_y[i],upper_z[i]));
      }

      /*! Return bounding box for time 1 */
      __forceinline BBox3fa bounds1(size_t i) const {
        return BBox3fa(Vec3fa(lower_x[i]+lower_dx[i],lower_y[i]+lower_dy[i],lower_z[i]+lower_dz[i]),
                      Vec3fa(upper_x[i]+upper_dx[i],upper_y[i]+upper_dy[i],upper_z[i]+upper_dz[i]));
      }

      /*! Returns extent of bounds of specified child. */
      __forceinline BBox3fa extend0(size_t i) const {
	return bounds0(i).size();
      }

      /*! Returns bounds of node. */
      __forceinline BBox3fa bounds() const {
        return BBox3fa(Vec3fa(reduce_min(min(lower_x,lower_x+lower_dx)),
                             reduce_min(min(lower_y,lower_y+lower_dy)),
                             reduce_min(min(lower_z,lower_z+lower_dz))),
                      Vec3fa(reduce_max(max(upper_x,upper_x+upper_dx)),
                             reduce_max(max(upper_y,upper_y+upper_dy)),
                             reduce_max(max(upper_z,upper_z+upper_dz))));
      }

      /*! swap two children of the node */
      __forceinline void swap(size_t i, size_t j)
      {
	assert(i<N && j<N);
	std::swap(children[i],children[j]);

	std::swap(lower_x[i],lower_x[j]);
	std::swap(upper_x[i],upper_x[j]);
	std::swap(lower_y[i],lower_y[j]);
	std::swap(upper_y[i],upper_y[j]);
	std::swap(lower_z[i],lower_z[j]);
	std::swap(upper_z[i],upper_z[j]);
	
	std::swap(lower_dx[i],lower_dx[j]);
	std::swap(upper_dx[i],upper_dx[j]);
	std::swap(lower_dy[i],lower_dy[j]);
	std::swap(upper_dy[i],upper_dy[j]);
	std::swap(lower_dz[i],lower_dz[j]);
	std::swap(upper_dz[i],upper_dz[j]);
      }

      /*! Returns reference to specified child */
      __forceinline       NodeRef& child(size_t i)       { assert(i<N); return children[i]; }
      __forceinline const NodeRef& child(size_t i) const { assert(i<N); return children[i]; }

      /*! intersection with single rays */
      __forceinline size_t intersect(size_t nearX, size_t nearY, size_t nearZ,
				     const sse3f& org, const sse3f& rdir, const sse3f& org_rdir, const ssef& tnear, const ssef& tfar, const float time,
				     ssef& dist) const
      {
	const size_t farX  = nearX ^ sizeof(ssef), farY  = nearY ^ sizeof(ssef), farZ  = nearZ ^ sizeof(ssef);
	const ssef* pNearX = (const ssef*)((const char*)&lower_x+nearX);
	const ssef* pNearY = (const ssef*)((const char*)&lower_x+nearY);
	const ssef* pNearZ = (const ssef*)((const char*)&lower_x+nearZ);
	const ssef tNearX = (ssef(pNearX[0]) + time*pNearX[6] - org.x) * rdir.x;
	const ssef tNearY = (ssef(pNearY[0]) + time*pNearY[6] - org.y) * rdir.y;
	const ssef tNearZ = (ssef(pNearZ[0]) + time*pNearZ[6] - org.z) * rdir.z;
	const ssef tNear = max(tnear,tNearX,tNearY,tNearZ);
	const ssef* pFarX = (const ssef*)((const char*)&lower_x+farX);
	const ssef* pFarY = (const ssef*)((const char*)&lower_x+farY);
	const ssef* pFarZ = (const ssef*)((const char*)&lower_x+farZ);
	const ssef tFarX = (ssef(pFarX[0]) + time*pFarX[6] - org.x) * rdir.x;
	const ssef tFarY = (ssef(pFarY[0]) + time*pFarY[6] - org.y) * rdir.y;
	const ssef tFarZ = (ssef(pFarZ[0]) + time*pFarZ[6] - org.z) * rdir.z;
	const ssef tFar = min(tfar,tFarX,tFarY,tFarZ);
	const size_t mask = movemask(tNear <= tFar);
	dist = tNear;
	return mask;
      }

      /*! intersection with ray packet of size 4 */
    __forceinline sseb intersect(const size_t i, const sse3f& org, const sse3f& rdir, const sse3f& org_rdir, const ssef& tnear, const ssef& tfar, const ssef& time, ssef& dist) const
    {
      const ssef vlower_x = ssef(lower_x[i]) + time * ssef(lower_dx[i]);
      const ssef vlower_y = ssef(lower_y[i]) + time * ssef(lower_dy[i]);
      const ssef vlower_z = ssef(lower_z[i]) + time * ssef(lower_dz[i]);
      const ssef vupper_x = ssef(upper_x[i]) + time * ssef(upper_dx[i]);
      const ssef vupper_y = ssef(upper_y[i]) + time * ssef(upper_dy[i]);
      const ssef vupper_z = ssef(upper_z[i]) + time * ssef(upper_dz[i]);
      
#if defined(__AVX2__)
      const ssef lclipMinX = msub(vlower_x,rdir.x,org_rdir.x);
      const ssef lclipMinY = msub(vlower_y,rdir.y,org_rdir.y);
      const ssef lclipMinZ = msub(vlower_z,rdir.z,org_rdir.z);
      const ssef lclipMaxX = msub(vupper_x,rdir.x,org_rdir.x);
      const ssef lclipMaxY = msub(vupper_y,rdir.y,org_rdir.y);
      const ssef lclipMaxZ = msub(vupper_z,rdir.z,org_rdir.z);
#else
      const ssef lclipMinX = (vlower_x - org.x) * rdir.x;
      const ssef lclipMinY = (vlower_y - org.y) * rdir.y;
      const ssef lclipMinZ = (vlower_z - org.z) * rdir.z;
      const ssef lclipMaxX = (vupper_x - org.x) * rdir.x;
      const ssef lclipMaxY = (vupper_y - org.y) * rdir.y;
      const ssef lclipMaxZ = (vupper_z - org.z) * rdir.z;
#endif

#if defined(__SSE4_1__)
      const ssef lnearP = maxi(maxi(mini(lclipMinX, lclipMaxX), mini(lclipMinY, lclipMaxY)), mini(lclipMinZ, lclipMaxZ));
      const ssef lfarP  = mini(mini(maxi(lclipMinX, lclipMaxX), maxi(lclipMinY, lclipMaxY)), maxi(lclipMinZ, lclipMaxZ));
      const sseb lhit   = maxi(lnearP,tnear) <= mini(lfarP,tfar);      
#else
      
      const ssef lnearP = max(max(min(lclipMinX, lclipMaxX), min(lclipMinY, lclipMaxY)), min(lclipMinZ, lclipMaxZ));
      const ssef lfarP  = min(min(max(lclipMinX, lclipMaxX), max(lclipMinY, lclipMaxY)), max(lclipMinZ, lclipMaxZ));
      const sseb lhit   = max(lnearP,tnear) <= min(lfarP,tfar);      
#endif
      dist = lnearP;
      return lhit;
    }

    /*! intersection with ray packet of size 8 */
#if defined(__AVX__)
    __forceinline avxb intersect(const size_t i, const avx3f& org, const avx3f& rdir, const avx3f& org_rdir, const avxf& tnear, const avxf& tfar, const avxf& time, avxf& dist) const
    {
      const avxf vlower_x = avxf(lower_x[i]) + time * avxf(lower_dx[i]);
      const avxf vlower_y = avxf(lower_y[i]) + time * avxf(lower_dy[i]);
      const avxf vlower_z = avxf(lower_z[i]) + time * avxf(lower_dz[i]);
      const avxf vupper_x = avxf(upper_x[i]) + time * avxf(upper_dx[i]);
      const avxf vupper_y = avxf(upper_y[i]) + time * avxf(upper_dy[i]);
      const avxf vupper_z = avxf(upper_z[i]) + time * avxf(upper_dz[i]);

#if defined(__AVX2__)
      const avxf lclipMinX = msub(vlower_x,rdir.x,org_rdir.x);
      const avxf lclipMinY = msub(vlower_y,rdir.y,org_rdir.y);
      const avxf lclipMinZ = msub(vlower_z,rdir.z,org_rdir.z);
      const avxf lclipMaxX = msub(vupper_x,rdir.x,org_rdir.x);
      const avxf lclipMaxY = msub(vupper_y,rdir.y,org_rdir.y);
      const avxf lclipMaxZ = msub(vupper_z,rdir.z,org_rdir.z);
      const avxf lnearP = maxi(maxi(mini(lclipMinX, lclipMaxX), mini(lclipMinY, lclipMaxY)), mini(lclipMinZ, lclipMaxZ));
      const avxf lfarP  = mini(mini(maxi(lclipMinX, lclipMaxX), maxi(lclipMinY, lclipMaxY)), maxi(lclipMinZ, lclipMaxZ));
      const avxb lhit   = maxi(lnearP,tnear) <= mini(lfarP,tfar);      
#else
      const avxf lclipMinX = (vlower_x - org.x) * rdir.x;
      const avxf lclipMinY = (vlower_y - org.y) * rdir.y;
      const avxf lclipMinZ = (vlower_z - org.z) * rdir.z;
      const avxf lclipMaxX = (vupper_x - org.x) * rdir.x;
      const avxf lclipMaxY = (vupper_y - org.y) * rdir.y;
      const avxf lclipMaxZ = (vupper_z - org.z) * rdir.z;
      const avxf lnearP = max(max(min(lclipMinX, lclipMaxX), min(lclipMinY, lclipMaxY)), min(lclipMinZ, lclipMaxZ));
      const avxf lfarP  = min(min(max(lclipMinX, lclipMaxX), max(lclipMinY, lclipMaxY)), max(lclipMinZ, lclipMaxZ));
      const avxb lhit   = max(lnearP,tnear) <= min(lfarP,tfar);      
#endif
      dist = lnearP;
      return lhit;
    }
#endif

    public:
      ssef lower_x;        //!< X dimension of lower bounds of all 4 children.
      ssef upper_x;        //!< X dimension of upper bounds of all 4 children.
      ssef lower_y;        //!< Y dimension of lower bounds of all 4 children.
      ssef upper_y;        //!< Y dimension of upper bounds of all 4 children.
      ssef lower_z;        //!< Z dimension of lower bounds of all 4 children.
      ssef upper_z;        //!< Z dimension of upper bounds of all 4 children.

      ssef lower_dx;        //!< X dimension of lower bounds of all 4 children.
      ssef upper_dx;        //!< X dimension of upper bounds of all 4 children.
      ssef lower_dy;        //!< Y dimension of lower bounds of all 4 children.
      ssef upper_dy;        //!< Y dimension of upper bounds of all 4 children.
      ssef lower_dz;        //!< Z dimension of lower bounds of all 4 children.
      ssef upper_dz;        //!< Z dimension of upper bounds of all 4 children.
    };

    /*! Node with unaligned bounds */
    struct UnalignedNode : public BaseNode
    {
      /*! Clears the node. */
      __forceinline void clear() 
      {
	naabb.l.vx = Vec3fa(nan);
	naabb.l.vy = Vec3fa(nan);
	naabb.l.vz = Vec3fa(nan);
	naabb.p    = Vec3fa(nan);
	BaseNode::clear();
      }

      /*! Sets bounding box. */
      __forceinline void set(size_t i, const NAABBox3fa& b) 
      {
        assert(i < N);

        AffineSpace3fa space = b.space;
        space.p -= b.bounds.lower;
        space = AffineSpace3fa::scale(1.0f/max(Vec3fa(1E-19f),b.bounds.upper-b.bounds.lower))*space;
        
        naabb.l.vx.x[i] = space.l.vx.x;
        naabb.l.vx.y[i] = space.l.vx.y;
        naabb.l.vx.z[i] = space.l.vx.z;

        naabb.l.vy.x[i] = space.l.vy.x;
        naabb.l.vy.y[i] = space.l.vy.y;
        naabb.l.vy.z[i] = space.l.vy.z;

        naabb.l.vz.x[i] = space.l.vz.x;
        naabb.l.vz.y[i] = space.l.vz.y;
        naabb.l.vz.z[i] = space.l.vz.z;

        naabb.p.x[i] = space.p.x;
        naabb.p.y[i] = space.p.y;
        naabb.p.z[i] = space.p.z;
      }

      /*! Sets ID of child. */
      __forceinline void set(size_t i, const NodeRef& childID) {
        //Node::set(i,childID);
	assert(i < N);
	children[i] = childID;
      }

      /*! Returns the extend of the bounds of the ith child */
      __forceinline Vec3fa extend(size_t i) const {
        assert(i<N);
        const Vec3fa vx(naabb.l.vx.x[i],naabb.l.vx.y[i],naabb.l.vx.z[i]);
        const Vec3fa vy(naabb.l.vy.x[i],naabb.l.vy.y[i],naabb.l.vy.z[i]);
        const Vec3fa vz(naabb.l.vz.x[i],naabb.l.vz.y[i],naabb.l.vz.z[i]);
        const Vec3fa p (naabb.p   .x[i],naabb.p   .y[i],naabb.p   .z[i]);
        return rsqrt(vx*vx + vy*vy + vz*vz);
      }

      /*! Returns reference to specified child */
      __forceinline       NodeRef& child(size_t i)       { assert(i<N); return children[i]; }
      __forceinline const NodeRef& child(size_t i) const { assert(i<N); return children[i]; }

      /*! intersect 4 OBBs with single ray */
      __forceinline size_t intersect(const sse3f& ray_org, const sse3f& ray_dir, 
				     const ssef& tnear, const ssef& tfar, ssef& dist)
      {
	const sse3f dir = xfmVector(naabb,ray_dir);
	//const sse3f nrdir = sse3f(ssef(-1.0f))/dir;
	const sse3f nrdir = sse3f(ssef(-1.0f))*rcp_safe(dir);
	const sse3f org = xfmPoint(naabb,ray_org);
	const sse3f tLowerXYZ = org * nrdir;     // (Vec3fa(zero) - org) * rdir;
	const sse3f tUpperXYZ = tLowerXYZ - nrdir; // (Vec3fa(one ) - org) * rdir;

#if defined(__SSE4_1__)
	const ssef tNearX = mini(tLowerXYZ.x,tUpperXYZ.x);
	const ssef tNearY = mini(tLowerXYZ.y,tUpperXYZ.y);
	const ssef tNearZ = mini(tLowerXYZ.z,tUpperXYZ.z);
	const ssef tFarX  = maxi(tLowerXYZ.x,tUpperXYZ.x);
	const ssef tFarY  = maxi(tLowerXYZ.y,tUpperXYZ.y);
	const ssef tFarZ  = maxi(tLowerXYZ.z,tUpperXYZ.z);
	const ssef tNear  = max(tnear, tNearX,tNearY,tNearZ);
	const ssef tFar   = min(tfar,  tFarX ,tFarY ,tFarZ );
	const sseb vmask = tNear <= tFar;
	dist = tNear;
	return movemask(vmask);
#else
	const ssef tNearX = min(tLowerXYZ.x,tUpperXYZ.x);
	const ssef tNearY = min(tLowerXYZ.y,tUpperXYZ.y);
	const ssef tNearZ = min(tLowerXYZ.z,tUpperXYZ.z);
	const ssef tFarX  = max(tLowerXYZ.x,tUpperXYZ.x);
	const ssef tFarY  = max(tLowerXYZ.y,tUpperXYZ.y);
	const ssef tFarZ  = max(tLowerXYZ.z,tUpperXYZ.z);
	const ssef tNear = max(tnear, tNearX,tNearY,tNearZ);
	const ssef tFar  = min(tfar,  tFarX ,tFarY ,tFarZ );
	const sseb vmask = tNear <= tFar;
	dist = tNear;
	return movemask(vmask);
#endif
      }

    public:
      AffineSpaceSSE3f naabb;   //!< non-axis aligned bounding boxes (bounds are [0,1] in specified space)
    };

    struct NodeSingleSpaceMB : public BaseNode
    {
      struct Precalculations {
	__forceinline Precalculations (const Ray& ray) {}
      };

      /*! Clears the node. */
      __forceinline void clear() 
      {
        space0 = one;
        //b0.lower = b0.upper = Vec3fa(nan);
        b1.lower = b1.upper = Vec3fa(nan);
        BaseNode::clear();
      }

      /*! Sets space and bounding boxes. */
      __forceinline void set(size_t i, const AffineSpace3fa& s0, const BBox3fa& a, const BBox3fa& c)
      {
        assert(i < N);

	AffineSpace3fa space = s0;
        space.p -= a.lower;
	Vec3fa scale = 1.0f/max(Vec3fa(1E-19f),a.upper-a.lower);
        space = AffineSpace3fa::scale(scale)*space;
	BBox3fa a1((a.lower-a.lower)*scale,(a.upper-a.lower)*scale);
	BBox3fa c1((c.lower-a.lower)*scale,(c.upper-a.lower)*scale);

        space0.l.vx.x[i] = space.l.vx.x; space0.l.vx.y[i] = space.l.vx.y; space0.l.vx.z[i] = space.l.vx.z; 
        space0.l.vy.x[i] = space.l.vy.x; space0.l.vy.y[i] = space.l.vy.y; space0.l.vy.z[i] = space.l.vy.z;
        space0.l.vz.x[i] = space.l.vz.x; space0.l.vz.y[i] = space.l.vz.y; space0.l.vz.z[i] = space.l.vz.z; 
        space0.p   .x[i] = space.p   .x; space0.p   .y[i] = space.p   .y; space0.p   .z[i] = space.p   .z; 

        /*b0.lower.x[i] = a1.lower.x; b0.lower.y[i] = a1.lower.y; b0.lower.z[i] = a1.lower.z;
	  b0.upper.x[i] = a1.upper.x; b0.upper.y[i] = a1.upper.y; b0.upper.z[i] = a1.upper.z;*/

        b1.lower.x[i] = c1.lower.x; b1.lower.y[i] = c1.lower.y; b1.lower.z[i] = c1.lower.z;
        b1.upper.x[i] = c1.upper.x; b1.upper.y[i] = c1.upper.y; b1.upper.z[i] = c1.upper.z;
      }

      /*! Sets ID of child. */
      __forceinline void set(size_t i, const NodeRef& childID) {
        //Node::set(i,childID);
	assert(i < N);
	children[i] = childID;
      }

      /*! Returns bounds of specified child. */
      __forceinline const BBox3fa bounds0(const size_t i) const { 
        assert(i < N);
        /*const Vec3fa lower(b0.lower.x[i],b0.lower.y[i],b0.lower.z[i]);
        const Vec3fa upper(b0.upper.x[i],b0.upper.y[i],b0.upper.z[i]);
        return BBox3fa(lower,upper);*/
	return empty; // FIXME: not yet implemented
      }

      /*! Returns the extend of the bounds of the ith child */
      __forceinline Vec3fa extend0(size_t i) const {
        assert(i < N);
        //return bounds0(i).size();
	return zero; // FIXME: no yet implemented
      }

      /*! intersect 4 OBBs with single ray */
      __forceinline size_t intersect(const Precalculations& pre,
				     const sse3f& ray_org, const sse3f& ray_dir, 
				     const ssef& tnear, const ssef& tfar, const float time, ssef& dist)
      {
	const ssef t0 = ssef(1.0f)-time, t1 = time;

	const AffineSpaceSSE3f xfm = space0;
	const sse3f b0_lower = zero;
	const sse3f b0_upper = one;
	const sse3f lower = t0*b0_lower + t1*b1.lower;
	const sse3f upper = t0*b0_upper + t1*b1.upper;
	
	const BBoxSSE3f bounds(lower,upper);
	const sse3f dir = xfmVector(xfm,ray_dir);
	const sse3f rdir = rcp_safe(dir); 
	const sse3f org = xfmPoint(xfm,ray_org);
	
	const sse3f tLowerXYZ = (bounds.lower - org) * rdir;
	const sse3f tUpperXYZ = (bounds.upper - org) * rdir;
	
#if defined(__SSE4_1__)
	const ssef tNearX = mini(tLowerXYZ.x,tUpperXYZ.x);
	const ssef tNearY = mini(tLowerXYZ.y,tUpperXYZ.y);
	const ssef tNearZ = mini(tLowerXYZ.z,tUpperXYZ.z);
	const ssef tFarX  = maxi(tLowerXYZ.x,tUpperXYZ.x);
	const ssef tFarY  = maxi(tLowerXYZ.y,tUpperXYZ.y);
	const ssef tFarZ  = maxi(tLowerXYZ.z,tUpperXYZ.z);
	const ssef tNear  = max(tnear, tNearX,tNearY,tNearZ);
	const ssef tFar   = min(tfar,  tFarX ,tFarY ,tFarZ );
	const sseb vmask = tNear <= tFar;
	dist = tNear;
	return movemask(vmask);
#else
	const ssef tNearX = min(tLowerXYZ.x,tUpperXYZ.x);
	const ssef tNearY = min(tLowerXYZ.y,tUpperXYZ.y);
	const ssef tNearZ = min(tLowerXYZ.z,tUpperXYZ.z);
	const ssef tFarX  = max(tLowerXYZ.x,tUpperXYZ.x);
	const ssef tFarY  = max(tLowerXYZ.y,tUpperXYZ.y);
	const ssef tFarZ  = max(tLowerXYZ.z,tUpperXYZ.z);
	const ssef tNear = max(tnear, tNearX,tNearY,tNearZ);
	const ssef tFar  = min(tfar,  tFarX ,tFarY ,tFarZ );
	const sseb vmask = tNear <= tFar;
	dist = tNear;
	return movemask(vmask);
#endif
      }

      __forceinline size_t intersect(const sse3f& ray_org, const sse3f& ray_dir, 
				     const ssef& tnear, const ssef& tfar, const float time, ssef& dist) { return 0; }

    public:
      AffineSpaceSSE3f space0;   
      //BBoxSSE3f b0;
      BBoxSSE3f b1;
    };

    struct NodeDualSpaceMB : public BaseNode
    {
      struct Precalculations {
	__forceinline Precalculations (const Ray& ray) {}
      };

      /*! Clears the node. */
      __forceinline void clear() 
      {
        space0 = space1 = one;
        //t0s0.lower = t0s0.upper = Vec3fa(nan);
        t1s0_t0s1.lower = t1s0_t0s1.upper = Vec3fa(nan);
        //t1s1.lower = t1s1.upper = Vec3fa(nan);
        BaseNode::clear();
      }

      static __forceinline AffineSpace3fa normalizeSpace(const AffineSpace3fa& other, const BBox3fa& bounds)
      {
	AffineSpace3fa space = other;
        space.p -= bounds.lower;
        space = AffineSpace3fa::scale(1.0f/max(Vec3fa(1E-19f),bounds.upper-bounds.lower))*space;
	return space;
      }

      /*! Sets spaces. */
      __forceinline void set(size_t i, const AffineSpace3fa& s0, const AffineSpace3fa& s1) 
      {
        assert(i < N);

        space0.l.vx.x[i] = s0.l.vx.x; space0.l.vx.y[i] = s0.l.vx.y; space0.l.vx.z[i] = s0.l.vx.z; 
        space0.l.vy.x[i] = s0.l.vy.x; space0.l.vy.y[i] = s0.l.vy.y; space0.l.vy.z[i] = s0.l.vy.z;
        space0.l.vz.x[i] = s0.l.vz.x; space0.l.vz.y[i] = s0.l.vz.y; space0.l.vz.z[i] = s0.l.vz.z; 
        space0.p   .x[i] = s0.p   .x; space0.p   .y[i] = s0.p   .y; space0.p   .z[i] = s0.p   .z; 

        space1.l.vx.x[i] = s1.l.vx.x; space1.l.vx.y[i] = s1.l.vx.y; space1.l.vx.z[i] = s1.l.vx.z;
        space1.l.vy.x[i] = s1.l.vy.x; space1.l.vy.y[i] = s1.l.vy.y; space1.l.vy.z[i] = s1.l.vy.z;
        space1.l.vz.x[i] = s1.l.vz.x; space1.l.vz.y[i] = s1.l.vz.y; space1.l.vz.z[i] = s1.l.vz.z;
        space1.p   .x[i] = s1.p   .x; space1.p   .y[i] = s1.p   .y; space1.p   .z[i] = s1.p   .z; 
      }

      /*! Sets bounding boxes. */
      __forceinline void set(size_t i, const BBox3fa& a, const BBox3fa& b, const BBox3fa& c)
      {
        assert(i < N);

        //t0s0.lower.x[i] = a.lower.x; t0s0.lower.y[i] = a.lower.y; t0s0.lower.z[i] = a.lower.z;
        //t0s0.upper.x[i] = a.upper.x; t0s0.upper.y[i] = a.upper.y; t0s0.upper.z[i] = a.upper.z;

        t1s0_t0s1.lower.x[i] = b.lower.x; t1s0_t0s1.lower.y[i] = b.lower.y; t1s0_t0s1.lower.z[i] = b.lower.z;
        t1s0_t0s1.upper.x[i] = b.upper.x; t1s0_t0s1.upper.y[i] = b.upper.y; t1s0_t0s1.upper.z[i] = b.upper.z;

        //t1s1.lower.x[i] = c.lower.x; t1s1.lower.y[i] = c.lower.y; t1s1.lower.z[i] = c.lower.z;
        //t1s1.upper.x[i] = c.upper.x; t1s1.upper.y[i] = c.upper.y; t1s1.upper.z[i] = c.upper.z;
      }

      /*! Sets ID of child. */
      __forceinline void set(size_t i, const NodeRef& childID) {
        //Node::set(i,childID);
	assert(i < N);
	children[i] = childID;
      }

      /*! Returns bounds of specified child. */
      __forceinline const BBox3fa bounds0(const size_t i) const { 
        assert(i < N);
        /*const Vec3fa lower(t0s0.lower.x[i],t0s0.lower.y[i],t0s0.lower.z[i]);
        const Vec3fa upper(t0s0.upper.x[i],t0s0.upper.y[i],t0s0.upper.z[i]);
        return BBox3fa(lower,upper);*/
	return empty; // FIXME: not implemented yet
      }

      /*! Returns the extend of the bounds of the ith child */
      __forceinline Vec3fa extend0(size_t i) const {
        assert(i < N);
        //return bounds0(i).size(); // FIXME: not implemented yet
	return zero;
      }

      /*! Returns reference to specified child */
      __forceinline       NodeRef& child(size_t i)       { assert(i<N); return children[i]; }
      __forceinline const NodeRef& child(size_t i) const { assert(i<N); return children[i]; }

      /*! intersect 4 OBBs with single ray */
      __forceinline size_t intersect(const Precalculations& pre,
				     const sse3f& ray_org, const sse3f& ray_dir, 
				     const ssef& tnear, const ssef& tfar, const float time, ssef& dist)
      {
	const ssef t0 = ssef(1.0f)-time, t1 = time;
	const sse3f t0s0_lower = zero;
	const sse3f t0s0_upper = one;
	const sse3f t1s1_lower = zero;
	const sse3f t1s1_upper = one;

	const AffineSpaceSSE3f xfm = t0*space0 + t1*space1;
	const sse3f lower = t0*t0*t0s0_lower + t0*t1*t1s0_t0s1.lower + t1*t1*t1s1_lower;
	const sse3f upper = t0*t0*t0s0_upper + t0*t1*t1s0_t0s1.upper + t1*t1*t1s1_upper;

	const BBoxSSE3f bounds(lower,upper);
	const sse3f dir = xfmVector(xfm,ray_dir);
	const sse3f rdir = rcp_safe(dir); 
	const sse3f org = xfmPoint(xfm,ray_org);
	
	const sse3f tLowerXYZ = (bounds.lower - org) * rdir;
	const sse3f tUpperXYZ = (bounds.upper - org) * rdir;
	
#if defined(__SSE4_1__)
	const ssef tNearX = mini(tLowerXYZ.x,tUpperXYZ.x);
	const ssef tNearY = mini(tLowerXYZ.y,tUpperXYZ.y);
	const ssef tNearZ = mini(tLowerXYZ.z,tUpperXYZ.z);
	const ssef tFarX  = maxi(tLowerXYZ.x,tUpperXYZ.x);
	const ssef tFarY  = maxi(tLowerXYZ.y,tUpperXYZ.y);
	const ssef tFarZ  = maxi(tLowerXYZ.z,tUpperXYZ.z);
	const ssef tNear  = max(tnear, tNearX,tNearY,tNearZ);
	const ssef tFar   = min(tfar,  tFarX ,tFarY ,tFarZ );
	const sseb vmask = tNear <= tFar;
	dist = tNear;
	return movemask(vmask);
#else
	const ssef tNearX = min(tLowerXYZ.x,tUpperXYZ.x);
	const ssef tNearY = min(tLowerXYZ.y,tUpperXYZ.y);
	const ssef tNearZ = min(tLowerXYZ.z,tUpperXYZ.z);
	const ssef tFarX  = max(tLowerXYZ.x,tUpperXYZ.x);
	const ssef tFarY  = max(tLowerXYZ.y,tUpperXYZ.y);
	const ssef tFarZ  = max(tLowerXYZ.z,tUpperXYZ.z);
	const ssef tNear = max(tnear, tNearX,tNearY,tNearZ);
	const ssef tFar  = min(tfar,  tFarX ,tFarY ,tFarZ );
	const sseb vmask = tNear <= tFar;
	dist = tNear;
	return movemask(vmask);
#endif
      }

      __forceinline size_t intersect(const sse3f& ray_org, const sse3f& ray_dir, 
				     const ssef& tnear, const ssef& tfar, const float time, ssef& dist) { return 0; }

    public:
      AffineSpaceSSE3f space0;   
      AffineSpaceSSE3f space1;   
      //BBoxSSE3f t0s0;
      BBoxSSE3f t1s0_t0s1;
      //BBoxSSE3f t1s1;
    };  

    
    struct NodeConeMB : public BaseNode
    {
      struct Precalculations 
      {
	__forceinline Precalculations (const Ray& ray)
	  : depth_scale(rsqrt(dot(ray.dir,ray.dir))), ray_space(frame(depth_scale*ray.dir).transposed()) {}
	
	float depth_scale;
	LinearSpace3fa ray_space;
      };

      /*! Clears the node. */
      __forceinline void clear() 
      {
	v0t0 = v0t1 = ssef(nan);
	v1t0 = v1t1 = ssef(nan);
	rt0 =  rt1 = float(nan);
        BaseNode::clear();
      }

      /*! Sets bounding boxes. */
      __forceinline void set(size_t i, const Vec3fa& v0t0, const Vec3fa& v0t1, const Vec3fa& v1t0, const Vec3fa& v1t1, const float rt0, const float rt1)
      {
        assert(i < N);
	this->v0t0.x[i] = v0t0.x; this->v0t0.y[i] = v0t0.y; this->v0t0.z[i] = v0t0.z;
	this->v0t1.x[i] = v0t1.x; this->v0t1.y[i] = v0t1.y; this->v0t1.z[i] = v0t1.z;
	this->v1t0.x[i] = v1t0.x; this->v1t0.y[i] = v1t0.y; this->v1t0.z[i] = v1t0.z;
	this->v1t1.x[i] = v1t1.x; this->v1t1.y[i] = v1t1.y; this->v1t1.z[i] = v1t1.z;
	this->rt0[i] = rt0;
	this->rt1[i] = rt1;
      }

      /*! Sets ID of child. */
      __forceinline void set(size_t i, const NodeRef& childID) {
        //Node::set(i,childID);
	assert(i < N);
	children[i] = childID;
      }

      /*! Returns bounds of specified child. */
      __forceinline const BBox3fa bounds0(const size_t i) const { 
        assert(i < N);
        /*const Vec3fa lower(t0s0.lower.x[i],t0s0.lower.y[i],t0s0.lower.z[i]);
        const Vec3fa upper(t0s0.upper.x[i],t0s0.upper.y[i],t0s0.upper.z[i]);
        return BBox3fa(lower,upper);*/
	return empty; // FIXME: not implemented yet
      }

      /*! Returns the extend of the bounds of the ith child */
      __forceinline Vec3fa extend0(size_t i) const {
        assert(i < N);
        //return bounds0(i).size(); // FIXME: not implemented yet
	return zero;
      }

      /*! Returns reference to specified child */
      __forceinline       NodeRef& child(size_t i)       { assert(i<N); return children[i]; }
      __forceinline const NodeRef& child(size_t i) const { assert(i<N); return children[i]; }

      /*! intersect 4 OBBs with single ray */
      __forceinline size_t intersect(const Precalculations& pre, 
				     const sse3f& ray_org, const sse3f& ray_dir, 
				     const ssef& tnear, const ssef& tfar, const float time, ssef& dist)
      {
	const ssef t0 = ssef(1.0f)-time, t1 = time;
	const sse3f v0 = t0*v0t0 + t1*v0t1;
	const sse3f v1 = t0*v1t0 + t1*v1t1;
	const ssef  r  = t0*rt0 + t1*rt1;
	const sse3f p0 = xfmVector(LinearSpaceSSE3f(pre.ray_space),v0-ray_org); 
	const sse3f p1 = xfmVector(LinearSpaceSSE3f(pre.ray_space),v1-ray_org);
	const sse3f v = p1-p0;
	const sse3f w = -p0;
	const ssef d0 = w.x*v.x + w.y*v.y;
	const ssef d1 = v.x*v.x + v.y*v.y;
	const ssef u = clamp(d0*rcp(d1),ssef(zero),ssef(one));
	const sse3f p = p0 + u*v;
	//const ssef t = p.z*pre.depth_scale;
	const ssef d2 = p.x*p.x + p.y*p.y; 
	//const ssef r = p.w;
	const ssef r2 = r*r;
	sseb valid = d2 <= r2; // & avxf(ray.tnear) < t & t < avxf(ray.tfar);*/

	dist = 0.0f;
	return movemask(valid);
      }

      __forceinline size_t intersect(const sse3f& ray_org, const sse3f& ray_dir, 
				     const ssef& tnear, const ssef& tfar, const float time, ssef& dist) { return 0; }

    public:
      sse3f v0t0, v0t1;
      sse3f v1t0, v1t1;
      ssef rt0, rt1;
    };  

    /*! swap the children of two nodes */
    __forceinline static void swap(Node* a, size_t i, Node* b, size_t j)
    {
      assert(i<N && j<N);
      std::swap(a->children[i],b->children[j]);
      std::swap(a->lower_x[i],b->lower_x[j]);
      std::swap(a->lower_y[i],b->lower_y[j]);
      std::swap(a->lower_z[i],b->lower_z[j]);
      std::swap(a->upper_x[i],b->upper_x[j]);
      std::swap(a->upper_y[i],b->upper_y[j]);
      std::swap(a->upper_z[i],b->upper_z[j]);
    }

    /*! compacts a node (moves empty children to the end) */
    __forceinline static void compact(Node* a)
    {
      /* find right most filled node */
      ssize_t j=N;
      for (j=j-1; j>=0; j--)
        if (a->child(j) != emptyNode)
          break;

      /* replace empty nodes with filled nodes */
      for (ssize_t i=0; i<j; i++) {
        if (a->child(i) == emptyNode) {
          a->swap(i,j);
          for (j=j-1; j>i; j--)
            if (a->child(j) != emptyNode)
              break;
        }
      }
    }

    /*! compacts a node (moves empty children to the end) */
    __forceinline static void compact(NodeMB* a)
    {
      /* find right most filled node */
      ssize_t j=N;
      for (j=j-1; j>=0; j--)
        if (a->child(j) != emptyNode)
          break;

      /* replace empty nodes with filled nodes */
      for (ssize_t i=0; i<j; i++) {
        if (a->child(i) == emptyNode) {
          a->swap(i,j);
          for (j=j-1; j>i; j--)
            if (a->child(j) != emptyNode)
              break;
        }
      }
    }

  public:

    /*! BVH4 default constructor. */
    BVH4 (const PrimitiveType& primTy, Scene* scene, bool listMode);

    /*! BVH4 destruction */
    ~BVH4 ();

    /*! BVH4 instantiations */
    static Accel* BVH4Triangle1vMB(Scene* scene);
    static Accel* BVH4Triangle4vMB(Scene* scene);

    static Accel* BVH4Bezier1v(Scene* scene);
    static Accel* BVH4Bezier1i(Scene* scene);
    
    static Accel* BVH4OBBBezier1v(Scene* scene, bool highQuality);
    static Accel* BVH4OBBBezier1i(Scene* scene, bool highQuality);
    static Accel* BVH4OBBBezier1iMB(Scene* scene, bool highQuality);

    static Accel* BVH4Triangle1(Scene* scene);
    static Accel* BVH4Triangle4(Scene* scene);
    static Accel* BVH4Triangle8(Scene* scene);
    static Accel* BVH4Triangle1v(Scene* scene);
    static Accel* BVH4Triangle4v(Scene* scene);
    static Accel* BVH4Triangle4i(Scene* scene);
    static Accel* BVH4SubdivPatch1(Scene* scene);
    static Accel* BVH4SubdivPatch1Cached(Scene* scene);
    static Accel* BVH4SubdivGrid(Scene* scene);
    static Accel* BVH4SubdivGridEager(Scene* scene);
    static Accel* BVH4SubdivGridLazy(Scene* scene);
    static Accel* BVH4UserGeometry(Scene* scene);
    
    static Accel* BVH4BVH4Triangle1Morton(Scene* scene);
    static Accel* BVH4BVH4Triangle1ObjectSplit(Scene* scene);
    static Accel* BVH4BVH4Triangle4ObjectSplit(Scene* scene);
    static Accel* BVH4BVH4Triangle1vObjectSplit(Scene* scene);
    static Accel* BVH4BVH4Triangle4vObjectSplit(Scene* scene);
    static Accel* BVH4BVH4Triangle4iObjectSplit(Scene* scene);

    static Accel* BVH4Triangle1SpatialSplit(Scene* scene);
    static Accel* BVH4Triangle4SpatialSplit(Scene* scene);
    static Accel* BVH4Triangle8SpatialSplit(Scene* scene);
    static Accel* BVH4Triangle1ObjectSplit(Scene* scene);
    static Accel* BVH4Triangle4ObjectSplit(Scene* scene);
    static Accel* BVH4Triangle8ObjectSplit(Scene* scene);
    static Accel* BVH4Triangle1vObjectSplit(Scene* scene);
    static Accel* BVH4Triangle4vObjectSplit(Scene* scene);
    static Accel* BVH4Triangle4iObjectSplit(Scene* scene);

    static Accel* BVH4Triangle1ObjectSplit(TriangleMesh* mesh);
    static Accel* BVH4Triangle4ObjectSplit(TriangleMesh* mesh);
    static Accel* BVH4Triangle1vObjectSplit(TriangleMesh* mesh);
    static Accel* BVH4Triangle4vObjectSplit(TriangleMesh* mesh);
    static Accel* BVH4Triangle4Refit(TriangleMesh* mesh);

    /*! initializes the acceleration structure */
    void init (size_t nodeSize, size_t numPrimitives, size_t numThreads);

    /*! Clears the barrier bits of a subtree. */
    void clearBarrier(NodeRef& node);

    /*! Propagate bounds for time t0 and time t1 up the tree. */
    std::pair<BBox3fa,BBox3fa> refit(Scene* scene, NodeRef node);

    LinearAllocatorPerThread alloc;

    FastAllocator alloc2;

    void *data_mem; /* additional memory, currently used for subdivpatch1cached memory */
    size_t size_data_mem;

    __forceinline Node* allocNode(LinearAllocatorPerThread::ThreadAllocator& thread) {
      Node* node = (Node*) thread.malloc(sizeof(Node),1 << alignment); node->clear(); return node; // FIXME: why 16 bytes aligned and not 64 bytes?
    }

    __forceinline NodeMB* allocNodeMB(LinearAllocatorPerThread::ThreadAllocator& thread) {
      NodeMB* node = (NodeMB*) thread.malloc(sizeof(NodeMB),1 << alignment); node->clear(); return node;
    }

    /*! allocates a new unaligned node */
    __forceinline UnalignedNode* allocUnalignedNode(LinearAllocatorPerThread::ThreadAllocator& thread) {
      UnalignedNode* node = (UnalignedNode*) thread.malloc(sizeof(UnalignedNode),1 << alignment); node->clear(); return node;
    }

    /*! allocates a new unaligned node */
    __forceinline UnalignedNodeMB* allocUnalignedNodeMB(LinearAllocatorPerThread::ThreadAllocator& thread) {
      UnalignedNodeMB* node = (UnalignedNodeMB*) thread.malloc(sizeof(UnalignedNodeMB),1 << alignment); node->clear(); return node;
    }

    __forceinline char* allocPrimitiveBlocks(LinearAllocatorPerThread::ThreadAllocator& thread, size_t num) {
      return (char*) thread.malloc(num*primTy.bytes,1 << alignment);
    }

    /*! Encodes a node */
    static __forceinline NodeRef encodeNode2(Node* node) {  // FIXME: make all static
      assert(!((size_t)node & align_mask)); 
      return NodeRef((size_t) node);
    }

    /*! Encodes a node */
    static __forceinline NodeRef encodeNode(void* node) {  // FIXME: template these functions
      assert(!((size_t)node & align_mask)); 
      return NodeRef((size_t) node);
    }

    /*! Encodes a node */
    static __forceinline NodeRef encodeNode(NodeMB* node) { 
      assert(!((size_t)node & align_mask)); 
      return NodeRef((size_t) node | tyNodeMB);
    }

    /*! Encodes an unaligned node */
    static __forceinline NodeRef encodeNode(UnalignedNode* node) { 
      return NodeRef((size_t) node | tyUnalignedNode);
    }

    /*! Encodes an unaligned motion blur node */
    static __forceinline NodeRef encodeNode(UnalignedNodeMB* node) { 
      return NodeRef((size_t) node |  tyUnalignedNodeMB);
    }
    
    /*! Encodes a leaf */
    static __forceinline NodeRef encodeLeaf(void* tri, size_t num) {
      assert(!((size_t)tri & align_mask)); 
      return NodeRef((size_t)tri | (tyLeaf+min(num,(size_t)maxLeafBlocks)));
    }

    /*! Encodes a leaf */
    static __forceinline NodeRef encodeTypedLeaf(void* ptr, size_t ty) {
      assert(!((size_t)ptr & align_mask)); 
      return NodeRef((size_t)ptr | (tyLeaf+ty));
    }
    
  public:
    
    /*! calculates the amount of bytes allocated */
    size_t bytesAllocated() {
      return alloc.bytes();
    }

  public:
    const PrimitiveType& primTy;       //!< primitive type stored in the BVH
    Scene* scene;                      //!< scene pointer
    bool listMode;                     //!< true if number of leaf items not encoded in NodeRef
    NodeRef root;                      //!< Root node
    size_t numPrimitives;
    size_t numVertices;

    /*! data arrays for fast builders */
  public:
    std::vector<BVH4*> objects;
  };

   // FIXME: move the below code to somewhere else
  typedef void (*createTriangleMeshAccelTy)(TriangleMesh* mesh, BVH4*& accel, Builder*& builder); 
  typedef Builder* (*BVH4BuilderTopLevelFunc)(BVH4* accel, Scene* scene, const createTriangleMeshAccelTy createTriangleMeshAccel);

#define DECLARE_TOPLEVEL_BUILDER(symbol)                                         \
  namespace isa   { extern Builder* symbol(BVH4* accel, Scene* scene, const createTriangleMeshAccelTy createTriangleMeshAccel); } \
  namespace sse41 { extern Builder* symbol(BVH4* accel, Scene* scene, const createTriangleMeshAccelTy createTriangleMeshAccel); } \
  namespace avx   { extern Builder* symbol(BVH4* accel, Scene* scene, const createTriangleMeshAccelTy createTriangleMeshAccel); } \
  namespace avx2  { extern Builder* symbol(BVH4* accel, Scene* scene, const createTriangleMeshAccelTy createTriangleMeshAccel); } \
  BVH4BuilderTopLevelFunc symbol;
}
