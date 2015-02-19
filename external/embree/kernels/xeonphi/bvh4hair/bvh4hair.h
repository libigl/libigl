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

#include "common/alloc.h"
#include "common/accel.h"
#include "common/scene.h"
#include "geometry/primitive.h"

namespace embree
{
  class BVH4Hair : public AccelData 
  {
  public:
    /*! branching width of the tree */
    static const size_t N = 4;

    /*! Masks the bits that store the number of items per leaf. */
    static const unsigned int encodingBits     = 4;
    static const unsigned int offset_mask      = ((unsigned int)-1) << encodingBits;
    static const unsigned int leaf_shift       = 3;
    static const unsigned int leaf_mask        = 1<<leaf_shift;  
    static const unsigned int items_mask       = (1<<(leaf_shift-1))-1;
    static const unsigned int aux_flag_mask    = 1<<(leaf_shift-1);

    //static const unsigned int alignednode_mask = aux_flag_mask; 
    static const unsigned int alignednode_mask = 1 << (leaf_shift+1);

    
    /*! Maximal depth of the BVH. */
    static const size_t maxBuildDepth = 26;
    static const size_t maxBuildDepthLeaf = maxBuildDepth+6;
    static const size_t maxDepth = maxBuildDepth + maxBuildDepthLeaf;
    static const int travCost = 1;

    /*! Empty node */
    static const unsigned int emptyNode = leaf_mask;

    /*! Invalid node */
    static const unsigned int invalidNode = (unsigned int)-1; 

    struct NodeRef
    {
      /*! Default constructor */
      __forceinline NodeRef () {}

      /*! Construction from integer */
      __forceinline NodeRef (unsigned id) : _id(id) { }

      /*! Cast to unsigned */
      __forceinline operator unsigned int() const { return _id; }
     
      /*! checks if this is a leaf */
      __forceinline unsigned int isLeaf() const { return _id & leaf_mask; }

      /*! checks if this is a leaf */
      __forceinline unsigned int isLeaf(const unsigned int mask) const { return _id & mask; }
      
      /*! checks if this is a node */
      __forceinline unsigned isNode() const { return (_id & leaf_mask) == 0; }
      
      /*! returns node pointer */

      __forceinline       void* node(      void* base) const { return (      void*)((      char*)base + (size_t)_id); }
      __forceinline const void* node(const void* base) const { return (const void*)((const char*)base + (size_t)_id); }

      __forceinline unsigned int items() const {
	assert( (_id & items_mask)+1 <= 4);
        return (_id & items_mask)+1;
      }      
      
      /*! returns leaf pointer */
      template<unsigned int scale=4>
	__forceinline const char* leaf(const void* base, unsigned int& num) const {
        assert(isLeaf());
        num = items();
        return (const char*)base + (_id & offset_mask)*scale;
      }

      /*! returns leaf pointer */
	template<unsigned int scale=4>
	__forceinline const char* leaf(const void* base) const {
	  assert(isLeaf());
	  return (const char*)base + (_id & offset_mask)*scale;
	}

	  __forceinline unsigned int offset() const {
	    return _id & offset_mask;
	  }

      __forceinline unsigned int offsetIndex() const {
        return _id >> encodingBits;
      }

      
      __forceinline unsigned int &id() { return _id; }
    private:
      unsigned int _id;
    };


    struct __aligned(64) UnalignedNode
    {
      UnalignedNode()
	{
	  assert(sizeof(UnalignedNode) == 192);
	}

      static float identityMatrix[16];
      static float  invalidMatrix[16];


      struct NodeStruct {
        float x,y,z;           // x,y, and z coordinates of bounds
        NodeRef data;         
      } lower[4], upper[4];    

      char xfm[4][16];

      __forceinline mic_i getChildren() const
      {
	return load16i((int*)lower);
      }

      __forceinline mic_f getRow(size_t i) const
      {
	return load16f_int8(xfm[i]);
      }

      /*! Returns bounds of specified child. */
      __forceinline BBox3fa bounds(size_t i) const {
	assert( i < 4 );
        Vec3fa l = *(Vec3fa*)&lower[i];
        Vec3fa u = *(Vec3fa*)&upper[i];
        return BBox3fa(l,u);
      }

      __forceinline void setInvalid(size_t i) 
      {
	lower[i].x = pos_inf;
	lower[i].y = pos_inf;
	lower[i].z = pos_inf;
	lower[i].data = invalidNode;

	upper[i].x = neg_inf;
	upper[i].y = neg_inf;
	upper[i].z = neg_inf;
	upper[i].data = invalidNode;
		
      }

      __forceinline void setInvalid() 
      {
	for (size_t i=0;i<4;i++)
	  setInvalid(i);

	for (size_t y=0;y<4;y++)
	  for (size_t x=0;x<16;x++)
	    xfm[y][x] = 0;
      }
      
      template<int PFHINT>
	__forceinline void prefetchNode() const
	{
	  prefetch<PFHINT>((char*)this + 0 * 64);
	  prefetch<PFHINT>((char*)this + 1 * 64);
	  prefetch<PFHINT>((char*)this + 2 * 64);
	}


      __forceinline const char &c_matrix(const size_t row,
					 const size_t column,
					 const size_t matrixID) const
      {
	assert(matrixID < 4);
	assert(row < 3);
	assert(column < 3);
	return xfm[row][matrixID*4+column];
      } 

      __forceinline char &c_matrix(const size_t row,
				   const size_t column,
				   const size_t matrixID) 
      {
	assert(matrixID < 4);
	assert(row < 3);
	assert(column < 3);
	return xfm[row][matrixID*4+column];
      } 

      __forceinline float matrix(const size_t row,
				 const size_t column,
				 const size_t matrixID) const
      {
	assert(matrixID < 4);
	assert(row < 3);
	assert(column < 3);
	return ((float)c_matrix(row,column,matrixID)) * 1.0f/127.0f;
      } 

      

      __forceinline void setMatrix(const BBox3fa &b, const size_t m)
      {
	lower[m].x = b.lower.x;
	lower[m].y = b.lower.y;
	lower[m].z = b.lower.z;

	upper[m].x = b.upper.x;
	upper[m].y = b.upper.y;
	upper[m].z = b.upper.z;

	c_matrix(0,0,m) = 127;
	c_matrix(1,1,m) = 127;
	c_matrix(2,2,m) = 127;
      }
   
      __forceinline void setMatrix(const LinearSpace3fa &mat, BBox3fa &b, const size_t m)
      {
	lower[m].x = b.lower.x;
	lower[m].y = b.lower.y;
	lower[m].z = b.lower.z;

	upper[m].x = b.upper.x;
	upper[m].y = b.upper.y;
	upper[m].z = b.upper.z;

#if 0
	for (size_t i=0;i<3;i++)
	  {
	    c_matrix(0,i,m) = (char)(127.0f * mat.vx[i]);
	    c_matrix(1,i,m) = (char)(127.0f * mat.vy[i]);
	    c_matrix(2,i,m) = (char)(127.0f * mat.vz[i]);
	  }
#else
	const mic_f vx = broadcast4to16f(&mat.vx) * 127.0f;
	const mic_f vy = broadcast4to16f(&mat.vy) * 127.0f;
	const mic_f vz = broadcast4to16f(&mat.vz) * 127.0f;
	store4f_int8(&c_matrix(0,0,m),vx);
	store4f_int8(&c_matrix(1,0,m),vy);
	store4f_int8(&c_matrix(2,0,m),vz);

#endif
	

      }



      /*! Returns reference to specified child */
      __forceinline       NodeRef& child(size_t i)       { return lower[i].data; }
      __forceinline const NodeRef& child(size_t i) const { return lower[i].data; }

      __forceinline       NodeRef& child_ref(size_t i)       { return ((BVH4Hair::NodeRef*)lower)[i]; }
      __forceinline const NodeRef& child_ref(size_t i) const { return ((BVH4Hair::NodeRef*)lower)[i]; }
     
      __forceinline NodeRef *nodeRefPtr() const { return (NodeRef*)lower; }

    };

    
    struct __aligned(64) AlignedNode
    {
    public:

      AlignedNode() 
	{
	  assert( sizeof(AlignedNode) == 192 );
	}
      struct NodeStruct {
        float x,y,z;           // x,y, and z coordinates of bounds
        NodeRef data;          
      } lower[4], upper[4];    

      mic_i dummy;

      template<int PFHINT>
	__forceinline void prefetchNode() const
	{
	  prefetch<PFHINT>((char*)this + 0 * 64);
	  prefetch<PFHINT>((char*)this + 1 * 64);
	}

      /*! Returns bounds of specified child. */
      __forceinline BBox3fa bounds(size_t i) const {
	assert( i < 4 );
        Vec3fa l = *(Vec3fa*)&lower[i];
        Vec3fa u = *(Vec3fa*)&upper[i];
        return BBox3fa(l,u);
      }

      __forceinline mic_f lowerXYZ(size_t i) const {
	return broadcast4to16f(&lower[i]);
      }

      __forceinline mic_f upperXYZ(size_t i) const {
	return broadcast4to16f(&upper[i]);
      }

      __forceinline mic_i getChildren() const
      {
	return load16i((int*)lower);
      }

      __forceinline void setInvalid(size_t i)
      {
	lower[i].x = pos_inf;
	lower[i].y = pos_inf;
	lower[i].z = pos_inf;
	lower[i].data = (unsigned int)BVH4Hair::invalidNode;

	upper[i].x = neg_inf;
	upper[i].y = neg_inf;
	upper[i].z = neg_inf;
	upper[i].data = (unsigned int)BVH4Hair::invalidNode;
      }

      __forceinline void setInvalid()
      {
	for (size_t i=0;i<4;i++)
	  setInvalid(i);
      }


      __forceinline       NodeRef &child(size_t i)       { 
	return  lower[i].data;
      }
      __forceinline const NodeRef &child(size_t i) const { 
	return  lower[i].data;
      }

      __forceinline       NodeRef& child_ref(size_t i)       { return ((NodeRef*)lower)[i]; }
      __forceinline const NodeRef& child_ref(size_t i) const { return ((NodeRef*)lower)[i]; }

      __forceinline NodeRef *nodeRefPtr() const { return (NodeRef*)lower; }

      __forceinline void setMatrix(const BBox3fa &b, const size_t m)
      {
	lower[m].x = b.lower.x;
	lower[m].y = b.lower.y;
	lower[m].z = b.lower.z;

	upper[m].x = b.upper.x;
	upper[m].y = b.upper.y;
	upper[m].z = b.upper.z;

	lower[m].data = 0;
	upper[m].data = 0;

      }

      __forceinline void setMatrix(const LinearSpace3fa &mat, BBox3fa &b, const size_t m)
      {
	FATAL("not implemented");
      }

    };


    NodeRef root;                      //!< Root node (can also be a leaf).

    const PrimitiveType& primTy;       //!< primitive type stored in BVH
    void* geometry;                    //!< pointer to geometry for intersection
    UnalignedNode *unaligned_nodes;
    void *accel;
    size_t size_node;
    size_t size_accel;

    __forceinline       void* nodePtr()       { return unaligned_nodes; }
    __forceinline const void* nodePtr() const { return unaligned_nodes; }

    __forceinline       void* primitivesPtr()       { return accel; }
    __forceinline const void* primitivesPtr() const { return accel; }

    size_t bytes () const {
      return size_node+size_accel;
    }


  BVH4Hair(const PrimitiveType& primTy, void* geometry = NULL) : primTy(primTy), 
      geometry(geometry), 
      root(emptyNode), 
      accel(NULL),
      size_node(0),
      size_accel(0)
      {	
	assert( sizeof(UnalignedNode) == 192 );
	unaligned_nodes = NULL;
      }

    
    static Accel* BVH4HairBinnedSAH(Scene* scene);
    


  };

  __forceinline std::ostream& operator<<(std::ostream &o, const BVH4Hair::AlignedNode &n)
    {
      for (size_t i=0;i<4;i++)
	{
	  o << "lower: [" << n.lower[i].x << "," << n.lower[i].y << "," << n.lower[i].z << "," << n.lower[i].data << "] ";
	  o << "upper: [" << n.upper[i].x << "," << n.upper[i].y << "," << n.upper[i].z << "," << n.upper[i].data << "] ";
	  o << std::endl;
	  }
      return o;
    }


  __forceinline std::ostream &operator<<(std::ostream &o, const BVH4Hair::UnalignedNode &n)
    {
      o << "ptr " << (void*)&n << std::endl;
      for (size_t m=0;m<4;m++)
	{
	  o << "matrix " << m << ": " << std::endl;
	  for (size_t y=0;y<4;y++)
	    {
	      for (size_t x=0;x<3;x++)
		o << n.matrix(y,x,m) << " ";
	      o << std::endl;
	    }
	}
      o << "children: ";
      for (size_t m=0;m<4;m++)
	{
	  o << n.child(m) << " ";
	  if (n.child(m).isLeaf())
	    o << "(LEAF: index " << n.child(m).offsetIndex() << " items " << n.child(m).items() << ") ";
	  else
	    o << "(NODE) ";
	}
      o << std::endl;

      return o;
    } 

};
