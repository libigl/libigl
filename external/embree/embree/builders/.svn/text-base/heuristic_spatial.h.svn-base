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

#ifndef __EMBREE_HEURISTIC_SPATIAL_H__
#define __EMBREE_HEURISTIC_SPATIAL_H__

#include "../geometry/geometry.h"
#include "../geometry/triangle_mesh.h"
#include "primref.h"

namespace embree
{
  /* Single threaded combined object and spatial binner. Performs the
   * same object binning procedure as the HeuristicBinning class. In
   * addition tries spatial splits, by splitting each dimension
   * spatially in the center. In contrast to object binning, spatial
   * binning potentially cuts triangles along the splitting plane and
   * sortes the corresponding parts into the left and right set. */
  
  template<int logBlockSize>
    class HeuristicSpatial
  {
    /*! Maximal number of bins for object binning. */
    static const size_t maxBins = 32;

    /*! Maximal number of times a spatial split is allowed to fail. */
    static const size_t maxFailed = 2;

    /*! Maximal amount of triangle duplications allowed */
    static const size_t duplicationPercentage = 100;  //!< up to 100% additional triangles allowed
    
    /*! Primitives whose right bound is smaller go to the left size. */
    static const float leftSplitPos;
    
    /*! Primitives whose left bound is larger go to the right size. */
    static const float rightSplitPos;

    /*! Clipping code. Splits a triangle with a clipping plane into two halves. */
    static std::pair<PrimRef,PrimRef> splitPrimRef(const PrimRef& ref, int dim, float pos, const TriangleMesh* geom);

  public:

    /*! We build the tree breadth first to better distribute triangle duplications over scene. */
    static const bool depthFirst = false;
    
    /*! stores bounding information for a set of primitives */
    class PrimInfo
    {
    public:
      __forceinline PrimInfo () 
        : num(0), numFailed(0), geomBounds(empty), centBounds(empty), duplications(NULL) {}

      __forceinline PrimInfo (size_t num, const BBox3f& geomBounds) 
        : num(num), numFailed(0), geomBounds(geomBounds), centBounds(geomBounds) 
      {
        duplications = new Atomic(size_t(duplicationPercentage*double(num)/100.0));
      }
      
      __forceinline PrimInfo (size_t num, const BBox3f& geomBounds, const BBox3f& centBounds)
        : num(num), numFailed(0), geomBounds(geomBounds), centBounds(centBounds) 
      {
        duplications = new Atomic(size_t(duplicationPercentage*double(num)/100.0));
      }

      __forceinline PrimInfo (size_t num, int numFailed, const BBox3f& geomBounds, const BBox3f& centBounds, Atomic* duplications) 
        : num(num), numFailed(numFailed), geomBounds(geomBounds), centBounds(centBounds), duplications(duplications) {}
      
      /*! returns the number of primitives */
      __forceinline size_t size() const { 
        return num; 
      }

      /*! return the surface area heuristic when creating a leaf */
      __forceinline float sah () const { 
        return halfArea(geomBounds)*blocks(num); 
      }

      /*! clears global data */
      __forceinline void clear () { 
        delete duplications; duplications = NULL;
      }

      /*! stream output */
      friend std::ostream& operator<<(std::ostream& cout, const PrimInfo& pinfo) {
        return cout << "PrimInfo { num = " << pinfo.num << ", failed = " << pinfo.numFailed << 
          ", geomBounds = " << pinfo.geomBounds << ", centBounds = " << pinfo.centBounds << "}";
      }
      
    public:
      size_t num;         //!< number of primitives
      int    numFailed;   //!< number of times a spatial split failed
      BBox3f geomBounds;  //!< geometry bounds of primitives
      BBox3f centBounds;  //!< centroid bounds of primitives
      Atomic* duplications; //!< maximal number of duplications allowed
    };
    
    /*! mapping from bounding boxes to bins and spatial mapping */
    class Mapping
    {
    public:
      
      /*! default constructor */
      __forceinline Mapping () {}
      
      /*! construct from centroid bounds and number of primitives */
      __forceinline Mapping (const PrimInfo& pinfo)
      {
        num   = min(maxBins,size_t(4.0f + 0.05f*pinfo.num));
        ofs   = pinfo.centBounds.lower;
        scale = rcp(pinfo.centBounds.upper - pinfo.centBounds.lower) * Vector3f((float)num);

        center = 0.5f*center2(pinfo.geomBounds);
        bleft  = pinfo.geomBounds.lower + leftSplitPos *pinfo.geomBounds.size();
        bright = pinfo.geomBounds.lower + rightSplitPos*pinfo.geomBounds.size();
      }
      
      /*! returns number of bins */
      __forceinline size_t size() const { return num; }
      
      /*! Computes the bin numbers for each dimension for a box. */
      __forceinline Vector3i bin_unsafe(const BBox3f& box) const {
        return Vector3i((center2(box) - ofs)*scale-Vector3f(0.5f));
      }
      
      /*! Computes the bin numbers for each dimension for a box. */
      __forceinline Vector3i bin(const BBox3f& box) const {
        return clamp(bin_unsafe(box),Vector3i(0),Vector3i(int(num-1)));
      }

      /*! Computes if object is on the left of a spatial split. */
      __forceinline Vector3b left(const BBox3f& box) const {
        return le_mask(box.upper,bleft);
      }

      /*! Computes if object is on the right of a spatial split. */
      __forceinline Vector3b right(const BBox3f& box) const {
        return ge_mask(box.lower,bright);
      }

      /*! for spatial split */
    public:
      Vector3f center;     //!< center for spatial splits
      Vector3f bleft;      //!< if upper bound of a primitive is smaller it goes to the left
      Vector3f bright;     //!< if lower bound of a primitive is larger it goes to the left

      /*! for object splits */
    public:
      size_t num;       //!< number of bins to use
      Vector3f ofs;        //!< offset to compute bin
      Vector3f scale;      //!< scaling factor to compute bin
    };
    
    /*! Stores information about an object split. */
    class Split
    {
    public:
      
      /*! create an invalid split by default */
      __forceinline Split () 
        : objectSAH(inf), dim(0), pos(0), spatialSAH(inf), sdim(0), numFailed(0), geom(NULL) {}
      
      /*! return SAH cost of performing a split */
      __forceinline float sah() const { return min(objectSAH,spatialSAH); } 
      
      /*! return true if this is a spatial split that introduces duplications */
      __forceinline bool spatial() const { 
        return spatialSAH < objectSAH;
      }

      /*! tests if a primitive belongs to the left side of the split */
      __forceinline int left(const PrimRef& prim) const {
        assert(objectSAH <= spatialSAH);
        return mapping.bin_unsafe(prim.bounds())[dim] < pos;
      }

      /*! splitting the primitive if required */
      __forceinline void split(const PrimRef& prim, PrimRef& lprim_o, PrimRef& rprim_o) const 
      {
        const BBox3f bounds = prim.bounds();
        const bool left  = mapping.left(bounds )[sdim];
        const bool right = mapping.right(bounds)[sdim];

        /*! Choose better side for primitives that can be put left or right. */
        if (left && right) {
          if (prim.upper[sdim]-mapping.bright[sdim] > mapping.bleft[sdim]-prim.lower[sdim]) {
            lprim_o = PrimRef(empty,size_t(-1)); rprim_o = prim;
          } else {
            lprim_o = prim; rprim_o = PrimRef(empty,size_t(-1));
          }
        }
        /*! These definitely go to the left. */
        else if (left) {
          lprim_o = prim; rprim_o = PrimRef(empty,size_t(-1));
        }
        /*! These definitely go to the right. */
        else if (right) {
          lprim_o = PrimRef(empty,size_t(-1)); rprim_o = prim;
        }
        /*! These have to get split and put to left and right. */
        else {
          std::pair<PrimRef,PrimRef> pair = splitPrimRef(prim,sdim,mapping.center[sdim],geom);
          lprim_o = pair.first; rprim_o = pair.second;
        }
      }
      
      /*! stream output */
      friend std::ostream& operator<<(std::ostream& cout, const Split& split) {
        return cout << "Split { objectSAH = " << split.objectSAH << 
          ", dim = " << split.dim << 
          ", pos = " << split.pos << 
          ", spatialSAH = " << split.spatialSAH << 
          ", sdim = " << split.sdim << 
          ", linfo = " << split.linfo <<
          ", rinfo = " << split.rinfo << "}";
      }
      
      /*! best object split */
    public:
      float objectSAH;    //!< SAH cost of performing best object split
      int dim;            //!< best object split dimension
      int pos;            //!< best object split position

      /*! best spatial split */
    public:
      float spatialSAH;    //!< SAH cost of performing best spatial split 
      int sdim;            //!< best spatial split dimension
      int numFailed;       //!< number of times a spatial split failed
      
    public:
      Mapping mapping;    //!< Mapping to bins
      PrimInfo linfo;     //!< Left geometry information
      PrimInfo rinfo;     //!< Right geometry information

    public:
      const TriangleMesh* geom;       //!< geometry
    };
    
    /*! default constructor */
    __forceinline HeuristicSpatial() {}
    
    /*! construction from geometry info */
    __forceinline HeuristicSpatial(const PrimInfo& pinfo, const RTCGeometry* geom)
      : pinfo(pinfo), mapping(pinfo), geom((const TriangleMesh*)geom) { clear(); }

    /*! clear the binner */
    __forceinline void clear()
    {
      for (size_t i=0; i<mapping.size(); i++) {
        counts[i] = 0;
        geomBounds[i][0] = geomBounds[i][1] = geomBounds[i][2] = empty;
        centBounds[i][0] = centBounds[i][1] = centBounds[i][2] = empty;
      }
      
      lcounts = rcounts = 0;
      lgeomBounds[0] = lgeomBounds[1] = lgeomBounds[2] = empty;
      rgeomBounds[0] = rgeomBounds[1] = rgeomBounds[2] = empty;
      lcentBounds[0] = lcentBounds[1] = lcentBounds[2] = empty;
      rcentBounds[0] = rcentBounds[1] = rcentBounds[2] = empty;
    }
    
    /*! bin an array of primitives */
    void bin(const PrimRef* prim, size_t num);
    
    /*! merge multiple binners into one */
    static void reduce(const HeuristicSpatial binners[], size_t num, HeuristicSpatial& binner_o);
    
    /*! calculate the best possible split */
    void best(Split& split_o);

  private:

    /*! calculate the best object split */
    void best_object_split(Split& split_o);

    /*! calculate the best  spatial split */
    void best_spatial_split(Split& split_o);

    /*! update info for left and right primitives */
    void update_object_split(Split& split_o);

    /*! update info for left and right primitives */
    bool update_spatial_split(Split& split_o);
    
    /*! Compute the number of blocks occupied for each dimension. */
    __forceinline static Vector3i blocks(const Vector3i& a) { return (a+Vector3i((1 << logBlockSize)-1)) >> logBlockSize; }
    
    /*! Compute the number of blocks occupied in one dimension. */
    __forceinline static int  blocks(size_t a) { return (int)((a+((1LL << logBlockSize)-1)) >> logBlockSize); }
    
  public:
    PrimInfo pinfo;                 //!< bounding information of geometry
    Mapping mapping;                //!< mapping from geometry to the bins
    
    /* counter and bounds for object binning */
  public:
    Vector3i   counts    [maxBins];    //< number of primitives mapped to bin
    BBox3f geomBounds[maxBins][4]; //< bounds for every bin in every dimension
    BBox3f centBounds[maxBins][4]; //< centroid bounds for every bin in every dimension

    /* counter and bounds for spatial binning */
  public:
    Vector3i lcounts, rcounts;        //< left and right count for spatial split
    BBox3f lgeomBounds[4];        //< left geometry bounds for spatial split
    BBox3f rgeomBounds[4];        //< right geometry bounds for spatial split
    BBox3f lcentBounds[4];        //< left centroid bounds for spatial split
    BBox3f rcentBounds[4];        //< right centroid bounds for spatial split
    
  public:
    const TriangleMesh* geom;       //!< geometry
  };
}

#endif
