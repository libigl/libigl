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

#include "heuristic_spatial.h"

namespace embree
{
  template<int logBlockSize>
  const size_t HeuristicSpatial<logBlockSize>::maxBins;

  template<int logBlockSize>
  const float HeuristicSpatial<logBlockSize>::leftSplitPos = 0.51f;
  
  template<int logBlockSize>
  const float HeuristicSpatial<logBlockSize>::rightSplitPos = 0.49f;
  
  template<int logBlockSize>
  __forceinline std::pair<PrimRef,PrimRef> HeuristicSpatial<logBlockSize>::splitPrimRef(const PrimRef& prim, int dim, float pos,
                                                                                        const TriangleMesh* mesh)
  {
    std::pair<BBox3f,BBox3f> pair(empty,empty);
    const TriangleMesh::Triangle& tri = mesh->triangle(prim.id());
    const Vector3f p0 = mesh->vertex(tri.v0);
    const Vector3f p1 = mesh->vertex(tri.v1);
    const Vector3f p2 = mesh->vertex(tri.v2);
    const Vector3f v[3] = { p0,p1,p2 };

    /* clip triangle to left and right box by processing all edges */
    Vector3f v1 = v[2];
    for (size_t i=0; i<3; i++)
    {
      Vector3f v0 = v1; v1 = v[i];
      float v0d = v0[dim], v1d = v1[dim];
      
      if (v0d <= pos) pair.first .grow(v0); // this point is on left side
      if (v0d >= pos) pair.second.grow(v0); // this point is on right side

      if ((v0d < pos && pos < v1d) || (v1d < pos && pos < v0d)) // the edge crosses the splitting location
      {
        assert((v1d-v0d) != 0.0f);
        float t = clamp((pos-v0d)/(v1d-v0d),0.0f,1.0f);
        Vector3f c = v0 + t*(v1-v0);
        pair.first .grow(c);
        pair.second.grow(c);
      }
    }
    assert(!pair.first .empty()); // happens if split does not hit triangle
    assert(!pair.second.empty()); // happens if split does not hit triangle

    return std::pair<PrimRef,PrimRef>(PrimRef(intersect(pair.first ,prim.bounds()), prim.id()),
                                      PrimRef(intersect(pair.second,prim.bounds()), prim.id()));
  }

  template<int logBlockSize>
  void HeuristicSpatial<logBlockSize>::bin(const PrimRef* prims, size_t num)
  {
    for (size_t i=0; i<num; i++)
    {
      
      /*! map even and odd primitive to bin */
      const BBox3f prim = prims[i].bounds(); const Vector3i bin = mapping.bin(prim); const Vector3f center = Vector3f(center2(prim));
      
      /*! increase bounds for bins for even primitive */
      const int b0 = bin.x; counts[b0][0]++; geomBounds[b0][0].grow(prim); centBounds[b0][0].grow(center);
      const int b1 = bin.y; counts[b1][1]++; geomBounds[b1][1].grow(prim); centBounds[b1][1].grow(center);
      const int b2 = bin.z; counts[b2][2]++; geomBounds[b2][2].grow(prim); centBounds[b2][2].grow(center);
      
      /*! calculate for each dimension if primitive goes to the left and right for spatial splits */
      const Vector3b left  = mapping.left (prim);
      const Vector3b right = mapping.right(prim);

      /*! Test spatial split in center of each dimension. */
      for (int maxDim=0; maxDim<3; maxDim++)
      {
        /*! Choose better side for primitives that can be put left or right. */
        if (left[maxDim] && right[maxDim]) {
          if (prim.upper[maxDim]-mapping.bright[maxDim] > mapping.bleft[maxDim]-prim.lower[maxDim]) {
            rgeomBounds[maxDim].grow(prim); rcentBounds[maxDim].grow(center); rcounts[maxDim]++;
          } else {
            lgeomBounds[maxDim].grow(prim); lcentBounds[maxDim].grow(center); lcounts[maxDim]++;
          }
        }
        /*! These definitely go to the left. */
        else if (left[maxDim]) {
          lgeomBounds[maxDim].grow(prim); lcentBounds[maxDim].grow(center); lcounts[maxDim]++;
        }
        /*! These definitely go to the right. */
        else if (right[maxDim]) {
          rgeomBounds[maxDim].grow(prim); rcentBounds[maxDim].grow(center); rcounts[maxDim]++;
        }
        /*! These have to get split and put to left and right. */
        else {
          std::pair<PrimRef,PrimRef> pair = splitPrimRef(prims[i],maxDim,mapping.center[maxDim],geom);
          lgeomBounds[maxDim].grow(pair.first. bounds()); lcentBounds[maxDim].grow(center2(pair.first. bounds())); lcounts[maxDim]++;
          rgeomBounds[maxDim].grow(pair.second.bounds()); rcentBounds[maxDim].grow(center2(pair.second.bounds())); rcounts[maxDim]++;
        }
      }
    }
  }

  template<int logBlockSize>
  void HeuristicSpatial<logBlockSize>::best_spatial_split(Split& split)
  {
    /* find best spatial split */
    for (int i=0; i<3; i++) 
    {
      size_t numFailed = pinfo.numFailed + (size_t(lcounts[i]) == pinfo.num && size_t(rcounts[i]) == pinfo.num);
      if (unlikely(numFailed > maxFailed)) continue;
      float u = float(ulp)*max(abs(pinfo.geomBounds.lower[i]),abs(pinfo.geomBounds.upper[i]));
      bool flat = pinfo.geomBounds.upper[i] - pinfo.geomBounds.lower[i] <= 128.0f*u;
      if (unlikely(flat)) continue;
      float sah = halfArea(lgeomBounds[i])*blocks(lcounts[i]) + halfArea(rgeomBounds[i])*blocks(rcounts[i]);
      if (unlikely(sah >= split.spatialSAH)) continue;
      split.spatialSAH = sah;
      split.sdim = i;
      split.numFailed = (int)numFailed;
    }
  }

  template<int logBlockSize>
  bool HeuristicSpatial<logBlockSize>::update_spatial_split(Split& split) 
  {
    size_t duplicatedTriangles = lcounts[split.sdim]+rcounts[split.sdim]-pinfo.size();
    atomic_t remainingDuplications = pinfo.duplications->sub(duplicatedTriangles);
    if (remainingDuplications < 0) {
      pinfo.duplications->add(duplicatedTriangles);
      return false;
    }

    new (&split.linfo) PrimInfo(lcounts[split.sdim],split.numFailed,lgeomBounds[split.sdim],lcentBounds[split.sdim],pinfo.duplications);
    new (&split.rinfo) PrimInfo(rcounts[split.sdim],split.numFailed,rgeomBounds[split.sdim],rcentBounds[split.sdim],pinfo.duplications);
    return true;
  }

  template<int logBlockSize>
  void HeuristicSpatial<logBlockSize>::best_object_split(Split& split)
  {
    /* Compute SAH for best object split. */
    Vector3f rAreas [maxBins];      //!< area of bounds of primitives on the right
    Vector3f rCounts[maxBins];      //!< blocks of primitives on the right
    
    /* sweep from right to left and compute parallel prefix of merged bounds */
    assert(mapping.size() > 0);
    Vector3i count = 0; BBox3f bx = empty; BBox3f by = empty; BBox3f bz = empty;
    for (size_t i=mapping.size()-1; i>0; i--)
    {
      count += counts[i];
      rCounts[i] = Vector3f(blocks(count));
      bx = merge(bx,geomBounds[i][0]); rAreas[i][0] = halfArea(bx);
      by = merge(by,geomBounds[i][1]); rAreas[i][1] = halfArea(by);
      bz = merge(bz,geomBounds[i][2]); rAreas[i][2] = halfArea(bz);
    }
    
    /* sweep from left to right and compute SAH */
    Vector3i ii = 1; Vector3f bestSAH = pos_inf; Vector3i bestPos = 0;
    count = 0; bx = empty; by = empty; bz = empty;
    for (size_t i=1; i<mapping.size(); i++, ii+=1)
    {
      count += counts[i-1];
      bx = merge(bx,geomBounds[i-1][0]); float Ax = halfArea(bx);
      by = merge(by,geomBounds[i-1][1]); float Ay = halfArea(by);
      bz = merge(bz,geomBounds[i-1][2]); float Az = halfArea(bz);
      const Vector3f lCount = Vector3f(blocks(count));
      const Vector3f lArea = Vector3f(Ax,Ay,Az);
      const Vector3f sah = lArea*lCount + rAreas[i]*rCounts[i];
      bestPos = select(lt_mask(sah,bestSAH),ii ,bestPos);
      bestSAH = select(lt_mask(sah,bestSAH),sah,bestSAH);
    }
    
    /* find best dimension */
    for (int i=0; i<3; i++) 
    {
      if (unlikely(pinfo.centBounds.lower[i] >= pinfo.centBounds.upper[i]))
        continue;
      
      if (bestSAH[i] < split.objectSAH && bestPos[i] != 0) {
        split.dim = i;
        split.pos = bestPos[i];
        split.objectSAH = bestSAH[i];
      }
    }
    split.mapping = mapping;
  }

  template<int logBlockSize>
  void HeuristicSpatial<logBlockSize>::update_object_split(Split& split)
  {
    /* early exit for invalid splits */
    const size_t pos = split.pos;
    const size_t dim = split.dim;
    if (unlikely(pos == 0)) return;
    
    /* calculate geometry info from binning data */
    size_t numLeft = 0, numRight = 0;
    BBox3f lcentBounds = empty, rcentBounds = empty;
    BBox3f lgeomBounds = empty, rgeomBounds = empty;
    for (size_t i=0; i<pos; i++) {
      numLeft += counts[i][dim];
      lcentBounds.grow(centBounds[i][dim]);
      lgeomBounds.grow(geomBounds[i][dim]);
    }
    for (size_t i=pos; i<mapping.size(); i++) {
      numRight += counts[i][dim];
      rcentBounds.grow(centBounds[i][dim]);
      rgeomBounds.grow(geomBounds[i][dim]);
    }
    assert(numLeft + numRight == pinfo.size());
    new (&split.linfo) PrimInfo(numLeft ,pinfo.numFailed,lgeomBounds,lcentBounds,pinfo.duplications);
    new (&split.rinfo) PrimInfo(numRight,pinfo.numFailed,rgeomBounds,rcentBounds,pinfo.duplications);
  }
  
  template<int logBlockSize>
  void HeuristicSpatial<logBlockSize>::best(Split& split)
  {
    /* find best spatial and object split */
    split.mapping = mapping;
    split.geom = geom;
    best_spatial_split(split);
    best_object_split (split);

    /* update prim info of left and right geometry */
    if (split.spatialSAH < split.objectSAH) {
      if (update_spatial_split(split)) return;
      split.spatialSAH = inf;
    }

    /* we also perform object split of spatial split failed */
    update_object_split(split);
  }
  
  template<int logBlockSize>
  void HeuristicSpatial<logBlockSize>::reduce(const HeuristicSpatial binners[], size_t num, HeuristicSpatial& binner_o)
  {
    binner_o = binners[0];
    for (size_t tid=1; tid<num; tid++) 
    {
      const HeuristicSpatial& binner = binners[tid];

      for (size_t bin=0; bin<binner.mapping.size(); bin++) 
      {
        for (size_t dim=0; dim<3; dim++) {
          binner_o.counts    [bin][dim] += binner.counts[bin][dim];
          binner_o.geomBounds[bin][dim].grow(binner.geomBounds[bin][dim]);
          binner_o.centBounds[bin][dim].grow(binner.centBounds[bin][dim]);
        }
      }

      binner_o.lcounts += binner.lcounts;
      binner_o.rcounts += binner.rcounts;
      for (size_t dim=0; dim<3; dim++) 
      {
        binner_o.lgeomBounds[dim].grow(binner.lgeomBounds[dim]);
        binner_o.rgeomBounds[dim].grow(binner.rgeomBounds[dim]);
        binner_o.lcentBounds[dim].grow(binner.lcentBounds[dim]);
        binner_o.rcentBounds[dim].grow(binner.rcentBounds[dim]);
      }
    }
  }
  
  /*! explicit template instantiations */
  template class HeuristicSpatial<0>;
  template class HeuristicSpatial<2>;
  template class HeuristicSpatial<3>;
}
