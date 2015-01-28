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
  void HeuristicSpatial<logBlockSize>::bin(const PrimRef* prims, size_t num)
  {
    for (size_t i=0; i<num; i++)
    {
      /*! map even and odd primitive to bin */
      const BBox3f prim = prims[i].bounds(); const Vec3ia bin = mapping.bin(prim); const Vec3fa center = Vec3fa(center2(prim));
      
      /*! increase bounds for bins for even primitive */
      const int b0 = bin.x; counts[b0][0]++; geomBounds[b0][0].extend(prim); centBounds[b0][0].extend(center);
      const int b1 = bin.y; counts[b1][1]++; geomBounds[b1][1].extend(prim); centBounds[b1][1].extend(center);
      const int b2 = bin.z; counts[b2][2]++; geomBounds[b2][2].extend(prim); centBounds[b2][2].extend(center);
      
      /*! calculate for each dimension if primitive goes to the left and right for spatial splits */
      const Vec3ba left  = mapping.left (prim);
      const Vec3ba right = mapping.right(prim);

      /*! Test spatial split in center of each dimension. */
      for (int maxDim=0; maxDim<3; maxDim++)
      {
        /*! Choose better side for primitives that can be put left or right. */
        if (left[maxDim] && right[maxDim]) {
          if (prim.upper[maxDim]-mapping.bright[maxDim] > mapping.bleft[maxDim]-prim.lower[maxDim]) {
            rgeomBounds[maxDim].extend(prim); rcentBounds[maxDim].extend(center); rcounts[maxDim]++;
          } else {
            lgeomBounds[maxDim].extend(prim); lcentBounds[maxDim].extend(center); lcounts[maxDim]++;
          }
        }
        /*! These definitely go to the left. */
        else if (left[maxDim]) {
          lgeomBounds[maxDim].extend(prim); lcentBounds[maxDim].extend(center); lcounts[maxDim]++;
        }
        /*! These definitely go to the right. */
        else if (right[maxDim]) {
          rgeomBounds[maxDim].extend(prim); rcentBounds[maxDim].extend(center); rcounts[maxDim]++;
        }
        /*! These have to get split and put to left and right. */
        else {
          PrimRef left,right; geom->split(prims[i],maxDim,mapping.center[maxDim],left,right);
          lgeomBounds[maxDim].extend(left. bounds()); lcentBounds[maxDim].extend(center2(left. bounds())); lcounts[maxDim]++;
          rgeomBounds[maxDim].extend(right.bounds()); rcentBounds[maxDim].extend(center2(right.bounds())); rcounts[maxDim]++;
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
    Vec3fa rAreas [maxBins];      //!< area of bounds of primitives on the right
    Vec3fa rCounts[maxBins];      //!< blocks of primitives on the right
    
    /* sweep from right to left and compute parallel prefix of merged bounds */
    assert(mapping.size() > 0);
    Vec3ia count = 0; BBox3f bx = empty; BBox3f by = empty; BBox3f bz = empty;
    for (size_t i=mapping.size()-1; i>0; i--)
    {
      count += counts[i];
      rCounts[i] = Vec3fa(blocks(count));
      bx = merge(bx,geomBounds[i][0]); rAreas[i][0] = halfArea(bx);
      by = merge(by,geomBounds[i][1]); rAreas[i][1] = halfArea(by);
      bz = merge(bz,geomBounds[i][2]); rAreas[i][2] = halfArea(bz);
    }
    
    /* sweep from left to right and compute SAH */
    Vec3ia ii = 1; Vec3fa bestSAH = pos_inf; Vec3ia bestPos = 0;
    count = 0; bx = empty; by = empty; bz = empty;
    for (size_t i=1; i<mapping.size(); i++, ii+=1)
    {
      count += counts[i-1];
      bx = merge(bx,geomBounds[i-1][0]); float Ax = halfArea(bx);
      by = merge(by,geomBounds[i-1][1]); float Ay = halfArea(by);
      bz = merge(bz,geomBounds[i-1][2]); float Az = halfArea(bz);
      const Vec3fa lCount = Vec3fa(blocks(count));
      const Vec3fa lArea = Vec3fa(Ax,Ay,Az);
      const Vec3fa sah = lArea*lCount + rAreas[i]*rCounts[i];
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
      lcentBounds.extend(centBounds[i][dim]);
      lgeomBounds.extend(geomBounds[i][dim]);
    }
    for (size_t i=pos; i<mapping.size(); i++) {
      numRight += counts[i][dim];
      rcentBounds.extend(centBounds[i][dim]);
      rgeomBounds.extend(geomBounds[i][dim]);
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
          binner_o.geomBounds[bin][dim].extend(binner.geomBounds[bin][dim]);
          binner_o.centBounds[bin][dim].extend(binner.centBounds[bin][dim]);
        }
      }

      binner_o.lcounts += binner.lcounts;
      binner_o.rcounts += binner.rcounts;
      for (size_t dim=0; dim<3; dim++) 
      {
        binner_o.lgeomBounds[dim].extend(binner.lgeomBounds[dim]);
        binner_o.rgeomBounds[dim].extend(binner.rgeomBounds[dim]);
        binner_o.lcentBounds[dim].extend(binner.lcentBounds[dim]);
        binner_o.rcentBounds[dim].extend(binner.rcentBounds[dim]);
      }
    }
  }
  
  /*! explicit template instantiations */
  template class HeuristicSpatial<0>;
  template class HeuristicSpatial<2>;
  template class HeuristicSpatial<3>;
}
