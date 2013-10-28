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

#ifndef __EMBREE_DISTRIBUTION_2D_H__
#define __EMBREE_DISTRIBUTION_2D_H__

#include "distribution1d.h"

namespace embree
{
  /*! 2D probability distribution. The probability distribution
   *  function (PDF) can be initialized with arbitrary data and be
   *  sampled. */
  class Distribution2D : public RefCount
  {
  public:

    /*! Default construction. */
    Distribution2D();

    /*! Construction from 2D distribution array f. */
    Distribution2D(const float** f, const size_t width, const size_t height);

    /*! Destruction. */
    ~Distribution2D();

    /*! Initialized the PDF and CDF arrays. */
    void init(const float** f, const size_t width, const size_t height);

    /*! read back the y-cdf and y-pdf data (for SPMD device
        offload). cdf_y is of size 'height+1', pdf_y of size 'height'  */
    void readback(float *cdf_y, float *pdf_y)
    { yDist.readback(cdf_y,pdf_y); }
    
    /*! read back the x-cdf and x-pdf data for given y (for SPMD
        device offload). cdf_x is of size 'width+1', pdf_x of size
        'width' */
    void readback(int y, float *cdf_x, float *pdf_x) 
    { xDists[y].readback(cdf_x,pdf_x); }

  public:

    /*! Draws a sample from the distribution. \param u is a pair of
     *  random numbers to use for sampling */
    Sample2f sample(const Vec2f& u) const;

    /*! Returns the probability density a sample would be drawn from
     *  location p. */
    float pdf(const Vec2f& p) const;

  private:
    size_t width;            //!< Number of elements in x direction
    size_t height;           //!< Number of elements in y direction
    Distribution1D yDist;    //!< Distribution to select between rows
    Distribution1D* xDists;  //!< One 1D Distribution per row
  };
}

#endif
