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

#include "distribution2d.isph"

/*! Construction of 1D distribution from distribution array f. */
void Distribution1D__Create(const uniform float* uniform f, 
                            const uniform uint size, 
                            uniform float* uniform pdf, 
                            uniform float* uniform cdf)
{
  /*! accumulate the function f */
  cdf[0] = 0.0f;
  for (uniform uint i=1; i<size+1; i++)
    cdf[i] = cdf[i-1] + f[i-1];
  
  /*! compute reciprocal sum */
  uniform float rcpSum = cdf[size] == 0.0f ? 0.0f : rcp(cdf[size]);

  /*! normalize the probability distribution and cumulative distribution */
  for (uniform uint i=1; i<=size; i++) {
    pdf[i-1] = f[i-1] * rcpSum * size;
    cdf[i] *= rcpSum;
  }
  cdf[size] = 1.0f;
}

//! sample a 1d cdf from given array and index range 
inline Sample1f Distribution1D__sample(const float u, const int _first, const int _len,
                                       const float* uniform cdf, const float* uniform pdf)
{
  // -------------------------------------------------------
  // upper_bound
  // -------------------------------------------------------
  const float val = u;
  int first = _first;
  int len = _len;

  while (len > 0) {
    const int __half = len >> 1;
    const int __middle = first + __half;
    if (val < cdf[__middle]) {
      len = __half;
    } else {
      first = __middle+1;
      len = len - __half - 1;
    }
  };
  
  // -------------------------------------------------------
  // sample
  // -------------------------------------------------------
  const int index = clamp(first-1,_first,_first+_len-1);
  const float fraction = (u - cdf[index]) * rcp(cdf[index+1] - cdf[index]);
  return make_Sample1f(index-_first+fraction,pdf[index]);
}

Sample2f Distribution2D__sample(const uniform Distribution2D* uniform this, const vec2f &u) 
{
  /*! use u.y to sample a row */
  const Sample1f sy = Distribution1D__sample(u.y,0,this->size.y,this->cdf_y,this->pdf_y);
  const int y = clamp((int)(sy.v),0,(int)(this->size.y)-1);
  
  /*! use u.x to sample inside the row */
  const int x0 = y*(this->size.x+1);
  const Sample1f sx = Distribution1D__sample(u.x,x0,this->size.x,this->cdf_x,this->pdf_x);
  return make_Sample2f(make_vec2f(sx.v,sy.v),sx.pdf*sy.pdf);
}

void Distribution2D__Destructor(uniform RefCount* uniform _this)
{ 
  uniform Distribution2D* uniform this = (uniform Distribution2D* uniform) _this;
  delete[] this->cdf_x;
  delete[] this->cdf_y;
  delete[] this->pdf_x;
  delete[] this->pdf_y;
  RefCount__Destructor(_this);
}

void Distribution2D__Constructor(uniform Distribution2D* uniform this,
                                 const uniform float* uniform f,
                                 const uniform vec2ui size) 
{
  RefCount__Constructor(&this->base,Distribution2D__Destructor);

  /*! create pdf and cdf for each row */
  uniform float* uniform pdf_x = uniform new uniform float[size.y*(size.x+1)]; // could be size.y*size.x, but this makes addressing later easier
  uniform float* uniform cdf_x = uniform new uniform float[size.y*(size.x+1)];

  /*! create pdf and cdf for sampling rows */
  uniform float* uniform fy    = uniform new uniform float[size.y];
  uniform float* uniform pdf_y = uniform new uniform float[size.y];
  uniform float* uniform cdf_y = uniform new uniform float[size.y+1];
  
  /*! compute y distribution and initialize row distributions */
  for (uniform uint y=0; y<size.y; y++)
  {
    /*! accumulate row to compute y distribution */
    fy[y] = 0.0f;
    for (uniform uint x=0; x<size.x; x++) {
      fy[y] += f[y*size.x+x];
    }
    
    /*! initialize distribution for current row */
    Distribution1D__Create(f+y*size.x, size.x, pdf_x+y*(size.x+1), cdf_x+y*(size.x+1));
  }
  
  /*! initializes the y distribution */
  Distribution1D__Create(fy,size.y,pdf_y,cdf_y);
  delete[] fy;

  /*! set data */
  this->size  = size;
  this->pdf_x = pdf_x;
  this->pdf_y = pdf_y;
  this->cdf_x = cdf_x;
  this->cdf_y = cdf_y;
}

uniform Distribution2D* uniform Distribution2D__new(const uniform float* uniform f,
                                                    const uniform vec2ui size)
{
  uniform Distribution2D* uniform this = uniform new uniform Distribution2D;
  Distribution2D__Constructor(this,f,size);
  return this;
}

