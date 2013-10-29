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

#include "accubuffer.isph"

void AccuBuffer__Destructor(uniform RefCount* uniform _this)
{
  uniform AccuBuffer* uniform this = (uniform AccuBuffer* uniform) _this;
  delete[] this->ptr; this->ptr = NULL;
  RefCount__Destructor(_this);
}

void AccuBuffer__Constructor(uniform AccuBuffer* uniform this, const uniform uint width, const uniform uint height)
{
  RefCount__Constructor(&this->base,AccuBuffer__Destructor);
  this->size.x = width;
  this->size.y = height;
  this->ptr = uniform new uniform vec4f[width*height];
  for (uniform int i=0; i<width*height; i++)
    this->ptr[i] = make_vec4f(0.0f,0.0f,0.0f,0.0f);
}

uniform AccuBuffer* uniform AccuBuffer__new(const uniform uint width, const uniform uint height)
{
  uniform AccuBuffer* uniform this = uniform new uniform AccuBuffer;
  AccuBuffer__Constructor(this,width,height);
  return this;
}

