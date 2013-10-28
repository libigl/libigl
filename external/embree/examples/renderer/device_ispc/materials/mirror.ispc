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

#include "material.isph"
#include "brdfs/reflection.isph"

struct Mirror
{
  Material base;
  Reflection reflection;
};

///////////////////////////////////////////////////////////////////////////////
// Implementation

void Mirror__shade(const uniform Material *uniform _this,
                  const Ray&                  ray,            
                  const Medium&               currentMedium,  
                  const DifferentialGeometry& dg,             
                  uniform CompositedBRDF&     brdfs)          
{
  uniform const Mirror *uniform this = (uniform const Mirror *uniform)_this;
  CompositedBRDF__add(&brdfs,&this->reflection.base);
}

inline void Mirror__Constructor(uniform Mirror *uniform this,
                               const uniform vec3f& reflectance)
{
  Material__Constructor(&this->base,Material__Destructor,Mirror__shade,NULL);
  Reflection__Constructor(&this->reflection, reflectance);
}

///////////////////////////////////////////////////////////////////////////////
// External API

export void* uniform Mirror__new(const uniform vec3f& reflectance)
{
  uniform Mirror *uniform m = uniform new uniform Mirror;
  Mirror__Constructor(m,reflectance);
  return m;
}
