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

void Material__selectNextMedium(uniform const Material *uniform _THIS,
                                varying Medium& currentMedium) 
{ 
}

void Material__Destructor(uniform RefCount *uniform this)
{
  LOG(print("Material__Destructor\n"));
  RefCount__Destructor(this);
}

void Material__Constructor(uniform Material *uniform this,
                           uniform DestructorFunc destructor,
                           uniform Material__shade shade,
                           uniform Material__selectNextMediumFunc selectNextMedium)
{
  LOG(print("Material__Constructor\n"));
  RefCount__Constructor(&this->base,destructor);
  this->shade = shade;
  if (selectNextMedium) this->selectNextMedium = selectNextMedium; 
  else                  this->selectNextMedium = Material__selectNextMedium;
}
