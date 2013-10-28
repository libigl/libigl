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

#include "ref.isph"

void RefCount__Destructor(uniform RefCount* uniform this) {
}

void RefCount__Constructor(uniform RefCount* uniform this,
                           uniform DestructorFunc destructor)
{
  this->counter = 0;
  this->destructor = destructor;
}

void RefCount__IncRef(uniform RefCount* uniform this) {
  if (this) this->counter++;
}

uniform bool RefCount__DecRef(uniform RefCount* uniform this) 
{
  if (this && --this->counter == 0) {
    if (this->destructor)
      this->destructor(this);
    delete this;
    return true;
  }
  return false;
}

void RefCount__Destroy(uniform RefCount* uniform this)
{
  if (this->counter == 0) {
    if (this->destructor)
      this->destructor(this);
    delete this;
  }
}

export void RefCount__incref(void* uniform _this) {
  RefCount__IncRef((uniform RefCount* uniform) _this);
}

export uniform bool RefCount__decref(void* uniform _this) {
  return RefCount__DecRef((uniform RefCount* uniform) _this);
}
