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

#include "framebuffer.isph"

void FrameBuffer__Destructor(uniform FrameBuffer* uniform this) {
  LOG(print("FrameBuffer__Destructor\n"));
  RefCount__Destructor(&this->base);
}

void FrameBuffer__Constructor(uniform FrameBuffer* uniform this,
                              const uniform uint width, 
                              const uniform uint height,
                              const uniform DestructorFunc destructor,
                              const uniform FrameBuffer__set set,
                              const uniform FrameBuffer__map map)
{
  LOG(print("FrameBuffer__Constructor\n"));
  RefCount__Constructor(&this->base,destructor);
  this->size.x = width;
  this->size.y = height;
  this->invSize.x = rcp(width);
  this->invSize.y = rcp(height);
  this->set = set;
  this->map = map;
}
