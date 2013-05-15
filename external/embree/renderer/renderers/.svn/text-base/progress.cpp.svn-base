// ======================================================================== //
// Copyright 2009-2012 Intel Corporation                                    //
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

#include "progress.h"
#include "sys/sysinfo.h"

namespace embree
{
  Progress::Progress (size_t num) 
    : curElement(0), numElements(num), numDrawn(0), terminalWidth(getTerminalWidth()) {}

  void Progress::start() {
    drawEmptyBar();
  }

  void Progress::drawEmptyBar() 
  {
    std::cout << "\r[" << std::flush;
    size_t width = max(ssize_t(2),terminalWidth-2);
    for (size_t i=0; i<width; i++) 
      std::cout << " ";
    std::cout << "]\r";
    std::cout << "[" << std::flush;
    numDrawn = 0;
  }

  void Progress::next() 
  {
    Lock<MutexSys> lock(mutex);
    size_t cur = curElement++;
    size_t width = max(ssize_t(2),terminalWidth-2);
    size_t numToDraw = cur*width/max((size_t)1,numElements);
    for (size_t i=numDrawn; i<numToDraw; i++) {
      std::cout << "+" << std::flush;
    }
    numDrawn = numToDraw;
  }

  void Progress::end() 
  {
    size_t width = max(ssize_t(2),terminalWidth-2);
    for (size_t i=numDrawn; i<width; i++) {
      std::cout << "+" << std::flush;
    }
    std::cout << "]" << std::endl;
  }
}
