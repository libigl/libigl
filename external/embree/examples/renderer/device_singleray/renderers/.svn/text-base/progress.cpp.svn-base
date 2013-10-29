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
    size_t width = max(ssize_t(2),ssize_t(terminalWidth-2));
    for (size_t i=0; i<width; i++) 
      std::cout << " ";
    std::cout << "]\r";
    std::cout << "[" << std::flush;
    numDrawn = 0;
  }

  void Progress::next() 
  {
    Lock<MutexSys> lock(mutex);
    ssize_t curTerminalWidth = getTerminalWidth();
    if (terminalWidth != size_t(curTerminalWidth)) {
      drawEmptyBar();
      terminalWidth = curTerminalWidth;
    }

    size_t cur = curElement++;
    size_t width = max(ssize_t(2),ssize_t(curTerminalWidth-2));
    size_t numToDraw = cur*width/max(ssize_t(1),ssize_t(numElements-1));
    for (size_t i=numDrawn; i<numToDraw; i++) {
      std::cout << "+" << std::flush;
    }
    numDrawn = numToDraw;
  }

  void Progress::end() {
    std::cout << "]" << std::endl;
  }
}
