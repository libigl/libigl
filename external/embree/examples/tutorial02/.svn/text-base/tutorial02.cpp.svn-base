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

#include "sys/platform.h"
#include "sys/ref.h"
#include "lexers/streamfilters.h"
#include "lexers/parsestream.h"
#include "../tutorials/glutdisplay.h"
#include "../tutorial_host/tutorials_host.h"

namespace embree
{
  /* name of the tutorial */
  const char* tutorialName = "tutorial02";

  /* flags */
  static int g_verbose = 0;

  /* output settings */
  static size_t g_width = 512;
  static size_t g_height = 512;
  static bool g_fullscreen = false;

  static void parseCommandLine(Ref<ParseStream> cin, const FileName& path)
  {
    while (true)
    {
      std::string tag = cin->getString();
      if (tag == "") return;

      /* parse command line parameters from a file */
      if (tag == "-c") {
        FileName file = path + cin->getFileName();
        parseCommandLine(new ParseStream(new LineCommentFilter(file, "#")), file.path());
      }

      /* parse camera parameters */
      else if (tag == "-vp") g_camera.from = cin->getVector3f();
      else if (tag == "-vi") g_camera.to = cin->getVector3f();
      else if (tag == "-vd") g_camera.to = g_camera.from + cin->getVector3f();
      else if (tag == "-vu") g_camera.up = cin->getVector3f();
      else if (tag == "-fov") g_camera.fov = cin->getFloat();

      /* frame buffer size */
      else if (tag == "-size") {
        g_width = cin->getInt();
        g_height = cin->getInt();
      }

      /* full screen mode */
      else if (tag == "-fullscreen") 
        g_fullscreen = true;
      
      /*! enable verbose output mode */
      else if (tag == "-verbose") {
        g_verbose = 1;
      }

      /* skip unknown command line parameter */
      else {
        std::cerr << "unknown command line parameter: " << tag << " ";
        while (cin->peek() != "" && cin->peek()[0] != '-') std::cerr << cin->getString() << " ";
        std::cerr << std::endl;
      }
    }
  }

  /* main function in embree namespace */
  int main(int argc, char** argv) 
  {
    /* create stream for parsing */
    Ref<ParseStream> stream = new ParseStream(new CommandLineStream(argc, argv));

    /* parse command line */  
    parseCommandLine(stream, FileName());

    /* initialize ISPC */
    init(g_verbose);

    /* initialize GLUT */
    initGlut("tutorial02",g_width,g_height,g_fullscreen,true);
    
    return 0;
  }
}

int main(int argc, char** argv) {
  return embree::main(argc, argv);
}
