// ======================================================================== //
// Copyright 2009-2014 Intel Corporation                                    //
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

#include "tasklogger.h"
#include "sys/sysinfo.h"
#include "math/math.h"
#include "math/vec2.h"
#include "math/bbox.h"
#include <fstream>
#include "string.h"

namespace embree
{
  bool TaskLogger::active = false;
  int64 TaskLogger::startCycle = 0;
  std::vector<TaskLogger*> TaskLogger::threads;

  bool TaskLogger::init (size_t numThreads)
  {
#if defined(RTCORE_TASKLOGGER)
    if (threads.size() != numThreads) {
      threads.clear();
      for (size_t i=0; i<numThreads; i++)
        threads.push_back(new TaskLogger((int)i));
    }
#endif
    return true;
  }
  
  void TaskLogger::start() 
  {
#if defined(RTCORE_TASKLOGGER)
    for (size_t i=0; i<threads.size(); i++)
      threads[i]->reset();
    
    startCycle = rdtsc();
    active = true;
#endif
  }

  void TaskLogger::stop() {
#if defined(RTCORE_TASKLOGGER)
    active = false;
#endif
  }

  namespace DRAW
  {
    int Black = 0;
    int Blue = 1;
    int Green = 2;
    int Cyan = 3;
    int Red = 4;
    int Magenta = 5;
    int Yellow = 6;
    int White = 7;
    int DarkBlue = 8;
    
    class Drawing 
    {
    public:

      Drawing(const char* file) 
      {
        id = 0;
        fout.open (file);
        fout << "#FIG 3.2" << "\n";
        fout << "# Created by Embree" << "\n";
        fout << "# METADATA <version>1.0</version>" << "\n";
        fout << "Portrait" << "\n";
        fout << "Flush left" << "\n";
        fout << "Metric" << "\n";
        fout << "A4" << "\n";
        fout << "100.00" << "\n";
        fout << "Single" << "\n";
        fout << "-2" << "\n";
        fout << "1200 2" << "\n";

        ofsX = 0;
        ofsY = 0;
        scaleX = 500;
        scaleY = 500;
      }

      void drawNamedBox(const char* name, Vec2f a, Vec2f b, int fillColor, int thickness = 0, int borderColor = Black) 
      {
        fout << "# METADATA <id>" << id++ << "</id>" << "\n";

        if (name) 
          fout << "# " << name << "\n";

        const int fillBrightness = 20;
        fout << "2 2 0 " << thickness << " " << borderColor << " " << fillColor << " 50 0 " << fillBrightness << " 4.000 2 0 7 0 0 5" << "\n";

        const int x0 = int(scaleX * a.x + ofsX);
        const int y0 = int(scaleY * a.y + ofsY);
        const int x1 = int(scaleX * b.x + ofsX);
        const int y1 = int(scaleY * b.y + ofsY);

        fout << "     ";
        fout << x0 << " " << y0 << " ";
        fout << x1 << " " << y0 << " ";
        fout << x1 << " " << y1 << " ";
        fout << x0 << " " << y1 << " ";
        fout << x0 << " " << y0 << "\n";
      }

      void drawBox(Vec2f a, Vec2f b, int fillColor, int thickness = 0, int borderColor = Black) {
        drawNamedBox(NULL, a, b, fillColor, thickness, borderColor) ;
      }

      void drawText(Vec2f p, const char* text, int size, int color) 
      {
         fout << "# METADATA <id>" << id++ << "</id>" << "\n";
         const int x0 = int(scaleX * p.x + ofsX);
         const int y0 = int(scaleY * p.y + ofsY);
         fout << "4 0 " << color << " 50 0 0 " << size << " 0.0000 3 873 150 " << x0 << " " << y0 << " " << text << "\\001" << "\n";
      }

      void drawPolyLine(Vec2f* points, int numPoints, int size, int color) 
      {
          fout << "# METADATA <id>" << id++ << "</id>" << "\n";
          fout << "2 1 0 2 0 7 50 0 -1 4.000 2 0 0 0 0 " << numPoints << "\n";
          for (int i=0; i<numPoints; i++) {
            const int x0 = int(scaleX * points[i].x + ofsX);
            const int y0 = int(scaleY * points[i].y + ofsY);
            fout << x0 << " " << y0 << " ";
          }
          fout << "\n";
      }

      void drawArrow(Vec2f a, Vec2f b, int color) 
      {
        fout << "# METADATA <id>" << id++ << "</id>" << "\n";
        fout << "2 1 0 2 0 7 50 0 -1 4.000 2 0 0 1 0 2" << "\n";
        fout << "1 1 2.00 120.00 240.00" << "\n";
        const int x0 = int(scaleX * a.x + ofsX);
        const int y0 = int(scaleY * a.y + ofsY);
        const int x1 = int(scaleX * b.x + ofsX);
        const int y1 = int(scaleY * b.y + ofsY);
        fout << x0 << " " << y0 << " ";
        fout << x1 << " " << y1 << "\n";
      }

      void drawLine(Vec2f a, Vec2f b, int color) 
      {
        fout << "# METADATA <id>" << id++ << "</id>" << "\n";
        fout << "2 1 0 2 0 7 50 0 -1 4.000 2 0 0 0 0 2" << "\n";
        const int x0 = int(scaleX * a.x + ofsX);
        const int y0 = int(scaleY * a.y + ofsY);
        const int x1 = int(scaleX * b.x + ofsX);
        const int y1 = int(scaleY * b.y + ofsY);
        fout << x0 << " " << y0 << " ";
        fout << x1 << " " << y1 << "\n";
      }

      ~Drawing() {
        fout.close();
      }

    public:
       int id;
       float ofsX, ofsY;
       float scaleX, scaleY;
       std::ofstream fout;
    };
  }

  /** store all logged data into FIG file */
  void TaskLogger::store(const char* fname)
  {
#if defined(RTCORE_TASKLOGGER)

    /** generate xfig drawing */
    const int64 xAxisStepSize = 1000000000;
    const char* xAxisUnit = "B";
		const int textSize = 8;
		const int lineSize = 12;

    BBox2f box;
    DRAW::Drawing sheet(fname); // DIN A4: width 21cm, height 29.7cm

    /** find start and end cycle */
    int64 minCycle = ((int64) 1) << 62;
    int64 maxCycle = 0;
    for (size_t i=0; i<threads.size(); i++) {
      for (size_t j=0; j<threads[i]->curTask; j++) {
        minCycle = min(minCycle,threads[i]->counters[j].start);
        maxCycle = max(maxCycle,threads[i]->counters[j].stop);
      }
    }
    minCycle = 0;
    
    //int64 dCycle = 900000000; 
    int64 dCycle = max(maxCycle - minCycle, 1LL); // to avoid div by 0
    
    /** compute usage and bandwidth statistics */
    const unsigned numUsage = 1024;
    int64 usedCycles[numUsage];
    for (size_t i=0; i<numUsage; i++) usedCycles[i] = 0;
    
    for (size_t i=0; i<threads.size(); i++) {
      for (size_t j=0; j<threads[i]->curTask; j++) {
        int64 startCycle = threads[i]->counters[j].start;
        int64 endCycle = threads[i]->counters[j].stop;
        for (unsigned b=0; b<numUsage; b++) {
          int64 clipMin = max(startCycle,minCycle+b*dCycle/numUsage);
          int64 clipMax = min(endCycle,minCycle+(b+1)*dCycle/numUsage);
          int64 overlap = clipMax-clipMin;
          if (overlap > 0) usedCycles[b] += overlap;
        }
      }
    }
    
    /** compute time scale */
    box.lower = Vec2f(3.0f,1.0f);
    box.upper = Vec2f(18.0f,17.0f);
    float scaleX = (box.upper.x-box.lower.x) / dCycle;
    float ofsX = box.lower.x - scaleX * minCycle;
    size_t numThreads = threads.size();

    /** draw task scheduling */
    float dy = min(1.0f,(box.upper.y-box.lower.y) / threads.size());
    for (size_t tid=0; tid<numThreads; tid++) 
    {
      const size_t i = tid;
      float top = box.lower.y + i * dy;
      float bottom = box.lower.y + (i+1) * dy - 0.2f;
      char str[256];
      sprintf(str,"thread%i",(int)tid);
      sheet.drawText(Vec2f(box.lower.x-2.0f,top+0.4f),str,textSize,DRAW::Black);
      
      TaskLogger* counters = threads[tid];
      for (size_t j=0; j<counters->curTask; j++) 
      {
        float t0 = scaleX * counters->counters[j].start + ofsX;
        float t1 = scaleX * counters->counters[j].stop + ofsX;
        
        int fillColor = DRAW::White;
        
        const char* name = counters->counters[j].name;
        if (name == NULL) name = "NULL";
        
        /*else if (!strcmp(name,"build::primrefgen")) fillColor = DRAW::Blue;
        else if (!strcmp(name,"build::full")) fillColor = DRAW::Blue;
        else if (!strcmp(name,"build::split")) fillColor = DRAW::Yellow;
        else if (!strcmp(name,"build::parsplit")) fillColor = DRAW::Yellow;*/

        if (!strcmp(name,"parallel_binning")) fillColor = DRAW::Blue;
        else if (!strcmp(name,"parallel_partition")) fillColor = DRAW::Blue;

        else if (!strcmp(name,"toplevel_open_parallel")) fillColor = DRAW::Blue;
        else if (!strcmp(name,"toplevel_build_subtrees")) fillColor = DRAW::Blue;
        else if (!strcmp(name,"toplevel_build_parallel")) fillColor = DRAW::Blue;

        else if (!strcmp(name,"build_parallel_morton")) fillColor = DRAW::Red;

        else if (!strcmp(name,"scene_build")) fillColor = DRAW::Black;

        else if (!strcmp(name,"BVH4Refit::sequential")) fillColor = DRAW::Blue;
        else if (!strcmp(name,"BVH4Refit::parallel")) fillColor = DRAW::DarkBlue;

        else if (!strcmp(name,"ISPCTask")) fillColor = DRAW::Green;

        else {
          fillColor = DRAW::Black;
        }
        
        if (counters->counters[j].elt == size_t(-1)) fillColor = DRAW::Black;
        
        char taskName[256];
        sprintf(taskName,"%s : %i",name,(int)counters->counters[j].elt);
        sheet.drawNamedBox(taskName,Vec2f(t0,top),Vec2f(t1,bottom),fillColor); //,1);
      }
    }
    sheet.drawLine(Vec2f(box.lower.x,box.lower.y),Vec2f(box.lower.x,box.upper.y),DRAW::Black);
    
    /** draw usage diagram */
    box.lower = Vec2f(3.0f,22.0f);
    box.upper = Vec2f(18.0f,18.0f);
    
    Vec2f points[numUsage];
    for (unsigned i=0; i<numUsage; i++) {
      float usage = float(numUsage)*usedCycles[i]/(threads.size()*dCycle);
      points[i].x = box.lower.x + i*(box.upper.x-box.lower.x)/float(numUsage);
      points[i].y = box.lower.y + usage*(box.upper.y-box.lower.y);
    }
    
    // x-axis
    sheet.drawArrow(Vec2f(box.lower.x,box.lower.y),Vec2f(box.upper.x,box.lower.y),DRAW::Black);
    int n=0;
    for (int64 i=minCycle; i<maxCycle; i+=xAxisStepSize, n++) 
    {
      char txt[256];
      sprintf(txt,"%i%s",n,xAxisUnit);
      float x = box.lower.x + (i-minCycle)*(box.upper.x-box.lower.x)/dCycle;
      sheet.drawLine(Vec2f(x,box.lower.y),Vec2f(x,box.lower.y+0.1f),DRAW::Black);
      sheet.drawText(Vec2f(x,box.lower.y+0.5f),txt,textSize,DRAW::Black);
    }
    
    // y-axis
    sheet.drawArrow(Vec2f(box.lower.x,box.lower.y),Vec2f(box.lower.x,box.upper.y),DRAW::Black);
    sheet.drawLine(Vec2f(box.lower.x,box.lower.y),Vec2f(box.lower.x-0.1f,box.lower.y),DRAW::Black);
    sheet.drawText(Vec2f(box.lower.x-0.7f,box.lower.y),"0%",textSize,DRAW::Black);
    sheet.drawText(Vec2f(box.lower.x-1.2f,box.upper.y),"100%",textSize,DRAW::Black);
    
    sheet.drawPolyLine(points,numUsage,lineSize,DRAW::Blue);
    sheet.drawText(Vec2f(0.5f+box.lower.x,0.5f+box.upper.y),"Usage [Percent]",textSize,DRAW::Black);
#endif
  }
}
