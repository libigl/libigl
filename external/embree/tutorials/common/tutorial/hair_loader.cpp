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

#include "hair_loader.h"

#define CONVERT_TO_BINARY 0

namespace embree
{
  float g_reduce_hair_segment_error = 0.0f;
  const int hair_bin_magick = 0x12EF3F90;

  int loadHairASCII(const FileName& fileName, OBJScene::HairSet* hairset, Vec3fa& offset)
  {  
    /* open hair file */
    FILE* f = fopen(fileName.c_str(),"r");
    if (!f) THROW_RUNTIME_ERROR("could not open " + fileName.str());

    char line[10000];
    fgets(line,10000,f);
    int numCurves = 0;
    
    while (fgets(line,10000,f) && !feof(f))
    {
      /* comment */
      if (line[0] == '#')
	continue;
      
      if (!strncmp(line,"Curve:",strlen("Curve:")))
      {
        char name[1000];
        unsigned int tracks, points;
        sscanf(line,"Curve: %s %d Tracks %d Points",name,&tracks,&points);

        /* skip Tracks line */
        fgets(line,10000,f);
        
        const int vertex_start_id = hairset->v.size();
        
        unsigned int id = 0;
        for (int i=0; i<points; i++)
        {
          fgets(line,10000,f);

          /* comment */
          if (line[0] == '#' || !strncmp(line," Tracks:",strlen(" Tracks:")))
            continue;

          Vec3fa v;
          if (i == 0) sscanf(line,"%d : Bezier %f %f %f %f",&id,&v.x,&v.y,&v.z,&v.w);
          else        sscanf(line,"%d : %f %f %f %f",&id,&v.x,&v.y,&v.z,&v.w);
          //printf("%d %d : %f %f %f %f \n",id,vertex_start_id+id,v.x,v.y,v.z,v.w);		
          v.x-=offset.x;
          v.y-=offset.y;
          v.z-=offset.z;
          hairset->v.push_back(v);
        }
        
        /* add indices to hair starts */
        for (int i=0; i<points-1; i+=3)
          hairset->hairs.push_back(OBJScene::Hair(vertex_start_id + i,numCurves));
	
        if (id != points-1) 
          THROW_RUNTIME_ERROR("hair parsing error");

        numCurves++;
      }
    }
    fclose(f);
    return numCurves;
  }

  int loadHairBin(const FileName& fileName, OBJScene::HairSet* hairset, Vec3fa& offset)
  {  
    FILE* fin = fopen(fileName.c_str(),"rb");
    if (!fin) THROW_RUNTIME_ERROR("could not open " + fileName.str());
    int magick; fread(&magick,sizeof(int),1,fin);
    if (magick != hair_bin_magick)
      THROW_RUNTIME_ERROR("invalid binary hair file " + fileName.str());
    int numHairs; fread(&numHairs,sizeof(int),1,fin);
    int numPoints; fread(&numPoints,sizeof(int),1,fin);
    int numSegments; fread(&numSegments,sizeof(int),1,fin);
    hairset->v.resize(numPoints);
    hairset->hairs.resize(numSegments);
    if (numPoints) fread(&hairset->v[0],sizeof(Vec3fa),numPoints,fin);
    if (numSegments) fread(&hairset->hairs[0],sizeof(OBJScene::Hair),numSegments,fin);
    fclose(fin);

    for (size_t i=0; i<numPoints; i++) {
      hairset->v[i].x-=offset.x;
      hairset->v[i].y-=offset.y;
      hairset->v[i].z-=offset.z;
    }
    return numHairs;
  }

  OBJScene::HairSet* reduce_hairs(OBJScene::HairSet* in, float relative_error)
  {
    OBJScene::HairSet* out = new OBJScene::HairSet;
    int numSegments = in->hairs.size();
    
    for (size_t i=0; i<numSegments; i++)
    {
      const OBJScene::Hair h0 = in->hairs[i+0];
      const Vec3fa a0 = in->v[h0.vertex+0];
      const Vec3fa a1 = in->v[h0.vertex+1];
      const Vec3fa a2 = in->v[h0.vertex+2];
      const Vec3fa a3 = in->v[h0.vertex+3];
      
      if (i+1 >= numSegments) 
      {
        if (out->v.size() == 0 || out->v.back() != a0) 
        out->v.push_back(a0);
        out->v.push_back(a1);
        out->v.push_back(a2);
        out->v.push_back(a3);
        out->hairs.push_back(OBJScene::Hair(out->v.size()-4,h0.id));
        continue;
      }
      
      const OBJScene::Hair h1 = in->hairs[i+1];
      const Vec3fa b0 = in->v[h1.vertex+0];
      const Vec3fa b1 = in->v[h1.vertex+1];
      const Vec3fa b2 = in->v[h1.vertex+2];
      const Vec3fa b3 = in->v[h1.vertex+3];
      
      if (h0.vertex+3 != h1.vertex || h0.id != h1.id) 
      {
        if (out->v.size() == 0 || out->v.back() != a0) 
        out->v.push_back(a0);
        out->v.push_back(a1);
        out->v.push_back(a2);
        out->v.push_back(a3);
        out->hairs.push_back(OBJScene::Hair(out->v.size()-4,h0.id));
        continue;
      }
      
      const Vec3fa c0 = a0;
      const Vec3fa c1 = 2*a1-a0;
      const Vec3fa c2 = 2*b2-b3;
      const Vec3fa c3 = b3;

      const Vec3fa p00 = c0;
      const Vec3fa p01 = c1;
      const Vec3fa p02 = c2;
      const Vec3fa p03 = c3;
      const Vec3fa p10 = 0.5f*(p00 + p01);
      const Vec3fa p11 = 0.5f*(p01 + p02);
      const Vec3fa p12 = 0.5f*(p02 + p03);
      const Vec3fa p20 = 0.5f*(p10 + p11);
      const Vec3fa p21 = 0.5f*(p11 + p12);
      const Vec3fa p30 = 0.5f*(p20 + p21);

      const float len = length(a3-a0) + length(b3-b0);
      const float err = length(p30-a3);

      if (err > relative_error*len)
      {
        if (out->v.size() == 0 || out->v.back() != a0) 
        out->v.push_back(a0);
        out->v.push_back(a1);
        out->v.push_back(a2);
        out->v.push_back(a3);
        out->hairs.push_back(OBJScene::Hair(out->v.size()-4,h0.id));
      }
      else
      {
        if (out->v.size() == 0 || out->v.back() != c0) 
        out->v.push_back(c0);
        out->v.push_back(c1);
        out->v.push_back(c2);
        out->v.push_back(c3);
        out->hairs.push_back(OBJScene::Hair(out->v.size()-4,h0.id));
        i++;
      }
    }

    delete in; 
    return out;
  }

  void loadHair(const FileName& fileName, OBJScene& scene, Vec3fa& offset)
  {
    /* add new hair set to scene */
    OBJScene::HairSet* hairset = new OBJScene::HairSet; 
#if CONVERT_TO_BINARY   
    offset = Vec3fa(zero);
#endif

    int numHairs = 0;
    if (fileName.ext() == "txt")
      numHairs = loadHairASCII(fileName,hairset,offset);
    else
      numHairs = loadHairBin(fileName,hairset,offset);
    
    /* reduce number of hairs */
    if (g_reduce_hair_segment_error != 0.0f)
    {
      std::ios_base :: fmtflags   flag = std::cout.flags();
      std           :: streamsize prec = std::cout.precision();
      std::cout << "reducing number of hair segments ... " << std::flush;
      std::cout.precision(3);
      std::cout << 1E-6*float(hairset->hairs.size()) << "M" << std::flush;
      for (size_t i=0; i<10; i++) {
        hairset = reduce_hairs(hairset,g_reduce_hair_segment_error);
        std::cout << " " << 1E-6*float(hairset->hairs.size()) << "M" << std::flush;
      }
      std::cout << " [DONE]" << std::endl;
      std::cout.flags    (flag);
      std::cout.precision(prec);
    }

    /* add hairset to scene */
    scene.hairsets.push_back(hairset);

    int numPoints = hairset->v.size();
    int numSegments = hairset->hairs.size();

#if CONVERT_TO_BINARY
    FILE* fout = fopen(fileName.setExt(".bin").c_str(),"wb");
    if (!fout) THROW_RUNTIME_ERROR("could not open " + fileName.str());
    fwrite(&hair_bin_magick,sizeof(int),1,fout);
    fwrite(&numHairs,sizeof(int),1,fout);
    fwrite(&numPoints,sizeof(int),1,fout);
    fwrite(&numSegments,sizeof(int),1,fout);
    if (numPoints) fwrite(&hairset->v[0],sizeof(Vec3fa),numPoints,fout);
    if (numSegments) fwrite(&hairset->hairs[0],sizeof(OBJScene::Hair),numSegments,fout);
    fclose(fout);
#endif
  }
}

  
