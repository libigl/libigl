// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "readWRL.h"
#include <iostream>

template <typename Scalar, typename Index>
IGL_INLINE bool igl::readWRL(
  const std::string wrl_file_name,
  std::vector<std::vector<Scalar > > & V,
  std::vector<std::vector<Index > > & F)
{
  std::vector<std::vector<Scalar> > VC;
  return readWRL(wrl_file_name,V,F,VC);
}

template <typename Scalar, typename Index>
IGL_INLINE bool igl::readWRL(
  const std::string wrl_file_name,
  std::vector<std::vector<Scalar > > & V,
  std::vector<std::vector<Index > > & F,
  std::vector<std::vector<Scalar > > & VC)
{
  using namespace std;
  // for using fseek/ftell, open with binary mode.
  FILE * wrl_file = fopen(wrl_file_name.c_str(),"rb");
  if(NULL==wrl_file)
  {
    printf("IOError: %s could not be opened...",wrl_file_name.c_str());
    return false;
  }
  return readWRL(wrl_file,V,F,VC);
}

template <typename Scalar, typename Index>
IGL_INLINE bool igl::readWRL(
  FILE * wrl_file,
  std::vector<std::vector<Scalar > > & V,
  std::vector<std::vector<Index > > & F,
  std::vector<std::vector<Scalar > >& VC)
{
  using namespace std;

  char line[1000];
  // Read lines until seeing "point ["
  // treat other lines in file as "comments"
  bool still_comments = true;
  string needle("point [");
  string haystack;
  while(still_comments)
  {
    long thisLine = ftell(wrl_file);
    if(fgets(line,1000,wrl_file) == NULL)
    {
      std::cerr<<"readWRL, reached EOF without finding \"point [\""<<std::endl;
      fclose(wrl_file);
      return false;
    }
    haystack = string(line);
    still_comments = string::npos == haystack.find(needle);

    // in case first vertex position is written in the same line.
    // i.e. "point[ %lf %lf %lf,"
    if(!still_comments)
    {
      fseek(wrl_file, thisLine, SEEK_SET);
      while(fgetc(wrl_file) != '[')
      {
      }
    }
  }

  // read points in sets of 3
  int floats_read = 3;
  double x,y,z;
  while(floats_read == 3)
  {
    floats_read = fscanf(wrl_file," %lf %lf %lf,",&x,&y,&z);
    if(floats_read == 3)
    {
      vector<Scalar > point;
      point.resize(3);
      point[0] = x;
      point[1] = y;
      point[2] = z;
      V.push_back(point);
      //printf("(%g, %g, %g)\n",x,y,z);
    }else if(floats_read != 0)
    {
      printf("ERROR: unrecognized format...\n");
      return false;
    }
  }

  // there are no guarantee of the order of nodes.
  fseek(wrl_file, 0, SEEK_SET);
  bool colorPerVertex = true;
  still_comments = true;
  needle = string("colorPerVertex TRUE");
  while(still_comments)
  {
    long thisLine = ftell(wrl_file);
    if(fgets(line,1000,wrl_file) == NULL)
    {
      colorPerVertex = false;
      break;
    }
    haystack = string(line);
    still_comments = string::npos == haystack.find(needle);
  }
  if(colorPerVertex)
  {
    // we only support colorPerVertex
    // we do NOT support colorIndex

    VC.reserve(V.size());
    fseek(wrl_file, 0, SEEK_SET);
    // Read lines until seeing "color ["
    // treat other lines in file as "comments"
    still_comments = true;
    needle = string("color [");
    while(still_comments)
    {
      long thisLine = ftell(wrl_file);
      if(fgets(line,1000,wrl_file) == NULL)
      {
        std::cerr<<"readWRL, reached EOF without finding \"color [\""<<std::endl;
        fclose(wrl_file);
        return false;
      }
      haystack = string(line);
      still_comments = string::npos == haystack.find(needle);

      // in case first vertex position is written in the same line.
      // i.e. "color[ %lf %lf %lf,"
      if(!still_comments)
      {
        fseek(wrl_file, thisLine, SEEK_SET);
        while(fgetc(wrl_file) != '[')
        {
        }
      }
    }
    // read points in sets of 3
    int floats_read = 3;
    double r,g,b;
    while(floats_read == 3)
    {
      floats_read = fscanf(wrl_file," %lf %lf %lf,",&r,&g,&b);
      if(floats_read == 3)
      {
        vector<Scalar > color;
        color.resize(3);
        color[0] = r;
        color[1] = g;
        color[2] = b;
        VC.push_back(color);
        //printf("(%g, %g, %g)\n",x,y,z);
      }else if(floats_read != 0)
      {
        printf("ERROR: unrecognized format...\n");
        return false;
      }
    }

    if(V.size() != VC.size())
    {
      std::cerr<<"readWRL, we do NOT support colorIndex."<<std::endl;
      VC.clear();
    }
  }

  // there are no guarantee of the order of nodes.
  fseek(wrl_file, 0, SEEK_SET);

  // Read lines until seeing "coordIndex ["
  // treat other lines in file as "comments"
  still_comments = true;
  needle = string("coordIndex [");
  while(still_comments)
  {
    long thisLine = ftell(wrl_file);
    if(fgets(line,1000,wrl_file) == NULL)
    {
      std::cerr<<"readWRL, reached EOF without finding \"coordIndex [\""<<std::endl;
      fclose(wrl_file);
      return false;
    }
    haystack = string(line);
    still_comments = string::npos == haystack.find(needle);

    // in case first vertex position is written in the same line.
    // i.e. "coordIndex [ %d %d %d -1,"
    if(!still_comments)
    {
      fseek(wrl_file, thisLine, SEEK_SET);
      while(fgetc(wrl_file) != '[')
      {
      }
    }
  }
  // read F
  int ints_read = 1;
  while(ints_read > 0)
  {
    // read new face indices (until hit -1)
    vector<Index > face;
    while(true)
    {
      // indices are 0-indexed
      int i;
      ints_read = fscanf(wrl_file," %d,",&i);
      if(ints_read > 0)
      {
        if(i>=0)
        {
          face.push_back(i);
        }else
        {
          F.push_back(face);
          break;
        }
      }else
      {
        // for the last entry, there are not '-1'
        if(fgetc(wrl_file) == ']')
        {
          F.push_back(face);
        }
        break;
      }
    }
  }

  fclose(wrl_file);
  return true;
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template bool igl::readWRL<double, int>(std::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >&, std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > >&, std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >&);
#endif
