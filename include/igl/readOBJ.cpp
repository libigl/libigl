// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "readOBJ.h"

#include "list_to_matrix.h"

#include <iostream>
#include <cstdio>
#include <fstream>

template <typename Scalar, typename Index>
IGL_INLINE bool igl::readOBJ(
  const std::string obj_file_name,
  std::vector<std::vector<Scalar > > & V,
  std::vector<std::vector<Scalar > > & TC,
  std::vector<std::vector<Scalar > > & N,
  std::vector<std::vector<Index > > & F,
  std::vector<std::vector<Index > > & FTC,
  std::vector<std::vector<Index > > & FN)
{
  // Open file, and check for error
  FILE * obj_file = fopen(obj_file_name.c_str(),"r");
  if(NULL==obj_file)
  {
    fprintf(stderr,"IOError: %s could not be opened...\n",
            obj_file_name.c_str());
    return false;
  }
  // File open was succesfull so clear outputs
  V.clear();
  TC.clear();
  N.clear();
  F.clear();
  FTC.clear();
  FN.clear();

  // variables an constants to assist parsing the .obj file
  // Constant strings to compare against
  std::string v("v");
  std::string vn("vn");
  std::string vt("vt");
  std::string f("f");
  std::string tic_tac_toe("#");
#ifndef IGL_LINE_MAX
#  define IGL_LINE_MAX 2048
#endif

  char line[IGL_LINE_MAX];
  int line_no = 1;
  while (fgets(line, IGL_LINE_MAX, obj_file) != NULL)
  {
    char type[IGL_LINE_MAX];
    // Read first word containing type
    if(sscanf(line, "%s",type) == 1)
    {
      // Get pointer to rest of line right after type
      char * l = &line[strlen(type)];
      if(type == v)
      {
        double x[4];
        int count =
        sscanf(l,"%lf %lf %lf %lf\n",&x[0],&x[1],&x[2],&x[3]);
        if(count != 3 && count != 4)
        {
          fprintf(stderr,
                  "Error: readOBJ() vertex on line %d should have 3 or 4 coordinates",
                  line_no);
          fclose(obj_file);
          return false;
        }
        std::vector<Scalar > vertex(count);
        for(int i = 0;i<count;i++)
        {
          vertex[i] = x[i];
        }
        V.push_back(vertex);
      }else if(type == vn)
      {
        double x[3];
        int count =
        sscanf(l,"%lf %lf %lf\n",&x[0],&x[1],&x[2]);
        if(count != 3)
        {
          fprintf(stderr,
                  "Error: readOBJ() normal on line %d should have 3 coordinates",
                  line_no);
          fclose(obj_file);
          return false;
        }
        std::vector<Scalar > normal(count);
        for(int i = 0;i<count;i++)
        {
          normal[i] = x[i];
        }
        N.push_back(normal);
      }else if(type == vt)
      {
        double x[3];
        int count =
        sscanf(l,"%lf %lf %lf\n",&x[0],&x[1],&x[2]);
        if(count != 2 && count != 3)
        {
          fprintf(stderr,
                  "Error: readOBJ() texture coords on line %d should have 2 "
                  "or 3 coordinates (%d)",
                  line_no,count);
          fclose(obj_file);
          return false;
        }
        std::vector<Scalar > tex(count);
        for(int i = 0;i<count;i++)
        {
          tex[i] = x[i];
        }
        TC.push_back(tex);
      }else if(type == f)
      {
        std::vector<Index > f;
        std::vector<Index > ftc;
        std::vector<Index > fn;
        // Read each "word" after type
        char word[IGL_LINE_MAX];
        int offset;
        while(sscanf(l,"%s%n",word,&offset) == 1)
        {
          // adjust offset
          l += offset;
          // Process word
          unsigned int i,it,in;
          if(sscanf(word,"%u/%u/%u",&i,&it,&in) == 3)
          {
            f.push_back(i-1);
            ftc.push_back(it-1);
            fn.push_back(in-1);
          }else if(sscanf(word,"%u/%u",&i,&it) == 2)
          {
            f.push_back(i-1);
            ftc.push_back(it-1);
          }else if(sscanf(word,"%u//%u",&i,&in) == 2)
          {
            f.push_back(i-1);
            fn.push_back(in-1);
          }else if(sscanf(word,"%u",&i) == 1)
          {
            f.push_back(i-1);
          }else
          {
            fprintf(stderr,
                    "Error: readOBJ() face on line %d has invalid element format\n",
                    line_no);
            fclose(obj_file);
            return false;
          }
        }
        if(
           (f.size()>0 && fn.size() == 0 && ftc.size() == 0) ||
           (f.size()>0 && fn.size() == f.size() && ftc.size() == 0) ||
           (f.size()>0 && fn.size() == 0 && ftc.size() == f.size()) ||
           (f.size()>0 && fn.size() == f.size() && ftc.size() == f.size()))
        {
          // No matter what add each type to lists so that lists are the
          // correct lengths
          F.push_back(f);
          FTC.push_back(ftc);
          FN.push_back(fn);
        }else
        {
          fprintf(stderr,
                  "Error: readOBJ() face on line %d has invalid format\n", line_no);
          fclose(obj_file);
          return false;
        }
      }else if(strlen(type) >= 1 && (type[0] == '#' ||
            type[0] == 'g'  ||
            type[0] == 's'  ||
            strcmp("usemtl",type)==0 ||
            strcmp("mtllib",type)==0))
      {
        //ignore comments or other shit
      }else
      {
        //ignore any other lines
        fprintf(stderr,
                "Warning: readOBJ() ignored non-comment line %d:\n  %s",
                line_no,
                line);
      }
    }else
    {
      // ignore empty line
    }
    line_no++;
  }
  fclose(obj_file);

  assert(F.size() == FN.size());
  assert(F.size() == FTC.size());

  return true;
}

template <typename Scalar, typename Index>
IGL_INLINE bool igl::readOBJ(
  const std::string obj_file_name,
  std::vector<std::vector<Scalar > > & V,
  std::vector<std::vector<Index > > & F)
{
  std::vector<std::vector<Scalar > > TC,N;
  std::vector<std::vector<Index > > FTC,FN;
  return readOBJ(obj_file_name,V,TC,N,F,FTC,FN);
}

template <typename DerivedV, typename DerivedF, typename DerivedT>
IGL_INLINE bool igl::readOBJ(
  const std::string str,
    Eigen::PlainObjectBase<DerivedV>& V,
    Eigen::PlainObjectBase<DerivedT>& TC,
    Eigen::PlainObjectBase<DerivedV>& CN,
    Eigen::PlainObjectBase<DerivedF>& F,
    Eigen::PlainObjectBase<DerivedF>& FTC,
    Eigen::PlainObjectBase<DerivedF>& FN)
{
  std::vector<std::vector<double> > vV,vTC,vN;
  std::vector<std::vector<int> > vF,vFTC,vFN;
  bool success = igl::readOBJ(str,vV,vTC,vN,vF,vFTC,vFN);
  if(!success)
  {
    // readOBJ(str,vV,vTC,vN,vF,vFTC,vFN) should have already printed an error
    // message to stderr
    return false;
  }
  bool V_rect = igl::list_to_matrix(vV,V);
  if(!V_rect)
  {
    // igl::list_to_matrix(vV,V) already printed error message to std err
    return false;
  }
  bool F_rect = igl::list_to_matrix(vF,F);
  if(!F_rect)
  {
    // igl::list_to_matrix(vF,F) already printed error message to std err
    return false;
  }
  if(!vN.empty())
  {
    bool VN_rect = igl::list_to_matrix(vN,CN);
    if(!VN_rect)
    {
      // igl::list_to_matrix(vV,V) already printed error message to std err
      return false;
    }
  }

  if(!vFN.empty() && !vFN[0].empty())
  {
    bool FN_rect = igl::list_to_matrix(vFN,FN);
    if(!FN_rect)
    {
      // igl::list_to_matrix(vV,V) already printed error message to std err
      return false;
    }
  }

  if(!vTC.empty())
  {

    bool T_rect = igl::list_to_matrix(vTC,TC);
    if(!T_rect)
    {
      // igl::list_to_matrix(vTC,T) already printed error message to std err
      return false;
    }
  }
  if(!vFTC.empty()&& !vFTC[0].empty())
  {

    bool FTC_rect = igl::list_to_matrix(vFTC,FTC);
    if(!FTC_rect)
    {
      // igl::list_to_matrix(vF,F) already printed error message to std err
      return false;
    }
  }
  return true;
}

template <typename DerivedV, typename DerivedF>
IGL_INLINE bool igl::readOBJ(
  const std::string str,
  Eigen::PlainObjectBase<DerivedV>& V,
  Eigen::PlainObjectBase<DerivedF>& F)
{
  std::vector<std::vector<double> > vV,vTC,vN;
  std::vector<std::vector<int> > vF,vFTC,vFN;
  bool success = igl::readOBJ(str,vV,vTC,vN,vF,vFTC,vFN);
  if(!success)
  {
    // readOBJ(str,vV,vTC,vN,vF,vFTC,vFN) should have already printed an error
    // message to stderr
    return false;
  }
  bool V_rect = igl::list_to_matrix(vV,V);
  if(!V_rect)
  {
    // igl::list_to_matrix(vV,V) already printed error message to std err
    return false;
  }
  bool F_rect = igl::list_to_matrix(vF,F);
  if(!F_rect)
  {
    // igl::list_to_matrix(vF,F) already printed error message to std err
    return false;
  }
  return true;
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
// generated by autoexplicit.sh
template bool igl::readOBJ<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(std::basic_string<char, std::char_traits<char>, std::allocator<char> >, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template bool igl::readOBJ<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3> >(std::string, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&);
template bool igl::readOBJ<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(std::basic_string<char, std::char_traits<char>, std::allocator<char> >, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
#endif
