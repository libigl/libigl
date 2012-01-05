//
//  IGL Lib - Simple C++ mesh library 
//
//  Copyright 2011, Daniele Panozzo. All rights reserved.

// History:
//  return type changed from void to bool  Alec 18 Sept 2011
//  added pure vector of vectors version that has much more support Alec 31 Oct
//    2011

#ifndef IGL_READOBJ_H
#define IGL_READOBJ_H

#include <Eigen/Core>
#include <string>
#include <vector>

namespace igl 
{
  // Read a mesh from an ascii obj file, filling in vertex positions, normals
  // and texture coordinates. Mesh may have faces of any number of degree
  //
  // Templates:
  //   Scalar  type for positions and vectors (will be read as double and cast
  //     to Scalar)
  //   Index  type for indices (will be read as int and cast to Index)
  // Inputs:
  //  str  path to .obj file
  // Outputs:
  //   V  double matrix of vertex positions  #V by 3
  //   F  #F list of face indices into vertex positions
  //   TC  double matrix of texture coordinats #TC by 2
  //   FTC  #F list of face indices into vertex texture coordinates
  //   N  double matrix of corner normals #N by 3
  //   FN  #F list of face indices into vertex normals
  // Returns true on success, false on errors
  template <typename Scalar, typename Index>
  inline bool readOBJ(
    const std::string obj_file_name, 
    std::vector<std::vector<Scalar > > & V,
    std::vector<std::vector<Scalar > > & TC,
    std::vector<std::vector<Scalar > > & N,
    std::vector<std::vector<Index > > & F,
    std::vector<std::vector<Index > > & FTC,
    std::vector<std::vector<Index > > & FN);

  //! Read a mesh from an ascii obj file
  // Inputs:
  //   str  path to .obj file
  // Outputs:
  //   V  eigen double matrix #V by 3
  //   F  eigen int matrix #F by 3
  //
  // KNOWN BUG: This only knows how to read *triangle* meshes. It will probably
  // crash or give garbage on anything else.
  //
  // KNOWN BUG: This only knows how to face lines without normal or texture
  // indices. It will probably crash or give garbage on anything else.
  inline bool readOBJ(const std::string str, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
}

// Implementation
#include "list_to_matrix.h"

#include <iostream>
#include <fstream>

template <typename Scalar, typename Index>
inline bool igl::readOBJ(
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
  // flag for whether vertex texture coordinates exist in file
  bool has_texture = false;
  // flag for whether vertex normals exist in file
  bool has_normals = false;
  // Constant strings to compare against
  std::string v("v");
  std::string vn("vn");
  std::string vt("vt");
  std::string f("f");
  std::string tic_tac_toe("#");
#ifndef LINE_MAX
#  define LINE_MAX 2048
#endif

  char line[LINE_MAX];
  int line_no = 1;
  while (fgets(line, LINE_MAX, obj_file) != NULL) 
  {
    char type[LINE_MAX];
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
        has_normals = true;
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
        has_texture = true;
        double x[3];
        int count = 
          sscanf(l,"%lf %lf %lf\n",&x[0],&x[1],&x[2]);
        if(count != 2 && count != 3)
        {
          fprintf(stderr, 
            "Error: readOBJ() vertex on line %d should have 2 or 3 coordinates", 
            line_no);
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
        char word[LINE_MAX];
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
      }else if(strlen(type) >= 1 && type[0] == '#')
      {
        //ignore comments
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

inline bool igl::readOBJ(const std::string str, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
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
  // Legacy
  if(F.cols() != 3)
  {
    fprintf(stderr,
      "Error: readOBJ(filename,V,F) is meant for reading triangle-only"
      " meshes. This mesh has faces all with size %d. See readOBJ.h for other"
      " options.\n",
      (int)F.cols());
    return false;
  }
  return true;
}
#endif
