// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "readSTL.h"
#include "list_to_matrix.h"

#include <iostream>
template <typename DerivedV, typename DerivedF, typename DerivedN>
IGL_INLINE bool igl::readSTL(
  const std::string & filename,
  Eigen::PlainObjectBase<DerivedV> & V,
  Eigen::PlainObjectBase<DerivedF> & F,
  Eigen::PlainObjectBase<DerivedN> & N)
{
  using namespace std;
  vector<vector<typename DerivedV::Scalar> > vV;
  vector<vector<typename DerivedN::Scalar> > vN;
  vector<vector<typename DerivedF::Scalar> > vF;
  if(!readSTL(filename,vV,vF,vN))
  {
    return false;
  }

  if(!list_to_matrix(vV,V))
  {
    return false;
  }

  if(!list_to_matrix(vF,F))
  {
    return false;
  }

  if(!list_to_matrix(vN,N))
  {
    return false;
  }
  return true;
}

template <typename TypeV, typename TypeF, typename TypeN>
IGL_INLINE bool igl::readSTL(
  const std::string & filename,
  std::vector<std::vector<TypeV> > & V,
  std::vector<std::vector<TypeF> > & F,
  std::vector<std::vector<TypeN> > & N)
{
  using namespace std;
  // Should test for ascii

  // Open file, and check for error
  FILE * stl_file = fopen(filename.c_str(),"r");
  if(NULL==stl_file)
  {
    fprintf(stderr,"IOError: %s could not be opened...\n",
            filename.c_str());
    return false;
  }
  V.clear();
  F.clear();
  N.clear();
#ifndef IGL_LINE_MAX
#  define IGL_LINE_MAX 2048
#endif
  char solid[IGL_LINE_MAX];
  if(fscanf(stl_file,"%s",solid)!=1)
  {
    // file too short
    cerr<<"IOError: "<<filename<<" too short."<<endl;
    goto close_false;
  }
  if(string("solid") == solid)
  {
    // Eat file name
    char name[IGL_LINE_MAX];
    if(NULL==fgets(name,IGL_LINE_MAX,stl_file))
    {
      cerr<<"IOError: "<<filename<<" too short."<<endl;
      goto close_false;
    }
    // ascii
    while(true)
    {
      int ret;
      char facet[IGL_LINE_MAX],normal[IGL_LINE_MAX];
      vector<TypeN > n(3);
      double nd[3];
      ret = fscanf(stl_file,"%s %s %lg %lg %lg",facet,normal,nd,nd+1,nd+2);
      if(string("endsolid") == facet)
      {
        break;
      }
      if(ret != 5 || string("facet") != facet || string("normal") != normal)
      {
        cerr<<"IOError: "<<filename<<" bad format (1)."<<endl;
        goto close_false;
      }
      // copy casts to Type
      n[0] = nd[0]; n[1] = nd[1]; n[2] = nd[2];
      N.push_back(n);
      char outer[IGL_LINE_MAX], loop[IGL_LINE_MAX];
      ret = fscanf(stl_file,"%s %s",outer,loop);
      if(ret != 2 || string("outer") != outer || string("loop") != loop)
      {
        cerr<<"IOError: "<<filename<<" bad format (2)."<<endl;
        goto close_false;
      }
      vector<TypeF> f;
      while(true)
      {
        char word[IGL_LINE_MAX];
        int ret = fscanf(stl_file,"%s",word);
        if(ret == 1 && string("endloop") == word)
        {
          break;
        }else if(ret == 1 && string("vertex") == word)
        {
          vector<TypeV> v(3);
          double vd[3];
          int ret = fscanf(stl_file,"%lg %lg %lg",vd,vd+1,vd+2);
          if(ret != 3)
          {
            cerr<<"IOError: "<<filename<<" bad format (3)."<<endl;
            goto close_false;
          }
          f.push_back(V.size());
          // copy casts to Type
          v[0] = vd[0]; v[1] = vd[1]; v[2] = vd[2];
          V.push_back(v);
        }else
        {
          cerr<<"IOError: "<<filename<<" bad format (4)."<<endl;
          goto close_false;
        }
      }
      F.push_back(f);
      char endfacet[IGL_LINE_MAX];
      ret = fscanf(stl_file,"%s",endfacet);
      if(ret != 1 || string("endfacet") != endfacet)
      {
        cerr<<"IOError: "<<filename<<" bad format (5)."<<endl;
        goto close_false;
      }
    }
    // read endfacet
    fclose(stl_file);
    goto close_true;
  }else
  {
    // Binary
    fclose(stl_file);
    stl_file = fopen(filename.c_str(),"rb");
    if(NULL==stl_file)
    {
      fprintf(stderr,"IOError: %s could not be opened...\n",
              filename.c_str());
      return false;
    }
    // Read 80 header
    char header[80];
    if(fread(header,sizeof(char),80,stl_file)!=80)
    {
      cerr<<"IOError: "<<filename<<" bad format (1)."<<endl;
      goto close_false;
    }
    // Read number of triangles
    unsigned int num_tri;
    if(fread(&num_tri,sizeof(unsigned int),1,stl_file)!=1)
    {
      cerr<<"IOError: "<<filename<<" bad format (2)."<<endl;
      goto close_false;
    }
    V.resize(num_tri*3,vector<TypeV >(3,0));
    N.resize(num_tri,vector<TypeN >(3,0));
    F.resize(num_tri,vector<TypeF >(3,0));
    for(int t = 0;t<(int)num_tri;t++)
    {
      // Read normal
      float n[3];
      if(fread(n,sizeof(float),3,stl_file)!=3)
      {
        cerr<<"IOError: "<<filename<<" bad format (3)."<<endl;
        goto close_false;
      }
      // Read each vertex
      for(int c = 0;c<3;c++)
      {
        F[t][c] = 3*t+c;
        N[t][c] = n[c];
        float v[3];
        if(fread(v,sizeof(float),3,stl_file)!=3)
        {
          cerr<<"IOError: "<<filename<<" bad format (4)."<<endl;
          goto close_false;
        }
        V[3*t+c][0] = v[0];
        V[3*t+c][1] = v[1];
        V[3*t+c][2] = v[2];
      }
      // Read attribute size
      unsigned short att_count;
      if(fread(&att_count,sizeof(unsigned short),1,stl_file)!=1)
      {
        cerr<<"IOError: "<<filename<<" bad format (5)."<<endl;
        goto close_false;
      }
    }
    goto close_true;
  }
close_false:
  fclose(stl_file);
  return false;
close_true:
  fclose(stl_file);
  return true;
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
// generated by autoexplicit.sh
template bool igl::readSTL<Eigen::Matrix<float, -1, 3, 1, -1, 3>, Eigen::Matrix<unsigned int, -1, 3, 1, -1, 3>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(std::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, Eigen::PlainObjectBase<Eigen::Matrix<float, -1, 3, 1, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<unsigned int, -1, 3, 1, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
// generated by autoexplicit.sh
template bool igl::readSTL<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(std::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
// generated by autoexplicit.sh
template bool igl::readSTL<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(std::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template bool igl::readSTL<Eigen::Matrix<float, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<float, -1, -1, 0, -1, -1> >(std::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, Eigen::PlainObjectBase<Eigen::Matrix<float, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<float, -1, -1, 0, -1, -1> >&);
#endif
