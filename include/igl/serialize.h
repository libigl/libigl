// 
// Copyright (C) 2014 Christian Schüller <schuellchr@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
// -----------------------------------------------------------------------------
// Functions to save and load a serialization of fundamental c++ data types to
// and from a binary file. STL containers, Eigen matrix types and nested data
// structures are also supported. To serialize a user defined class implement
// the interface Serializable.
//
// See also: xml/serialize_xml.h
// -----------------------------------------------------------------------------
// TODOs:
// * loops of pointers and from multiple objects
// * cross-platform compatibility (big-, little-endian)
// -----------------------------------------------------------------------------

#ifndef IGL_SERIALIZE_H
#define IGL_SERIALIZE_H

#include <type_traits>
#include <iostream>
#include <numeric>
#include <vector>
#include <set>
#include <map>

#include <Eigen/Dense>
#include <Eigen/Sparse>

//#define SERIALIZE(x) igl::serialize(x,#x,buffer);
//#define DESERIALIZE(x) igl::deserialize(x,#x,buffer);

namespace igl
{

  // serializes the given object either to a file or to a provided buffer
  // Templates:
  //   T  type of the object to serialize
  // Inputs:
  //   obj        object to serialize
  //   objectName unique object name,used for the identification
  //   filename   name of the file containing the serialization
  // Outputs: 
  //   buffer     binary serialization
  //
  template <typename T>
  void serialize(const T& obj,const std::string& filename);
  template <typename T>
  void serialize(const T& obj,const std::string& objectName,std::vector<char>& buffer);

  // deserializes the given data from a file or buffer back to the provided object
  //
  // Templates:
  //   T  type of the object to serialize
  // Inputs:
  //   buffer     binary serialization
  //   objectName unique object name, used for the identification
  //   filename   name of the file containing the serialization
  // Outputs: 
  //   obj        object to load back serialization to 
  //
  template <typename T>
  void deserialize(T& obj,const std::string& filename);
  template <typename T>
  void deserialize(T& obj,const std::string& objectName,const std::vector<char>& buffer);

  // interface for user defined types
  struct Serializable
  {
    virtual void Serialize(std::vector<char>& buffer) const = 0;
    virtual void Deserialize(const std::vector<char>& buffer) = 0;
  };
  // example:
  //
  // class Test : public igl::Serializable {
  //   
  //   int var;
  // 
  //   void Serialize(std::vector<char>& buffer) {
  //     serialize(var,"var1",buffer);
  //   }
  //   void Deserialize(const std::vector<char>& buffer) {
  //     deserialize(var,"var1",buffer);
  //   }
  // }

  // internal functions
  namespace detail
  {
    // fundamental types
    template <typename T>
    std::enable_if_t<std::is_fundamental<T>::value,size_t> getByteSize(const T& obj);
    template <typename T>
    std::enable_if_t<std::is_fundamental<T>::value> serialize(const T& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter);
    template <typename T>
    std::enable_if_t<std::is_fundamental<T>::value> deserialize(T& obj,std::vector<char>::const_iterator& iter);

    // std::string
    size_t getByteSize(const std::string& obj);
    void serialize(const std::string& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter);
    void deserialize(std::string& obj,std::vector<char>::iterator& iter);

    // Serializable
    template <typename T>
    std::enable_if_t<std::is_base_of<Serializable,T>::value,size_t> getByteSize(const T& obj);
    template <typename T>
    std::enable_if_t<std::is_base_of<Serializable,T>::value> serialize(const T& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter);
    template <typename T>
    std::enable_if_t<std::is_base_of<Serializable,T>::value> deserialize(T& obj,std::vector<char>::const_iterator& iter);

    // stl containers
    // std::pair
    template <typename T1,typename T2>
    size_t getByteSize(const std::pair<T1,T2>& obj);
    template <typename T1,typename T2>
    void serialize(const std::pair<T1,T2>& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter);
    template <typename T1,typename T2>
    void deserialize(std::pair<T1,T2>& obj,std::vector<char>::const_iterator& iter);

    // std::vector
    template <typename T1,typename T2>
    size_t getByteSize(const std::vector<T1,T2>& obj);
    template <typename T1,typename T2>
    void serialize(const std::vector<T1,T2>& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter);
    template <typename T1,typename T2>
    void deserialize(std::vector<T1,T2>& obj,std::vector<char>::const_iterator& iter);

    // std::set
    template <typename T>
    size_t getByteSize(const std::set<T>& obj);
    template <typename T>
    void serialize(const std::set<T>& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter);
    template <typename T>
    void deserialize(std::set<T>& obj,std::vector<char>::const_iterator& iter);

    // std::map
    template <typename T1,typename T2>
    size_t getByteSize(const std::map<T1,T2>& obj);
    template <typename T1,typename T2>
    void serialize(const std::map<T1,T2>& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter);
    template <typename T1,typename T2>
    void deserialize(std::map<T1,T2>& obj,std::vector<char>::const_iterator& iter);

    // Eigen types
    template<typename T,int R,int C,int P,int MR,int MC>
    size_t getByteSize(const Eigen::Matrix<T,R,C,P,MR,MC>& obj);
    template<typename T,int R,int C,int P,int MR,int MC>
    void serialize(const Eigen::Matrix<T,R,C,P,MR,MC>& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter);
    template<typename T,int R,int C,int P,int MR,int MC>
    void deserialize(Eigen::Matrix<T,R,C,P,MR,MC>& obj,std::vector<char>::const_iterator& iter);

    template<typename T,int P,typename I>
    size_t getByteSize(const Eigen::SparseMatrix<T,P,I>& obj);
    template<typename T,int P,typename I>
    void serialize(const Eigen::SparseMatrix<T,P,I>& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter);
    template<typename T,int P,typename I>
    void deserialize(Eigen::SparseMatrix<T,P,I>& obj,std::vector<char>::const_iterator& iter);

    // pointers
    template <typename T>
    std::enable_if_t<std::is_pointer<T>::value,size_t> getByteSize(const T& obj);
    template <typename T>
    std::enable_if_t<std::is_pointer<T>::value> serialize(const T& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter);
    template <typename T>
    std::enable_if_t<std::is_pointer<T>::value> deserialize(T& obj,std::vector<char>::const_iterator& iter);

    // compile time type serializable check
    template <typename T>
    struct is_stl_container { static const bool value = false; };
    template <typename T1,typename T2>
    struct is_stl_container<std::pair<T1,T2> > { static const bool value = true; };
    template <typename T1,typename T2>
    struct is_stl_container<std::vector<T1,T2> > { static const bool value = true; };
    template <typename T>
    struct is_stl_container<std::set<T> > { static const bool value = true; };
    template <typename T1,typename T2>
    struct is_stl_container<std::map<T1,T2> > { static const bool value = true; };

    template <typename T>
    struct is_eigen_type { static const bool value = false; };
    template <typename T,int R,int C,int P,int MR,int MC>
    struct is_eigen_type<Eigen::Matrix<T,R,C,P,MR,MC> > { static const bool value = true; };
    template <typename T,int P,typename I>
    struct is_eigen_type<Eigen::SparseMatrix<T,P,I> > { static const bool value = true; };

    template <typename T>
    struct is_serializable {
      using T0 = typename  std::remove_pointer<T>::type;
      static const bool value = std::is_fundamental<T0>::value || std::is_same<std::string,T0>::value || std::is_base_of<Serializable,T0>::value
        || is_stl_container<T0>::value || is_eigen_type<T0>::value;
    };
  }
}

#include "serialize.cpp"

#endif