// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
// ---------------------------------------------------------------------------
// serialize.h
// author: Christian Schüller <schuellchr@gmail.com>
// -----------------------------------------------------------------------------
// Functions to save and load a serialization of fundamental c++ data types to
// and from a xml file. STL containers, Eigen matrix types and nested data
// structures are also supported. To serialize a user defined class implement
// the interface XMLSerializable.
//
// See also: serialize.h
// -----------------------------------------------------------------------------

#ifndef IGL_SERIALIZABLE_XML_H
#define IGL_SERIALIZABLE_XML_H

#include <type_traits>
#include <iostream>
#include <vector>
#include <set>
#include <map>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <igl/serialize.h>

#include "tinyxml2.h"

//#define SERIALIZE_XML(x) igl::serialize_xml(x,#x,doc,element);
//#define DESERIALIZE_XML(x) igl::deserialize_xml(x,#x,,doc,element);

namespace igl
{
  // serializes the given object either to a xml file or to the provided doc data
  //
  // Templates:
  //   T  type of the object to serialize
  // Inputs:
  //   obj        object to serialize
  //   objectName unique object name,used for the identification
  //   filename   name of the file containing the serialization
  //   binary     set to true to serialize the object in binary format (faster for big data)
  //   overwrite  set to true to update the serialization in an existing xml file
  //   element    tinyxml2 virtual representation of the current xml node
  // Outputs: 
  //   doc        contains current tinyxml2 virtual representation of the xml data
  //
  template <typename T>
  void serialize_xml(const T& obj,const std::string& filename);
  template <typename T>
  void serialize_xml(const T& obj,const std::string& objectName,const std::string& filename, bool binary = false,bool overwrite = false);
  template <typename T>
  void serialize_xml(const T& obj,const std::string& objectName,tinyxml2::XMLDocument* doc,tinyxml2::XMLElement* element,bool binary = false);

  // deserializes the given data from a xml file or doc data back to the provided object
  //
  // Templates:
  //   T  type of the object to serialize
  // Inputs:
  //
  //   objectName unique object name,used for the identification
  //   filename   name of the file containing the serialization
  //   binary     set to true to serialize the object in binary format (faster for big data)
  //   overwrite  set to true to update the serialization in an existing xml file
  //   doc        contains current tinyxml2 virtual representation of the xml data
  //   element    tinyxml2 virtual representation of the current xml node
  // Outputs: 
  //   obj        object to load back serialization to 
  //
  template <typename T>
  void deserialize_xml(T& obj,const std::string& filename);
  template <typename T>
  void deserialize_xml(T& obj,const std::string& objectName,const std::string& filename);
  template <typename T>
  void deserialize_xml(T& obj,const std::string& objectName,const tinyxml2::XMLDocument* doc,const tinyxml2::XMLElement* element);

  // interface for user defined types
  struct XMLSerializable : public Serializable
  {
    virtual void Serialize(std::vector<char>& buffer) const = 0;
    virtual void Deserialize(const std::vector<char>& buffer) = 0;
    virtual void Serialize(tinyxml2::XMLDocument* doc,tinyxml2::XMLElement* element) const = 0;
    virtual void Deserialize(const tinyxml2::XMLDocument* doc,const tinyxml2::XMLElement* element) = 0;
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
  //   void Serialize(tinyxml2::XMLDocument* doc,tinyxml2::XMLElement* element) const {
  //     serialize_xml(var,"var1",doc,element);
  //   }
  //   void Deserialize(const tinyxml2::XMLDocument* doc,const tinyxml2::XMLElement* element) {
  //     deserialize_xml(var,"var1",doc,element);
  //   }
  // }

  // internal functions
  namespace detail_xml
  {
    // fundamental types
    template <typename T>
    std::enable_if_t<std::is_fundamental<T>::value> serialize(const T& obj,tinyxml2::XMLDocument* doc,tinyxml2::XMLElement* element,const std::string& name);
    template <typename T>
    std::enable_if_t<std::is_fundamental<T>::value> deserialize(T& obj,const tinyxml2::XMLDocument* doc,const tinyxml2::XMLElement* element,const std::string& name);

    // std::string
    void serialize(const std::string& obj,tinyxml2::XMLDocument* doc,tinyxml2::XMLElement* element,const std::string& name);
    void deserialize(std::string& obj,const tinyxml2::XMLDocument* doc,const tinyxml2::XMLElement* element,const std::string& name);

    // Serializable
    template <typename T>
    std::enable_if_t<std::is_base_of<XMLSerializable,T>::value> serialize(const T& obj,tinyxml2::XMLDocument* doc,tinyxml2::XMLElement* element,const std::string& name);
    template <typename T>
    std::enable_if_t<std::is_base_of<XMLSerializable,T>::value> deserialize(T& obj,const tinyxml2::XMLDocument* doc,const tinyxml2::XMLElement* element,const std::string& name);

    // STL containers
    template <typename T1, typename T2>
    void serialize(const std::pair<T1,T2>& obj,tinyxml2::XMLDocument* doc,tinyxml2::XMLElement* element,const std::string& name);
    template <typename T1,typename T2>
    void deserialize(std::pair<T1,T2>& obj,const tinyxml2::XMLDocument* doc,const tinyxml2::XMLElement* element,const std::string& name);

    template <typename T1,typename T2>
    void serialize(const std::vector<T1,T2>& obj,tinyxml2::XMLDocument* doc,tinyxml2::XMLElement* element,const std::string& name);
    template <typename T1,typename T2>
    void deserialize(std::vector<T1,T2>& obj,const tinyxml2::XMLDocument* doc,const tinyxml2::XMLElement* element,const std::string& name);

    template <typename T>
    void serialize(const std::set<T>& obj,tinyxml2::XMLDocument* doc,tinyxml2::XMLElement* element,const std::string& name);
    template <typename T>
    void deserialize(std::set<T>& obj,const tinyxml2::XMLDocument* doc,const tinyxml2::XMLElement* element,const std::string& name);

    template <typename T1,typename T2>
    void serialize(const std::map<T1,T2>& obj,tinyxml2::XMLDocument* doc,tinyxml2::XMLElement* element,const std::string& name);
    template <typename T1,typename T2>
    void deserialize(std::map<T1,T2>& obj,const tinyxml2::XMLDocument* doc,const tinyxml2::XMLElement* element,const std::string& name);

    // Eigen types
    template<typename T,int R,int C,int P,int MR,int MC>
    void serialize(const Eigen::Matrix<T,R,C,P,MR,MC>& obj,tinyxml2::XMLDocument* doc,tinyxml2::XMLElement* element,const std::string& name);
    template<typename T,int R,int C,int P,int MR,int MC>
    void deserialize(Eigen::Matrix<T,R,C,P,MR,MC>& obj,const tinyxml2::XMLDocument* doc,const tinyxml2::XMLElement* element,const std::string& name);

    template<typename T,int P,typename I>
    void serialize(const Eigen::SparseMatrix<T,P,I>& obj,tinyxml2::XMLDocument* doc,tinyxml2::XMLElement* element,const std::string& name);
    template<typename T,int P,typename I>
    void deserialize(Eigen::SparseMatrix<T,P,I>& obj,const tinyxml2::XMLDocument* doc,const tinyxml2::XMLElement* element,const std::string& name);

    // pointers
    template <typename T>
    std::enable_if_t<std::is_pointer<T>::value> serialize(const T& obj,tinyxml2::XMLDocument* doc,tinyxml2::XMLElement* element,const std::string& name);
    template <typename T>
    std::enable_if_t<std::is_pointer<T>::value> deserialize(T& obj,const tinyxml2::XMLDocument* doc,const tinyxml2::XMLElement* element,const std::string& name);

    // helper functions
    tinyxml2::XMLElement* getElement(tinyxml2::XMLDocument* doc,tinyxml2::XMLElement* element,const std::string& name);
    void getAttribute(const char* src,bool& dest);
    void getAttribute(const char* scr,char& dest);
    void getAttribute(const char* src,std::string& dest);
    void getAttribute(const char* src,float& dest);
    void getAttribute(const char* src,double& dest);
    template<typename T>
    std::enable_if_t<std::is_integral<T>::value && std::is_unsigned<T>::value> getAttribute(const char* src,T& dest);
    template<typename T>
    std::enable_if_t<std::is_integral<T>::value && !std::is_unsigned<T>::value> getAttribute(const char* src,T& dest);
    void replaceSubString(std::string& str,const std::string& search,const std::string& replace);
    void encodeXMLElementName(std::string& name);
    void decodeXMLElementName(std::string& name);
    std::string base64_encode(unsigned char const* bytes_to_encode,unsigned int in_len);
    std::string base64_decode(std::string const& encoded_string);

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

#include "serialize_xml.cpp";

#endif