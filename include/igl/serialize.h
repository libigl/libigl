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
// the interface Serializable or SerializableBase.
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
#include <fstream>
#include <cstdint>
#include <numeric>
#include <vector>
#include <set>
#include <map>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "igl_inline.h"

//#define SERIALIZE(x) igl::serialize(x,#x,buffer);
//#define DESERIALIZE(x) igl::deserialize(x,#x,buffer);

namespace igl
{

  // Serializes the given object either to a file or to a provided buffer
  // Templates:
  //   T  type of the object to serialize
  // Inputs:
  //   obj        object to serialize
  //   objectName unique object name,used for the identification
  //   overwrite  set to true to overwrite an existing file
  //   filename   name of the file containing the serialization
  // Outputs:
  //   buffer     binary serialization
  //
  template <typename T>
  IGL_INLINE bool serialize(const T& obj,const std::string& filename);
  template <typename T>
  IGL_INLINE bool serialize(const T& obj,const std::string& objectName,const std::string& filename,bool overwrite = false);
  template <typename T>
  IGL_INLINE bool serialize(const T& obj,const std::string& objectName,std::vector<char>& buffer);

  // Deserializes the given data from a file or buffer back to the provided object
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
  IGL_INLINE bool deserialize(T& obj,const std::string& filename);
  template <typename T>
  IGL_INLINE bool deserialize(T& obj,const std::string& objectName,const std::string& filename);
  template <typename T>
  IGL_INLINE bool deserialize(T& obj,const std::string& objectName,const std::vector<char>& buffer);

  // User defined types have to derive from the class Serializable
  // and add their member variables in InitSerialization like the
  // following:
  //
  // class Test : public igl::Serializable {
  //
  //   int var;
  //
  //   void InitSerialization() {
  //     this->Add(var,"var");
  //   }
  // };

  // Base interface for user defined types
  struct SerializableBase
  {
    virtual void Serialize(std::vector<char>& buffer) const = 0;
    virtual void Deserialize(const std::vector<char>& buffer) = 0;
  };

  // Convenient interface for user defined types
  class Serializable: public SerializableBase
  {
  private:

    template <typename T>
    struct SerializationObject : public SerializableBase
    {
      bool Binary;
      std::string Name;
      T* Object;

      void Serialize(std::vector<char>& buffer) const {
        igl::serialize(*Object,Name,buffer);
      }

      void Deserialize(const std::vector<char>& buffer) {
        igl::deserialize(*Object,Name,buffer);
      }
    };

    mutable bool initialized;
    mutable std::vector<SerializableBase*> objects;

  public:

    // Override this function to add your member variables which should be serialized
    IGL_INLINE virtual void InitSerialization() = 0;

    // Following functions can be overridden to handle the specific events.
    // Return false to prevent the de-/serialization of an object.
    IGL_INLINE virtual bool PreSerialization() const;
    IGL_INLINE virtual void PostSerialization() const;
    IGL_INLINE virtual bool PreDeserialization();
    IGL_INLINE virtual void PostDeserialization();

    // Default implementation of SerializableBase interface
    IGL_INLINE void Serialize(std::vector<char>& buffer) const;
    IGL_INLINE void Deserialize(const std::vector<char>& buffer);

    // Default constructor, destructor, assignment and copy constructor
    IGL_INLINE Serializable();
    IGL_INLINE Serializable(const Serializable& obj);
    IGL_INLINE ~Serializable();
    IGL_INLINE Serializable& operator=(const Serializable& obj);

    // Use this function to add your variables which should be serialized
    template <typename T>
    IGL_INLINE void Add(T& obj,std::string name,bool binary = false);
  };

  // internal functions
  namespace detail
  {
    // fundamental types
    template <typename T>
    IGL_INLINE typename std::enable_if<std::is_fundamental<T>::value,size_t>::type getByteSize(const T& obj);
    template <typename T>
    IGL_INLINE typename std::enable_if<std::is_fundamental<T>::value>::type serialize(const T& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter);
    template <typename T>
    IGL_INLINE typename std::enable_if<std::is_fundamental<T>::value>::type deserialize(T& obj,std::vector<char>::const_iterator& iter);

    // std::string
    IGL_INLINE size_t getByteSize(const std::string& obj);
    IGL_INLINE void serialize(const std::string& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter);
    IGL_INLINE void deserialize(std::string& obj,std::vector<char>::const_iterator& iter);

    // SerializableBase
    template <typename T>
    IGL_INLINE typename std::enable_if<std::is_base_of<SerializableBase,T>::value,size_t>::type getByteSize(const T& obj);
    template <typename T>
    IGL_INLINE typename std::enable_if<std::is_base_of<SerializableBase,T>::value>::type serialize(const T& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter);
    template <typename T>
    IGL_INLINE typename std::enable_if<std::is_base_of<SerializableBase,T>::value>::type deserialize(T& obj,std::vector<char>::const_iterator& iter);

    // stl containers
    // std::pair
    template <typename T1,typename T2>
    IGL_INLINE size_t getByteSize(const std::pair<T1,T2>& obj);
    template <typename T1,typename T2>
    IGL_INLINE void serialize(const std::pair<T1,T2>& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter);
    template <typename T1,typename T2>
    IGL_INLINE void deserialize(std::pair<T1,T2>& obj,std::vector<char>::const_iterator& iter);

    // std::vector
    template <typename T1,typename T2>
    IGL_INLINE size_t getByteSize(const std::vector<T1,T2>& obj);
    template <typename T1,typename T2>
    IGL_INLINE void serialize(const std::vector<T1,T2>& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter);
    template <typename T1,typename T2>
    IGL_INLINE void deserialize(std::vector<T1,T2>& obj,std::vector<char>::const_iterator& iter);

    // std::set
    template <typename T>
    IGL_INLINE size_t getByteSize(const std::set<T>& obj);
    template <typename T>
    IGL_INLINE void serialize(const std::set<T>& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter);
    template <typename T>
    IGL_INLINE void deserialize(std::set<T>& obj,std::vector<char>::const_iterator& iter);

    // std::map
    template <typename T1,typename T2>
    IGL_INLINE size_t getByteSize(const std::map<T1,T2>& obj);
    template <typename T1,typename T2>
    IGL_INLINE void serialize(const std::map<T1,T2>& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter);
    template <typename T1,typename T2>
    IGL_INLINE void deserialize(std::map<T1,T2>& obj,std::vector<char>::const_iterator& iter);

    // Eigen types
    template<typename T,int R,int C,int P,int MR,int MC>
    IGL_INLINE size_t getByteSize(const Eigen::Matrix<T,R,C,P,MR,MC>& obj);
    template<typename T,int R,int C,int P,int MR,int MC>
    IGL_INLINE void serialize(const Eigen::Matrix<T,R,C,P,MR,MC>& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter);
    template<typename T,int R,int C,int P,int MR,int MC>
    IGL_INLINE void deserialize(Eigen::Matrix<T,R,C,P,MR,MC>& obj,std::vector<char>::const_iterator& iter);

    template<typename T,int P,typename I>
    IGL_INLINE size_t getByteSize(const Eigen::SparseMatrix<T,P,I>& obj);
    template<typename T,int P,typename I>
    IGL_INLINE void serialize(const Eigen::SparseMatrix<T,P,I>& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter);
    template<typename T,int P,typename I>
    IGL_INLINE void deserialize(Eigen::SparseMatrix<T,P,I>& obj,std::vector<char>::const_iterator& iter);

    // pointers
    template <typename T>
    IGL_INLINE typename std::enable_if<std::is_pointer<T>::value,size_t>::type getByteSize(const T& obj);
    template <typename T>
    IGL_INLINE typename std::enable_if<std::is_pointer<T>::value>::type serialize(const T& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter);
    template <typename T>
    IGL_INLINE typename std::enable_if<std::is_pointer<T>::value>::type deserialize(T& obj,std::vector<char>::const_iterator& iter);

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
      static const bool value = std::is_fundamental<T0>::value || std::is_same<std::string,T0>::value || std::is_base_of<SerializableBase,T0>::value
        || is_stl_container<T0>::value || is_eigen_type<T0>::value;
    };
  }
}

#ifndef IGL_STATIC_LIBRARY
  #include "serialize.cpp"
#endif

#endif
