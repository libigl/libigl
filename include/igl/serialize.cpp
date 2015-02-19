//
// Copyright (C) 2014 Christian Schüller <schuellchr@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "serialize.h"

namespace igl
{
  template <typename T>
  IGL_INLINE bool serialize(const T& obj,const std::string& filename)
  {
    return serialize(obj,"obj",filename,true);
  }

  template <typename T>
  IGL_INLINE bool serialize(const T& obj,const std::string& objectName,const std::string& filename,bool overwrite)
  {
    bool success = false;

    std::vector<char> buffer;

    std::ios_base::openmode mode = std::ios::out | std::ios::binary;

    if(overwrite)
      mode |= std::ios::trunc;
    else
      mode |= std::ios::app;

    std::ofstream file(filename.c_str(),mode);

    if(file.is_open())
    {
      serialize(obj,objectName,buffer);

      file.write(&buffer[0],buffer.size());

      file.close();

      success = true;
    }
    else
    {
      std::cerr << "saving binary serialization failed!" << std::endl;
    }

    return success;
  }

  template <typename T>
  IGL_INLINE bool serialize(const T& obj,const std::string& objectName,std::vector<char>& buffer)
  {
    // serialize object data
    size_t size = serialization::getByteSize(obj);
    std::vector<char> tmp(size);
    auto it = tmp.begin();
    serialization::serialize(obj,tmp,it);

    std::string objectType(typeid(obj).name());
    size_t newObjectSize = tmp.size();
    size_t newHeaderSize = serialization::getByteSize(objectName) + serialization::getByteSize(objectType) + sizeof(size_t);
    size_t curSize = buffer.size();
    size_t newSize = curSize + newHeaderSize + newObjectSize;

    buffer.resize(newSize);

    std::vector<char>::iterator iter = buffer.begin()+curSize;

    // serialize object header (name/type/size)
    serialization::serialize(objectName,buffer,iter);
    serialization::serialize(objectType,buffer,iter);
    serialization::serialize(newObjectSize,buffer,iter);

    // copy serialized data to buffer
    iter = std::copy(tmp.begin(),tmp.end(),iter);

    return true;
  }

  template <typename T>
  IGL_INLINE bool deserialize(T& obj,const std::string& filename)
  {
    return deserialize(obj,"obj",filename);
  }

  template <typename T>
  IGL_INLINE bool deserialize(T& obj,const std::string& objectName,const std::string& filename)
  {
    bool success = false;

    std::ifstream file(filename.c_str(),std::ios::binary);

    if(file.is_open())
    {
      file.seekg(0,std::ios::end);
      std::streamoff size = file.tellg();
      file.seekg(0,std::ios::beg);

      std::vector<char> buffer(size);
      file.read(&buffer[0],size);

      deserialize(obj,objectName,buffer);
      file.close();

      success = true;
    }
    else
    {
      std::cerr << "Loading binary serialization failed!" << std::endl;
    }

    return success;
  }

  template <typename T>
  IGL_INLINE bool deserialize(T& obj,const std::string& objectName,const std::vector<char>& buffer)
  {
    bool success = false;

    // find suitable object header
    auto objectIter = buffer.cend();
    auto iter = buffer.cbegin();
    while(iter != buffer.end())
    {
      std::string name;
      std::string type;
      size_t size;
      serialization::deserialize(name,iter);
      serialization::deserialize(type,iter);
      serialization::deserialize(size,iter);

      if(name == objectName && type == typeid(obj).name())
      {
        objectIter = iter;
        //break; // find first suitable object header
      }

      iter+=size;
    }

    if(objectIter != buffer.end())
    {
      serialization::deserialize(obj,objectIter);
      success = true;
    }
    else
    {
      obj = T();
    }

    return success;
  }

  // Wrapper function which combines both, de- and serialization
  template <typename T>
  IGL_INLINE bool serializer(bool s,T& obj,std::string& filename)
  {
    return s ? serialize(obj,filename) : deserialize(obj,filename);
  }

  template <typename T>
  IGL_INLINE bool serializer(bool s,T& obj,std::string& objectName,const std::string& filename,bool overwrite)
  {
    return s ? serialize(obj,objectName,filename,overwrite) : deserialize(obj,objectName,filename);
  }

  template <typename T>
  IGL_INLINE bool serializer(bool s,T& obj,std::string& objectName,std::vector<char>& buffer)
  {
    return s ? serialize(obj,objectName,buffer) : deserialize(obj,objectName,buffer);
  }

  IGL_INLINE bool Serializable::PreSerialization() const
  {
    return true;
  }

  IGL_INLINE void Serializable::PostSerialization() const
  {
  }

  IGL_INLINE bool Serializable::PreDeserialization()
  {
    return true;
  }

  IGL_INLINE void Serializable::PostDeserialization()
  {
  }

  IGL_INLINE void Serializable::Serialize(std::vector<char>& buffer) const
  {
    if(this->PreSerialization())
    {
      if(initialized == false)
      {
        objects.clear();
        (const_cast<Serializable*>(this))->InitSerialization();
        initialized = true;
      }

      for(const auto& v : objects)
      {
        v->Serialize(buffer);
      }

      this->PostSerialization();
    }
  }

  IGL_INLINE void Serializable::Deserialize(const std::vector<char>& buffer)
  {
    if(this->PreDeserialization())
    {
      if(initialized == false)
      {
        objects.clear();
        (const_cast<Serializable*>(this))->InitSerialization();
        initialized = true;
      }

      for(auto& v : objects)
      {
        v->Deserialize(buffer);
      }

      this->PostDeserialization();
    }
  }

  IGL_INLINE Serializable::Serializable()
  {
    initialized = false;
  }

  IGL_INLINE Serializable::Serializable(const Serializable& obj)
  {
    initialized = false;
    objects.clear();
  }

  IGL_INLINE Serializable::~Serializable()
  {
    initialized = false;
    objects.clear();
  }

  IGL_INLINE Serializable& Serializable::operator=(const Serializable& obj)
  {
    if(this != &obj)
    {
      if(initialized)
      {
        initialized = false;
        objects.clear();
      }
    }
    return *this;
  }

  template <typename T>
  IGL_INLINE void Serializable::Add(T& obj,std::string name,bool binary)
  {
    auto object = new SerializationObject<T>();
    object->Binary = binary;
    object->Name = name;
    object->Object = std::unique_ptr<T>(&obj);

    objects.push_back(object);
  }

  namespace serialization
  {
    // not serializable
    template <typename T>
    IGL_INLINE typename std::enable_if<!is_serializable<T>::value,size_t>::type getByteSize(const T& obj)
    {
      return sizeof(std::vector<char>::size_type);
    }

    template <typename T>
    IGL_INLINE typename std::enable_if<!is_serializable<T>::value>::type serialize(const T& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter)
    {
      // data
      std::vector<char> tmp;
      serialize(obj,tmp);

      // size
      size_t size = buffer.size();
      serialization::serialize(tmp.size(),buffer,iter);
      size_t cur = iter - buffer.begin();

      buffer.resize(size+tmp.size());
      iter = buffer.begin()+cur;
      iter = std::copy(tmp.begin(),tmp.end(),iter);
    }

    template <typename T>
    IGL_INLINE typename std::enable_if<!is_serializable<T>::value>::type deserialize(T& obj,std::vector<char>::const_iterator& iter)
    {
      std::vector<char>::size_type size;
      serialization::deserialize(size,iter);

      std::vector<char> tmp;
      tmp.resize(size);
      std::copy(iter,iter+size,tmp.begin());

      deserialize(obj,tmp);
      iter += size;
    }

    // fundamental types

    template <typename T>
    IGL_INLINE typename std::enable_if<std::is_fundamental<T>::value,size_t>::type getByteSize(const T& obj)
    {
      return sizeof(T);
    }

    template <typename T>
    IGL_INLINE typename std::enable_if<std::is_fundamental<T>::value>::type serialize(const T& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter)
    {
      //serialization::updateMemoryMap(obj,sizeof(T),memoryMap);
      const uint8_t* ptr = reinterpret_cast<const uint8_t*>(&obj);
      iter = std::copy(ptr,ptr+sizeof(T),iter);
    }

    template <typename T>
    IGL_INLINE typename std::enable_if<std::is_fundamental<T>::value>::type deserialize(T& obj,std::vector<char>::const_iterator& iter)
    {
      uint8_t* ptr = reinterpret_cast<uint8_t*>(&obj);
      std::copy(iter,iter+sizeof(T),ptr);
      iter += sizeof(T);
    }

    // std::string

    IGL_INLINE size_t getByteSize(const std::string& obj)
    {
      return getByteSize(obj.length())+obj.length()*sizeof(uint8_t);
    }

    IGL_INLINE void serialize(const std::string& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter)
    {
      serialization::serialize(obj.length(),buffer,iter);
      for(const auto& cur : obj)
      {
        serialization::serialize(cur,buffer,iter);
      }
    }

    IGL_INLINE void deserialize(std::string& obj,std::vector<char>::const_iterator& iter)
    {
      size_t size;
      serialization::deserialize(size,iter);

      std::string str(size,'\0');
      for(size_t i=0; i<size; ++i)
      {
        serialization::deserialize(str.at(i),iter);
      }

      obj = str;
    }

    // SerializableBase

    template <typename T>
    IGL_INLINE typename std::enable_if<std::is_base_of<SerializableBase,T>::value,size_t>::type getByteSize(const T& obj)
    {
      return sizeof(std::vector<char>::size_type);
    }

    template <typename T>
    IGL_INLINE typename std::enable_if<std::is_base_of<SerializableBase,T>::value>::type serialize(const T& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter)
    {
      // data
      std::vector<char> tmp;
      obj.Serialize(tmp);

      // size
      size_t size = buffer.size();
      serialization::serialize(tmp.size(),buffer,iter);
      size_t cur = iter - buffer.begin();

      buffer.resize(size+tmp.size());
      iter = buffer.begin()+cur;
      iter = std::copy(tmp.begin(),tmp.end(),iter);
    }

    template <typename T>
    IGL_INLINE typename std::enable_if<std::is_base_of<SerializableBase,T>::value>::type deserialize(T& obj,std::vector<char>::const_iterator& iter)
    {
      std::vector<char>::size_type size;
      serialization::deserialize(size,iter);

      std::vector<char> tmp;
      tmp.resize(size);
      std::copy(iter,iter+size,tmp.begin());

      obj.Deserialize(tmp);
      iter += size;
    }

    // STL containers

    // std::pair

    template <typename T1,typename T2>
    IGL_INLINE size_t getByteSize(const std::pair<T1,T2>& obj)
    {
      return getByteSize(obj.first)+getByteSize(obj.second);
    }

    template <typename T1,typename T2>
    IGL_INLINE void serialize(const std::pair<T1,T2>& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter)
    {
      serialization::serialize(obj.first,buffer,iter);
      serialization::serialize(obj.second,buffer,iter);
    }

    template <typename T1,typename T2>
    IGL_INLINE void deserialize(std::pair<T1,T2>& obj,std::vector<char>::const_iterator& iter)
    {
      serialization::deserialize(obj.first,iter);
      serialization::deserialize(obj.second,iter);
    }

    // std::vector

    template <typename T1,typename T2>
    IGL_INLINE size_t getByteSize(const std::vector<T1,T2>& obj)
    {
      return std::accumulate(obj.begin(),obj.end(),sizeof(size_t),[](const size_t& acc,const T1& cur) { return acc+getByteSize(cur); });
    }

    template <typename T1,typename T2>
    IGL_INLINE void serialize(const std::vector<T1,T2>& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter)
    {
      size_t size = obj.size();
      serialization::serialize(size,buffer,iter);
      for(const auto& cur : obj)
      {
        serialization::serialize(cur,buffer,iter);
      }
    }

    template <typename T1,typename T2>
    IGL_INLINE void deserialize(std::vector<T1,T2>& obj,std::vector<char>::const_iterator& iter)
    {
      size_t size;
      serialization::deserialize(size,iter);

      obj.resize(size);
      for(auto& v : obj)
      {
        serialization::deserialize(v,iter);
      }
    }

    //std::set

    template <typename T>
    IGL_INLINE size_t getByteSize(const std::set<T>& obj)
    {
      return std::accumulate(obj.begin(),obj.end(),getByteSize(obj.size()),[](const size_t& acc,const T& cur) { return acc+getByteSize(cur); });
    }

    template <typename T>
    IGL_INLINE void serialize(const std::set<T>& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter)
    {
      serialization::serialize(obj.size(),buffer,iter);
      for(const auto& cur : obj)
      {
        serialization::serialize(cur,buffer,iter);
      }
    }

    template <typename T>
    IGL_INLINE void deserialize(std::set<T>& obj,std::vector<char>::const_iterator& iter)
    {
      size_t size;
      serialization::deserialize(size,iter);

      obj.clear();
      for(size_t i=0; i<size; ++i)
      {
        T val;
        serialization::deserialize(val,iter);
        obj.insert(val);
      }
    }

    // std::map

    template <typename T1,typename T2>
    IGL_INLINE size_t getByteSize(const std::map<T1,T2>& obj)
    {
      return std::accumulate(obj.begin(),obj.end(),sizeof(size_t),[](const size_t& acc,const std::pair<T1,T2>& cur) { return acc+getByteSize(cur); });
    }

    template <typename T1,typename T2>
    IGL_INLINE void serialize(const std::map<T1,T2>& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter)
    {
      serialization::serialize(obj.size(),buffer,iter);
      for(const auto& cur : obj)
      {
        serialization::serialize(cur,buffer,iter);
      }
    }

    template <typename T1,typename T2>
    IGL_INLINE void deserialize(std::map<T1,T2>& obj,std::vector<char>::const_iterator& iter)
    {
      size_t size;
      serialization::deserialize(size,iter);

      obj.clear();
      for(size_t i=0; i<size; ++i)
      {
        std::pair<T1,T2> pair;
        serialization::deserialize(pair,iter);
        obj.insert(pair);
      }
    }

    // Eigen types
    template<typename T,int R,int C,int P,int MR,int MC>
    IGL_INLINE size_t getByteSize(const Eigen::Matrix<T,R,C,P,MR,MC>& obj)
    {
      // space for numbers of rows,cols and data
      return 2*sizeof(typename Eigen::Matrix<T,R,C,P,MR,MC>::Index)+sizeof(T)*obj.rows()*obj.cols();
    }

    template<typename T,int R,int C,int P,int MR,int MC>
    IGL_INLINE void serialize(const Eigen::Matrix<T,R,C,P,MR,MC>& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter)
    {
      serialization::serialize(obj.rows(),buffer,iter);
      serialization::serialize(obj.cols(),buffer,iter);
      size_t size = sizeof(T)*obj.rows()*obj.cols();
      auto ptr = reinterpret_cast<const uint8_t*>(obj.data());
      iter = std::copy(ptr,ptr+size,iter);
    }

    template<typename T,int R,int C,int P,int MR,int MC>
    IGL_INLINE void deserialize(Eigen::Matrix<T,R,C,P,MR,MC>& obj,std::vector<char>::const_iterator& iter)
    {
      typename Eigen::Matrix<T,R,C,P,MR,MC>::Index rows,cols;
      serialization::deserialize(rows,iter);
      serialization::deserialize(cols,iter);
      size_t size = sizeof(T)*rows*cols;
      obj.resize(rows,cols);
      auto ptr = reinterpret_cast<uint8_t*>(obj.data());
      std::copy(iter,iter+size,ptr);
      iter+=size;
    }

    template<typename T,int P,typename I>
    IGL_INLINE size_t getByteSize(const Eigen::SparseMatrix<T,P,I>& obj)
    {
      // space for numbers of rows,cols,nonZeros and tripplets with data (rowIdx,colIdx,value)
      size_t size = sizeof(Eigen::SparseMatrix<T,P,I>::Index);
      return 3*size+(sizeof(T)+2*size)*obj.nonZeros();
    }

    template<typename T,int P,typename I>
    IGL_INLINE void serialize(const Eigen::SparseMatrix<T,P,I>& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter)
    {
      serialization::serialize(obj.rows(),buffer,iter);
      serialization::serialize(obj.cols(),buffer,iter);
      serialization::serialize(obj.nonZeros(),buffer,iter);

      for(int k=0;k<obj.outerSize();++k)
      {
        for(typename Eigen::SparseMatrix<T,P,I>::InnerIterator it(obj,k);it;++it)
        {
          serialization::serialize(it.row(),buffer,iter);
          serialization::serialize(it.col(),buffer,iter);
          serialization::serialize(it.value(),buffer,iter);
        }
      }
    }

    template<typename T,int P,typename I>
    IGL_INLINE void deserialize(Eigen::SparseMatrix<T,P,I>& obj,std::vector<char>::const_iterator& iter)
    {
      typename Eigen::SparseMatrix<T,P,I>::Index rows,cols,nonZeros;
      serialization::deserialize(rows,iter);
      serialization::deserialize(cols,iter);
      serialization::deserialize(nonZeros,iter);

      obj.resize(rows,cols);
      obj.setZero();

      std::vector<Eigen::Triplet<T,I> > triplets;
      for(int i=0;i<nonZeros;i++)
      {
        typename Eigen::SparseMatrix<T,P,I>::Index rowId,colId;
        serialization::deserialize(rowId,iter);
        serialization::deserialize(colId,iter);
        T value;
        serialization::deserialize(value,iter);
        triplets.push_back(Eigen::Triplet<T,I>(rowId,colId,value));
      }
      obj.setFromTriplets(triplets.begin(),triplets.end());
    }

    // pointers

    template <typename T>
    IGL_INLINE typename std::enable_if<std::is_pointer<T>::value,size_t>::type getByteSize(const T& obj)
    {
      size_t size = sizeof(bool);

      if(obj)
        size += getByteSize(*obj);

      return size;
    }

    template <typename T>
    IGL_INLINE typename std::enable_if<std::is_pointer<T>::value>::type serialize(const T& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter)
    {
      serialization::serialize(obj == nullptr,buffer,iter);

      if(obj)
        serialization::serialize(*obj,buffer,iter);
    }

    template <typename T>
    IGL_INLINE typename std::enable_if<std::is_pointer<T>::value>::type deserialize(T& obj,std::vector<char>::const_iterator& iter)
    {
      bool isNullPtr;
      serialization::deserialize(isNullPtr,iter);

      if(isNullPtr)
      {
        if(obj)
        {
          std::cout << "deserialization: possible memory leak for '" << typeid(obj).name() << "'" << std::endl;
          obj = nullptr;
        }
      }
      else
      {
        if(obj)
          std::cout << "deserialization: possible memory leak for '" << typeid(obj).name() << "'" << std::endl;

        obj = new typename std::remove_pointer<T>::type();
        serialization::deserialize(*obj,iter);
      }
    }

    // std::shared_ptr and std::unique_ptr

    /*template <typename T>
    IGL_INLINE typename std::enable_if<serialization::is_smart_ptr<T>::value,size_t>::type getByteSize(const T& obj)
    {
      return getByteSize(obj.get());
    }

    template <typename T>
    IGL_INLINE typename std::enable_if<serialization::is_smart_ptr<T>::value>::type serialize(const T& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter)
    {
      serialize(obj.get(),buffer,iter);
    }

    template <template<typename> class T0,typename T1>
    IGL_INLINE typename std::enable_if<serialization::is_smart_ptr<T0<T1> >::value>::type deserialize(T0<T1>& obj,std::vector<char>::const_iterator& iter)
    {
      bool isNullPtr;
      serialization::deserialize(isNullPtr,iter);

      if(isNullPtr)
      {
        obj.reset();
      }
      else
      {
        obj = T0<T1>(new T1());
        serialization::deserialize(*obj,iter);
      }
    }

    // std::weak_ptr

    template <typename T>
    IGL_INLINE size_t getByteSize(const std::weak_ptr<T>& obj)
    {
      return sizeof(size_t);
    }

    template <typename T>
    IGL_INLINE void serialize(const std::weak_ptr<T>& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter)
    {

    }

    template <typename T>
    IGL_INLINE void deserialize(std::weak_ptr<T>& obj,std::vector<char>::const_iterator& iter)
    {

    }*/

    // functions to overload for non-intrusive serialization
    template <typename T>
    IGL_INLINE void serialize(const T& obj,std::vector<char>& buffer)
    {
      std::cerr << typeid(obj).name() << " is not serializable: derive from igl::Serializable or overload the function igl::serialization::serialize(const T& obj,std::vector<char>& buffer)" << std::endl;
    }

    template <typename T>
    IGL_INLINE void deserialize(T& obj,const std::vector<char>& buffer)
    {
      std::cerr << typeid(obj).name() << " is not drserializable: derive from igl::Serializable or overload the function igl::serialization::deserialize(T& obj, const std::vector<char>& buffer)" << std::endl;
    }

    // helper functions

    template <typename T>
    IGL_INLINE void updateMemoryMap(T& obj,size_t size,std::map<std::uintptr_t,IndexedPointerBase*>& memoryMap)
    {
      // check if object is already serialized
      auto startPtr = new IndexedPointer<T>();
      startPtr->Object = &obj;
      auto startBasePtr = static_cast<IndexedPointerBase*>(startPtr);
      startBasePtr->Type = IndexedPointerBase::BEGIN;
      auto startAddress = reinterpret_cast<std::uintptr_t>(&obj);
      auto p = std::pair<std::uintptr_t,IndexedPointerBase*>(startAddress,startBasePtr);

      auto el = memoryMap.insert(p);
      auto iter = ++el.first; // next elememt
      if(el.second && (iter == memoryMap.end() || iter->second->Type != IndexedPointerBase::END))
      {
        // not yet serialized
        auto endPtr = new IndexedPointer<T>();
        auto endBasePtr = static_cast<IndexedPointerBase*>(endPtr);
        endBasePtr->Type = IndexedPointerBase::END;
        auto endAddress = reinterpret_cast<std::uintptr_t>(&obj) + size - 1;
        auto p = std::pair<std::uintptr_t,IndexedPointerBase*>(endAddress,endBasePtr);

        // insert end address
        memoryMap.insert(el.first,p);
      }
      else
      {
        // already serialized

        // remove inserted address
        memoryMap.erase(el.first);
      }
    }
  }
}
