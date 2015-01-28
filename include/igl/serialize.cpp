//
// Copyright (C) 2014 Christian Sch�ller <schuellchr@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "serialize.h"

namespace igl
{
  template <typename T>
  IGL_INLINE void serialize(const T& obj,const std::string& filename)
  {
    std::ofstream fout(filename.c_str(),std::ios::out | std::ios::binary);

    if(fout.is_open())
    {
      std::vector<char> buffer;
      serialize(obj,"obj",buffer);

      fout.write(&buffer[0],buffer.size());

      fout.close();
    }
    else
    {
      std::cerr << "Saving binary serialization failed!" << std::endl;
    }
  }

  template <typename T>
  IGL_INLINE void serialize(const T& obj,const std::string& objectName,std::vector<char>& buffer)
  {
    static_assert(detail::is_serializable<T>::value,"'igl::serialize': type is not serializable");

    // serialize object data
    size_t size = detail::getByteSize(obj);
    std::vector<char> tmp(size);
    detail::serialize(obj,tmp,tmp.begin());

    std::string objectType(typeid(obj).name());

    size_t curSize = buffer.size();
    size_t objSize = tmp.size();
    size_t newSize = curSize + detail::getByteSize(objectName) + detail::getByteSize(objectType) + sizeof(size_t)+objSize;

    buffer.resize(newSize);

    std::vector<char>::iterator iter = buffer.begin()+curSize;

    // serialize object header (name/type/size)
    detail::serialize(objectName,buffer,iter);
    detail::serialize(objectType,buffer,iter);
    detail::serialize(objSize,buffer,iter);

    // copy serilized data to buffer
    iter = std::copy(tmp.begin(),tmp.end(),iter);
  }

  template <typename T>
  IGL_INLINE void deserialize(T& obj,const std::string& filename)
  {
    std::ifstream file(filename.c_str(),std::ios::binary);

    if(file.is_open())
    {
      file.seekg(0,std::ios::end);
      int size = file.tellg();
      file.seekg(0,std::ios::beg);

      std::vector<char> buffer(size);
      file.read(&buffer[0],size);

      deserialize(obj,"obj",buffer);
      file.close();
    }
    else
    {
      std::cerr << "Loading binary serialization failed!" << std::endl;
    }
  }

  template <typename T>
  IGL_INLINE void deserialize(T& obj,const std::string& objectName,const std::vector<char>& buffer)
  {
    static_assert(detail::is_serializable<T>::value,"'igl::deserialize': type is not deserializable");

    // find suitable object header
    std::vector<char>::const_iterator iter = buffer.begin();
    while(iter != buffer.end())
    {
      std::string name;
      std::string type;
      size_t size;
      detail::deserialize(name,iter);
      detail::deserialize(type,iter);
      detail::deserialize(size,iter);

      if(name == objectName && type == typeid(obj).name())
        break;
      else
        iter+=size;
    }

    if(iter != buffer.end())
      detail::deserialize(obj,iter);
    else
      obj = T();
  }

  namespace detail
  {
    // fundamental types

    template <typename T>
    IGL_INLINE typename std::enable_if<std::is_fundamental<T>::value,size_t>::type getByteSize(const T& obj)
    {
      return sizeof(T);
    }

    template <typename T>
    IGL_INLINE typename std::enable_if<std::is_fundamental<T>::value>::type serialize(const T& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter)
    {
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
      detail::serialize(obj.length(),buffer,iter);
      for(const auto& cur : obj) { detail::serialize(cur,buffer,iter); }
    }

    IGL_INLINE void deserialize(std::string& obj,std::vector<char>::const_iterator& iter)
    {
      size_t size;
      detail::deserialize(size,iter);

      std::string str(size,'\0');
      for(size_t i=0; i<size; ++i)
      {
        detail::deserialize(str.at(i),iter);
      }

      obj = str;
    }

    // Serializable

    template <typename T>
    IGL_INLINE typename std::enable_if<std::is_base_of<Serializable,T>::value,size_t>::type getByteSize(const T& obj)
    {
      return sizeof(std::vector<char>::size_type);
    }

    template <typename T>
    IGL_INLINE typename std::enable_if<std::is_base_of<Serializable,T>::value>::type serialize(const T& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter)
    {
      // data
      std::vector<char> tmp;
      obj.Serialize(tmp);

      // size
      size_t size = buffer.size();
      detail::serialize(tmp.size(),buffer,iter);
      size_t cur = iter - buffer.begin();

      buffer.resize(size+tmp.size());
      iter = buffer.begin()+cur;
      iter = std::copy(tmp.begin(),tmp.end(),iter);
    }

    template <typename T>
    IGL_INLINE typename std::enable_if<std::is_base_of<Serializable,T>::value>::type deserialize(T& obj,std::vector<char>::const_iterator& iter)
    {
      std::vector<char>::size_type size;
      detail::deserialize(size,iter);

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
      detail::serialize(obj.first,buffer,iter);
      detail::serialize(obj.second,buffer,iter);
    }

    template <typename T1,typename T2>
    IGL_INLINE void deserialize(std::pair<T1,T2>& obj,std::vector<char>::const_iterator& iter)
    {
      detail::deserialize(obj.first,iter);
      detail::deserialize(obj.second,iter);
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
      detail::serialize(size,buffer,iter);
      for(const auto& cur : obj)
      {
        detail::serialize(cur,buffer,iter);
      }
    }

    template <typename T1,typename T2>
    IGL_INLINE void deserialize(std::vector<T1,T2>& obj,std::vector<char>::const_iterator& iter)
    {
      size_t size;
      detail::deserialize(size,iter);

      obj.resize(size);
      for(size_t i=0; i<size; ++i)
      {
        detail::deserialize(obj[i],iter);
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
      detail::serialize(obj.size(),buffer,iter);
      for(const auto& cur : obj) { detail::serialize(cur,buffer,iter); }
    }

    template <typename T>
    IGL_INLINE void deserialize(std::set<T>& obj,std::vector<char>::const_iterator& iter)
    {
      size_t size;
      detail::deserialize(size,iter);

      obj.clear();
      for(size_t i=0; i<size; ++i)
      {
        T val;
        detail::deserialize(val,iter);
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
      detail::serialize(obj.size(),buffer,iter);
      for(const auto& cur : obj) { detail::serialize(cur,buffer,iter); }
    }

    template <typename T1,typename T2>
    IGL_INLINE void deserialize(std::map<T1,T2>& obj,std::vector<char>::const_iterator& iter)
    {
      size_t size;
      detail::deserialize(size,iter);

      obj.clear();
      for(size_t i=0; i<size; ++i)
      {
        std::pair<T1,T2> pair;
        detail::deserialize(pair,iter);
        obj.insert(pair);
      }
    }

    // Eigen types
    template<typename T,int R,int C,int P,int MR,int MC>
    IGL_INLINE size_t getByteSize(const Eigen::Matrix<T,R,C,P,MR,MC>& obj)
    {
      // space for numbers of rows,cols and data
      return 2*sizeof(Eigen::Matrix<T,R,C,P,MR,MC>::Index)+sizeof(T)*obj.rows()*obj.cols();
    }

    template<typename T,int R,int C,int P,int MR,int MC>
    IGL_INLINE void serialize(const Eigen::Matrix<T,R,C,P,MR,MC>& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter)
    {
      detail::serialize(obj.rows(),buffer,iter);
      detail::serialize(obj.cols(),buffer,iter);
      size_t size = sizeof(T)*obj.rows()*obj.cols();
      const uint8_t* ptr = reinterpret_cast<const uint8_t*>(obj.data());
      iter = std::copy(ptr,ptr+size,iter);
    }

    template<typename T,int R,int C,int P,int MR,int MC>
    IGL_INLINE void deserialize(Eigen::Matrix<T,R,C,P,MR,MC>& obj,std::vector<char>::const_iterator& iter)
    {
      typename Eigen::Matrix<T,R,C,P,MR,MC>::Index rows,cols;
      detail::deserialize(rows,iter);
      detail::deserialize(cols,iter);
      size_t size = sizeof(T)*rows*cols;
      obj.resize(rows,cols);
      uint8_t* ptr = reinterpret_cast<uint8_t*>(obj.data());
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
      detail::serialize(obj.rows(),buffer,iter);
      detail::serialize(obj.cols(),buffer,iter);
      detail::serialize(obj.nonZeros(),buffer,iter);

      for(int k=0;k<obj.outerSize();++k)
      {
        for(typename Eigen::SparseMatrix<T,P,I>::InnerIterator it(obj,k);it;++it)
        {
          detail::serialize(it.row(),buffer,iter);
          detail::serialize(it.col(),buffer,iter);
          detail::serialize(it.value(),buffer,iter);
        }
      }
    }

    template<typename T,int P,typename I>
    IGL_INLINE void deserialize(Eigen::SparseMatrix<T,P,I>& obj,std::vector<char>::const_iterator& iter)
    {
      typename Eigen::SparseMatrix<T,P,I>::Index rows,cols,nonZeros;
      detail::deserialize(rows,iter);
      detail::deserialize(cols,iter);
      detail::deserialize(nonZeros,iter);

      obj.resize(rows,cols);
      obj.setZero();

      std::vector<Eigen::Triplet<T,I> > triplets;
      for(int i=0;i<nonZeros;i++)
      {
        typename Eigen::SparseMatrix<T,P,I>::Index rowId,colId;
        detail::deserialize(rowId,iter);
        detail::deserialize(colId,iter);
        T value;
        detail::deserialize(value,iter);
        triplets.push_back(Eigen::Triplet<T,I>(rowId,colId,value));
      }
      obj.setFromTriplets(triplets.begin(),triplets.end());
    }

    // pointers

    template <typename T>
    IGL_INLINE typename std::enable_if<std::is_pointer<T>::value,size_t>::type getByteSize(const T& obj)
    {
      bool isNullPtr = (obj == NULL);

      size_t size = sizeof(bool);

      if(isNullPtr == false)
        size += getByteSize(*obj);

      return size;
    }

    template <typename T>
    IGL_INLINE typename std::enable_if<std::is_pointer<T>::value>::type serialize(const T& obj,std::vector<char>& buffer,std::vector<char>::iterator& iter)
    {
      bool isNullPtr = (obj == NULL);

      detail::serialize(isNullPtr,buffer,iter);

      if(isNullPtr == false)
        detail::serialize(*obj,buffer,iter);
    }

    template <typename T>
    IGL_INLINE typename std::enable_if<std::is_pointer<T>::value>::type deserialize(T& obj,std::vector<char>::const_iterator& iter)
    {
      bool isNullPtr;
      detail::deserialize(isNullPtr,iter);

      if(isNullPtr)
      {
        if(obj != NULL)
        {
          std::cout << "deserialization: possible memory leak for '" << typeid(obj).name() << "'" << std::endl;
          obj = NULL;
        }
      }
      else
      {
        if(obj != NULL)
          std::cout << "deserialization: possible memory leak for '" << typeid(obj).name() << "'" << std::endl;

        obj = new typename std::remove_pointer<T>::type();

        detail::deserialize(*obj,iter);
      }
    }
  }
}
