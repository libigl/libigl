// 
// Copyright (C) 2014 Christian Schüller <schuellchr@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
/* -----------------------------------------------------------------------------
 
 Header library which allows to save and load a serialization of basic c++ data
 types like char, char*, std::string, bool, uint, int, float, double to and
 from a xml file.  Containers like std::vector, std::std::pair, Eigen dense and
 sparse matrices are supported as well as combinations of them (like
 vector<pair<string,bool>> or vector<vector<int>>).
 
 To serialize an arbitrary object use the XMLSerializable interface or even
 simpler the XMLSerialization class.
 
 The serialized objects are organised in groups in the xml file and have
 their own names which must be unique within one group.
 
 You can find examples how to use it in the test case class XMLSerializerTest.
 
 -----------------------------------------------------------------------------
TODOs:
* handle memory leak when deserializing to pointers
* loops of pointers and from multiple objects
* NULL pointers
* Versioning
 -----------------------------------------------------------------------------
Bugs:
* Doesn't handle RowMajor Eigen matrices
 ----------------------------------------------------------------------------- */
#ifndef IGL_XML_SERIALIZER_H
#define IGL_XML_SERIALIZER_H

#include <igl/igl_inline.h>
#include <igl/xml/old_version/XMLSerializable.h>

#include <iostream>
//#include <array>
#include <vector>
#include <set>
#include <map>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "tinyxml2.h"

namespace igl
{
  namespace
  {
    // Utility functions
    void EncodeXMLElementName(std::string& name);
    void DecodeXMLElementName(std::string& name);
    void ReplaceSubString(std::string& str, const std::string& search, const std::string& replace);
    
    // Forward declaration
    class XMLSerializer;
    
    /**
     * class XMLSerializableObject
     * internal usage
     */
    class XMLSerializableObject : public ::igl::XMLSerializable
    {
    public:
      
      XMLSerializableObject(const std::string& name, const std::string& group);
      virtual ~XMLSerializableObject();
      
      // set attribute conversion functions
      void SetAttribute(tinyxml2::XMLElement* element, const char* name, char& dest);
      void SetAttribute(tinyxml2::XMLElement* element, const char* name, char*& dest);
      void SetAttribute(tinyxml2::XMLElement* element, const char* name, std::string& dest);
      void SetAttribute(tinyxml2::XMLElement* element, const char* name, bool& dest);
      void SetAttribute(tinyxml2::XMLElement* element, const char* name, unsigned int& dest);
      void SetAttribute(tinyxml2::XMLElement* element, const char* name, int& dest);
      void SetAttribute(tinyxml2::XMLElement* element, const char* name, float& dest);
      void SetAttribute(tinyxml2::XMLElement* element, const char* name, double& dest);
      
      // get attribute conversion functions
      void GetAttribute(const char* src, char& dest);
      void GetAttribute(const char* src, char*& dest);
      void GetAttribute(const char* src, std::string& dest);
      void GetAttribute(const char* src, bool& dest);
      void GetAttribute(const char* src, unsigned int& dest);
      void GetAttribute(const char* src, int& dest);
      void GetAttribute(const char* src, float& dest);
      void GetAttribute(const char* src, double& dest);
      
      // Initialization
      
      // Basic data types
      using XMLSerializable::Init;
      void Init(char& val);
      void Init(char*& val);
      void Init(std::string& val);
      void Init(bool& val);
      void Init(unsigned int& val);
      void Init(int& val);
      void Init(float& val);
      void Init(double& val);
      
      // XMLSerializable*
      template<typename T>
      void Init(T& obj);
      template<typename T>
      void Init(T*& obj);
      
      // STL containers
      /*template<typename T, int S>
       void Init(std::array<T,S>& obj);*/
      template<typename T0, typename T1>
      void Init(std::pair<T0,T1>& obj);
      template<typename T>
      void Init(std::vector<T>& obj);
      template<typename T>
      void Init(std::set<T>& obj);
      template<typename T0, typename T1>
      void Init(std::map<T0,T1>& obj);
      
      // Eigen types
      template<typename T, int R, int C>
      void Init(Eigen::Matrix<T,R,C>& obj);
      template<typename T>
      void Init(Eigen::SparseMatrix<T>& obj);
      
      // Serialization
      
      // Basic data types
      using XMLSerializable::Serialize;
      bool Serialize(char& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      bool Serialize(char*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      bool Serialize(std::string& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      bool Serialize(std::string*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      bool Serialize(bool obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      bool Serialize(bool*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      bool Serialize(unsigned int& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      bool Serialize(unsigned int*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      bool Serialize(int& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      bool Serialize(int*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      bool Serialize(float& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      bool Serialize(float*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      bool Serialize(double& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      bool Serialize(double*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      
      // XMLSerializable*
      template<typename T>
      bool Serialize(T& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      template<typename T>
      bool Serialize(T*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      
      // STL containers
      
      /*template<typename T, size_t S>
       bool Serialize(std::array<T,S>& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
       template<typename T, size_t S>
       bool Serialize(std::array<T,S>*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);*/
      
      template<typename T0, typename T1>
      bool Serialize(std::pair<T0,T1>& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      template<typename T0, typename T1>
      bool Serialize(std::pair<T0,T1>*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      
      template<typename T>
      bool Serialize(std::vector<T>& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      template<typename T>
      bool Serialize(std::vector<T>*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      bool Serialize(std::vector<bool>& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      bool Serialize(std::vector<bool>*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      bool Serialize(std::vector<int>& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      bool Serialize(std::vector<int>*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);

      template<typename T>
      bool Serialize(std::set<T>& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      template<typename T>
      bool Serialize(std::set<T>*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);

      template<typename T0, typename T1>
      bool Serialize(std::map<T0,T1>& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      template<typename T0, typename T1>
      bool Serialize(std::map<T0,T1>*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      bool Serialize(std::map<int,int>& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      bool Serialize(std::map<int,int>*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      
      // Eigen types
      template<typename T, int R, int C>
      bool Serialize(Eigen::Matrix<T,R,C>& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      template<typename T, int R, int C>
      bool Serialize(Eigen::Matrix<T,R,C>*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      
      template<typename T>
      bool Serialize(Eigen::SparseMatrix<T>& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      template<typename T>
      bool Serialize(Eigen::SparseMatrix<T>*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      
      // Deserialization
      
      // Basic data types
      using XMLSerializable::Deserialize;
      bool Deserialize(char& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      bool Deserialize(char*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      bool Deserialize(std::string& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      bool Deserialize(std::string*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      bool Deserialize(bool& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      bool Deserialize(bool*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      bool Deserialize(unsigned int& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      bool Deserialize(unsigned int*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      bool Deserialize(int& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      bool Deserialize(int*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      bool Deserialize(float& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      bool Deserialize(float*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      bool Deserialize(double& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      bool Deserialize(double*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      
      // XMLSerializable*
      template<typename T>
      bool Deserialize(T& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      template<typename T>
      bool Deserialize(T*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      
      // STL containers
      
      /*template<typename T, size_t S>
       bool Deserialize(std::array<T,S>& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
       template<typename T, size_t S>
       bool Deserialize(std::array<T,S>*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);*/
      
      template<typename T0, typename T1>
      bool Deserialize(std::pair<T0,T1>& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      template<typename T0, typename T1>
      bool Deserialize(std::pair<T0,T1>*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      
      template<typename T>
      bool Deserialize(std::vector<T>& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      template<typename T>
      bool Deserialize(std::vector<T>*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      bool Deserialize(std::vector<bool>& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      bool Deserialize(std::vector<bool>*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      bool Deserialize(std::vector<int>& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      bool Deserialize(std::vector<int>*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);

      template<typename T>
      bool Deserialize(std::set<T>& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      template<typename T>
      bool Deserialize(std::set<T>*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);

      template<typename T0, typename T1>
      bool Deserialize(std::map<T0,T1>& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      template<typename T0, typename T1>
      bool Deserialize(std::map<T0,T1>*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      bool Deserialize(std::map<int,int>& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      bool Deserialize(std::map<int,int>*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      
      // Eigen types
      template<typename T, int R, int C>
      bool Deserialize(Eigen::Matrix<T,R,C>& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      template<typename T, int R, int C>
      bool Deserialize(Eigen::Matrix<T,R,C>*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      
      template<typename T>
      bool Deserialize(Eigen::SparseMatrix<T>& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      template<typename T>
      bool Deserialize(Eigen::SparseMatrix<T>*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
      
    private:
      
      template<typename T>
      bool setElementAttribute(T& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      template<typename T>
      bool getElementAttribute(T& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
    };
    
    /**
     * class XMLSerializableInstance
     * internal usage
     */
    template<typename T>
    class XMLSerializableInstance : public XMLSerializableObject
    {
    public:
      
      T& Object;
      T DefaultValue;
      
      XMLSerializableInstance(T& obj, const std::string& name, const std::string group);
      XMLSerializableInstance(T& obj, const std::string& name, const std::string group, T defaultValue);
      ~XMLSerializableInstance();
      
      // XMLSerializable interface implementation
      void Init();
      bool Serialize(tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element);
      bool Deserialize(tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element);
    };
    
    /**
     * struct XMLSerializerGroup
     * internal usage
     */
    struct XMLSerializerGroup
    {
      std::string Name;
      std::vector<XMLSerializable*>* Objects;
    };
    
    /**
     * class XMLSerializer
     * This is the core class which takes care of saving and loading of serialization of object structures.
     */
    class XMLSerializer
    {
    public:
      
      /**
       * Serializes an object to a file
       */
      template<typename T>
      static bool SaveObject(T& object, const char* filename);
      template<typename T>
      static bool SaveObject(T& object, const std::string& objectName, const std::string& groupName, const char* filename, bool overwrite);
      
      /**
       * Loads the serialization of an object from a file.
       */
      template<typename T>
      static bool LoadObject(T& object, const char* filename);
      template<typename T>
      static bool LoadObject(T& object, const std::string& objectName, const std::string& groupName, const char* filename);
      
      /**
       * Constructor which sets the default group
       */
      XMLSerializer(const std::string& defaultGroup);
      ~XMLSerializer();
      
      /**
       * Save the serialization of all groups to file.
       * Parameter overwrite specifies if file gets overwritten or updated
       */
      bool Save(const char* filename, bool overwrite);
      bool Save(const std::string& groupName, const char* filename, bool overwrite);
      bool Save(const std::string& objectName, const std::string& groupName, const char* filename, bool overwrite);
      
      /**
       * Save the serialization of all groups to a XMLDocument instance.
       */
      bool SaveToXMLDoc(tinyxml2::XMLDocument* doc);
      bool SaveToXMLDoc(const std::string& groupName, tinyxml2::XMLDocument* doc);
      bool SaveToXMLDoc(const std::string& objectName, const std::string& groupName, tinyxml2::XMLDocument* doc);
      
      /**
       * Save the serialization of a group with a new provided name to given XMLElement instance.
       */
      bool SaveGroupToXMLElement(tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      bool SaveGroupToXMLElement(const std::string& groupName, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
      
      /**
       * Load the serialization from a file.
       */
      bool Load(const char* filename);
      bool Load(const std::string& groupName, const char* filename);
      bool Load(const std::string& objectName, const std::string& groupName, const char* filename);
      
      /**
       * Load the serialization from an XMLDocument instance.
       */
      bool LoadFromXMLDoc(tinyxml2::XMLDocument* doc);
      bool LoadFromXMLDoc(const std::string& groupName, tinyxml2::XMLDocument* doc);
      bool LoadFromXMLDoc(const std::string& objectName, const std::string& groupName, tinyxml2::XMLDocument* doc);
      
      /**
       * Load the serialization from a XMLElement instance to given group.
       */
      bool LoadGroupFromXMLElement(const std::string& groupName, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element);
      
      /**
       * Set/Get current group. Every object which is added afterwards will be in this group, except it specifies another group.
       */
      void SetCurrentGroup(const std::string& group);
      std::string GetCurrentGroup();
      
      /**
       * Add an object to the serializer. Can be simple types like char, char*, string, unsigned int, int, float, double or containers like std::array, std::pair, std::vector.
       * Also Eigen dense or sparse matrices are supported and all objects of type Serializable* and combinations of thoses types like vector<vector>, vector<pair> or even vector<pair<vector,Serializable*>>>.
       * Also pointers to those objects can be used (for instance like vector<vector<pair<int,float>*>*>).
       * char* is also possible as base type and represents a array of chars, but be carefull that the pointer is not just a copy but a valid instance in the current programm scope.
       */
      
      // Basic types
      bool Add(char& obj, const std::string& name);
      bool Add(char*& obj, const std::string& name);
      bool Add(std::string& obj, const std::string& name);
      bool Add(bool& obj, const std::string& name);
      bool Add(unsigned int& obj, const std::string& name);
      bool Add(int& obj, const std::string& name);
      bool Add(float& obj, const std::string& name);
      bool Add(double& obj, const std::string& name);
      
      bool Add(char& obj, const std::string& name, char defaultValue);
      bool Add(char*& obj, const std::string& name, char* defaultValue);
      bool Add(std::string& obj, const std::string& name, std::string defaultValue);
      bool Add(bool& obj, const std::string& name, bool defaultValue);
      bool Add(unsigned int& obj, const std::string& name, unsigned int defaultValue);
      bool Add(int& obj, const std::string& name, int defaultValue);
      bool Add(float& obj, const std::string& name, float defaultValue);
      bool Add(double& obj, const std::string& name, double defaultValue);
      
      // XMLSerializable*
      template<typename T>
      bool Add(T& object, const std::string& name);
      template<typename T>
      bool Add(T& object, const std::string& name, T defaultValue);
      
      // STL containers
      /*template<typename T, size_t S>
       bool Add(std::array<T,S>& obj, const std::string& name);*/
      template<typename T0, typename T1>
      bool Add(std::pair<T0,T1>& obj, const std::string& name);
      template<typename T>
      bool Add(std::vector<T>& obj, const std::string& name);
      template<typename T>
      bool Add(std::set<T>& obj, const std::string& name);
      template<typename T0, typename T1>
      bool Add(std::map<T0,T1>& obj, const std::string& name);
      
      // Eigen types
      template<typename T, int R, int C>
      bool Add(Eigen::Matrix<T,R,C>& obj, const std::string& name);
      template<typename T>
      bool Add(Eigen::SparseMatrix<T>& obj, const std::string& name);
      
    private:
      
      std::map<std::string,XMLSerializerGroup*>::iterator currentGroup;
      std::map<std::string,XMLSerializerGroup*> groups;
      
      template<typename T>
      bool add(T& object, const std::string& name);
      template<typename T>
      bool add(T& object, const std::string& name, T defaultValue);
      bool addObjectToGroup(XMLSerializable* object, const std::string& group);
      bool addObjectToGroup(XMLSerializable* object, std::map<std::string,XMLSerializerGroup*>::iterator it);
      std::map<std::string,XMLSerializerGroup*>::iterator setGetGroup(const std::string& group);
      tinyxml2::XMLDocument* openDoc(const char* filename);
      tinyxml2::XMLElement* findAddGroup(tinyxml2::XMLDocument* doc, const char* groupName);
    };
    
    int numForbiddenChars = 8;
    char forbiddenChars[] = {' ','/','~','#','&','>','<','='};
    
    void ReplaceSubString(std::string& str, const std::string& search, const std::string& replace)
    {
      size_t pos = 0;
      while ((pos = str.find(search, pos)) != std::string::npos)
      {
        str.replace(pos, search.length(), replace);
        pos += replace.length();
      }
    }
    
    void EncodeXMLElementName(std::string& name)
    {
      // must not start with a digit
      if(isdigit(*name.begin()))
      {
        name = ":::" + name;
      }
      
      std::stringstream stream;
      for(int i=0;i<numForbiddenChars;i++)
      {
        std::string search;
        search = forbiddenChars[i];
        std::stringstream replaces;
        replaces << ":" << (int)forbiddenChars[i];
        std::string replace = replaces.str();
        
        ReplaceSubString(name,search,replace);
      }
    }
    
    void DecodeXMLElementName(std::string& name)
    {
      if(name.find("::", 0) == 0)
        name.replace(0,3,"");

      for(auto chr : forbiddenChars)
      {
        std::stringstream searchs;
        searchs << ":" << (int)chr;
        std::string search = searchs.str();
        std::string replace;
        replace = chr;
        
        ReplaceSubString(name,search,replace);
      }
    }
    
    XMLSerializableObject::XMLSerializableObject(const std::string& name, const std::string& group)
    {
      std::string groupName = group;
      std::string objectName = name;
      
      EncodeXMLElementName(groupName);
      EncodeXMLElementName(objectName);
      
      Name = objectName;
    }
    
    XMLSerializableObject::~XMLSerializableObject()
    {
    }
    
    // set attribute conversion functions
    void XMLSerializableObject::SetAttribute(tinyxml2::XMLElement* element, const char* name, char& dest)
    {
      element->SetAttribute(name,dest);
    }
    
    void XMLSerializableObject::SetAttribute(tinyxml2::XMLElement* element, const char* name, char*& dest)
    {
      element->SetAttribute(name,const_cast<const char*>(dest));
    }
    
    void XMLSerializableObject::SetAttribute(tinyxml2::XMLElement* element, const char* name, std::string& dest)
    {
      element->SetAttribute(name,dest.c_str());
    }
    
    void XMLSerializableObject::SetAttribute(tinyxml2::XMLElement* element, const char* name, bool& dest)
    {
      element->SetAttribute(name,dest);
    }
    
    void XMLSerializableObject::SetAttribute(tinyxml2::XMLElement* element, const char* name, unsigned int& dest)
    {
      element->SetAttribute(name,dest);
    }
    
    void XMLSerializableObject::SetAttribute(tinyxml2::XMLElement* element, const char* name, int& dest)
    {
      element->SetAttribute(name,dest);
    }
    
    void XMLSerializableObject::SetAttribute(tinyxml2::XMLElement* element, const char* name, float& dest)
    {
      element->SetAttribute(name,dest);
    }
    
    void XMLSerializableObject::SetAttribute(tinyxml2::XMLElement* element, const char* name, double& dest)
    {
      element->SetAttribute(name,dest);
    }
    
    // get attribute conversion functions
    void XMLSerializableObject::GetAttribute(const char* src, char& dest)
    {
      dest = (char)atoi(src);
    }
    
    void XMLSerializableObject::GetAttribute(const char* src, char*& dest)
    {
      unsigned int length = strlen(src)+1;
      dest = new char[length];
      strcpy(dest, src);
    }
    
    void XMLSerializableObject::GetAttribute(const char* src, std::string& dest)
    {
      dest = src;
    }
    
    void XMLSerializableObject::GetAttribute(const char* src, bool& dest)
    {
      tinyxml2::XMLUtil::ToBool(src,&dest);
    }
    
    void XMLSerializableObject::GetAttribute(const char* src, unsigned int& dest)
    {
      tinyxml2::XMLUtil::ToUnsigned(src,&dest);
    }
    
    void XMLSerializableObject::GetAttribute(const char* src, int& dest)
    {
      tinyxml2::XMLUtil::ToInt(src,&dest);
    }
    
    void XMLSerializableObject::GetAttribute(const char* src, float& dest)
    {
      tinyxml2::XMLUtil::ToFloat(src,&dest);
    }
    
    void XMLSerializableObject::GetAttribute(const char* src, double& dest)
    {
      tinyxml2::XMLUtil::ToDouble(src,&dest);
    }
    
    // specify default value of types
    void XMLSerializableObject::Init(char& val)
    {
      val = '0';
    }
    
    void XMLSerializableObject::Init(char*& val)
    {
      val = NULL;
    }
    
    void XMLSerializableObject::Init(std::string& val)
    {
      val = "";
    }
    
    void XMLSerializableObject::Init(bool& val)
    {
      val = false;
    }
    
    void XMLSerializableObject::Init(unsigned int& val)
    {
      val = 0;
    }
    
    void XMLSerializableObject::Init(int& val)
    {
      val = 0;
    }
    
    void XMLSerializableObject::Init(float& val)
    {
      val = 0.0f;
    }
    
    void XMLSerializableObject::Init(double& val)
    {
      val = 0.000000000000000;
    }

    template<typename T>
    void XMLSerializableObject::Init(T& obj)
    {
      XMLSerializable* object = static_cast<XMLSerializable*>(&obj);
      object->Init();
    }
    
    template<typename T>
    void XMLSerializableObject::Init(T*& obj)
    {
      XMLSerializable* object = static_cast<XMLSerializable*>(obj);
      object->Init();
    }
    
    /*template<typename T, int S>
     void XMLSerializableObject::Init(std::array<T,S>& obj)
     {
     for(unsigned int i=0;i<obj.size();i++)
     Init(obj[i]);
     }*/
    
    template<typename T0, typename T1>
    void XMLSerializableObject::Init(std::pair<T0,T1>& obj)
    {
      Init(obj.first);
      Init(obj.second);
    }
    
    template<typename T>
    void XMLSerializableObject::Init(std::vector<T>& obj)
    {
      obj.clear();
    }

    template<typename T>
    void XMLSerializableObject::Init(std::set<T>& obj)
    {
      obj.clear();
    }

    template<typename T0, typename T1>
    void XMLSerializableObject::Init(std::map<T0,T1>& obj)
    {
      obj.clear();
    }
    
    template<typename T, int R, int C>
    void XMLSerializableObject::Init(Eigen::Matrix<T,R,C>& obj)
    {
      obj.setZero(obj.rows(),obj.cols());
    }
    
    template<typename T>
    void  XMLSerializableObject::Init(Eigen::SparseMatrix<T>& obj)
    {
      obj.setZero();
    }
    
    bool XMLSerializableObject::Serialize(char& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      return setElementAttribute(obj,doc,element,name);
    }
    
    // overload function for char*, it interpreted as char array and can be used to handle strings
    bool XMLSerializableObject::Serialize(char*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      return setElementAttribute(obj,doc,element,name);
    }
    
    bool XMLSerializableObject::Serialize(std::string& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      return setElementAttribute(obj,doc,element,name);
    }
    
    bool XMLSerializableObject::Serialize(std::string*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      return Serialize(*obj,doc,element,name);
    }
    
    bool XMLSerializableObject::Serialize(bool obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      return setElementAttribute(obj,doc,element,name);
    }
    
    bool XMLSerializableObject::Serialize(bool*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      return Serialize(*obj,doc,element,name);
    }
    
    bool XMLSerializableObject::Serialize(unsigned int& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      return setElementAttribute(obj,doc,element,name);
    }
    
    bool XMLSerializableObject::Serialize(unsigned int*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      return Serialize(*obj,doc,element,name);
    }
    
    bool XMLSerializableObject::Serialize(int& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      return setElementAttribute(obj,doc,element,name);
    }
    
    bool XMLSerializableObject::Serialize(int*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      return Serialize(*obj,doc,element,name);
    }
    
    bool XMLSerializableObject::Serialize(float& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      return setElementAttribute(obj,doc,element,name);
    }
    
    bool XMLSerializableObject::Serialize(float*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      return Serialize(*obj,doc,element,name);
    }
    
    bool XMLSerializableObject::Serialize(double& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      return setElementAttribute(obj,doc,element,name);
    }
    
    bool XMLSerializableObject::Serialize(double*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      return Serialize(*obj,doc,element,name);
    }
    
    template<typename T>
    bool XMLSerializableObject::Serialize(T& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      // Serialize object implementing XMLSerializable interface
      XMLSerializable* object = static_cast<XMLSerializable*>(&obj);
      
      tinyxml2::XMLElement* child = doc->NewElement(name.c_str());
      element->InsertEndChild(child);
      
      return object->Serialize(doc,child);
    }
    
    template<typename T>
    bool XMLSerializableObject::Serialize(T*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      // Serialize object implementing XMLSerializable interface
      XMLSerializable* object = static_cast<XMLSerializable*>(obj);
      
      tinyxml2::XMLElement* child = doc->NewElement(name.c_str());
      element->InsertEndChild(child);
      
      if(object != NULL)
        return object->Serialize(doc,child);
      
      return false;
    }
    
    bool XMLSerializableObject::Deserialize(char& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      return getElementAttribute(obj,doc,element,name);
    }
    
    // template specialisation for char*, it interpreted as char array and can be used to handle strings
    bool XMLSerializableObject::Deserialize(char*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      return getElementAttribute(obj,doc,element,name);
    }
    
    bool XMLSerializableObject::Deserialize(std::string& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      return getElementAttribute(obj,doc,element,name);
    }
    
    bool XMLSerializableObject::Deserialize(std::string*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      return Deserialize(*obj,doc,element,name);
    }
    
    bool XMLSerializableObject::Deserialize(bool& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      return getElementAttribute(obj,doc,element,name);
    }
    
    bool XMLSerializableObject::Deserialize(bool*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      return Deserialize(*obj,doc,element,name);
    }
    
    bool XMLSerializableObject::Deserialize(unsigned int& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      return getElementAttribute(obj,doc,element,name);
    }
    
    bool XMLSerializableObject::Deserialize(unsigned int*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      return Deserialize(*obj,doc,element,name);
    }
    
    bool XMLSerializableObject::Deserialize(int& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      return getElementAttribute(obj,doc,element,name);
    }
    
    bool XMLSerializableObject::Deserialize(int*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      return Deserialize(*obj,doc,element,name);
    }
    
    bool XMLSerializableObject::Deserialize(float& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      return getElementAttribute(obj,doc,element,name);
    }
    
    bool XMLSerializableObject::Deserialize(float*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      return Deserialize(*obj,doc,element,name);
    }
    
    bool XMLSerializableObject::Deserialize(double& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      return getElementAttribute(obj,doc,element,name);
    }
    
    bool XMLSerializableObject::Deserialize(double*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      return Deserialize(*obj,doc,element,name);
    }
    
    template<typename T>
    bool XMLSerializableObject::Deserialize(T& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      obj = T();
      XMLSerializable* object = static_cast<XMLSerializable*>(&obj);
      
      const tinyxml2::XMLElement* child = element->FirstChildElement(name.c_str());
      
      if(child != NULL)
      {
        object->Name = child->FirstChild()->Value();
        obj.Deserialize(doc,child);
      }
      else
      {
        obj.Init();
        return false;
      }
      
      return true;
    }
    
    template<typename T>
    bool XMLSerializableObject::Deserialize(T*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      obj = NULL;
      const tinyxml2::XMLElement* child = element->FirstChildElement(name.c_str());
      
      if(child != NULL)
      {
        const tinyxml2::XMLNode* node = child->FirstChild();

        if(node != NULL)
        {
          obj = new T();
          static_cast<XMLSerializable*>(obj)->Name = node->Value();
          obj->Deserialize(doc,child);
          
          return true;
        }
      }
      else
      {
        obj = new T();
        obj->Init();
      }
      
      return false;
    }
    
    /*
     template<typename T, size_t S>
     bool XMLSerializableObject::Serialize(std::array<T,S>& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
     {
     tinyxml2::XMLElement* ar = doc->NewElement(name.c_str());
     element->InsertEndChild(ar);
     
     ar->SetAttribute("size",(unsigned int)obj.size());
     
     std::stringstream num;
     for(unsigned int i=0;i<obj.size();i++)
     {
     num.str("");
     num << "value" << i;
     Serialize(obj[i],doc,ar,num.str());
     }
     
     return true;
     }
     
     template<typename T, size_t S>
     bool XMLSerializableObject::Serialize(std::array<T,S>*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
     {
     return Serialize(*obj,doc,element,name);
     }
     
     template<typename T, size_t S>
     bool XMLSerializableObject::Deserialize(std::array<T,S>& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
     {
     bool res = true;
     
     const tinyxml2::XMLElement* child = element->FirstChildElement(name.c_str());
     if(child != NULL)
     {
     int size = child->UnsignedAttribute("size");
     size = S < size ? S : size;
     
     std::stringstream num;
     const tinyxml2::XMLAttribute* attribute = NULL;
     for(unsigned int i=0;i<size;i++)
     {
     num.str("");
     num << "value" << i;
     
     res &= Deserialize(obj[i],doc,child,num.str());
     }
     }
     else
     return false;
     
     return res;
     }
     
     template<typename T, size_t S>
     bool XMLSerializableObject::Deserialize(std::array<T,S>*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
     {
     obj = new std::array<T,S>();
     return Deserialize(*obj,doc,element,name);
     }
     */
    template<typename T0, typename T1>
    bool XMLSerializableObject::Serialize(std::pair<T0,T1>& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      bool res = true;
      
      tinyxml2::XMLElement* pair = doc->NewElement(name.c_str());
      element->InsertEndChild(pair);
      
      res &= Serialize(obj.first,doc,pair,"first");
      res &= Serialize(obj.second,doc,pair,"second");
      
      return res;
    }
    
    template<typename T0, typename T1>
    bool XMLSerializableObject::Serialize(std::pair<T0,T1>*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      return Serialize(*obj,doc,element,name);
    }
    
    template<typename T0, typename T1>
    bool XMLSerializableObject::Deserialize(std::pair<T0,T1>& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      bool res = true;
      
      const tinyxml2::XMLElement* child = element->FirstChildElement(name.c_str());
      if(child != NULL)
      {
        res &= Deserialize(obj.first,doc,child,"first");
        res &= Deserialize(obj.second,doc,child,"second");
      }
      else
      {
        Init(obj.first);
        Init(obj.second);
        return false;
      }
      
      return res;
    }
    
    template<typename T0, typename T1>
    bool XMLSerializableObject::Deserialize(std::pair<T0,T1>*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      obj = new std::pair<T0,T1>();
      return Deserialize(*obj,doc,element,name);
    }
    
    template<typename T>
    bool XMLSerializableObject::Serialize(std::vector<T>& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      tinyxml2::XMLElement* vector = doc->NewElement(name.c_str());
      element->InsertEndChild(vector);
      
      vector->SetAttribute("size",(unsigned int)obj.size());
      
      std::stringstream num;
      for(unsigned int i=0;i<obj.size();i++)
      {
        num.str("");
        num << "value" << i;
        Serialize(obj[i],doc,vector,num.str());
      }
      
      return true;
    }

    template<typename T>
    bool XMLSerializableObject::Serialize(std::vector<T>*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      return Serialize(*obj,doc,element,name);
    }

    bool XMLSerializableObject::Serialize(std::vector<bool>& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      tinyxml2::XMLElement* vector = doc->NewElement(name.c_str());
      element->InsertEndChild(vector);
      
      const unsigned int size = obj.size();
      
      vector->SetAttribute("size",size);
      
      std::stringstream ms;
      ms << "\n";
      for(unsigned int i=0;i<size;i++)
        ms << obj[i] << ",";
      ms << "\n";
      
      std::string mString = ms.str();
      if(mString.size() > 1)
        mString[mString.size()-2] = '\0';
      
      vector->SetAttribute("vector_bool",mString.c_str());
      
      return true;
    }

    bool XMLSerializableObject::Serialize(std::vector<bool>*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      return Serialize(*obj,doc,element,name);
    }

    bool XMLSerializableObject::Serialize(std::vector<int>& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      tinyxml2::XMLElement* vector = doc->NewElement(name.c_str());
      element->InsertEndChild(vector);
      
      const unsigned int size = obj.size();
      
      vector->SetAttribute("size",size);
      
      std::stringstream ms;
      ms << "\n";
      for(unsigned int i=0;i<size;i++)
        ms << obj[i] << ",";
      ms << "\n";
      
      std::string mString = ms.str();
      if(mString.size() > 1)
        mString[mString.size()-2] = '\0';
      
      vector->SetAttribute("vector_int",mString.c_str());
      
      return true;
    }

    bool XMLSerializableObject::Serialize(std::vector<int>*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      return Serialize(*obj,doc,element,name);
    }
    
    template<typename T>
    bool XMLSerializableObject::Deserialize(std::vector<T>& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      bool res = true;
      obj.clear();
      
      const tinyxml2::XMLElement* child = element->FirstChildElement(name.c_str());
      if(child != NULL)
      {
        unsigned int size = child->UnsignedAttribute("size");
        obj.resize(size);
        
        std::stringstream num;
        for(unsigned int i=0;i<size;i++)
        {
          num.str("");
          num << "value" << i;
          
          res &= Deserialize(obj[i],doc,child,num.str());
        }
      }
      else
        return false;
      
      return res;
    }

    template<typename T>
    bool XMLSerializableObject::Deserialize(std::vector<T>*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      obj = new std::vector<T>();
      return Deserialize(*obj,doc,element,name);
    }

    bool XMLSerializableObject::Deserialize(std::vector<bool>& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      const tinyxml2::XMLElement* child = element->FirstChildElement(name.c_str());
      if(child != NULL)
      {
        const unsigned int size = child->UnsignedAttribute("size");
        
        obj.resize(size);
        
        const tinyxml2::XMLAttribute* attribute = child->FindAttribute("vector_bool");
        if(attribute == NULL)
        {
          Init(obj);
          return false;
        }
        
        char* matTemp;
        GetAttribute(attribute->Value(),matTemp);
        
        std::string line, ssize;
        std::stringstream mats;
        mats.str(matTemp);
        
        int r=0;
        std::string val;
        // for each line
        getline(mats,line); // matrix starts with an empty line
        while(getline(mats,line))
        {
          // get current line
          std::stringstream liness(line);
          
          for(unsigned int c=0;c<size-1;c++)
          {
            // split line
            getline(liness, val, ',');
            
            // push pack the data if any
            if(!val.empty())
            {
              bool t;
              GetAttribute(val.c_str(),t);
              obj[c] = t;
            }
          }
          
          getline(liness, val);
          bool t;
          GetAttribute(val.c_str(),t);
          obj[size-1] = t;
          
          r++;
        }
      }
      else
      {
        Init(obj);
        return false;
      }
      
      return true;
    }

    bool XMLSerializableObject::Deserialize(std::vector<bool>*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      obj = new std::vector<bool>();
      return Deserialize(*obj,doc,element,name);
    }

    bool XMLSerializableObject::Deserialize(std::vector<int>& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      const tinyxml2::XMLElement* child = element->FirstChildElement(name.c_str());
      if(child != NULL)
      {
        const unsigned int size = child->UnsignedAttribute("size");
        
        obj.resize(size);
        
        const tinyxml2::XMLAttribute* attribute = child->FindAttribute("vector_int");
        if(attribute == NULL)
        {
          Init(obj);
          return false;
        }
        
        char* matTemp;
        GetAttribute(attribute->Value(),matTemp);
        
        std::string line, ssize;
        std::stringstream mats;
        mats.str(matTemp);
        
        int r=0;
        std::string val;
        // for each line
        getline(mats,line); // matrix starts with an empty line
        while(getline(mats,line))
        {
          // get current line
          std::stringstream liness(line);
          
          for(unsigned int c=0;c<size-1;c++)
          {
            // split line
            getline(liness, val, ',');
            
            // push pack the data if any
            if(!val.empty())
            {
              int t;
              GetAttribute(val.c_str(),t);
              obj[c] = t;
            }
          }
          
          getline(liness, val);
          int t;
          GetAttribute(val.c_str(),t);
          obj[size-1] = t;
          
          r++;
        }
      }
      else
      {
        Init(obj);
        return false;
      }
      
      return true;
    }

    bool XMLSerializableObject::Deserialize(std::vector<int>*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      obj = new std::vector<int>();
      return Deserialize(*obj,doc,element,name);
    }

    template<typename T>
    bool XMLSerializableObject::Serialize(std::set<T>& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      tinyxml2::XMLElement* set = doc->NewElement(name.c_str());
      element->InsertEndChild(set);
      
      set->SetAttribute("size",(unsigned int)obj.size());
      
      std::stringstream num;
      typename std::set<T>::iterator iter = obj.begin();
      for(int i=0;iter!=obj.end();iter++,i++)
      {
        num.str("");
        num << "value" << i;
        Serialize((T)*iter,doc,set,num.str());
      }
      
      return true;
    }

    template<typename T>
    bool XMLSerializableObject::Serialize(std::set<T>*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      return Serialize(*obj,doc,element,name);
    }

    template<typename T>
    bool XMLSerializableObject::Deserialize(std::set<T>& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      bool res = true;
      obj.clear();
      
      const tinyxml2::XMLElement* child = element->FirstChildElement(name.c_str());
      if(child != NULL)
      {
        unsigned int size = child->UnsignedAttribute("size");
        
        std::stringstream num;
        typename std::set<T>::iterator iter = obj.begin();
        for(int i=0;i<size;i++)
        {
          num.str("");
          num << "value" << i;
          
          T val;
          res &= Deserialize(val,doc,child,num.str());
          obj.insert(val);
        }
      }
      else
        return false;
      
      return res;
    }

    template<typename T>
    bool XMLSerializableObject::Deserialize(std::set<T>*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      obj = new std::set<T>();
      return Deserialize(*obj,doc,element,name);
    }

    template<typename T0, typename T1>
    bool XMLSerializableObject::Serialize(std::map<T0,T1>& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      tinyxml2::XMLElement* map = doc->NewElement(name.c_str());
      element->InsertEndChild(map);
      
      map->SetAttribute("size",(unsigned int)obj.size());
      
      std::stringstream num;
      typename std::map<T0,T1>::iterator iter = obj.begin();
      for(int i=0;iter!=obj.end();iter++,i++)
      {
        num.str("");
        num << "value" << i;
        Serialize((std::pair<T0,T1>)*iter,doc,map,num.str());
      }
      
      return true;
    }

    template<typename T0, typename T1>
    bool XMLSerializableObject::Serialize(std::map<T0,T1>*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      return Serialize(*obj,doc,element,name);
    }

    bool XMLSerializableObject::Serialize(std::map<int,int>& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      tinyxml2::XMLElement* map = doc->NewElement(name.c_str());
      element->InsertEndChild(map);
      
      const unsigned int size = obj.size();
      map->SetAttribute("size",size);
      
      std::stringstream ms;
      ms << "\n";
      for(std::map<int,int>::iterator it=obj.begin();it!=obj.end();it++)
        ms << it->first << "," << it->second << ";";
      ms << "\n";

      std::string mString = ms.str();
      if(mString.size() > 1)
        mString[mString.size()-2] = '\0';
      
      map->SetAttribute("map_int_int",mString.c_str());

      return true;
    }

    bool XMLSerializableObject::Serialize(std::map<int,int>*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      return Serialize(*obj,doc,element,name);
    }

    template<typename T0, typename T1>
    bool XMLSerializableObject::Deserialize(std::map<T0,T1>& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      bool res = true;
      obj.clear();
      
      const tinyxml2::XMLElement* child = element->FirstChildElement(name.c_str());
      if(child != NULL)
      {
        unsigned int size = child->UnsignedAttribute("size");
        
        std::stringstream num;
        typename std::map<T0,T1>::iterator iter = obj.begin();
        for(int i=0;i<size;i++)
        {
          num.str("");
          num << "value" << i;
          
          std::pair<T0,T1> pair;
          res &= Deserialize(pair,doc,child,num.str());
          obj.insert(pair);
        }
      }
      else
        return false;
      
      return res;
    }

    template<typename T0, typename T1>
    bool XMLSerializableObject::Deserialize(std::map<T0,T1>*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      obj = new std::map<T0,T1>();
      return Deserialize(*obj,doc,element,name);
    }

    bool XMLSerializableObject::Deserialize(std::map<int,int>& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      bool res = true;
      
      const tinyxml2::XMLElement* child = element->FirstChildElement(name.c_str());
      if(child != NULL)
      {
        unsigned int size = child->UnsignedAttribute("size");

        obj.clear();
        
        const tinyxml2::XMLAttribute* attribute = child->FindAttribute("map_int_int");
        if(attribute == NULL)
        {
          Init(obj);
          return false;
        }

        char* matTemp;
        GetAttribute(attribute->Value(),matTemp);
        
        std::string line, ssize;
        std::stringstream mats;
        mats.str(matTemp);
        
        std::string val;
        // for each line
        getline(mats,line); // starts with an empty line
        while(getline(mats,line))
        {
          // get current line
          std::stringstream liness(line);
          
          for(unsigned int c=0;c<size;c++)
          {
            std::pair<int,int> pair;
            
            // split line
            getline(liness, val, ',');
            
            // push pack the data if any
            if(!val.empty())
              GetAttribute(val.c_str(),pair.first);

            // split line
            getline(liness, val, ';');

            // push pack the data if any
            if(!val.empty())
              GetAttribute(val.c_str(),pair.second);

            obj.insert(pair);
          }
        }
      }
      else
      {
        Init(obj);
        return false;
      }
      
      return res;
    }

    bool XMLSerializableObject::Deserialize(std::map<int,int>*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      obj = new std::map<int,int>();
      return Deserialize(*obj,doc,element,name);
    }

    template<typename T, int R, int C>
    bool XMLSerializableObject::Serialize(Eigen::Matrix<T,R,C>& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      tinyxml2::XMLElement* matrix = doc->NewElement(name.c_str());
      element->InsertEndChild(matrix);
      
      const unsigned int rows = obj.rows();
      const unsigned int cols = obj.cols();
      
      matrix->SetAttribute("rows",rows);
      matrix->SetAttribute("cols",cols);
      
      char buffer[200];
      std::stringstream ms;
      ms << "\n";
      for(unsigned int r=0;r<rows;r++)
      {
        for(unsigned int c=0;c<cols;c++)
        {
          tinyxml2::XMLUtil::ToStr(obj(r,c),buffer,200);
          ms << buffer << ",";
        }
        ms << "\n";
      }
      
      std::string mString = ms.str();
      if(mString.size() > 1)
        mString[mString.size()-2] = '\0';
      
      matrix->SetAttribute("matrix",mString.c_str());
      
      return true;
    }
    
    template<typename T, int R, int C>
    bool XMLSerializableObject::Serialize(Eigen::Matrix<T,R,C>*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      return Serialize(*obj,doc,element,name);
    }
    
    template<typename T, int R, int C>
    bool XMLSerializableObject::Deserialize(Eigen::Matrix<T,R,C>& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      const tinyxml2::XMLElement* child = element->FirstChildElement(name.c_str());
      bool initialized = false;
      if(child != NULL)
      {
        const unsigned int rows = child->UnsignedAttribute("rows");
        const unsigned int cols = child->UnsignedAttribute("cols");
        
        if(rows > 0 && cols > 0)
        {
          obj.resize(rows,cols);
        
          const tinyxml2::XMLAttribute* attribute = child->FindAttribute("matrix");
          if(attribute == NULL)
          {
            Init(obj);
            return false;
          }
        
          char* matTemp;
          GetAttribute(attribute->Value(),matTemp);
        
          std::string line, srows, scols;
          std::stringstream mats;
          mats.str(matTemp);
        
          int r=0;
          std::string val;
          // for each line
          getline(mats,line);
          while(getline(mats,line))
          {
            // get current line
            std::stringstream liness(line);
          
            for(unsigned int c=0;c<cols-1;c++)
            {
              // split line
              getline(liness, val, ',');
            
              // push pack the data if any
              if(!val.empty())
                GetAttribute(val.c_str(),obj.coeffRef(r,c));
            }
          
            getline(liness, val);
            GetAttribute(val.c_str(),obj.coeffRef(r,cols-1));
          
            r++;
          }
          initialized = true;
        }
      }
      
      if(!initialized)
      {
        Init(obj);
        return false;
      }
      
      return true;
    }
    
    template<typename T, int R, int C>
    bool XMLSerializableObject::Deserialize(Eigen::Matrix<T,R,C>*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      obj = new Eigen::PlainObjectBase<T>();
      return Deserialize(*obj,doc,element,name);
    }
    
    template<typename T>
    bool XMLSerializableObject::Serialize(Eigen::SparseMatrix<T>& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      tinyxml2::XMLElement* matrix = doc->NewElement(name.c_str());
      element->InsertEndChild(matrix);
      
      const unsigned int rows = obj.rows();
      const unsigned int cols = obj.cols();
      
      matrix->SetAttribute("rows",rows);
      matrix->SetAttribute("cols",cols);
      
      char buffer[200];
      std::stringstream ms;
      ms << "\n";
      for (int k=0;k<obj.outerSize();++k)
      {
        for (typename Eigen::SparseMatrix<T>::InnerIterator it(obj,k);it;++it)
        {
          tinyxml2::XMLUtil::ToStr(it.value(),buffer,200);
          ms << it.row() << "," << it.col() << "," << buffer << "\n";
        }
      }
      
      std::string mString = ms.str();
      if(mString.size() > 0)
        mString[mString.size()-1] = '\0';
      
      matrix->SetAttribute("matrix",mString.c_str());
      
      return true;
    }
    
    template<typename T>
    bool XMLSerializableObject::Serialize(Eigen::SparseMatrix<T>*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      return Serialize(*obj,doc,element,name);
    }
    
    template<typename T>
    bool XMLSerializableObject::Deserialize(Eigen::SparseMatrix<T>& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      const tinyxml2::XMLElement* child = element->FirstChildElement(name.c_str());
      bool initialized = false;
      if(child != NULL)
      {
        const unsigned int rows = child->UnsignedAttribute("rows");
        const unsigned int cols = child->UnsignedAttribute("cols");
        
        if(rows > 0 && cols > 0)
        {
          obj.resize(rows,cols);
          obj.setZero();
        
          const tinyxml2::XMLAttribute* attribute = child->FindAttribute("matrix");
          if(attribute == NULL)
          {
            Init(obj);
            return false;
          }
        
          char* matTemp;
          GetAttribute(attribute->Value(),matTemp);
        
          std::string line, srows, scols;
          std::stringstream mats;
          mats.str(matTemp);
        
          std::vector<Eigen::Triplet<T> > triplets;
          int r=0;
          std::string val;
        
          // for each line
          getline(mats,line);
          while(getline(mats,line))
          {
            // get current line
            std::stringstream liness(line);
          
            // row
            getline(liness, val, ',');
            int row = atoi(val.c_str());
            // col
            getline(liness, val, ',');
            int col = atoi(val.c_str());
            // val
            getline(liness, val);
            T value;
            GetAttribute(val.c_str(),value);
          
            triplets.push_back(Eigen::Triplet<T>(row,col,value));
          
            r++;
          }
        
          obj.setFromTriplets(triplets.begin(),triplets.end());
          initialized = true;
        }
      }
      
      if(!initialized)
      {
        Init(obj);
        return false;
      }
      
      return true;
    }
    
    template<typename T>
    bool XMLSerializableObject::Deserialize(Eigen::SparseMatrix<T>*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      obj = new Eigen::SparseMatrix<T>();
      return Deserialize(*obj,doc,element,name);
    }
    
    template<typename T>
    bool XMLSerializableObject::setElementAttribute(T& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      tinyxml2::XMLElement* child = doc->NewElement(name.c_str());
      element->InsertEndChild(child);
      SetAttribute(child,"val",obj);
      return true;
    }
    
    template<typename T>
    bool XMLSerializableObject::getElementAttribute(T& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
    {
      // basic data type
      const tinyxml2::XMLElement* child = element->FirstChildElement(name.c_str());
      if(child != NULL)
      {
        XMLSerializableObject::GetAttribute(child->Attribute("val"),obj);
        return true;
      }
      else
      {
        Init(obj);
        return false;
      }
    }
    
    template<typename T>
    XMLSerializableInstance<T>::XMLSerializableInstance(T& obj, const std::string& name, const std::string group)
    : XMLSerializableObject(name, group), Object(obj)
    {
      XMLSerializableObject::Init(DefaultValue);
    }
    
    template<typename T>
    XMLSerializableInstance<T>::XMLSerializableInstance(T& obj, const std::string& name, const std::string group, T defaultValue)
    : XMLSerializableObject(name, group), Object(obj), DefaultValue(defaultValue)
    {
    }
    
    template<typename T>
    XMLSerializableInstance<T>::~XMLSerializableInstance()
    {
    }
    
    template<typename T>
    void XMLSerializableInstance<T>::Init()
    {
      XMLSerializableObject::Init(DefaultValue);
    }
    
    template<typename T>
    bool XMLSerializableInstance<T>::Serialize(tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element)
    {
      return XMLSerializableObject::Serialize(Object,doc,element,Name);
    }
    
    template<typename T>
    bool XMLSerializableInstance<T>::Deserialize(tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element)
    {
      return XMLSerializableObject::Deserialize(Object,doc,element,Name);
    }
    
    template<typename T>
    bool XMLSerializer::SaveObject(T& object, const char* filename)
    {
      return SaveObject(object,"Object","Serialization",filename,true);
    }
    
    template<typename T>
    bool XMLSerializer::SaveObject(T& object, const std::string& objectName, const std::string& groupName, const char* filename, bool overwrite)
    {
      bool result = true;
      XMLSerializer* serializer = new XMLSerializer(groupName);
      result &= serializer->Add(object,objectName);
      result &= serializer->Save(objectName,groupName,filename,overwrite);
      delete serializer;
      return result;
    }
    
    template<typename T>
    bool XMLSerializer::LoadObject(T& object, const char* filename)
    {
      return LoadObject(object,"Object","Serialization",filename);
    }
    
    template<typename T>
    bool XMLSerializer::LoadObject(T& object, const std::string& objectName, const std::string& groupName, const char* filename)
    {
      bool result = true;
      XMLSerializer* serializer = new XMLSerializer(groupName);
      result &= serializer->Add(object,objectName);
      result &= serializer->Load(objectName,groupName,filename);
      delete serializer;
      return result;
    }
    
    XMLSerializer::XMLSerializer(const std::string& defaultGroup)
    {
      SetCurrentGroup(defaultGroup);
    }
    
    XMLSerializer::~XMLSerializer()
    {
      std::map<std::string,XMLSerializerGroup*>::iterator it;
      for (it=groups.begin();it!=groups.end();it++)
      {
        delete it->second->Objects;
        delete it->second;
      }
    }
    
    bool XMLSerializer::Save(const char* filename, bool overwrite)
    {
      tinyxml2::XMLDocument* doc = new tinyxml2::XMLDocument();
      
      if(overwrite == false)
      {
        // Check if file exists
        tinyxml2::XMLError error = doc->LoadFile(filename);
        if(error != tinyxml2::XML_NO_ERROR)
          doc->Clear();
      }
      
      if(SaveToXMLDoc(doc) == false)
        return false;
      
      // Save
      tinyxml2::XMLError error = doc->SaveFile(filename);
      if(error != tinyxml2::XML_NO_ERROR)
      {
        doc->PrintError();
        return false;
      }
      
      delete doc;
      
      return true;
    }
    
    bool XMLSerializer::SaveToXMLDoc(tinyxml2::XMLDocument* doc)
    {
      std::map<std::string,XMLSerializerGroup*>::iterator it;
      for (it=groups.begin();it!=groups.end();it++)
      {
        // Update group
        tinyxml2::XMLElement* element = doc->FirstChildElement(it->first.c_str());
        if(element != NULL)
        {
          element->DeleteChildren();
        }
        else
        {
          element = doc->NewElement(it->first.c_str());
          doc->InsertEndChild(element);
        }
        
        std::vector<XMLSerializable*>* group = it->second->Objects;
        for(unsigned  int i=0;i<group->size();i++)
        {
          if((*group)[i]->Serialize(doc,element) == false)
            return false;
        }
      }
      
      return true;
    }
    
    bool XMLSerializer::Save(const std::string& groupName, const char* filename, bool overwrite)
    {
      tinyxml2::XMLDocument* doc = new tinyxml2::XMLDocument();
      
      if(overwrite == false)
      {
        // Check if file exists
        tinyxml2::XMLError error = doc->LoadFile(filename);
        if(error != tinyxml2::XML_NO_ERROR)
          doc->Clear();
      }
      
      if(SaveToXMLDoc(groupName, doc) == false)
        return false;
      
      // Save
      tinyxml2::XMLError error = doc->SaveFile(filename);
      if(error != tinyxml2::XML_NO_ERROR)
      {
        doc->PrintError();
        return false;
      }
      
      delete doc;
      
      return true;
    }
    
    bool XMLSerializer::SaveToXMLDoc(const std::string& groupName, tinyxml2::XMLDocument* doc)
    {
      std::string gn = groupName;
      EncodeXMLElementName(gn);
      
      std::map<std::string,XMLSerializerGroup*>::iterator it = groups.find(gn);
      if(it == groups.end())
        return false;
      
      // Update group
      tinyxml2::XMLElement* element = doc->FirstChildElement(it->first.c_str());
      if(element != NULL)
      {
        element->DeleteChildren();
      }
      else
      {
        element = doc->NewElement(it->first.c_str());
        doc->InsertEndChild(element);
      }
      
      std::vector<XMLSerializable*>* groups = it->second->Objects;
      for(unsigned int i=0;i<groups->size();i++)
      {
        if((*groups)[i]->Serialize(doc,element) == false)
          return false;
      }
      
      return true;
    }
    
    bool XMLSerializer::Save(const std::string& objectName, const std::string& groupName, const char* filename, bool overwrite)
    {
      tinyxml2::XMLDocument* doc = new tinyxml2::XMLDocument();
      
      if(overwrite == false)
      {
        // Check if file exists
        tinyxml2::XMLError error = doc->LoadFile(filename);
        if(error != tinyxml2::XML_NO_ERROR)
          doc->Clear();
      }
      
      if(SaveToXMLDoc(objectName, groupName, doc) == false)
        return false;
      
      // Save
      tinyxml2::XMLError error = doc->SaveFile(filename);
      if(error != tinyxml2::XML_NO_ERROR)
      {
        doc->PrintError();
        return false;
      }
      
      delete doc;
      
      return true;
    }
    
    bool XMLSerializer::SaveToXMLDoc(const std::string& objectName, const std::string& groupName, tinyxml2::XMLDocument* doc)
    {
      std::string gn = groupName;
      EncodeXMLElementName(gn);
      
      std::string on = objectName;
      EncodeXMLElementName(on);
      
      std::map<std::string,XMLSerializerGroup*>::iterator it = groups.find(gn);
      if(it == groups.end())
        return false;
      
      // Get/Add group
      tinyxml2::XMLElement* element = findAddGroup(doc, it->first.c_str());
      
      // Serialize
      std::vector<XMLSerializable*>* groups = it->second->Objects;
      bool found = false;
      for(unsigned int i=0;i<groups->size();i++)
      {
        if((*groups)[i]->Name == on)
        {
          found = true;
          
          tinyxml2::XMLElement* child = element->FirstChildElement(on.c_str());
          if(child != NULL)
          {
            element->DeleteChild(child);
          }
          
          if((*groups)[i]->Serialize(doc,element) == false)
            return false;
        }
      }
      
      return found;
    }
    
    bool XMLSerializer::SaveGroupToXMLElement(tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      std::string test = currentGroup->first;
      return SaveGroupToXMLElement(currentGroup->first,doc,element,name);
    }
    
    bool XMLSerializer::SaveGroupToXMLElement(const std::string& groupName, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
    {
      std::string gn = groupName;
      EncodeXMLElementName(gn);
      
      std::map<std::string,XMLSerializerGroup*>::iterator it = groups.find(gn);
      if(it == groups.end())
        return false;
      
      // Add new group
      tinyxml2::XMLElement* group = doc->NewElement(name.c_str());
      element->InsertEndChild(group);
      
      std::vector<XMLSerializable*>* groups = it->second->Objects;
      for(unsigned int i=0;i<groups->size();i++)
      {
        if((*groups)[i]->Serialize(doc,group) == false)
          return false;
      }
      
      return true;
    }
    
    bool XMLSerializer::Load(const char* filename)
    {
      tinyxml2::XMLDocument* doc = openDoc(filename);
      if(doc == NULL)
        return false;
      
      if(LoadFromXMLDoc(doc) == false)
        return false;
      
      delete doc;
      
      return true;
    }
    
    bool XMLSerializer::LoadFromXMLDoc(tinyxml2::XMLDocument* doc)
    {
      std::map<std::string,XMLSerializerGroup*>::iterator it;
      for (it=groups.begin();it!=groups.end();it++)
      {
        tinyxml2::XMLElement* element = doc->FirstChildElement(it->first.c_str());
        /*if(element == NULL)
          return false;*/
        
        // Deserialize
        std::vector<XMLSerializable*>* group = it->second->Objects;
        for(unsigned int i=0;i<group->size();i++)
        {
          if(element == NULL || (*group)[i]->Deserialize(doc,element) == false)
            (*group)[i]->Init(); // Load default value;
        }
      }
      
      return true;
    }
    
    bool XMLSerializer::Load(const std::string& groupName, const char* filename)
    {
      tinyxml2::XMLDocument* doc = openDoc(filename);
      if(doc == NULL)
        return false;
      
      if(LoadFromXMLDoc(groupName, doc) == false)
        return false;
      
      delete doc;
      
      return true;
    }
    
    bool XMLSerializer::LoadFromXMLDoc(const std::string& groupName, tinyxml2::XMLDocument* doc)
    {
      std::string gn = groupName;
      EncodeXMLElementName(gn);
      
      std::map<std::string,XMLSerializerGroup*>::iterator it = groups.find(gn);
      if(it == groups.end())
        return false;
      
      tinyxml2::XMLElement* element = doc->FirstChildElement(it->first.c_str());
      /*if(element == NULL)
        return false;*/
      
      // Deserialize
      std::vector<XMLSerializable*>* groups = it->second->Objects;
      for(unsigned int i=0;i<groups->size();i++)
      {
        if(element == NULL || (*groups)[i]->Deserialize(doc,element) == false)
          (*groups)[i]->Init(); // Load default value;
      }
      
      return true;
    }
    
    bool XMLSerializer::Load(const std::string& objectName, const std::string& groupName, const char* filename)
    {
      tinyxml2::XMLDocument* doc = openDoc(filename);
      if(doc == NULL)
        return false;
      
      if(LoadFromXMLDoc(objectName,groupName,doc) == false)
        return false;
      
      delete doc;
      
      return true;
    }
    
    bool XMLSerializer::LoadFromXMLDoc(const std::string& objectName, const std::string& groupName, tinyxml2::XMLDocument* doc)
    {
      std::string gn = groupName;
      EncodeXMLElementName(gn);
      
      std::string on = objectName;
      EncodeXMLElementName(on);
      
      std::map<std::string,XMLSerializerGroup*>::iterator it = groups.find(gn);
      if(it == groups.end())
        return false;
      
      tinyxml2::XMLElement* element = doc->FirstChildElement(it->first.c_str());
      /*if(element == NULL)
        return false;*/
      
      // Deserialize
      std::vector<XMLSerializable*>* groups = it->second->Objects;
      bool found = false;
      for(unsigned int i=0;i<groups->size();i++)
      {
        if((*groups)[i]->Name == on)
        {
          found = true;
          if(element == NULL || (*groups)[i]->Deserialize(doc,element) == false)
            (*groups)[i]->Init(); // Load default value;
        }
      }
      
      return found;
    }
    
    bool XMLSerializer::LoadGroupFromXMLElement(const std::string& groupName, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element)
    {
      std::string gn = groupName;
      EncodeXMLElementName(gn);
      
      std::map<std::string,XMLSerializerGroup*>::iterator it = groups.find(gn);
      
      if(it == groups.end())
        return false;
      
      const tinyxml2::XMLElement* group = element->FirstChildElement(groupName.c_str());
      if(group == NULL)
        return false;
      
      // Deserialize
      std::vector<XMLSerializable*>* groups = it->second->Objects;
      for(unsigned int i=0;i<groups->size();i++)
      {
        if(element == NULL || (*groups)[i]->Deserialize(doc,group) == false)
          (*groups)[i]->Init(); // Load default value;
      }
      
      return true;
    }
    
    void XMLSerializer::SetCurrentGroup(const std::string& group)
    {
      currentGroup = setGetGroup(group);
    }
    
    std::string XMLSerializer::GetCurrentGroup()
    {
      return currentGroup->first;
    }
    
    template<typename T>
    bool XMLSerializer::Add(T& obj, const std::string& name)
    {
      XMLSerializable* object = dynamic_cast<XMLSerializable*>(&obj);
      
      object->Name = name;
      return addObjectToGroup(object,currentGroup);
    }
    
    bool XMLSerializer::Add(char& obj, const std::string& name)
    {
      return add(obj,name);
    }
    
    bool XMLSerializer::Add(char*& obj, const std::string& name)
    {
      return add(obj,name);
    }
    
    bool XMLSerializer::Add(std::string& obj, const std::string& name)
    {
      return add(obj,name);
    }
    
    bool XMLSerializer::Add(bool& obj, const std::string& name)
    {
      return add(obj,name);
    }
    
    bool XMLSerializer::Add(unsigned int& obj, const std::string& name)
    {
      return add(obj,name);
    }
    
    bool XMLSerializer::Add(int& obj, const std::string& name)
    {
      return add(obj,name);
    }
    
    bool XMLSerializer::Add(float& obj, const std::string& name)
    {
      return add(obj,name);
    }
    
    bool XMLSerializer::Add(double& obj, const std::string& name)
    {
      return add(obj,name);
    }
    
    /*template<typename T, size_t S>
     bool XMLSerializer::Add(std::array<T,S>& obj, const std::string& name)
     {
     return add(obj,name);
     }*/
    
    template<typename T0, typename T1>
    bool XMLSerializer::Add(std::pair<T0,T1>& obj, const std::string& name)
    {
      return add(obj,name);
    }
    
    template<typename T>
    bool XMLSerializer::Add(std::vector<T>& obj, const std::string& name)
    {
      return add(obj,name);
    }

    template<typename T>
    bool XMLSerializer::Add(std::set<T>& obj, const std::string& name)
    {
      return add(obj,name);
    }

    template<typename T0, typename T1>
    bool XMLSerializer::Add(std::map<T0,T1>& obj, const std::string& name)
    {
      return add(obj,name);
    }
    
    template<typename T, int R, int C>
    bool XMLSerializer::Add(Eigen::Matrix<T,R,C>& obj, const std::string& name)
    {
      return add(obj,name);
    }
    
    template<typename T>
    bool XMLSerializer::Add(Eigen::SparseMatrix<T>& obj, const std::string& name)
    {
      return add(obj,name);
    }
    
    template<typename T>
    bool XMLSerializer::Add(T& object, const std::string& name, T defaultValue)
    {
      return false;
    }
    
    bool XMLSerializer::Add(char& obj, const std::string& name, char defaultValue)
    {
      return add(obj,name,defaultValue);
    }
    
    bool XMLSerializer::Add(char*& obj, const std::string& name, char* defaultValue)
    {
      return add(obj,name,defaultValue);
    }
    
    bool XMLSerializer::Add(std::string& obj, const std::string& name, std::string defaultValue)
    {
      return add(obj,name,defaultValue);
    }
    
    bool XMLSerializer::Add(bool& obj, const std::string& name, bool defaultValue)
    {
      return add(obj,name,defaultValue);
    }
    
    bool XMLSerializer::Add(unsigned int& obj, const std::string& name, unsigned int defaultValue)
    {
      return add(obj,name,defaultValue);
    }
    
    bool XMLSerializer::Add(int& obj, const std::string& name, int defaultValue)
    {
      return add(obj,name,defaultValue);
    }
    
    bool XMLSerializer::Add(float& obj, const std::string& name, float defaultValue)
    {
      return add(obj,name,defaultValue);
    }
    
    bool XMLSerializer::Add(double& obj, const std::string& name, double defaultValue)
    {
      return add(obj,name,defaultValue);
    }
    
    template<typename T>
    bool XMLSerializer::add(T& obj, const std::string& name)
    {
      XMLSerializable* object = new XMLSerializableInstance<T>(obj,name,currentGroup->first);
      return addObjectToGroup(object,currentGroup);
    }
    
    template<typename T>
    bool XMLSerializer::add(T& obj, const std::string& name, T defaultValue)
    {
      XMLSerializable* object = new XMLSerializableInstance<T>(obj,name,currentGroup->first,defaultValue);
      return addObjectToGroup(object,currentGroup);
    }
    
    bool XMLSerializer::addObjectToGroup(XMLSerializable* obj, const std::string& group)
    {
      std::map<std::string,XMLSerializerGroup*>::iterator it = setGetGroup(group);
      return addObjectToGroup(obj, it);
    }
    
    bool XMLSerializer::addObjectToGroup(XMLSerializable* object, std::map<std::string,XMLSerializerGroup*>::iterator it)
    {
      std::vector<XMLSerializable*>* objects = it->second->Objects;
      for(unsigned int i=0;i<objects->size();i++)
      {
        if((*objects)[i]->Name == object->Name)
          return false;
      }
      
      objects->push_back(object);
      
      return true;
    }
    
    std::map<std::string,XMLSerializerGroup*>::iterator XMLSerializer::setGetGroup(const std::string& group)
    {
      std::string groupName = group;
      EncodeXMLElementName(groupName);
      
      std::map<std::string,XMLSerializerGroup*>::iterator it = groups.find(groupName);
      if(it == groups.end())
      {
        XMLSerializerGroup* newGroup = new XMLSerializerGroup();
        newGroup->Objects = new std::vector<XMLSerializable*>();
        groups[groupName] = newGroup;
        it = groups.find(groupName);
      }
      
      return it;
    }
    
    tinyxml2::XMLDocument* XMLSerializer::openDoc(const char* filename)
    {
      tinyxml2::XMLDocument* doc = new tinyxml2::XMLDocument();
      
      tinyxml2::XMLError error = doc->LoadFile(filename);
      if(error != tinyxml2::XML_NO_ERROR)
      {
        doc->PrintError();
        doc = NULL;
      }
      
      return doc;
    }
    
    tinyxml2::XMLElement* XMLSerializer::findAddGroup(tinyxml2::XMLDocument* doc, const char* groupName)
    {
      tinyxml2::XMLElement* group = doc->FirstChildElement(groupName);
      if(group == NULL)
      {
        group = doc->NewElement(groupName);
        doc->InsertEndChild(group);
      }
      return group;
    }
  }
}
#endif
