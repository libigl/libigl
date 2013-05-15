/* ---------------------------------------------------------------------------
// XMLSerializer.h
// Author: Christian Schüller on 08/05/13.
------------------------------------------------------------------------------

This class allows to save and load a serialization of basic c++ data types like
char, char*, std::string, bool, uint, int, float, double to and from a xml file.
Containers like std::vector, std::std::pair, Eigen dense and sparse matrices are supported
as well as combination of them (like vector<pair<string,bool>> or vector<vector<int>>).
To serialize an arbitary object use the XMLSerializable interface.
 
The serialized objects are organised in groups in the xml file and have
their own names which must be unique within one group.

You can find examples how to use it in the test case class XMLSerializerTest.

----------------------------------------------------------------------------*/
#ifndef XML_SERIALIZER_H
#define XML_SERIALIZER_H

#include <iostream>
#include <array>
#include <vector>
#include <map>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "tinyxml2.h"

namespace igl
{

  void EncodeXMLElementName(std::string& name);
  //void DecodeXMLElementName(std::string& name);
  void ReplaceSubString(std::string& str, const std::string& search, const std::string& replace);

  // Forward declaration
  class XMLSerializer;

  /**
     * interface XMLSerializable
     * Inherit from this interface to have full control over the serialization of you user defined class.
     */ 
  class XMLSerializable
  {
  public:
    std::string Name;

    /**
     * This function gets called if the objects were not found during deserialization.
     * Initialize your objects as you like. 
     */ 
    virtual void Init() = 0;
    /**
     * Serialize your stuff within this function.
     * It contains the current serialization xml file. You can use SaveToXMLDoc or SaveGroupToXMLElement to add your objects.
     */ 
    virtual bool Serialize(tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element) = 0;

    /**
     * Deserialize your stuff within this function.
     * It contains the current serialization xml file. You can use LoadFromXMLDoc or LoadGroupFromXMLElement to read out your objects.
     */ 
    virtual bool Deserialize(tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element) = 0;
  };

  /**
     * class XMLSerialization
     * Inherit from this class to have the easiest way to serialize your user defined class.
     */
  class XMLSerialization : public XMLSerializable
  {
  public:
    igl::XMLSerializer* xmlSerializer;

    /**
     * Default implementation of XMLSerializable interface
     */
    virtual void Init();
    virtual bool Serialize(tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element);
    virtual bool Deserialize(tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element);

    XMLSerialization(const std::string& name);
    ~XMLSerialization();

    /**
     * Following functions can be overwritten to handle the specific events.
     * Return false to prevent serialization of object.
     */
    virtual bool BeforeSerialization();
    virtual void AfterSerialization();
    virtual bool BeforeDeserialization();
    virtual void AfterDeserialization();
  };


  /**
     * class XMLSerializableObject
     * internal usage
     */
  class XMLSerializableObject : public XMLSerializable
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

    // initialize objects   
    template<typename T>
    void Init(T& obj);
    template<typename T>
    void Init(T*& obj);
    template<typename T, int S>
    void Init(std::array<T,S>& obj);
    template<typename T0, typename T1>
    void Init(std::pair<T0,T1>& obj);
    template<typename T>
    void Init(std::vector<T>& obj);
    template<typename T, int R, int C>
    void Init(Eigen::Matrix<T,R,C>& obj);
    template<typename T>
    void Init(Eigen::SparseMatrix<T>& obj);

    // base types and XMLSerializable objects
    template<typename T>
    bool Serialize(T& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
    template<typename T>
    bool Serialize(T*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
    template<typename T>
    bool Deserialize(T& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
    template<typename T>
    bool Deserialize(T*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);

    // std::array<T>
    template<typename T, int S>
    bool Serialize(std::array<T,S>& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
    template<typename T, int S>
    bool Serialize(std::array<T,S>*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
    template<typename T, int S>
    bool Deserialize(std::array<T,S>& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
    template<typename T, int S>
    bool Deserialize(std::array<T,S>*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);

    // std::pair<T0,T1>
    template<typename T0, typename T1>
    bool Serialize(std::pair<T0,T1>& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
    template<typename T0, typename T1>
    bool Serialize(std::pair<T0,T1>*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
    template<typename T0, typename T1>
    bool Deserialize(std::pair<T0,T1>& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
    template<typename T0, typename T1>
    bool Deserialize(std::pair<T0,T1>*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);

    // std::vector<T>
    template<typename T>
    bool Serialize(std::vector<T>& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
    template<typename T>
    bool Serialize(std::vector<T>*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
    template<typename T>
    bool Deserialize(std::vector<T>& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
    template<typename T>
    bool Deserialize(std::vector<T>*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
   
    // Eigen dense matrix
    template<typename T, int R, int C>
    bool Serialize(Eigen::Matrix<T,R,C>& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
    template<typename T, int R, int C>
    bool Serialize(Eigen::Matrix<T,R,C>*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
    template<typename T, int R, int C>
    bool Deserialize(Eigen::Matrix<T,R,C>& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);
    template<typename T, int R, int C>
    bool Deserialize(Eigen::Matrix<T,R,C>*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name);

    // Eigen sparse matrix
    template<typename T>
    bool Serialize(Eigen::SparseMatrix<T>& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
    template<typename T>
    bool Serialize(Eigen::SparseMatrix<T>*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name);
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
    template<typename T>
    bool Add(T& object, const std::string& name);
    template<typename T>
    bool Add(T& object, const std::string& name, T defaultValue);

    // stl containers
    template<typename T, int S>
    bool Add(std::array<T,S>& obj, const std::string& name);
    template<typename T0, typename T1>
    bool Add(std::pair<T0,T1>& obj, const std::string& name);
    template<typename T>
    bool Add(std::vector<T>& obj, const std::string& name);
    template<typename T, int R, int C>
    
    // eigen matrices
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
  
  /**
     * class XMLSerializerTest
     * Used to test the functionality of the library and also shows howto use it.
     */
  class XMLSerializerTest : public igl::XMLSerialization
  {
  public:

    int testInt;
    std::vector<float> testVector;

    XMLSerializerTest();

    bool Test();
  };
}

#ifdef IGL_HEADER_ONLY
#include "XMLSerializer.cpp"
#endif

#endif
