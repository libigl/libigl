#ifndef IGL_XML_SERIALIZATION_H
#define IGL_XML_SERIALIZATION_H

#include <igl/xml/serialize_xml.h>

namespace igl
{
  class XMLSerialization: public XMLSerializable
  {
  private:
    
    template <typename T>
    struct XMLSerializationObject: public XMLSerializable
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

      void Serialize(tinyxml2::XMLDocument* doc,tinyxml2::XMLElement* element) const {
        igl::serialize_xml(*Object,Name,doc,element,Binary);
      }

      void Deserialize(const tinyxml2::XMLDocument* doc,const tinyxml2::XMLElement* element) {
        igl::deserialize_xml(*Object,Name,doc,element);
      }
    };

    mutable bool initialized;
    mutable std::vector<XMLSerializable*> objects;
  
  public:

    /**
    * Override this function to add your member variables which should be
    * serialized:
    *
    * this->Add(var1);
    * this->Add(var2);
    * ...
    */
    virtual void InitSerialization() = 0;

    /**
    * Following functions can be overridden to handle the specific events.
    * Return false to prevent the de-/serialization of an object.
    */
    virtual bool BeforeSerialization() const { return true; }
    virtual void AfterSerialization() const {}
    virtual bool BeforeDeserialization() { return true; }
    virtual void AfterDeserialization() {}

    /**
    * Default implementation of XMLSerializable interface
    */
    void Serialize(std::vector<char>& buffer) const
    {
      if(this->BeforeSerialization())
      {
        if(initialized == false)
        {
          objects.clear();
          (const_cast<XMLSerialization*>(this))->InitSerialization();
          initialized = true;
        }

        for(int i=0;i<objects.size();i++)
          objects[i]->Serialize(buffer);

        this->AfterSerialization();
      }
    }

    void Deserialize(const std::vector<char>& buffer)
    {
      if(this->BeforeDeserialization())
      {
        if(initialized == false)
        {
          objects.clear();
          (const_cast<XMLSerialization*>(this))->InitSerialization();
          initialized = true;
        }

        for(int i=0;i<objects.size();i++)
          objects[i]->Deserialize(buffer);

        this->AfterDeserialization();
      }
    }

    void Serialize(tinyxml2::XMLDocument* doc,tinyxml2::XMLElement* element) const
    {
      if(this->BeforeSerialization())
      {
        if(initialized == false)
        {
          objects.clear();
          (const_cast<XMLSerialization*>(this))->InitSerialization();
          initialized = true;
        }

        for(int i=0;i<objects.size();i++)
          objects[i]->Serialize(doc,element);

        this->AfterSerialization();
      }
    }

    void Deserialize(const tinyxml2::XMLDocument* doc,const tinyxml2::XMLElement* element)
    {
      if(this->BeforeDeserialization())
      {
        if(initialized == false)
        {
          objects.clear();
          (const_cast<XMLSerialization*>(this))->InitSerialization();
          initialized = true;
        }

        for(int i=0;i<objects.size();i++)
          objects[i]->Deserialize(doc,element);

        this->AfterDeserialization();
      }
    }

    /**
    * Default constructor, destructor, assignment and copy constructor
    */

    XMLSerialization()
    {
      initialized = false;
    }
    
    XMLSerialization(const XMLSerialization& obj)
    {
      initialized = false;
      objects.clear();
    }

    ~XMLSerialization()
    {
      initialized = false;
      objects.clear();
    }


    XMLSerialization& operator=(const XMLSerialization& obj)
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

    /**
    * Use this function to add your variables which should be serialized.
    */
  
    template <typename T>
    void Add(T& obj,std::string name, bool binary = false)
    {
      XMLSerializationObject<T>* object = new XMLSerializationObject<T>();
      object->Binary = binary;
      object->Name = name;
      object->Object = &obj;

      objects.push_back(object);
    }
  };

}

#endif