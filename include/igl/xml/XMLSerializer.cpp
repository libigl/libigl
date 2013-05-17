#include "XMLSerializer.h"

//TODO: these should probably be static. Do they need to be in the igl
//namespace?
namespace igl
{
  int numForbiddenChars = 8;
  char forbiddenChars[] = {' ','/','~','#','&','>','<','='};
}

void igl::ReplaceSubString(std::string& str, const std::string& search, const std::string& replace)
{
  size_t pos = 0;
  while ((pos = str.find(search, pos)) != std::string::npos)
  {
    str.replace(pos, search.length(), replace);
    pos += replace.length();
  }
}

void igl::EncodeXMLElementName(std::string& name)
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

/*void igl::DecodeXMLElementName(std::string& name)
{
  if(name.find("::", 0) == 0)
    name.replace(0,3,"");

  std::stringstream stream;
  for(unsigned int i=0;i<numForbiddenChars;i++)
  {
    std::stringstream searchs;
    searchs << ":" << (int)forbiddenChars[i];
    std::string search = searchs.str();
    std::string replace;
    replace = forbiddenChars[i];
      
    ReplaceSubString(name,search,replace);
  }
}*/

void igl::XMLSerialization::Init()
{
}

bool igl::XMLSerialization::Serialize(tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element)
{
  bool serialized = false;
  
  if(BeforeSerialization())
  {
    serialized = xmlSerializer->SaveGroupToXMLElement(doc,element,Name);
    AfterSerialization();
  }

  return serialized;
}

bool igl::XMLSerialization::Deserialize(tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element)
{
  bool serialized = false;
  
  if(BeforeDeserialization())
  {
    serialized = xmlSerializer->LoadGroupFromXMLElement(Name,doc,element);
    AfterDeserialization();
  }
  
  return serialized;
}

igl::XMLSerialization::XMLSerialization(const std::string& name)
{
  Name = name;
  xmlSerializer = new igl::XMLSerializer(name);
}

igl::XMLSerialization::~XMLSerialization()
{
  delete xmlSerializer;
}

bool igl::XMLSerialization::BeforeSerialization()
{
  return true;
}

void igl::XMLSerialization::AfterSerialization()
{
}

bool igl::XMLSerialization::BeforeDeserialization()
{
  return true;
}

void igl::XMLSerialization::AfterDeserialization()
{
}

igl::XMLSerializableObject::XMLSerializableObject(const std::string& name, const std::string& group)
{
  std::string groupName = group;
  std::string objectName = name;
    
  igl::EncodeXMLElementName(groupName);
  igl::EncodeXMLElementName(objectName);
    
  Name = objectName;
}

igl::XMLSerializableObject::~XMLSerializableObject()
{
}

// set attribute conversion functions
void igl::XMLSerializableObject::SetAttribute(tinyxml2::XMLElement* element, const char* name, char& dest)
{ 
  element->SetAttribute(name,dest);
}
  
void igl::XMLSerializableObject::SetAttribute(tinyxml2::XMLElement* element, const char* name, char*& dest)
{ 
  element->SetAttribute(name,const_cast<const char*>(dest));
}

void igl::XMLSerializableObject::SetAttribute(tinyxml2::XMLElement* element, const char* name, std::string& dest)
{
  element->SetAttribute(name,dest.c_str());
}

void igl::XMLSerializableObject::SetAttribute(tinyxml2::XMLElement* element, const char* name, bool& dest)
{
  element->SetAttribute(name,dest);
}
  
void igl::XMLSerializableObject::SetAttribute(tinyxml2::XMLElement* element, const char* name, unsigned int& dest)
{
  element->SetAttribute(name,dest);
}

void igl::XMLSerializableObject::SetAttribute(tinyxml2::XMLElement* element, const char* name, int& dest)
{
  element->SetAttribute(name,dest);
}

void igl::XMLSerializableObject::SetAttribute(tinyxml2::XMLElement* element, const char* name, float& dest)
{
  element->SetAttribute(name,dest);
}

void igl::XMLSerializableObject::SetAttribute(tinyxml2::XMLElement* element, const char* name, double& dest)
{
  element->SetAttribute(name,dest);
}

// get attribute conversion functions
void igl::XMLSerializableObject::GetAttribute(const char* src, char& dest)
{ 
  dest = (char)atoi(src);
}
  
void igl::XMLSerializableObject::GetAttribute(const char* src, char*& dest)
{ 
  unsigned int length = strlen(src)+1;
  dest = new char[length];
  strcpy(dest, src);
}

void igl::XMLSerializableObject::GetAttribute(const char* src, std::string& dest)
{
  dest = src;
}
  
void igl::XMLSerializableObject::GetAttribute(const char* src, bool& dest)
{
  tinyxml2::XMLUtil::ToBool(src,&dest);
}
  
void igl::XMLSerializableObject::GetAttribute(const char* src, unsigned int& dest)
{
  tinyxml2::XMLUtil::ToUnsigned(src,&dest);
}

void igl::XMLSerializableObject::GetAttribute(const char* src, int& dest)
{
  tinyxml2::XMLUtil::ToInt(src,&dest);
}

void igl::XMLSerializableObject::GetAttribute(const char* src, float& dest)
{
  tinyxml2::XMLUtil::ToFloat(src,&dest);
}

void igl::XMLSerializableObject::GetAttribute(const char* src, double& dest)
{
  tinyxml2::XMLUtil::ToDouble(src,&dest);
}

// specify default value of types
void igl::XMLSerializableObject::Init(char& val)
{
  val = '0';
}

void igl::XMLSerializableObject::Init(char*& val)
{
  val = NULL;
}

void igl::XMLSerializableObject::Init(std::string& val)
{
  val = "";
}

void igl::XMLSerializableObject::Init(bool& val)
{
  val = false;
}

void igl::XMLSerializableObject::Init(unsigned int& val)
{
  val = 0;
}

void igl::XMLSerializableObject::Init(int& val)
{
  val = 0;
}

void igl::XMLSerializableObject::Init(float& val)
{
  val = 0.0f;
}

void igl::XMLSerializableObject::Init(double& val)
{
  val = 0.000000000000000;
}

template<typename T>
void igl::XMLSerializableObject::Init(T*& obj)
{
  igl::XMLSerializable* object = static_cast<igl::XMLSerializable*>(obj);
  object->Init();  
}

template<typename T, int S>
void igl::XMLSerializableObject::Init(std::array<T,S>& obj)
{
  for(unsigned int i=0;i<obj.size();i++)
    Init(obj[i]);  
}

template<typename T0, typename T1>
void igl::XMLSerializableObject::Init(std::pair<T0,T1>& obj)
{
  Init(obj.first);
  Init(obj.second);
}

template<typename T>
void igl::XMLSerializableObject::Init(std::vector<T>& obj)
{
  obj.clear();
}

template<typename T, int R, int C>
void igl::XMLSerializableObject::Init(Eigen::Matrix<T,R,C>& obj)
{
  obj.setZero(obj.rows(),obj.cols());
}

template<typename T>
void  igl::XMLSerializableObject::Init(Eigen::SparseMatrix<T>& obj)
{
  obj.setZero();
}

bool igl::XMLSerializableObject::Serialize(char& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
{
  return setElementAttribute(obj,doc,element,name);
}

// overload function for char*, it interpreted as char array and can be used to handle strings
bool igl::XMLSerializableObject::Serialize(char*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
{
  return setElementAttribute(obj,doc,element,name);
}

bool igl::XMLSerializableObject::Serialize(std::string& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
{
  return setElementAttribute(obj,doc,element,name);
}

bool igl::XMLSerializableObject::Serialize(std::string*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
{
  return Serialize(*obj,doc,element,name);
}

bool igl::XMLSerializableObject::Serialize(bool& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
{
  return setElementAttribute(obj,doc,element,name);
}

bool igl::XMLSerializableObject::Serialize(bool*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
{
  return Serialize(*obj,doc,element,name);
}

bool igl::XMLSerializableObject::Serialize(unsigned int& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
{
  return setElementAttribute(obj,doc,element,name);
}

bool igl::XMLSerializableObject::Serialize(unsigned int*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
{
  return Serialize(*obj,doc,element,name);
}

bool igl::XMLSerializableObject::Serialize(int& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
{
  return setElementAttribute(obj,doc,element,name);
}

bool igl::XMLSerializableObject::Serialize(int*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
{
  return Serialize(*obj,doc,element,name);
}

bool igl::XMLSerializableObject::Serialize(float& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
{
  return setElementAttribute(obj,doc,element,name);
}

bool igl::XMLSerializableObject::Serialize(float*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
{
  return Serialize(*obj,doc,element,name);
}

bool igl::XMLSerializableObject::Serialize(double& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
{
  return setElementAttribute(obj,doc,element,name);
}

bool igl::XMLSerializableObject::Serialize(double*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
{
  return Serialize(*obj,doc,element,name);
}

template<typename T>
bool igl::XMLSerializableObject::Serialize(T& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
{ 
  return false;
}

template<typename T>
bool igl::XMLSerializableObject::Serialize(T*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
{
  // Serialize object implementing XMLSerializable interface
  igl::XMLSerializable* object = static_cast<igl::XMLSerializable*>(obj);
  
  tinyxml2::XMLElement* child = doc->NewElement(name.c_str());
  element->InsertEndChild(child);
  
  return object->Serialize(doc,child);
}

bool igl::XMLSerializableObject::Deserialize(char& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
{
  return getElementAttribute(obj,doc,element,name);      
}

// template specialisation for char*, it interpreted as char array and can be used to handle strings
bool igl::XMLSerializableObject::Deserialize(char*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
{
  return getElementAttribute(obj,doc,element,name);    
}

bool igl::XMLSerializableObject::Deserialize(std::string& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
{
  return getElementAttribute(obj,doc,element,name);
}

bool igl::XMLSerializableObject::Deserialize(std::string*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
{
  return Deserialize(*obj,doc,element,name);
}

bool igl::XMLSerializableObject::Deserialize(bool& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
{
  return getElementAttribute(obj,doc,element,name);
}

bool igl::XMLSerializableObject::Deserialize(bool*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
{
  return Deserialize(*obj,doc,element,name);
}

bool igl::XMLSerializableObject::Deserialize(unsigned int& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
{
  return getElementAttribute(obj,doc,element,name);
}

bool igl::XMLSerializableObject::Deserialize(unsigned int*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
{
  return Deserialize(*obj,doc,element,name);
}

bool igl::XMLSerializableObject::Deserialize(int& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
{
  return getElementAttribute(obj,doc,element,name);
}

bool igl::XMLSerializableObject::Deserialize(int*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
{
  return Deserialize(*obj,doc,element,name);
}

bool igl::XMLSerializableObject::Deserialize(float& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
{
  return getElementAttribute(obj,doc,element,name);
}

bool igl::XMLSerializableObject::Deserialize(float*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
{
  return Deserialize(*obj,doc,element,name);
}

bool igl::XMLSerializableObject::Deserialize(double& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
{
  return getElementAttribute(obj,doc,element,name);
}

bool igl::XMLSerializableObject::Deserialize(double*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
{
  return Deserialize(*obj,doc,element,name);
}

template<typename T>
bool igl::XMLSerializableObject::Deserialize(T& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
{
  return false;  
}

template<typename T>
bool igl::XMLSerializableObject::Deserialize(T*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
{
  obj = new T();
  igl::XMLSerializable* object = static_cast<igl::XMLSerializable*>(obj);
  
  const tinyxml2::XMLElement* child = element->FirstChildElement(name.c_str());

  object->Name = child->FirstChild()->Value();

  if(child != NULL)
  {  
    obj->Deserialize(doc,child);
  }
  else
  {
    obj->Init();
    return false;
  }

  return true;
}

template<typename T, size_t S>
bool igl::XMLSerializableObject::Serialize(std::array<T,S>& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
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
bool igl::XMLSerializableObject::Serialize(std::array<T,S>*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
{
  return Serialize(*obj,doc,element,name);
}

template<typename T, size_t S>
bool igl::XMLSerializableObject::Deserialize(std::array<T,S>& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
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
bool igl::XMLSerializableObject::Deserialize(std::array<T,S>*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
{
  obj = new std::array<T,S>();
  return Deserialize(*obj,doc,element,name);
}

template<typename T0, typename T1>
bool igl::XMLSerializableObject::Serialize(std::pair<T0,T1>& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
{
  bool res = true;

  tinyxml2::XMLElement* pair = doc->NewElement(name.c_str());
  element->InsertEndChild(pair);

  res &= Serialize(obj.first,doc,pair,"first");
  res &= Serialize(obj.second,doc,pair,"second");

  return res;
}

template<typename T0, typename T1>
bool igl::XMLSerializableObject::Serialize(std::pair<T0,T1>*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
{
  return Serialize(*obj,doc,element,name);
}

template<typename T0, typename T1>
bool igl::XMLSerializableObject::Deserialize(std::pair<T0,T1>& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
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
bool igl::XMLSerializableObject::Deserialize(std::pair<T0,T1>*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
{
  obj = new std::pair<T0,T1>();
  return Deserialize(*obj,doc,element,name);
}

template<typename T>
bool igl::XMLSerializableObject::Serialize(std::vector<T>& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
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
bool igl::XMLSerializableObject::Deserialize(std::vector<T>& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
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
bool igl::XMLSerializableObject::Serialize(std::vector<T>*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
{
  return Serialize(*obj,doc,element,name);
}

template<typename T>
bool igl::XMLSerializableObject::Deserialize(std::vector<T>*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
{
  obj = new std::vector<T>();
  return Deserialize(*obj,doc,element,name);
}

template<typename T, int R, int C>
bool igl::XMLSerializableObject::Serialize(Eigen::Matrix<T,R,C>& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
{
  tinyxml2::XMLElement* matrix = doc->NewElement(name.c_str());
  element->InsertEndChild(matrix);
  
  const unsigned int rows = obj.rows();
  const unsigned int cols = obj.cols();

  matrix->SetAttribute("rows",rows);
  matrix->SetAttribute("cols",cols);

  std::stringstream ms;
  ms << "\n";
  for(unsigned int r=0;r<rows;r++)
  {
    for(unsigned int c=0;c<cols;c++)
    {
      ms << obj(r,c) << ",";
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
bool igl::XMLSerializableObject::Serialize(Eigen::Matrix<T,R,C>*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
{
  return Serialize(*obj,doc,element,name);
}

template<typename T, int R, int C>
bool igl::XMLSerializableObject::Deserialize(Eigen::Matrix<T,R,C>& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
{
  const tinyxml2::XMLElement* child = element->FirstChildElement(name.c_str());
  if(child != NULL)
  {
    const unsigned int rows = child->UnsignedAttribute("rows");
    const unsigned int cols = child->UnsignedAttribute("cols");
    
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
  }
  else
  {
    Init(obj);
    return false;
  }
    
  return true;
}

template<typename T, int R, int C>
bool igl::XMLSerializableObject::Deserialize(Eigen::Matrix<T,R,C>*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
{
  obj = new Eigen::PlainObjectBase<T>();
  return Deserialize(*obj,doc,element,name);
}

template<typename T>
bool igl::XMLSerializableObject::Serialize(Eigen::SparseMatrix<T>& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
{
  tinyxml2::XMLElement* matrix = doc->NewElement(name.c_str());
  element->InsertEndChild(matrix);

  const unsigned int rows = obj.rows();
  const unsigned int cols = obj.cols();

  matrix->SetAttribute("rows",rows);
  matrix->SetAttribute("cols",cols);
    
  std::stringstream ms;
  ms << "\n";
  for (int k=0;k<obj.outerSize();++k)
  {
    for (typename Eigen::SparseMatrix<T>::InnerIterator it(obj,k);it;++it)
    {
      ms << it.row() << "," << it.col() << "," << it.value() << "\n";
    }
  }

  std::string mString = ms.str();
  if(mString.size() > 0)
    mString[mString.size()-1] = '\0';

  matrix->SetAttribute("matrix",mString.c_str());

  return true;
}

template<typename T>
bool igl::XMLSerializableObject::Serialize(Eigen::SparseMatrix<T>*& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
{
  return Serialize(*obj,doc,element,name);
}

template<typename T>
bool igl::XMLSerializableObject::Deserialize(Eigen::SparseMatrix<T>& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
{
  const tinyxml2::XMLElement* child = element->FirstChildElement(name.c_str());
  if(child != NULL)
  {
    const unsigned int rows = child->UnsignedAttribute("rows");
    const unsigned int cols = child->UnsignedAttribute("cols");
    
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

    std::vector<Eigen::Triplet<T>> triplets;
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
  }
  else
  {
    Init(obj);
    return false;
  }
    
  return true;
}

template<typename T>
bool igl::XMLSerializableObject::Deserialize(Eigen::SparseMatrix<T>*& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
{
  obj = new Eigen::SparseMatrix<T>();
  return Deserialize(*obj,doc,element,name);
}

template<typename T>
bool igl::XMLSerializableObject::setElementAttribute(T& obj, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
{
  tinyxml2::XMLElement* child = doc->NewElement(name.c_str());
  element->InsertEndChild(child);
  SetAttribute(child,"val",obj);
  return true;
}

template<typename T>
bool igl::XMLSerializableObject::getElementAttribute(T& obj, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element, const std::string& name)
{
// basic data type
  const tinyxml2::XMLElement* child = element->FirstChildElement(name.c_str());
  if(child != NULL)
  {  
    igl::XMLSerializableObject::GetAttribute(child->Attribute("val"),obj);
    return true;
  }
  else
  {
    Init(obj);
    return false;
  }
}

template<typename T>
igl::XMLSerializableInstance<T>::XMLSerializableInstance(T& obj, const std::string& name, const std::string group)
  : XMLSerializableObject(name, group), Object(obj)
{
  igl::XMLSerializableObject::Init(DefaultValue);
}

template<typename T>
igl::XMLSerializableInstance<T>::XMLSerializableInstance(T& obj, const std::string& name, const std::string group, T defaultValue)
  : XMLSerializableObject(name, group), Object(obj), DefaultValue(defaultValue)
{
}

template<typename T>
igl::XMLSerializableInstance<T>::~XMLSerializableInstance()
{
}

template<typename T>
void igl::XMLSerializableInstance<T>::Init()
{
  igl::XMLSerializableObject::Init(DefaultValue);
}

template<typename T>
bool igl::XMLSerializableInstance<T>::Serialize(tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element)
{
  return igl::XMLSerializableObject::Serialize(Object,doc,element,Name);
}

template<typename T>
bool igl::XMLSerializableInstance<T>::Deserialize(tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element)
{
  return igl::XMLSerializableObject::Deserialize(Object,doc,element,Name);
}

template<typename T>
bool igl::XMLSerializer::SaveObject(T& object, const char* filename)
{
  return SaveObject(object,"Object","Serialization",filename,true);
}

template<typename T>
bool igl::XMLSerializer::SaveObject(T& object, const std::string& objectName, const std::string& groupName, const char* filename, bool overwrite)
{
  bool result = true;
  XMLSerializer* serializer = new XMLSerializer(groupName);
  result &= serializer->Add(object,objectName);
  result &= serializer->Save(objectName,groupName,filename,overwrite);
  delete serializer;
  return result;
}

template<typename T>
bool igl::XMLSerializer::LoadObject(T& object, const char* filename)
{
  return LoadObject(object,"Object","Serialization",filename);
}

template<typename T>
bool igl::XMLSerializer::LoadObject(T& object, const std::string& objectName, const std::string& groupName, const char* filename)
{
  bool result = true;
  XMLSerializer* serializer = new XMLSerializer(groupName);
  result &= serializer->Add(object,objectName);
  result &= serializer->Load(objectName,groupName,filename);
  delete serializer;
  return result;
}
    
igl::XMLSerializer::XMLSerializer(const std::string& defaultGroup) 
{
  SetCurrentGroup(defaultGroup);
}

igl::XMLSerializer::~XMLSerializer()
{
  std::map<std::string,XMLSerializerGroup*>::iterator it;
  for (it=groups.begin();it!=groups.end();it++)
  {
    delete it->second->Objects;
    delete it->second;
  }
}

bool igl::XMLSerializer::Save(const char* filename, bool overwrite)
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

bool igl::XMLSerializer::SaveToXMLDoc(tinyxml2::XMLDocument* doc)
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

bool igl::XMLSerializer::Save(const std::string& groupName, const char* filename, bool overwrite)
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

bool igl::XMLSerializer::SaveToXMLDoc(const std::string& groupName, tinyxml2::XMLDocument* doc)
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

bool igl::XMLSerializer::Save(const std::string& objectName, const std::string& groupName, const char* filename, bool overwrite)
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

bool igl::XMLSerializer::SaveToXMLDoc(const std::string& objectName, const std::string& groupName, tinyxml2::XMLDocument* doc)
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

bool igl::XMLSerializer::SaveGroupToXMLElement(tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
{
  return SaveGroupToXMLElement(currentGroup->first,doc,element,name);
}

bool igl::XMLSerializer::SaveGroupToXMLElement(const std::string& groupName, tinyxml2::XMLDocument* doc, tinyxml2::XMLElement* element, const std::string& name)
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

bool igl::XMLSerializer::Load(const char* filename)
{
  tinyxml2::XMLDocument* doc = openDoc(filename);
  if(doc == NULL)
    return false;

  if(LoadFromXMLDoc(doc) == false)
    return false;

  delete doc;

  return true;
}

bool igl::XMLSerializer::LoadFromXMLDoc(tinyxml2::XMLDocument* doc)
{
  std::map<std::string,XMLSerializerGroup*>::iterator it;
  for (it=groups.begin();it!=groups.end();it++)
  {
    tinyxml2::XMLElement* element = doc->FirstChildElement(it->first.c_str());
    if(element == NULL)
    return false;
      
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

bool igl::XMLSerializer::Load(const std::string& groupName, const char* filename)
{
  tinyxml2::XMLDocument* doc = openDoc(filename);
  if(doc == NULL)
    return false;

  if(LoadFromXMLDoc(groupName, doc) == false)
    return false;

  delete doc;
    
  return true;
}

bool igl::XMLSerializer::LoadFromXMLDoc(const std::string& groupName, tinyxml2::XMLDocument* doc)
{
  std::string gn = groupName;
  EncodeXMLElementName(gn);

  std::map<std::string,XMLSerializerGroup*>::iterator it = groups.find(gn);
  if(it == groups.end())
    return false;

  tinyxml2::XMLElement* element = doc->FirstChildElement(it->first.c_str());
  if(element == NULL)
    return false;

  // Deserialize
  std::vector<XMLSerializable*>* groups = it->second->Objects;
  for(unsigned int i=0;i<groups->size();i++)
  {
    if(element == NULL || (*groups)[i]->Deserialize(doc,element) == false)
      (*groups)[i]->Init(); // Load default value;
  }
    
  return true;
}

bool igl::XMLSerializer::Load(const std::string& objectName, const std::string& groupName, const char* filename)
{
  tinyxml2::XMLDocument* doc = openDoc(filename);
  if(doc == NULL)
    return false;

  if(LoadFromXMLDoc(objectName,groupName,doc) == false)
    return false;

  delete doc;
    
  return true;
}

bool igl::XMLSerializer::LoadFromXMLDoc(const std::string& objectName, const std::string& groupName, tinyxml2::XMLDocument* doc)
{
  std::string gn = groupName;
  EncodeXMLElementName(gn);
    
  std::string on = objectName;
  EncodeXMLElementName(on);
    
  std::map<std::string,XMLSerializerGroup*>::iterator it = groups.find(gn);
  if(it == groups.end())
    return false;

  tinyxml2::XMLElement* element = doc->FirstChildElement(it->first.c_str());
  if(element == NULL)
    return false;

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

bool igl::XMLSerializer::LoadGroupFromXMLElement(const std::string& groupName, tinyxml2::XMLDocument* doc, const tinyxml2::XMLElement* element)
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

void igl::XMLSerializer::SetCurrentGroup(const std::string& group)
{
  currentGroup = setGetGroup(group);
}

std::string igl::XMLSerializer::GetCurrentGroup()
{
  return currentGroup->first;
}
  
template<typename T>
bool igl::XMLSerializer::Add(T& obj, const std::string& name)
{
  igl::XMLSerializable* object = static_cast<igl::XMLSerializable*>(obj);
  
  object->Name = name;
  return addObjectToGroup(object,currentGroup);  
}

bool igl::XMLSerializer::Add(char& obj, const std::string& name)
{
  return add(obj,name);
}

bool igl::XMLSerializer::Add(char*& obj, const std::string& name)
{
  return add(obj,name);
}

bool igl::XMLSerializer::Add(std::string& obj, const std::string& name)
{
  return add(obj,name);
}

bool igl::XMLSerializer::Add(bool& obj, const std::string& name)
{
  return add(obj,name);
}

bool igl::XMLSerializer::Add(unsigned int& obj, const std::string& name)
{
  return add(obj,name);
}

bool igl::XMLSerializer::Add(int& obj, const std::string& name)
{
  return add(obj,name);
}

bool igl::XMLSerializer::Add(float& obj, const std::string& name)
{
  return add(obj,name);
}

bool igl::XMLSerializer::Add(double& obj, const std::string& name)
{
  return add(obj,name);
}

template<typename T, size_t S>
bool igl::XMLSerializer::Add(std::array<T,S>& obj, const std::string& name)
{
  return add(obj,name);
}

template<typename T0, typename T1>
bool igl::XMLSerializer::Add(std::pair<T0,T1>& obj, const std::string& name)
{
  return add(obj,name);
}

template<typename T>
bool igl::XMLSerializer::Add(std::vector<T>& obj, const std::string& name)
{
  return add(obj,name);
}

template<typename T, int R, int C>
bool igl::XMLSerializer::Add(Eigen::Matrix<T,R,C>& obj, const std::string& name)
{
  return add(obj,name);
}

template<typename T>
bool igl::XMLSerializer::Add(Eigen::SparseMatrix<T>& obj, const std::string& name)
{
  return add(obj,name);
}

template<typename T>
bool igl::XMLSerializer::Add(T& object, const std::string& name, T defaultValue)
{
  return false;    
}

bool igl::XMLSerializer::Add(char& obj, const std::string& name, char defaultValue)
{
  return add(obj,name,defaultValue);
}

bool igl::XMLSerializer::Add(char*& obj, const std::string& name, char* defaultValue)
{
  return add(obj,name,defaultValue);
}

bool igl::XMLSerializer::Add(std::string& obj, const std::string& name, std::string defaultValue)
{
  return add(obj,name,defaultValue);
}

bool igl::XMLSerializer::Add(bool& obj, const std::string& name, bool defaultValue)
{
  return add(obj,name,defaultValue);
}

bool igl::XMLSerializer::Add(unsigned int& obj, const std::string& name, unsigned int defaultValue)
{
  return add(obj,name,defaultValue);
}

bool igl::XMLSerializer::Add(int& obj, const std::string& name, int defaultValue)
{
  return add(obj,name,defaultValue);
}

bool igl::XMLSerializer::Add(float& obj, const std::string& name, float defaultValue)
{
  return add(obj,name,defaultValue);
}

bool igl::XMLSerializer::Add(double& obj, const std::string& name, double defaultValue)
{
  return add(obj,name,defaultValue);
}

template<typename T>
bool igl::XMLSerializer::add(T& obj, const std::string& name)
{
  XMLSerializable* object = new igl::XMLSerializableInstance<T>(obj,name,currentGroup->first);
  return addObjectToGroup(object,currentGroup);
}

template<typename T>
bool igl::XMLSerializer::add(T& obj, const std::string& name, T defaultValue)
{
  igl::XMLSerializable* object = new igl::XMLSerializableInstance<T>(obj,name,currentGroup->first,defaultValue);
  return addObjectToGroup(object,currentGroup);
}
  
bool igl::XMLSerializer::addObjectToGroup(XMLSerializable* obj, const std::string& group)
{
  std::map<std::string,XMLSerializerGroup*>::iterator it = setGetGroup(group);
  return addObjectToGroup(obj, it);
}
  
bool igl::XMLSerializer::addObjectToGroup(XMLSerializable* object, std::map<std::string,XMLSerializerGroup*>::iterator it)
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

std::map<std::string,igl::XMLSerializerGroup*>::iterator igl::XMLSerializer::setGetGroup(const std::string& group)
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

tinyxml2::XMLDocument* igl::XMLSerializer::openDoc(const char* filename)
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

tinyxml2::XMLElement* igl::XMLSerializer::findAddGroup(tinyxml2::XMLDocument* doc, const char* groupName)
{
  tinyxml2::XMLElement* group = doc->FirstChildElement(groupName);
  if(group == NULL)
  {
    group = doc->NewElement(groupName);
    doc->InsertEndChild(group);
  }
  return group;
}

igl::XMLSerializerTest::XMLSerializerTest()
  : XMLSerialization("testObject")
{
  xmlSerializer->Add(testInt,"testInt");
  xmlSerializer->Add(testVector,"testVector");

  testInt = 10;
  
  testVector.push_back(1.0001f);
  testVector.push_back(2.0001f);
  testVector.push_back(3.0001f);
}
/*
Use the following lines to run the XMLSerializer tests:

igl::XMLSerializerTest* test = new igl::XMLSerializerTest();
bool success = test->Test();
delete test;
*/
bool igl::XMLSerializerTest::Test()
{
  // test vars 0
  char testChar0 = 'a';
  char* testCharArray0 = "testCharArray0";
  std::string testString0 = "testString0";
  bool testBool0 = true;
  unsigned int testUnsigned0 = 13;
  int testInt0 = 1000;
  float testFloat0 = 0.00001f;
  double testDouble0 = 0.01000000005;

  std::array<float,2> testArray0 = {{0.001f,1.001f}};

  std::pair<int,bool> testPair0;
  testPair0.first = 5;
  testPair0.second = true;

  std::vector<int> testVector0;
  testVector0.push_back(1);
  testVector0.push_back(2);
  testVector0.push_back(3);

  Eigen::MatrixXi testDenseMatrix0;
  testDenseMatrix0 = Eigen::Matrix3i::Identity();

  Eigen::SparseMatrix<double> testSparseMatrix0(3,3);
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.push_back(Eigen::Triplet<double>(0,0,1));
  triplets.push_back(Eigen::Triplet<double>(1,1,2));
  triplets.push_back(Eigen::Triplet<double>(2,2,3));
  testSparseMatrix0.setFromTriplets(triplets.begin(),triplets.end());

  igl::XMLSerializerTest* testObject0 = new igl::XMLSerializerTest();

  std::vector<std::pair<int,bool>*> testComplex10;
  testComplex10.push_back(&testPair0);
  testComplex10.push_back(&testPair0);

  std::vector<igl::XMLSerializerTest*> testComplex20;
  testComplex20.push_back(testObject0);
  testComplex20.push_back(testObject0);
  
  // test vars 1
  char testChar1 = 'b';
  char* testCharArray1 = "testCharArray1";
  std::string testString1 = "testString1";
  bool testBool1 = false;
  unsigned int testUnsigned1 = 12;
  int testInt1 = -1000;
  float testFloat1 = -0.00001f;
  double testDouble1 = -0.000000000001;

  std::array<float,2> testArray1 = {{-0.001f,-1.001f}};

  std::pair<int,bool> testPair1;
  testPair1.first = -5;
  testPair1.second = false;

  std::vector<int> testVector1;
  testVector1.push_back(-1);
  testVector1.push_back(-2);
  testVector1.push_back(-3);

  Eigen::MatrixXi testDenseMatrix1;
  testDenseMatrix1 = Eigen::Matrix3i::Identity();
  testDenseMatrix1 *= 2;

  Eigen::SparseMatrix<double> testSparseMatrix1(3,3);
  triplets.clear();
  triplets.push_back(Eigen::Triplet<double>(0,0,-1));
  triplets.push_back(Eigen::Triplet<double>(1,1,-2));
  triplets.push_back(Eigen::Triplet<double>(2,2,-3));
  testSparseMatrix1.setFromTriplets(triplets.begin(),triplets.end());

  igl::XMLSerializerTest* testObject1 = new igl::XMLSerializerTest();
  testObject1->testInt = 3;
  testObject1->testVector.clear();

  std::vector<std::pair<int,bool>*> testComplex11;
  testComplex11.push_back(&testPair1);
  testComplex11.push_back(&testPair1);

  std::vector<igl::XMLSerializerTest*> testComplex21;
  testComplex21.push_back(testObject1);
  testComplex21.push_back(testObject1);
  
  // test Serializer
  bool testResult = true;
  
  // test static functions
  
  // check static updating
  int testInt = 5;
  testResult &= igl::XMLSerializer::SaveObject(testInt0,"test.xml");
  testResult &= igl::XMLSerializer::SaveObject(testChar0,"testChar","test","test.xml",false);
  testResult &= igl::XMLSerializer::LoadObject(testInt,"test.xml");
  testResult &= testInt0 == testInt;

  // check static overwriting
  testInt = 3;
  testResult &= igl::XMLSerializer::SaveObject(testInt0,"test.xml");
  testResult &= igl::XMLSerializer::SaveObject(testChar0,"testChar","test","test.xml",true);
  testResult &= igl::XMLSerializer::LoadObject(testInt,"test.xml") == false;
  testResult &= testInt0 != testInt;

  // check static based type saving
  testResult &= igl::XMLSerializer::SaveObject(testChar0,"testChar","test0","test.xml",true);
  testResult &= igl::XMLSerializer::SaveObject(testCharArray0,"testCharArray","test0","test.xml",false);
  testResult &= igl::XMLSerializer::SaveObject(testString0,"testString","test0","test.xml",false);
  testResult &= igl::XMLSerializer::SaveObject(testBool0,"testBool","test0","test.xml",false);
  testResult &= igl::XMLSerializer::SaveObject(testUnsigned0,"testUnsigned","test0","test.xml",false);
  testResult &= igl::XMLSerializer::SaveObject(testInt0,"testInt","test0","test.xml",false);
  testResult &= igl::XMLSerializer::SaveObject(testFloat0,"testFloat","test0","test.xml",false);
  testResult &= igl::XMLSerializer::SaveObject(testDouble0,"testDouble","test0","test.xml",false);
  testResult &= igl::XMLSerializer::SaveObject(testArray0,"testArray","test1","test.xml",false);
  testResult &= igl::XMLSerializer::SaveObject(testPair0,"testPair","test1","test.xml",false);
  testResult &= igl::XMLSerializer::SaveObject(testVector0,"testVector","test1","test.xml",false);
  testResult &= igl::XMLSerializer::SaveObject(testDenseMatrix0,"testDense","test1","test.xml",false);
  testResult &= igl::XMLSerializer::SaveObject(testSparseMatrix0,"testSparse","test1","test.xml",false);
  testResult &= igl::XMLSerializer::SaveObject(testObject0,"testObject","test1","test.xml",false);
  testResult &= igl::XMLSerializer::SaveObject(testComplex10,"testComplex1","test1","test.xml",false);
  testResult &= igl::XMLSerializer::SaveObject(testComplex20,"testComplex2","test1","test.xml",false);

  testResult &= igl::XMLSerializer::LoadObject(testChar1,"testChar","test0","test.xml");
  testResult &= igl::XMLSerializer::LoadObject(testCharArray1,"testCharArray","test0","test.xml");
  testResult &= igl::XMLSerializer::LoadObject(testString1,"testString","test0","test.xml");
  testResult &= igl::XMLSerializer::LoadObject(testBool1,"testBool","test0","test.xml");
  testResult &= igl::XMLSerializer::LoadObject(testUnsigned1,"testUnsigned","test0","test.xml");
  testResult &= igl::XMLSerializer::LoadObject(testInt1,"testInt","test0","test.xml");
  testResult &= igl::XMLSerializer::LoadObject(testFloat1,"testFloat","test0","test.xml");
  testResult &= igl::XMLSerializer::LoadObject(testDouble1,"testDouble","test0","test.xml");
  testResult &= igl::XMLSerializer::LoadObject(testArray1,"testArray","test1","test.xml");
  testResult &= igl::XMLSerializer::LoadObject(testPair1,"testPair","test1","test.xml");
  testResult &= igl::XMLSerializer::LoadObject(testVector1,"testVector","test1","test.xml");
  testResult &= igl::XMLSerializer::LoadObject(testDenseMatrix1,"testDense","test1","test.xml");
  testResult &= igl::XMLSerializer::LoadObject(testSparseMatrix1,"testSparse","test1","test.xml");
  testResult &= igl::XMLSerializer::LoadObject(testObject1,"testObject","test1","test.xml");
  testResult &= igl::XMLSerializer::LoadObject(testComplex11,"testComplex1","test1","test.xml");
  testResult &= igl::XMLSerializer::LoadObject(testComplex21,"testComplex2","test1","test.xml");

  testResult &= testChar0 == testChar1;
  testResult &= strcmp(testCharArray0,testCharArray1) == 0;
  testResult &= testString0 == testString1;
  testResult &= testBool0 == testBool1;
  testResult &= testUnsigned0 == testUnsigned1;
  testResult &= testInt0 == testInt1;
  testResult &= testFloat0 == testFloat1;
  testResult &= testDouble0 == testDouble1;
  for(unsigned int i=0;i<testArray0.size();i++)
    testResult &= testArray0[i] == testArray1[i];
  testResult &= testPair0.first == testPair1.first;
  testResult &= testPair0.second == testPair1.second;
  testResult &= testVector0.size() == testVector1.size();
  for(unsigned int i=0;i<testVector0.size();i++)
    testResult &= testVector0[i] == testVector1[i];
  testResult &= (testDenseMatrix0-testDenseMatrix1).sum() == 0;
  testResult &= (testSparseMatrix0-testSparseMatrix1).norm() == 0;
  testResult &= testObject0->testInt == testObject1->testInt;
  testResult &= testObject0->testVector.size() == testObject1->testVector.size();
  for(unsigned int i=0;i<testObject0->testVector.size();i++)
    testResult &= testObject0->testVector[i] == testObject1->testVector[i];
  testResult &= testComplex10.size() == testComplex11.size();
  for(unsigned int i=0;i<testComplex10.size();i++)
  {
    testResult &= testComplex10[i]->first == testComplex11[0]->first;
    testResult &= testComplex10[i]->second == testComplex11[0]->second;
  }
  testResult &= testComplex20.size() == testComplex21.size();
  for(unsigned int i=0;i<testComplex20.size();i++)
  {
    testResult &= ((XMLSerializerTest*)testComplex20[i])->testInt == ((XMLSerializerTest*)testComplex21[i])->testInt;
    testResult &= ((XMLSerializerTest*)testComplex20[i])->testVector.size() == ((XMLSerializerTest*)testComplex21[i])->testVector.size();
    for(unsigned int j=0;j<((XMLSerializerTest*)testComplex20[i])->testVector.size();j++)
       testResult &= ((XMLSerializerTest*)testComplex20[i])->testVector[j] == ((XMLSerializerTest*)testComplex21[i])->testVector[j];
  }

  igl::XMLSerializer* s = new igl::XMLSerializer("test0");

  // clear file
  testResult &= s->Save("test.xml",true);

  testResult &= s->Add(testChar0,"testChar");
  testResult &= s->Add(testCharArray0,"testCharArray");
  testResult &= s->Add(testString0,"testString");
  testResult &= s->Add(testBool0,"testBool");
  testResult &= s->Add(testUnsigned0,"testUnsigned");
  testResult &= s->Add(testInt0,"testInt");
  testResult &= s->Add(testFloat0,"testFloat");
  testResult &= s->Add(testDouble0,"testDouble");

  s->SetCurrentGroup("test1");

  testResult &= s->Add(testArray0,"testArray");
  testResult &= s->Add(testPair0,"testPair");
  testResult &= s->Add(testVector0,"testVector");
  testResult &= s->Add(testDenseMatrix0,"testDenseMatrix");
  testResult &= s->Add(testSparseMatrix0,"testSparseMatrix");
  testResult &= s->Add(testObject0,"testObject");
  testResult &= s->Add(testComplex10,"testComplex1");
  testResult &= s->Add(testComplex20,"testComplex2");
  
  // Test single attribute save load
  testResult &= s->Save("testComplex2","test1","test.xml",true);
  testResult &= s->Load("testComplex2","test1","test.xml");
  
  testResult &= testComplex20.size() == testComplex21.size();
  for(unsigned int i=0;i<testComplex20.size();i++)
  {
    testResult &= ((XMLSerializerTest*)testComplex20[i])->testInt == ((XMLSerializerTest*)testComplex21[i])->testInt;
    testResult &= ((XMLSerializerTest*)testComplex20[i])->testVector.size() == ((XMLSerializerTest*)testComplex21[i])->testVector.size();
    for(unsigned int j=0;j<((XMLSerializerTest*)testComplex20[i])->testVector.size();j++)
        testResult &= ((XMLSerializerTest*)testComplex20[i])->testVector[j] == ((XMLSerializerTest*)testComplex21[i])->testVector[j];
  }

  // Test group save load without overriding
  testResult &= s->Save("test0","test.xml",false);
  testResult &= s->Load("test0","test.xml");
  testResult &= s->Load("testComplex2","test1","test.xml");

  testResult &= testChar0 == testChar1;
  testResult &= strcmp(testCharArray0,testCharArray1) == 0;
  testResult &= testString0 == testString1;
  testResult &= testBool0 == testBool1;
  testResult &= testUnsigned0 == testUnsigned1;
  testResult &= testInt0 == testInt1;
  testResult &= testFloat0 == testFloat1;
  testResult &= testDouble0 == testDouble1;

  // Test full serialization save load
  testResult &= s->Save("test.xml",true);
  testResult &= s->Load("test.xml");
  
  testResult &= testChar0 == testChar1;
  testResult &= strcmp(testCharArray0,testCharArray1) == 0;
  testResult &= testString0 == testString1;
  testResult &= testBool0 == testBool1;
  testResult &= testUnsigned0 == testUnsigned1;
  testResult &= testInt0 == testInt1;
  testResult &= testFloat0 == testFloat1;
  testResult &= testDouble0 == testDouble1;
  for(unsigned int i=0;i<testArray0.size();i++)
    testResult &= testArray0[i] == testArray1[i];
  testResult &= testPair0.first == testPair1.first;
  testResult &= testPair0.second == testPair1.second;
  testResult &= testVector0.size() == testVector1.size();
  for(unsigned int i=0;i<testVector0.size();i++)
    testResult &= testVector0[i] == testVector1[i];
  testResult &= (testDenseMatrix0-testDenseMatrix1).sum() == 0;
  testResult &= (testSparseMatrix0-testSparseMatrix1).norm() == 0;
  testResult &= testObject0->testInt == testObject1->testInt;
  testResult &= testObject0->testVector.size() == testObject1->testVector.size();
  for(unsigned int i=0;i<testObject0->testVector.size();i++)
    testResult &= testObject0->testVector[i] == testObject1->testVector[i];
  testResult &= testComplex10.size() == testComplex11.size();
  for(unsigned int i=0;i<testComplex10.size();i++)
  {
    testResult &= testComplex10[i]->first == testComplex11[0]->first;
    testResult &= testComplex10[i]->second == testComplex11[0]->second;
  }
  testResult &= testComplex20.size() == testComplex21.size();
  for(unsigned int i=0;i<testComplex20.size();i++)
  {
    testResult &= ((XMLSerializerTest*)testComplex20[i])->testInt == ((XMLSerializerTest*)testComplex21[i])->testInt;
    testResult &= ((XMLSerializerTest*)testComplex20[i])->testVector.size() == ((XMLSerializerTest*)testComplex21[i])->testVector.size();
    for(unsigned int j=0;j<((XMLSerializerTest*)testComplex20[i])->testVector.size();j++)
        testResult &= ((XMLSerializerTest*)testComplex20[i])->testVector[j] == ((XMLSerializerTest*)testComplex21[i])->testVector[j];
  }

  delete s;

  return testResult;
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization

#endif
