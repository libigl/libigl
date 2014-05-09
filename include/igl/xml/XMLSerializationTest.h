// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
/* ---------------------------------------------------------------------------
// XMLSerializationTest.h
// Author: Christian Schüller <schuellchr@gmail.com>
------------------------------------------------------------------------------
 Used to demonstrates howto use the XMLSerialization class.
----------------------------------------------------------------------------*/
#ifndef IGL_XML_SERIALIZATION_TEST_H
#define IGL_XML_SERIALIZATION_TEST_H

#include <igl/xml/XMLSerialization.h>

namespace igl
{

  class XMLSerializationTest : public ::igl::XMLSerialization
  {
  public:
      
    int testInt;
    std::vector<float> testVector;
      
    XMLSerializationTest();
    void InitSerialization();
      
    bool Test();
  };

  XMLSerializerTest::XMLSerializerTest()
    : XMLSerialization("testObject")
  {
  }

  void XMLSerializerTest::InitSerialization()
  {
    xmlSerializer->Add(testInt,"testInt");
    xmlSerializer->Add(testVector,"testVector");
      
    testInt = 10;
      
    testVector.push_back(1.0001f);
    testVector.push_back(2.0001f);
    testVector.push_back(3.0001f);
  }
}

#endif
