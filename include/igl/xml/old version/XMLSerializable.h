// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
/* ---------------------------------------------------------------------------
 // XMLSerializable.h
 // Author: Christian Schüller <schuellchr@gmail.com>
 ------------------------------------------------------------------------------
 Inherit from this abstract class to have full control over the serialization
 of your user defined class.
 ----------------------------------------------------------------------------*/
#ifndef IGL_XML_SERIALIZABLE_H
#define IGL_XML_SERIALIZABLE_H

#include <iostream>
#include <tinyxml2.h>

namespace igl
{
  namespace
  {
    class XMLSerializable
    {
    public:
      /**
      * Default name of serializable object
      */
      std::string Name;

      /**
        * This function gets called if the objects were not found during deserialization.
        * Initialize your objects as you like.
        */
      virtual void Init() = 0;
      /**
        * Serialize your stuff within this function.
        * doc is the current serialization xml file. You can use SaveToXMLDoc to add your objects.
        */
      virtual void Serialize(tinyxml2::XMLDocument* doc,tinyxml2::XMLElement* element) = 0;

      /**
        * Deserialize your stuff within this function.
        * doc is the current serialization xml file. You can use LoadFromXMLDoc to read out your objects.
        */
      virtual void Deserialize(tinyxml2::XMLDocument* doc,const tinyxml2::XMLElement* element) = 0;
    };
  }
}

#endif
