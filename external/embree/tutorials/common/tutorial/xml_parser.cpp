// ======================================================================== //
// Copyright 2009-2013 Intel Corporation                                    //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#include "xml_parser.h"

#include <fstream>

namespace embree
{
  //////////////////////////////////////////////////////////////////////////////
  ///                           XML Input
  //////////////////////////////////////////////////////////////////////////////

  /*! parse a list of XML comments */
  void parseComments(Ref<Stream<Token> >& cin)
  {
    while (cin->peek() == Token::Sym("<!--")) {
      cin->drop();
      while (cin->peek() != Token::Sym("-->")) cin->drop();
      cin->drop();
    }
  }

  /*! parse XML parameter */
  void parseParm(Ref<Stream<Token> >& cin, std::map<std::string,std::string>& parms)
  {
    std::string name = cin->get().Identifier();
    if (cin->get() != Token::Sym("=")) THROW_RUNTIME_ERROR(cin->unget().Location().str()+": symbol \"=\" expected");
    parms[name] = cin->get().String();
  }

  /*! parse XML header */
  Ref<XML> parseHeader(Ref<Stream<Token> >& cin)
  {
    Ref<XML> xml = new XML;
    if (cin->get() != Token::Sym("<?")) THROW_RUNTIME_ERROR(cin->unget().Location().str()+": wrong XML header");
    xml->name = cin->get().Identifier();
    parseComments(cin);
    while (cin->peek() != Token::Sym("?>")) {
      parseParm(cin,xml->parms);
      parseComments(cin);
    }
    cin->drop();
    return xml;
  }

  /*! parse XML tag */
  Ref<XML> parseXML(Ref<Stream<Token> >& cin)
  {
    Ref<XML> xml = new XML;
    xml->loc = cin->peek().Location();

    /* parse tag opening */
    if (cin->get() != Token::Sym("<")) THROW_RUNTIME_ERROR(cin->unget().Location().str()+": tag expected");

    xml->name = cin->get().Identifier();
    parseComments(cin);
    while (cin->peek() != Token::Sym("/>") && cin->peek() != Token::Sym(">")) {
      parseParm(cin,xml->parms);
      parseComments(cin);
    }
    if (cin->peek() == Token::Sym("/>")) {
      cin->drop();
      return xml;
    }
    cin->drop();

    /* parse body token list */
    parseComments(cin);
    while (cin->peek() != Token::Sym("<") && cin->peek() != Token::Sym("</")) {
      xml->body.push_back(cin->get());
      parseComments(cin);
    }

    /* the body also contains children */
    if (cin->peek() == Token::Sym("<")) {
      while (cin->peek() != Token::Sym("</")) {
        xml->children.push_back(parseXML(cin));
        parseComments(cin);
      }
    }

    /* parse tag closing */
    if (cin->get() != Token::Sym("</")    ) THROW_RUNTIME_ERROR(cin->unget().Location().str()+": symbol \"</\" expected");
    if (cin->get() != Token::Id(xml->name)) THROW_RUNTIME_ERROR(cin->unget().Location().str()+": closing "+xml->name+" expected");
    if (cin->get() != Token::Sym(">")     ) THROW_RUNTIME_ERROR(cin->unget().Location().str()+": symbol \">\" expected");

    return xml;
  }

  /* load XML from token stream */
  Ref<XML> parseXML(Ref<Stream<int> > chars, bool hasHeader = true, bool hasTail = false)
  {
    /* create lexer for XML file */
    std::vector<std::string> symbols;
    symbols.push_back("<!--");
    symbols.push_back("-->");
    symbols.push_back("<?");
    symbols.push_back("?>");
    symbols.push_back("</");
    symbols.push_back("/>");
    symbols.push_back("<");
    symbols.push_back(">");
    symbols.push_back("=");
    Ref<Stream<Token> > cin = new TokenStream(chars,TokenStream::alpha + TokenStream::ALPHA + "_", TokenStream::separators, symbols);

    if (hasHeader) parseHeader(cin);
    parseComments(cin);
    Ref<XML> xml = parseXML(cin);
    parseComments(cin);

    if (!hasTail)
      if (cin->peek() != Token::Eof()) THROW_RUNTIME_ERROR(cin->peek().Location().str()+": end of file expected");

    return xml;
  }

  /*! load XML file from stream */
  std::istream& operator>>(std::istream& cin, Ref<XML>& xml) {
    xml = parseXML(new StdStream(cin),false,true);
    return cin;
  }

  /*! load XML file from disk */
  Ref<XML> parseXML(const FileName& fileName) {
    return parseXML(new FileStream(fileName),true,false);
  }


  //////////////////////////////////////////////////////////////////////////////
  ///                           XML Output
  //////////////////////////////////////////////////////////////////////////////


  /* indent to some hierarchy level using spaces */
  void indent(std::ostream& cout, size_t depth) {
    for (size_t i=0; i<2*depth; i++) cout << " ";
  }

  /* store XML to a stream */
  std::ostream& emitXML(std::ostream& cout, const Ref<XML>& xml, size_t depth = 0)
  {
    /* print header */
    if (depth == 0) cout << "<?xml version=\"1.0\"?>" << std::endl << std::endl;

    /* print tag opening */
    indent(cout,depth); cout << "<" << xml->name;
    for (std::map<std::string,std::string>::const_iterator i=xml->parms.begin(); i!=xml->parms.end(); i++)
      cout << " " << i->first << "=" << "\"" << i->second << "\"";
    if (xml->children.size() == 0 && xml->body.size() == 0) {
      cout << "/>" << std::endl;
      return cout;
    }
    cout << ">";

    bool compact = xml->body.size() < 16 && xml->children.size() == 0;
    if (!compact) cout << std::endl;

    /* print token list */
    if (xml->body.size()) {
      if (!compact) indent(cout,depth+1);
      for (size_t i=0; i<xml->body.size(); i++)
        cout << xml->body[i] << (i!=xml->body.size()-1?" ":"");
      if (!compact) cout << std::endl;
    }

    /* print children */
    for (size_t i=0; i<xml->children.size(); i++)
      emitXML(cout,xml->children[i],depth+1);

    /* print tag closing */
    if (!compact) indent(cout,depth);
    return cout << "</" << xml->name << ">" << std::endl;
  }

  /* store XML to stream */
  std::ostream& operator<<(std::ostream& cout, const Ref<XML>& xml) {
    return emitXML(cout,xml);
  }

  /*! store XML to disk */
  void emitXML(const FileName& fileName, const Ref<XML>& xml)
  {
    std::ofstream cout(fileName.c_str());
    if (!cout.is_open()) THROW_RUNTIME_ERROR("cannot open file " + fileName.str() + " for writing");
    emitXML(cout,xml);
    cout.close();
  }
}
