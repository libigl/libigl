// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2019 Jeremie Dumas <jeremie.dumas@ens-lyon.org>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_DEMANGLE_H
#define IGL_DEMANGLE_H

#include <string>
#ifdef _MSC_VER
#else // clang or gcc
#include <cxxabi.h>
#include <cstdlib>
#endif // clang or gcc branch of _MSC_VER

namespace igl
{

  inline void replaceAll(std::string& str, const std::string& from, const std::string& to) {
    if(from.empty())
        return;
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
    }
  }

  inline void replaceString(std::string& str) {
    replaceAll(str, "std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >", "std::string>");
    replaceAll(str, "std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >", "std::string");
  }

#ifdef _MSC_VER
  //! Demangles the type encoded in a string
  /*! @internal */
  inline std::string demangle( std::string const & name )
  { std::string str = name; replaceString(str); return str; }

  //! Gets the demangled name of a type
  /*! @internal */
  template <class T> inline
  std::string demangledName()
  { return typeid( T ).name(); }

#else // clang or gcc

  //! Demangles the type encoded in a string
  /*! @internal */
  inline std::string demangle(std::string mangledName)
  {
    int status = 0;
    char *demangledName = nullptr;
    std::size_t len;

    demangledName = abi::__cxa_demangle(mangledName.c_str(), 0, &len, &status);

    std::string retName(demangledName);
    free(demangledName);

    replaceString(retName);
    return retName;
  }

  //! Gets the demangled name of a type
  /*! @internal */
  template<class T> inline
  std::string demangledName()
  { return demangle(typeid(T).name()); }

#endif // clang or gcc branch of _MSC_VER

} // namespace libigl

#endif

