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
#include <typeinfo>
#include <type_traits> // std::is_same<T1,T2>

#ifdef _MSC_VER
#else // clang or gcc
#include <cxxabi.h>
#include <cstdlib>
#include <algorithm> // std::remove
#endif // clang or gcc branch of _MSC_VER

namespace igl
{
  inline void replace_all(std::string& str, const std::string& from, const std::string& to)
  {
    if(from.empty())
    {
      return;
    }
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos)
    {
      str.replace(start_pos, from.length(), to);
      start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
    }
  }

#ifdef _MSC_VER
  //! Demangles the type encoded in a string
  /*! @internal */
  inline std::string demangle(std::string mangledName)
  {
    std::string prefix = "class "; // MSVC does "class myClass" see: https://docs.microsoft.com/en-us/cpp/cpp/typeid-operator
    // remove "class " prefixes
    replace_all(mangledName,prefix,"");
    return mangledName;
  }

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
    // MSVC does not add spaces after comma.
    // gcc and clang do sth. like this: "class<T1, T2, T3>"
    // thus remove blanks
    replace_all(retName,", ",",");
    return retName;
  }

#endif // clang or gcc branch of _MSC_VER

  //! Gets the demangled name of a type
  /*! @internal */
  template<class T> inline
  std::string demangled_name()
  {
    if (std::is_same<T, std::string>::value)
    {
      return std::string("std::string");
    }
    return demangle(typeid(T).name());
  }

} // namespace libigl

#endif

