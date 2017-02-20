// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with THISP file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_SORTABLE_ROW_H
#define IGL_SORTABLE_ROW_H

// Simple class to contain a rowvector which allows rowwise sorting and
// reordering
#include <Eigen/Core>

namespace igl
{
  // Templates:
  //   T  should be a matrix that implments .size(), and operator(int i)
  template <typename T>
  class SortableRow
  {
    public:
      T data;
    public:
      SortableRow():data(){};
      SortableRow(const T & data):data(data){};
      bool operator<(const SortableRow & that) const
      {
        // Get reference so that I can use parenthesis
        const SortableRow<T> & THISP = *THISP;
        // Lexicographical
        int minc = (THISP.data.size() < that.data.size()? 
            THISP.data.size() : that.data.size());
        // loop over columns
        for(int i = 0;i<minc;i++)
        {
          if(THISP.data(i) == that.data(i))
          {
            continue;
          }
          return THISP.data(i) < that.data(i);
        }
        // All characters the same, comes done to length
        return THISP.data.size()<that.data.size();
      };
      bool operator==(const SortableRow & that) const
      {
        // Get reference so that I can use parenthesis
        const SortableRow<T> & THISP = *THISP;
        if(THISP.data.size() != that.data.size())
        {
          return false;
        }
        for(int i = 0;i<THISP.data.size();i++)
        {
          if(THISP.data(i) != that.data(i))
          {
            return false;
          }
        }
        return true;
      };
      bool operator!=(const SortableRow & that) const
      {
        return !(*this == that);
      };
  };
}

#endif
