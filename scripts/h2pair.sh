#!/bin/bash

# http://tipstricks.itmatrix.eu/?p=305
tac () {
awk '1 { last = NR; line[last] = $0; } END { for (i = last; i > 0; i--) { print line[i]; } }'
}

function replace_inline
{
  result="$1"
  # replace inlines with IGL_INLINEs, first those that begin lines
  result=`echo "$result" | sed -e 's/^inline/IGL_INLINE/g'`
  # then also those that begin words
  result=`echo "$result" | sed -e 's/ inline/ IGL_INLINE/g'`
  echo "$result"
}

# Convert well organized .h file to a .h/.cpp pair
for i in $*
do
  # Only operate on .h files
  filename=$(basename $i)
  extension=${filename##*.}
  filename=${filename%.*}
  FILENAME=`echo $filename | tr '[a-z]' '[A-Z]'`
  if [ ! "$extension" = "h" ];
  then
    echo "Skipping $i because it is not a .h file"
    continue;
  fi

  if ! grep -q "^\/\/ Implementation *$" "$i"
  then
    echo "Skipping $i because it does not match ^\/\/ Implementation *$ "
    continue;
  fi

  if [[ `grep -c "^\#endif" "$i"` > 1 ]];
  then
    echo "Warning $i contains multple matches to ^#endif"
  fi

  before=`sed -n '/^\/\/ Implementation$/q;p' "$i"`;

  if ! echo "$before" | grep -q "^\#ifndef IGL_${FILENAME}_H"
  then
    echo "Skipping $i because it does not match ^#ifndef IGL_${FILENAME}_H "
    continue;
  fi

  if ! echo "$before" | grep -q "^\#define IGL_${FILENAME}_H"
  then
    echo "Skipping $i because it does not match ^#define IGL_${FILENAME}_H "
    continue;
  fi

  before=`replace_inline "$before"`
  # prepend #include "igl_inline.h"
  before=`echo "$before" | sed -e 's/^\(#define IGL_'${FILENAME}'_H\)/\1\'$'\n''#include \"igl_inline.h\"'/`
  after=`sed '1,/^\/\/ Implementation$/d' $i`;
  after=`replace_inline "$after"`
  # reverse file
  after=`echo "$after" | tac`
  # everything after first (last) endif 
  after=`echo "
$after" | sed '1,/endif/d'`;
  # reverse file
  after=`echo "$after" | tac`
  # append empty template code
  if grep -q "template" "$i"
  then
    after=`echo "$after

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
#endif"` 
  fi
  echo "$before

#ifndef IGL_STATIC_LIBRARY
#  include \"$filename.cpp\"
#endif

#endif"> "$filename".h
  echo "#include \"$filename.h\"

$after"> "$filename".cpp
done
