#!/bin/bash

if [ "$#" -ne 0 ]
then
  USAGE="Usage:
  scripts/license.sh

Checks and appends MPL2 license to all files.
";
  echo "$USAGE"
  exit 1
fi

# By default use Alec's email address.
HEADER="\
This file is part of libigl, a simple c++ geometry processing library.

Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>

This Source Code Form is subject to the terms of the Mozilla Public License 
v. 2.0. If a copy of the MPL was not distributed with this file, You can 
obtain one at http://mozilla.org/MPL/2.0/.";


HEADER_TMPFILE="libigl.header.tmpfile";
CPPHEADER=`echo "$HEADER" | sed -e "s/^/\/\/ /"`
#echo -e "$(echo "$HEADER" | sed -e "s/^/\/\/ /")\n$(cat test.h)" > test.h

FIRST_LINE_OF_HEADER=`echo -e "$HEADER" | head -1`

# Files from third-party sources
EXCEPT="tga|MCTables|marching_cubes"

SRC_FILES=`find include/igl -name \*.h -print -o -name \*.cpp -print`
SRC_FILES=`echo "$SRC_FILES" | egrep -v $EXCEPT`
for s in $SRC_FILES
do
  if ! head -1 $s | grep -F "$FIRST_LINE_OF_HEADER" >/dev/null
  then
    echo "Appending license header to $s..."
    echo "$CPPHEADER" |cat - $s > /tmp/out && mv /tmp/out $s
  fi;
done
