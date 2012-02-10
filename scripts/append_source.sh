#!/bin/bash

# Replace #include "filename.cpp" with #include "../../source/filename.cpp"
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
  new=`cat "$filename.h" | sed -e "s/# *include \"$filename.cpp\"/#  include \"..\/..\/source\/$filename.cpp\"/"`
  echo "$new" > $i
done

