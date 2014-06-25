#!/bin/bash

if [ "$#" -lt 2 ]
then
  USAGE="Usage:
  scripts/name_change.sh old_name new_name";
  echo "$USAGE"
  exit 1
fi

old="$1"
new="$2"

old_cpp=`find . -name "$old.cpp"`
if [ -z "$old_cpp" ]
then
  echo "Error: $old.cpp does not exist."
  exit 1
fi

old_h=`find . -name "$old.h"`
if [ -z "$old_h" ]
then
  echo "Error: $old.h does not exist."
  exit 1
fi

grep -rI "$old" *

read -r -p "Are you sure? [y/N] " response
if [[ ! $response =~ ^([yY][eE][sS]|[yY])$ ]]
then
  exit 1
fi

