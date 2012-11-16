#!/bin/bash

for i in {1..20}
do
  make clean &>/dev/null
  printf "$i  "
  t=`(time make -j$i lib&>/dev/null) 2>&1 | grep real | sed -e "s/real[^0-9]*//g"`
  echo "$t"
done
