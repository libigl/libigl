#!/bin/bash

echo "<html>"
echo "  <head><title>IGL LIB documentation</title></head>"
echo "  <body>"
echo "    <h1>IGL LIB</h1>"
echo "    <h2>Headers:</h2>"
echo "    <table>"
echo "      <tr><th>.h file</th><th>Functions</th></tr>"
# loop over all headers
for h in include/igl/*.h;
do
  b=`basename $h`
  if [ -e $h ]
  then
    printf "      <tr id='$b'><td>$b</td>"
  fi
  # portion of file inside namespace 
  html_nsp=`cat $h | perl -ne 'BEGIN{$p = 0;$\="<br>";} $o=$p;$p ^= $_=~"[{}]";print if $o && $p;'`
  #html_nsp=`echo -e $nsp | sed -e "s/</\&lt;/g" | sed -e "s/>/\&gt;/g" | tr '\n' '<br>'`
  #html_nsp=`echo -e $nsp | sed -e "s/</\&lt;/g" | sed -e "s/>/\&gt;/g" | sed -e "s/%/%%/g"`
  printf "<td>$html_nsp</td>"
  # Try to find functions and corresponding comments
  echo "</tr>"
done

echo "  </body>"
echo "</html>"
