#!/bin/bash

# THIS would be much easier with a ruby script that parsed each .h/.cpp pair
# and created a <tr> row.
#

echo "<html>"
echo "  <head>"
echo "    <title>IGL LIB documentation</title>"
echo '    <link href="./style.css" rel="stylesheet" type="text/css">'
echo "  </head>"
echo "  <body>"
echo "    <h1>IGL LIB</h1>"
echo "    <h2>Headers:</h2>"
echo "    <table>"
echo "      <tr><th>.h file</th><th>Functions</th></tr>"
# loop over all headers
odd="0"
for h in include/igl/*.h;
do
  b=`basename $h`
  if [ -e $h ]
  then
    printf "      <tr id='$b' class=d$odd><td>$b</td>"
  fi
  # portion of file inside namespace 
  html_nsp=`cat $h | \
    perl -ne 'BEGIN{$p = 0} $o=$p;$p ^= $_=~"[{}]";print if $o && $p;' | \
    sed -e "s/</\&lt;/g" | sed -e "s/>/\&gt;/g" | sed -e "s/%/%%/g" | \
    sed -e "s/^\( *[^ \/].*\)$/<pre><code>\1<\/code><\/pre>/g" |  \
    sed -e ':a' -e 'N' -e '$!ba' -e 's/<\/code><\/pre>\n<pre><code>/\\\n/g' | \
    sed -e "s/^\(.*[^ ].*\)$/\1<br>/g"`;
  printf "<td>$html_nsp</td>"
  # Try to find functions and corresponding comments
  echo "</tr>"
  odd=`echo "($odd+1)%2" | bc`
done

echo "  </body>"
echo "</html>"
