#!/bin/bash

# THIS would be much easier with a ruby script that parsed each .h/.cpp pair
# and created a <tr> row.
#

echo "<html>"
echo "  <head>"
echo "    <title>libigl auto-documentation</title>"
echo '    <link href="./style.css" rel="stylesheet" type="text/css">'
echo "  </head>"
echo "  <body class=article_body>"
echo "  <div class=article>"
echo "    <a href=.><img src=libigl-logo.jpg alt='igl logo' class=center></a>"
echo "    <h1>libigl</h1>"
echo "    Automatically generated documentation for <a href=.>libigl</a>."
echo "    <h2>Headers:</h2>"
echo "    <table>"
echo "      <tr><th>.h file</th><th>Functions</th></tr>"
# loop over all headers
odd="0"
for h in include/igl/*.h;
do
  b=`basename $h`
  # only consider files that exist as proper .h/.cpp files (Those that don't
  # are mostly utilitarian or poorly written)
  if [ -e "${h%.h}.cpp" ]
  then
    printf "      <tr id='$b' class=d$odd><td>$b</td>"
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
  fi
done

echo "    </table>"
echo "    <p>See also: <a href=tutorial.html>tutorial</a>, <a href=style_guidelines.html>style guidelines</a>, <a href=file-formats/index.html>file formats</a></p>"
echo "    <p>Automatically generated on `date` by scripts/doc.sh.</p>"
echo "  </div>"
echo "  </body>"
echo "</html>"
