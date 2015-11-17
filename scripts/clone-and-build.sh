#!/bin/bash

# This script 
#  - clones the libigl github repo into a temporary directory
#  - updates submodules
#  - compiles static library
#  - compiles tutorial
# If any steps fail then an email is sent to a list of recipients (add yourself
# if you'd like to be alerted to failures, see below regarding spam filters).
#
# You can set up an automated nightly build-check at 2:30 AM using:
#
#     crontab -e
#
# and then add the line:
#
#     30 2 * * * source /Users/ajx/.profile; /usr/local/igl/libigl/scripts/clone-and-build.sh
#
# replace the path above with the **full path** to this file.
#



# Report that a command has failed and list its output
#
#     report_error command output
#
# Note: Most mail clients will mark these messages as spam even if the local
# sender (e.g. ajx@luftmatratze.local) is in your contact list. To force these
# message to be not marked as spam, create a filter.
#
# You can test your spam settings with
#
# echo "t3st 3mail" | mail -s "libigl test email" youremail@gmail.com
#
# This will almost certainly end up in your spam, but you can find your default
# email address as the sender and add a filter.
#
function report_error {
  subject="$(echo -e "libigl compilation failure: $1 \nContent-Type: text/html")" 
  recipients="alecjacobson@gmail.com"
  pre_style="style='background-color: #c3e0f0; overflow: auto; padding-left: 8px; padding-right: 8px; padding-top: 4px; padding-bottom: 4px; border: 1px solid #999;'";
  html_content="
<p>
The following command failed during the last nightly libigl build:
<pre $pre_style><code>$(echo -e "$1" | \
  sed \
  's/&/\&amp;/g; s/</\&lt;/g; s/>/\&gt;/g; s/"/\&quot;/g; s/'"'"'/\&#39;/g')
</code></pre>
</p>

<p>
The command produced the following standard output/error before failing:
<pre $pre_style><code>$(echo -e "$2" | \
  sed \
  's/&/\&amp;/g; s/</\&lt;/g; s/>/\&gt;/g; s/"/\&quot;/g; s/'"'"'/\&#39;/g')
</code></pre>
</p>

<p>
<a href=https://github.com/libigl/libigl/>libigl github</a> | <a href=https://github.com/libigl/libigl/commits/master>commits</a>
</p>
"
  echo -e "$html_content" | mail -s "$subject" $recipients
} 

# Runs the arguements as a command as usual, but if the command fails send an
# email using `report_error` and exit without continuing
#
#     guard command arg1 arg2 ...
#
function guard {
  command="$@"
  pwd
  echo "$command"
  if ! output=$($command 2>&1) ;
  then
    report_error "$command" "$output"
    echo "'$command' failed. Report sent."
    exit 1 
  fi
}

function grep_std_1
{
  (! grep -rI "std::__1" *)
}

set -o xtrace
# Clone repo
guard rm -rf /var/tmp/libigl
cd /var/tmp/
guard git clone --recursive git@github.com:libigl/libigl.git
cd libigl
# Build static library
mkdir lib
cd lib
# Redundant paths make it clear which command is failing
guard cmake -DCMAKE_BUILD_TYPE=Debug  ../optional/
guard make -C ../lib
# Build tutorial with default settings
mkdir ../tutorial/build
cd ../tutorial/build
guard cmake -DCMAKE_BUILD_TYPE=Debug  ../../tutorial/
guard make -C ../../tutorial/build
# Build tutorial with static library
cd ../
mkdir build-use-static
cd build-use-static
guard cmake -DCMAKE_BUILD_TYPE=Debug -DLIBIGL_USE_STATIC_LIBRARY=ON ../../tutorial/
guard make -C ../../tutorial/build-use-static
# check that no files contain `std::__1` often coming from templates copied
# from clang output. These will fail on gcc
cd ../../include/igl
guard grep_std_1
