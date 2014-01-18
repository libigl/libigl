#!/bin/bash
set -e
make clean
zip -9 -r --exclude=@exclude.lst  libigl.zip ../libigl
scp libigl.zip $WEBSORKINE:www/htdocs-igl/projects/libigl/
cp *.html ~/Documents/IGL-website/projects/libigl/
cp file-formats/* ~/Documents/IGL-website/projects/libigl/file-formats/
VERSION=`grep -v "^\#" VERSION.txt | tr -d '\n'`
cd ~/Documents/IGL-website/projects/libigl/
svn up
sed -ie "s/RELEASE_HISTORY\.txt>[0-9\.]*<\/a>/RELEASE_HISTORY\.txt>$VERSION<\/a>/" index.php
svn ci -m "update igl version: $VERSION"
ssh websorkine@web-login.inf.ethz.ch "svn update www/htdocs-igl/"
