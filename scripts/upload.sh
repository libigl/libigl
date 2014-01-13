#!/bin/bash
make clean
zip -9 -r --exclude=@exclude.lst  libigl.zip ../libigl
scp libigl.zip $WEBSORKINE:www/htdocs-igl/projects/libigl/
cp *.html ~/Documents/IGL-website/projects/libigl/
cp file-formats/* ~/Documents/IGL-website/projects/libigl/file-formats/
