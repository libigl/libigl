#!/bin/sh

echo "Installing Embree v2.4.0 ... "

INCLUDE_INSTALL_DIR=/usr/local/include
LIBRARY_INSTALL_DIR=/usr/local/lib

mkdir -p $INCLUDE_INSTALL_DIR
mkdir -p $LIBRARY_INSTALL_DIR

echo "  installing include files in" $INCLUDE_INSTALL_DIR
cp -r include/embree2 $INCLUDE_INSTALL_DIR

echo "  installing library files in" $LIBRARY_INSTALL_DIR
cp -r lib/x64/libembree.2.4.0.dylib $LIBRARY_INSTALL_DIR
ln -sf $LIBRARY_INSTALL_DIR/libembree.2.4.0.dylib $LIBRARY_INSTALL_DIR/libembree.2.dylib
ln -sf $LIBRARY_INSTALL_DIR/libembree.2.dylib     $LIBRARY_INSTALL_DIR/libembree.dylib
