#!/bin/bash

echo "Installing Embree v2.4.0 ... "

INCLUDE_INSTALL_DIR=/usr/local/include
LIBRARY_INSTALL_DIR=/usr/local/lib

mkdir -p $INCLUDE_INSTALL_DIR
mkdir -p $LIBRARY_INSTALL_DIR

echo "  installing include files in" $INCLUDE_INSTALL_DIR
cp -r include/embree2 $INCLUDE_INSTALL_DIR

echo "  installing Xeon library files in" $LIBRARY_INSTALL_DIR
cp -r lib/x64/libembree.so.2.4.0 $LIBRARY_INSTALL_DIR
ln -sf $LIBRARY_INSTALL_DIR/libembree.so.2.4.0 $LIBRARY_INSTALL_DIR/libembree.so.2
ln -sf $LIBRARY_INSTALL_DIR/libembree.so.2     $LIBRARY_INSTALL_DIR/libembree.so

echo "  installing Xeon Phi library files in" $LIBRARY_INSTALL_DIR
cp -r lib/x64/libembree_xeonphi.so.2.4.0 $LIBRARY_INSTALL_DIR
ln -sf $LIBRARY_INSTALL_DIR/libembree_xeonphi.so.2.4.0 $LIBRARY_INSTALL_DIR/libembree_xeonphi.so.2
ln -sf $LIBRARY_INSTALL_DIR/libembree_xeonphi.so.2     $LIBRARY_INSTALL_DIR/libembree_xeonphi.so
