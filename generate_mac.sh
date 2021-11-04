#!/bin/bash

BUILD_DIR="build"
if [ -d ${BUILD_DIR} ] ;then
	rm -rf ${BUILD_DIR}
fi

mkdir ${BUILD_DIR}
pushd ${BUILD_DIR} 

cmake -G Xcode -DCMAKE_CXX_FLAGS="-std=c++11"   ..

popd ${BUILD_DIR} 
