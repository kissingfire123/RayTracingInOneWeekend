#!/bin/bash

BUILD_DIR="build"
if [ -d ${BUILD_DIR} ] ;then
	rm -rf ${BUILD_DIR}
fi

mkdir ${BUILD_DIR}
pushd ${BUILD_DIR} 

cmake -G Xcode \
 -DCMAKE_BUILD_TYPE=Release \
 ..

popd ${BUILD_DIR} 
