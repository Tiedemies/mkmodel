#! /bin/bash
rm -rf Debug/
mkdir Debug
cmake -DCMAKE_BUILD_TYPE=Debug
make clean
make
make install
