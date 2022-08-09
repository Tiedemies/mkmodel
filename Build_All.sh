#! /bin/bash
rm -rf Debug/
mkdir Debug
cmake -DCMAKE_BUILD_TYPE=Debug
make clean
make
make install
rm -rf Release/
mkdir Release
cmake -DCMAKE_BUILD_TYPE=Release
make clean
make
make install
