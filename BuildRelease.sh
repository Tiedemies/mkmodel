#! /bin/bash
rm -rf Release/
mkdir Release
cmake -DCMAKE_BUILD_TYPE=Release
make clean
make
make install
