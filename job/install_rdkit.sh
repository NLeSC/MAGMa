#!/bin/sh

wget https://github.com/rdkit/rdkit/archive/Release_2014_09_2.tar.gz
tar -zxf Release_2014_09_2.tar.gz
cd rdkit-Release_2014_09_2/
cd External/INCHI-API/
./download-inchi.sh
cd ../..
mkdir build
cd build
cmake -DRDK_BUILD_INCHI_SUPPORT=ON ..
make -j 4
make install
