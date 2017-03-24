cd build_release
cmake -DPAR_VERSION=OFF -DCMAKE_BUILD_TYPE=Release .
make 
cd ..
cd build_debug
cmake -DPAR_VERSION=OFF -DCMAKE_BUILD_TYPE=Debug .
make
cd ..
cd doc
doxygen kadath.dox
cd ..
