cd build_release
cmake -DPAR_VERSION=ON -DCMAKE_BUILD_TYPE=Release .
make -j7
cd ..
cd build_debug
cmake -DPAR_VERSION=ON -DCMAKE_BUILD_TYPE=Debug .
make -j7
cd ..

