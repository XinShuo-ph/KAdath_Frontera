mkdir build_seq_release
cd build_seq_release
cmake -DPAR_VERSION=OFF -DCMAKE_BUILD_TYPE=Release ..
make 
cd ..
mkdir build_seq_debug
cd build_seq_debug
cmake -DPAR_VERSION=OFF -DCMAKE_BUILD_TYPE=Debug ..
make
cd ..
cd doc
doxygen kadath.dox
cd ..
