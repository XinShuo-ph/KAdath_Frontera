mkdir build_mpi_release
cd build_mpi_release
cmake -DPAR_VERSION=ON -DCMAKE_BUILD_TYPE=Release ..
make 
cd ..
cd build_mpi_debug
cmake -DPAR_VERSION=ON -DCMAKE_BUILD_TYPE=Debug ..
make
cd ..
cd doc
doxygen kadath.dox
cd ..
