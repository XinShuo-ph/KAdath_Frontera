echo "The sequential version is not compatible with this branch (yet), please fix src_seq/do_newton.cpp first"
#cd build_release
#cmake -DPAR_VERSION=OFF -DCMAKE_BUILD_TYPE=Release .
#make 
#cd ..
#cd build_debug
#cmake -DPAR_VERSION=OFF -DCMAKE_BUILD_TYPE=Debug .
#make
#cd ..
#cd doc
#doxygen kadath.dox
#cd ..