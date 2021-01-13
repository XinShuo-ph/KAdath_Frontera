
echo

show_help()
{
  echo -e "Usage\n"
  echo "installer.sh [options]"
  echo
  echo "Configure, build and compile Kadath in the desired directory."
  echo
  echo "Options"
  echo "  --help, -h               = Print usage information and exit."
  echo "  -j <number of jobs>      = Set the number of concurrent jobs for the compilation process."
  echo "  --verbose,-v             = Toggles the build in verbose mode."
  echo "  --build-directory <path> = Sets the build directory for the Kadath library cmake project."
  echo "  --debug-mode, -g         = Compile the library in debug mode"
  echo "  --release-mode, -r       = Compile the library in release mode with level 3 compiler "
  echo "                             optimizations."
  echo "  --mpi-parallel           = Build the MPI-parallel version of the library."
  echo "  --sequential             = Build the sequential version of the library (default)."
  echo "  --with-mkl               = Build the MKL version of the library."
  echo "  --without-mkl            = Build the ScaLAPACK-based version of the library (default)."
  echo "  --arch <architecture>    = Preconfigure for a specific platform. Current possible values"
  echo "                             are \"ubuntu\" and \"mesopsl\"."
  echo "  --with-doxydoc           = Enables the creation of the Doxygen documentation (disabled by "
  echo "                             default)."
}

use_default_build_dir="true"
build_dir=build
build_mode="Release"
build_mode_d="release"
build_doc="false"
parallel_build="Off"
parallel_build_d="seq"
verbosity="Off"
compile_nb_proc=4
with_mkl="Off"
preconf_arch="none"
opt_cmake_build_type=""
opt_cmake_par_version=""
opt_cmake_mkl_version=""
opt_cmake_mpi_c_compiler=""
opt_cmake_mpi_cxx_compiler=""

while [ -n "$1" ]; do # while loop starts
	case "$1" in
	--help)
	  show_help
	  exit
	  ;;
	-help)
	  show_help
	  exit
	  ;;
	-h)
	  show_help
	  exit
	  ;;
	--arch)
	  preconf_arch="$2"
#	  echo "Preconfiguration for $preconf_arch"
	  shift
	  ;;
	--with-mkl)
	  with_mkl="On"
#	  echo "MKL version enabled."
	  ;;
	--without-mkl)
	  with_mkl="Off"
#	  echo "MKL version disabled."
	  ;;
	--mpi-parallel)
	  parallel_build="On"
	  parallel_build_d="mpi"
#	  echo "Parallel version"
	  ;;
	--sequential)
	  parallel_build="Off"
	  parallel_build_d="seq"
#	  echo "Sequential version"
	  ;;
	-j)
		compile_nb_proc="$2"
#		echo "Number of proc for compilation : $param"
		shift
		;;
	-v)
		verbosity="On"
#		echo "Verbose build enabled"
		;;
	--verbose)
		verbosity="On"
#		echo "Verbose build enabled"
		;;
	-d)
		build_dir="$2"
		use_default_build_dir="false"
#		echo "Build directory : $build_dir"
		shift
		;;
	--build-directory)
		build_dir="$2"
		use_default_build_dir="false"
#		echo "Build directory : $build_dir"
		shift
		;;
  -g)
    build_mode="Debug"
    build_mode_d="debug"
#    echo "Build mode : $build_mode"
    ;;
  --debug-mode)
    build_mode="Debug"
    build_mode_d="debug"
#    echo "Build mode : $build_mode"
    ;;
  -r)
    build_mode="Release"
    build_mode_d="release"
#    echo "Build mode : $build_mode"
    ;;
  --release-mode)
    build_mode="Release"
    build_mode_d="release"
#    echo "Build mode : $build_mode"
		;;
  --with-doxydoc)
    build_doc="true"
#    echo "Doxygen documentation enabled."
    ;;
	*) echo "Option $1 not recognized" ;;
	esac
	shift
done

cmake_local_file=Cmake/CMakeLocal.cmake
cmake_local_backup_file=Cmake/CMakeLocal.cmake.backup
if [ "$preconf_arch" != "none" ]; then
  if [ -f "$cmake_local_file" ]; then
    echo "$cmake_local_file file already exists."
    echo "This file will be overwritten for the ${preconf_arch} configuration."
    echo "A backup file will be available under $cmake_local_backup_file in the Cmake directory."
    if [ -f "$cmake_local_backup_file" ]; then
      echo
      echo "WARNING: A file with $cmake_local_backup_file has been found and will be overwritten."
      echo -n "         Continue ? Y/N... "
      read continue_preconf
      if [ "$continue_preconf" != "y" ] && [ "$continue_preconf" != "Y" ]; then
        return
      fi
    fi
  fi
  echo
fi

case $preconf_arch in
  "mesopsl")
    mv "$cmake_local_file" "$cmake_local_backup_file"
    cp "$cmake_local_file".mesopsl "$cmake_local_file"
    parallel_build=On
    with_mkl=On
    opt_cmake_mpi_c_compiler="-DMPI_C_COMPILER=/shared/apps/intel/compilers_and_libraries_2019/linux/mpi/intel64/bin/mpiicc"
    opt_cmake_mpi_cxx_compiler="-DMPI_CXX_COMPILER=/shared/apps/intel/compilers_and_libraries_2019/linux/mpi/intel64/bin/mpiicpc"
    ;;
  "ubuntu")
    mv "$cmake_local_file" "$cmake_local_backup_file"
    cp "$cmake_local_file".ubuntu14 "$cmake_local_file"
    ;;
  "none")
  ;;
  *)
    echo "Unknown preconfigured architecture : $preconf_arch"
    ;;
esac

opt_cmake_par_version="-DPAR_VERSION=${parallel_build}"
opt_cmake_mkl_version="-DMKL_VERSION=${with_mkl}"
opt_cmake_build_type="-DCMAKE_BUILD_TYPE=${build_mode}"

cmake_options="${opt_cmake_build_type} ${opt_cmake_par_version} ${opt_cmake_mkl_version} ${opt_cmake_mpi_c_compiler} ${opt_cmake_mpi_cxx_compiler}"

if [ $use_default_build_dir = "true" ];
then
  if [ "$preconf_arch" = "none" ]
  then
    build_dir="${build_dir}-${parallel_build_d}-${build_mode_d}"
  else
    build_dir="${build_dir}-${preconf_arch}-${build_mode_d}"
  fi
fi

echo "Kadath installer configuration :"
echo "  - build directory ..................: $build_dir"
echo "  - build mode .......................: $build_mode"
echo "  - mpi-parallel version .............: $parallel_build"
if [ $parallel_build = "On" ]
then
  if [ $with_mkl = "On" ]
  then
    echo "  - parallel linear algebra library ..: Intel Math Kernel Library"
  else
    echo "  - parallel linear algebra library ..: ScaLAPACK"
  fi
fi
echo "  - verbose build ....................: $verbosity"
echo "  - number of process for compilation : ${compile_nb_proc}"
echo
echo -n "Proceed ? (Y/N)... "
read proceed

if [ "$proceed" = "Y" ] || [ "$proceed" = "y" ]
then
  echo "mkdir $build_dir"
  mkdir "$build_dir"
  echo "cd $build_dir"
  cd "$build_dir" || return
  echo "cmake ${cmake_options} .."
  cmake "${cmake_options}" ..
  if [ "$verbosity" = "On" ]; then
    echo "make -j${compile_nb_proc} VERBOSE=1"
    make -j"${compile_nb_proc}" VERBOSE=1
  else
    echo "make -j${compile_nb_proc}"
    make -j"${compile_nb_proc}"
  fi
  if [ $build_doc = "true" ]
  then
    echo "make doc"
    make doc
  fi
  echo "cd.."
  cd ..
  echo
  echo "Kadath installation done."
elif [ "$proceed" = "N" ] || [ "$proceed" = "n" ]
then
  echo
  echo "Installation stopped by user, aborting... (use the -h option for any help)."
  return
else
  echo
  echo "Bad answer. Aborting... "
  return
fi
