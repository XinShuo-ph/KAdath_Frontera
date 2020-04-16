/*
    Copyright 2017 Philippe Grandclement

    This file is part of Kadath.

    Kadath is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Kadath is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Kadath.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __HEADCPP_HPP_
#define __HEADCPP_HPP_

#include <iomanip>
#include <iostream>
#include <fstream>
#include <typeinfo>
#include <cmath>
#include <cassert>
#include "stdlib.h"
#include "config.h"
#include "memory.hpp"

using namespace std ;

constexpr double PRECISION = 1e-14;

namespace Kadath {

  using std::cos ;
  using std::sin ;
  using std::sqrt ;
  using std::pow ;

  using std::log10 ;
  using std::log ;
  using std::exp ;
  using std::acos ;
  using std::asin ;

  using std::abs ;
  using std::fabs ;
  using std::max ;
  using std::min ;

  using std::cos;
  using std::sin;
  using std::tan ;
  using std::atan ;
  using std::atanh ; 
  using std::cosh;
  using std::sinh;

  /**
   * Set of enumerators used to select the computational model to use for the matrix-related computations (matrix
   * coefficient calculation, linear system solve, etc.).
   */
  enum class Computational_model
  {
      //! value for fully sequential computations.
      sequential,
      //! Value for sequantial matrix comptuation with gpu accelerated linear solver.
      gpu_sequential,
      //! value for fully MPI parallel computations.
      mpi_parallel,
      //! value for hybrid MPI / GPU computations.
      gpu_mpi_parallel
  };

  inline constexpr const char* computational_model_name(Computational_model c)
  {
      return c == Computational_model::sequential ? "sequential" :
                (c == Computational_model::mpi_parallel ? "mpi_parallel" : "gpu_mpi_parallel");
  }



#ifdef PAR_VERSION
#ifdef ENABLE_GPU_USE
    constexpr Computational_model default_computational_model = Computational_model::gpu_mpi_parallel;
#else // #ifdef ENABLE_GPU_USE
  constexpr Computational_model default_computational_model = Computational_model::mpi_parallel;
#endif//#ifdef ENABLE_GPU_USE
#else //#ifdef PAR_VERSION
#ifdef ENABLE_GPU_USE
  constexpr Computational_model default_computational_model = Computational_model::gpu_sequential;
#else //#ifdef ENABLE_GPU_USE
  constexpr Computational_model default_computational_model = Computational_model::sequential;
#endif //#ifdef ENABLE_GPU_USE
#endif //#ifdef PAR_VERSION
}
#endif
