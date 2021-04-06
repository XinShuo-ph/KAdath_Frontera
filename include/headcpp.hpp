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

#include "config.h"
#include <type_traits>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <typeinfo>
#include <cmath>

#ifdef USE_CXX_STANDARD_17_OR_HIGHER
#define CXX_17_ATTRIBUTES(...) [[__VA_ARGS__]]
#else
#define CXX_17_ATTRIBUTES(...)
#endif

#ifdef ENABLE_ASSERTS
#include <cassert>
#else
#undef assert
#define assert(c) ((void)0)
#endif

#include "stdlib.h"
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

    /**
     * Template alias to select the (usually) optimal access type, either to pass an argument
     * or to return an object. \c cutoff_factor allows to let bigger objects be passed or
     * returned by value instead of reference to const.
     */
    template<typename T,std::size_t cutoff_factor=1> using optimal_access_type =
        typename std::conditional<(sizeof(T) <= cutoff_factor*sizeof(T*)),T,T const &>::type;


#ifdef PAR_VERSION
  constexpr Computational_model default_computational_model = Computational_model::mpi_parallel;
#else //#ifdef PAR_VERSION
  constexpr Computational_model default_computational_model = Computational_model::sequential;
#endif //#ifdef PAR_VERSION
}
#endif
