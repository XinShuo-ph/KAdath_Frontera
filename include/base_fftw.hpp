/*
    Copyright 2018 Ludwig Jens Papenfort

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

#ifndef __BASE_FFTW_HPP_
#define __BASE_FFTW_HPP_

#include <fftw3.h>
#include "profiled_object.hpp"
#include <unordered_map>

namespace Kadath {
// buffer and plan, keep them to save time
struct fftw_precomp_t
{
    int size;
    double* buffer;
    fftw_plan plan;

  fftw_precomp_t(int const n, fftw_r2r_kind FFTW_TRANSFORM) :
                                                              size{n},
                                                            buffer(fftw_alloc_real(n)),
                                                              plan(fftw_plan_r2r_1d(n, buffer, buffer, FFTW_TRANSFORM, FFTW_MEASURE)) {}
  ~fftw_precomp_t() {
    fftw_destroy_plan(plan);
    fftw_free(buffer);
  }

    inline void execute() {fftw_execute(plan);}
};

}
#endif
