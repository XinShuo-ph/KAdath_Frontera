/*
    Copyright 2020 sauliac

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

#include "solvers.hpp"

init_param_setter_buffer(error,Tolerance,1.e-12);
init_param_setter_buffer(nb_iteration,NbIter,-1);
init_param_setter_buffer(duration,ElapsedTime,0.);
init_param_setter_buffer(error_decrease,MinImprovement,0.);
init_boolean_setter_buffer(output, Output,true);

namespace Kadath {
    void Setter<output__>::apply(Solver * pthis, output__ const & x) {
        pthis->set_output (x.value);
    }
}