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
#include "kerr.hpp"

using namespace Kadath ;

int main(int argc, char** argv) {
    Arguments_parser arg_parser{argc,argv};
    auto max_iterations = arg_parser.get_option_value<int>("-iter");
    auto max_nb_omega = arg_parser.get_option_value<int>("-nomega");
    auto nb_points = arg_parser.get_option_value<int>("-npts");
    if(!nb_points.second) nb_points.first = 11;

    Kerr_init kerr_init{nb_points.first};
    if(max_iterations.second) kerr_init.newton_max_iterations = max_iterations.first;

    // build all internal data.
    kerr_init.build_space_and_system();

    kerr_init.do_newton();

    Kerr kerr{kerr_init};
    // this sequential demo is slow, so by default, we'll just make one iteration
    kerr.nbr_max_omega_val = max_nb_omega.second ? max_nb_omega.first : 1;
    // Computes initial guess from kerr_init data (the space is moved from kerr_init):
    kerr.reset_initial_guess();

    while(kerr.increment_omega())
    {
        //re-build the system with the new value of omega
        kerr.reset_system();
        // perform Newton-Raphson method for this value
        kerr.do_newton();
    }
    return EXIT_SUCCESS ;
}

