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
#include "codes_utilities.hpp"
#include "schwarz.hpp"
#include "mpi.h"
#include "magma_interface.hpp"


using namespace Kadath ;

int main(int argc,char** argv) {
    int rc = MPI_Init (&argc, &argv) ;
    int rank = 0 ;
    MPI_Comm_rank (MPI_COMM_WORLD, &rank) ;
#ifdef ENABLE_GPU_USE
    if(rank==0)
	{
		TESTING_CHECK(magma_init());
		magma_print_environment();
	}
#endif
    Arguments_parser arg_parser{argc,argv};
    auto max_iterations = arg_parser.get_option_value<int>("-iter");
    auto nb_points = arg_parser.get_option_value<int>("-npts");
    if(!nb_points.second) nb_points.first = 13;

    Schwarz schwarz_solver{nb_points.first};
    if(max_iterations.second) schwarz_solver.newton_max_iterations = max_iterations.first;

    schwarz_solver.build_space_and_system();

    schwarz_solver.do_newton();

    schwarz_solver.check_solution();

    if(rank==0) {
        cout << "Error max : " << schwarz_solver.get_error_l_infinity() << endl ;
        cout << "Error L2  : " << schwarz_solver.get_error_l_2() << endl;
    }

    schwarz_solver.finalize();

    if(rank==0) schwarz_solver.profiling_log(std::cout);
#ifdef ENABLE_GPU_USE
    if(rank==0)
	{
		TESTING_CHECK(magma_finalize());
	}
#endif
    MPI_Finalize() ;
    return EXIT_SUCCESS ;
}

