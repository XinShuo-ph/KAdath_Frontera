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
    auto max_iterations = arg_parser.get_option_value<int>("-niter","Sets the maximum number of iteration for Newton-Raphson method. A null or negative sets this limit to infinity.",0);
    auto nb_points = arg_parser.get_option_value<int>("-npts","Sets the number of collocation points (note that this value is constraint by the spectral method).",13);
    auto verbosity_level = arg_parser.get_option_value<int>("-v","Sets the verbosity level",1);
    auto tolerance = arg_parser.get_option_value<double>("-tol","Sets the tolerance level for the approximation success checking",1.e-8);
    bool const show_help {arg_parser.find_option("-h","Display this help message.")};
    if(show_help) {
        if(rank == 0) arg_parser.display(std::cout);
        MPI_Barrier(MPI_COMM_WORLD);
        return 0;
    }

    Schwarz schwarz_solver{nb_points.first};
    schwarz_solver.tolerance = tolerance.first;
    schwarz_solver.set_verbosity(verbosity_level.first);
    if(max_iterations.second) schwarz_solver.newton_max_iterations = max_iterations.first;

    schwarz_solver.initialize();

    schwarz_solver.do_newton();

    schwarz_solver.check_solution();

    if(rank==0 && verbosity_level.first > 0) {
        double const err_linfty {schwarz_solver.get_error_l_infinity()};
        double const err_l2 {schwarz_solver.get_error_l_2()};
        bool const success_linfty{err_linfty <= tolerance.first};
        bool const success_l2{err_l2 <= tolerance.first};
        std::cout << "Error max :  " << err_linfty << (success_linfty ? " <= " : " > ")
                  << tolerance.first << (success_linfty ? " ==> OK " : "  :( ") << std::endl;
        cout << "Error L2  :  " << err_l2 << (success_l2 ? " <= " : " > ")
                  << tolerance.first << (success_l2 ? " ==> OK " : "  :( ") << std::endl;
        if(success_l2 && success_linfty) {
            std::cout << "==> SUCCESS !" << std::endl;
        } else if(success_l2 xor success_linfty) {
            std::cout << "==> half success..." << std::endl;
        }
        else std::cout << "==> bad approximation, maybe try a finer discretization or more iterations ?" << std::endl;
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

