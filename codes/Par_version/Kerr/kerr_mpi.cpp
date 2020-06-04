#include "kerr.hpp"
#include "mpi.h"
#include "magma_interface.hpp"


int main(int argc, char** argv) {
		
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
    auto max_nb_omega = arg_parser.get_option_value<int>("-nomega");
    auto nb_points = arg_parser.get_option_value<int>("-npts");
    if(!nb_points.second) nb_points.first = 17;

    Kerr_init kerr_init{nb_points.first};
    kerr_init.set_mpi_rank(rank);
    if(max_iterations.second) kerr_init.set_newton_max_iterations(max_iterations.first);

    // build all internal data.
    kerr_init.build_space_and_system();

    kerr_init.do_newton();

    kerr_init.finalize();
    kerr_init.profiling_log(std::cout);

    Kerr kerr{kerr_init};
    kerr.set_mpi_rank(rank);

    kerr.set_nbr_max_omega_val(max_nb_omega.second ? max_nb_omega.first : 40);
    kerr.reset_initial_guess();

    while(kerr.increment_omega())
    {
        //re-build the system with the new value of omega
        kerr.reset_system();
        // perform Newton-Raphson method for this value
        kerr.do_newton();
    }
    kerr.finalize();
    kerr.profiling_log(std::cout);
#ifdef ENABLE_GPU_USE
	if(rank==0)
	{
		TESTING_CHECK(magma_finalize());
	}
#endif
	MPI_Finalize() ;
	return EXIT_SUCCESS ;
}

