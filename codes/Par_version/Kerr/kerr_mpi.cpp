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
    auto block_size = arg_parser.get_option_value<int>("-bsize","Sets block size of ScaLAPack cyclic decomposition",64);
    auto max_iterations = arg_parser.get_option_value<int>("-niter","Sets the maximum number of iteration for Newton-Raphson method. A null or negative sets this limit to infinity.",0);
    auto max_nb_omega = arg_parser.get_option_value<int>("-nomega","Sets the number of increments toward Omega to perform.", 40);
    auto nb_points = arg_parser.get_option_value<int>("-npts","Sets the number of collocation points (note that this value is constraint by the spectral method).",17);
    auto verbosity_level = arg_parser.get_option_value<int>("-v","Sets the verbosity level",1);
    auto gpu_enabled = arg_parser.get_option_value<bool>("-gpu","Enables GPU acceleration for linear solvers",false);
    auto tolerance = arg_parser.get_option_value<double>("-tol","Sets the tolerance value for Newton-Raphson's algorithms",default_tolerance);
    auto max_time = arg_parser.get_option_value<double>("-time","Sets a time limit for each Newton-Raphsons to perform (in seconds, negative values disable this chronometer). ",-1.);
    bool const show_help {arg_parser.find_option("-h","Display this help message.")};
    if(show_help) {
        if(rank == 0) arg_parser.display(std::cout);
        return 0;
    }

    Kerr_init kerr_init{nb_points.first};
    kerr_init.set_block_size(block_size.first);
    kerr_init.mpi_rank = rank;
    kerr_init.set_verbosity(verbosity_level.first);
    if(max_iterations.second) kerr_init.get_newton_solver().set_target_nb_iteration(max_iterations.first);
    if(gpu_enabled.second) kerr_init.get_newton_solver().set_enable_gpu(gpu_enabled.first);
    if(tolerance.second) kerr_init.get_newton_solver().set_target_error(tolerance.first);
    if(max_time.second) kerr_init.get_newton_solver().set_target_duration(max_time.first);

    // build all internal data.
    kerr_init.build_space_and_system();

    kerr_init.do_newton();

    kerr_init.finalize();
    kerr_init.profiling_log(std::cout);

    Kerr kerr{kerr_init};
    kerr.mpi_rank = rank;

    kerr.nbr_max_omega_val = max_nb_omega.second ? max_nb_omega.first : 40;
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

