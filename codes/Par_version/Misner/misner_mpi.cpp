#include "kadath_bispheric.hpp"
#include "codes_utilities.hpp"
#include "mpi.h"
#include "magma_interface.hpp"

using namespace Kadath ;
int main(int argc,char ** argv) {

    int rc = MPI_Init(&argc, &argv) ;
    assert(rc == MPI_SUCCESS);
    int rank = 0 ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;
#ifdef ENABLE_GPU_USE
    if(rank==0)
	{
		TESTING_CHECK(magma_init());
		magma_print_environment();
	}
#endif
    Arguments_parser arg_parser{argc,argv};
    auto nb_points = arg_parser.get_option_value<int>("-npts","Sets the number of collocation points (note that this value is constraint by the spectral method).",7);
    auto verbosity_level = arg_parser.get_option_value<int>("-v","Sets the verbosity level",1);
    bool const show_help {arg_parser.find_option("-h","Display this help message.")};
    if(show_help) {
        if(rank == 0) arg_parser.display(std::cout);
        return 0;
    }


    // Number of points
    int nbr = nb_points.first ;

    double scale = 1. ;

    // Parameters of the bispherical :
    double ecart = 10. * scale;
    double r1 = 1. * scale ;
    double r2 = 1. * scale ;
    double rext = ecart*1.5 * scale ;

    // Espace :
    int type_base = CHEB_TYPE ;
    Space_bispheric space(type_base, ecart, r1, r2, rext, nbr) ;

    // Legendre or Chebyshev as one wishes...
    Scalar conf(space) ;
    conf = 1 ;
    conf.std_base() ;

    // System (outside the spheres)
    System_of_eqs syst (space, 2, 7) ;
    syst.add_var ("P", conf) ;
    syst.add_cst ("a", r1) ;
    syst.add_cst ("b", r2) ;

    // Equations :
    space.add_bc_sphere_one (syst, "dn(P) + 0.5 / a * P = 0") ;
    space.add_bc_sphere_two (syst, "dn(P) + 0.5 / b * P = 0") ;

    /// Bispherical part
    for (int d=2 ; d<=6 ; d++)
        syst.add_eq_inside (d, "lap(P) =0") ;
    space.add_matching (syst, "P") ;
    space.add_matching (syst, "dn(P)") ;

    // matching with outer domain
    for (int d=2 ; d<=6 ; d++)
        syst.add_eq_matching_import (d, OUTER_BC, "dn(P) = import(dn(P))") ;
    syst.add_eq_matching_import(7, INNER_BC, "P=import(P)") ;
    syst.add_eq_inside (7, "lap(P)=0") ;

    // Outer BC
    syst.add_eq_bc (7, OUTER_BC, "P=1") ;

    double conv ;
    bool endloop = false ;
    while (!endloop) {
        endloop = syst.do_newton(1e-8, conv);
    }

    syst.finalize_profiling();
    if(rank==0 && verbosity_level.first > 1) profiling_report(syst,std::cout);

#ifdef ENABLE_GPU_USE
    if(rank==0)
	{
		TESTING_CHECK(magma_finalize());
	}
#endif
    MPI_Finalize() ;
    return EXIT_SUCCESS ;
}