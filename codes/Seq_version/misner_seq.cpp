#include "kadath_bispheric.hpp"

using namespace Kadath ;
int main() {

	// Number of points
	int nbr = 7 ;

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
	int ite = 1 ;
	while (!endloop) {
		endloop = syst.do_newton(1e-8, conv) ;
		cout << "Newton iteration " << ite << " " << conv << endl ;
		ite++ ;
	}

    syst.finalize_profiling();
    profiling_report(syst,std::cout);
	return EXIT_SUCCESS ;
}
