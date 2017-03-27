#include "kadath_bispheric.hpp"

using namespace Kadath ;
int main() {

	// Number of points
	int nbr = 7 ;

	// Parameters of the bispherical :
	double ecart = 10. ;
	double r1 = 1. ;
 	double r2 = 1. ;
	double rext = ecart*1.5 ;
	
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
	space.add_eq (syst, "lap(P) =0", "P", "dn(P)") ;
	space.add_bc_outer (syst, "P=1") ;

	double conv ;
	bool endloop = false ;
	int ite = 1 ;
	while (!endloop) {
		endloop = syst.do_newton(1e-8, conv) ;
		cout << "Newton iteration " << ite << " " << conv << endl ;
		ite++ ;
	}
	return EXIT_SUCCESS ;
}
