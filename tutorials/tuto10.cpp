#include "kadath.hpp"

using namespace Kadath ;
int main() {

	// Example of a spherical space 
	// Consists of several 3D spherical domains

	// 3D :
	int dim = 3 ;

	// Number of points
	Dim_array res (dim) ;
	res.set(0) = 13 ; // Points in r 
	res.set(1) = 9 ; // Points in theta
	res.set(2) = 8 ; // Points in phi

	// Absolute coordinates of the center
	Point center (dim) ;
	for (int i=1 ; i<=dim ; i++)
		center.set(i) = 0 ;

	// Number of domains and boundaries :
	int ndom = 4 ;
	Array<double> bounds (ndom-1) ;
	bounds.set(0) = 1 ; bounds.set(1) = 2 ; bounds.set(2) = 4 ;

	// Chebyshev or Legendre :
	int type_coloc = CHEB_TYPE ;

	// Spherical space constructor
	Space_spheric space(type_coloc, center, res, bounds) ;
    
	// The tensorial basis
	Base_tensor basis (space, CARTESIAN_BASIS) ; 
	Metric_tensor gmet (space, COV, basis) ;

	// Affectation to flat metric everywhere
	for (int i=1 ; i<=3 ; i++)  {
     		gmet.set(i,i) = 1. ;
    		for (int j=i+1 ; j<=3 ; j++)
        		gmet.set(i,j) = 0 ;
	}

	// Non trivial stuff in domain 1
	for (int i=1 ; i<=3 ; i++)
     		for (int j=i ; j<=3 ; j++)
           		gmet.set(i,j).set_domain(1) = 
              			gmet(i,j)(1) + space.get_domain(1)->get_cart(i)*space.get_domain(1)->get_cart(j) ;

	// Standard base in that case
	gmet.std_base() ;
	//cout << gmet << endl ;

	// Construction of a constant metric
	// The numerical value is contained in a Metric_tensor that can be covariant or contravariant.
	Metric_const met (gmet) ;

	// Passing it to a system 
	System_of_eqs syst (space, 1) ;
	met.set_system (syst, "g") ;

	// One can get the connection coefficients (reserved word Gam)
	syst.add_def ("Christo_ij^k = Gam_ij^k") ;
	//cout << syst.give_val_def("Christo") << endl ;
	
	// One can also get the Ricci tensor (reserved word R)
	syst.add_def ("Ricci_ij = R_ij") ;
	//cout << syst.give_val_def("Ricci") << endl ;

	// Construction of the fields for inner and outer BC
	Metric_tensor ginner (gmet) ; // Contains the initial value of gmet
	// No need to call std_base 
	Metric_tensor gouter (space, COV, basis) ;
	for (int i=1 ; i<=3 ; i++) {
		gouter.set(i,i) = 1. ;
		for (int j=i+1 ; j<=3 ; j++)
			gouter.set(i,j) = 0 ;
	}
	gouter.std_base() ;
	
	// Construction of a unknown metric
	// The initial value numerical value is contained in a Metric_tensor.
	Metric_general metgen (gmet) ;

	// Passing it to a new system 
	System_of_eqs systtwo  (space, 1) ;
	metgen.set_system (systtwo, "g") ;

	// Pass the fields ginner and gouter as constants
	systtwo.add_cst ("gin", ginner) ;
	systtwo.add_cst ("gout", gouter) ;

	// Random equations
	systtwo.add_eq_bc (1, INNER_BC, "g_ij = gin_ij") ;
	systtwo.add_eq_inside (1, "lap(g_ij) = 0") ;
	systtwo.add_eq_bc (1, OUTER_BC, "g_ij = gout_ij") ;

	// Newton-Raphson solver
	double prec = 1e-8 ;
	double error ;
	bool endloop = false ;
	while (!endloop) {
           endloop = systtwo.do_newton(prec, error) ;
	}
	
	// The solution is contained in gmet
	// One can check that it is different from ginner that contains the initial value of gmet
	for (int i=1 ; i<=3 ; i++)
		for (int j=i ; j<=3 ; j++)
			cout << "Diff in comp " << i << ", " << j << " : " << diffmax(gmet(i,j)(1), ginner(i,j)(1)) << endl ;

	return EXIT_SUCCESS ;
}

