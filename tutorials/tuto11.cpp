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
	
	// Fields for inner and outer BC
	Metric_tensor ginner (gmet) ; // Contains the initial value of gmet
	// No need to call std_base 
	Metric_tensor gouter (space, COV, basis) ;
	for (int i=1 ; i<=3 ; i++) {
		gouter.set(i,i) = 1. ;
		for (int j=i+1 ; j<=3 ; j++)
			gouter.set(i,j) = 0 ;
	}
	gouter.std_base() ;

	// General metric
	Metric_general met (gmet) ;

	// The system
	System_of_eqs syst (space, 1) ;
	met.set_system (syst, "g") ;
	
	// One can print the number of unknowns
	// This corresponds roughly to the number of coefficients of 6 components.
	cout << "Number of unknowns : " << syst.get_nbr_unknowns() << endl ;

	// Pass the constants
	syst.add_cst ("gin", ginner) ;
	syst.add_cst ("gout", gouter) ;

	// Now the equations
	syst.add_eq_bc (1, INNER_BC, "g_ij = gin_ij") ;
	syst.add_eq_inside (1, "lap(g_ij) = g_ik * R_j^k + g_jk * R_i^k") ;
	syst.add_eq_bc (1, OUTER_BC, "g_ij = gout_ij") ;

	// Compute the initial errors
	// This ensure that the code properly determines the number of conditions
	Array<double> errors (syst.sec_member()) ;
	
	// One can print the number of conditions arising from the equations
	cout << "Wrong number of conditions : " << syst.get_nbr_conditions() << endl ;

	// Constructor with two arguments.
	// First : the number of components to be considered.
	// Second : the rank (i.e. the number of indices)
	List_comp used (6, 2) ; // Here 6 components of a rank-2 tensor.

	// Passing of the components by hand
	used.set(0)->set(0) = 1 ; used.set(0)->set(1) = 1 ; // Component xx
	used.set(1)->set(0) = 1 ; used.set(1)->set(1) = 2 ; // Component xy
	used.set(2)->set(0) = 1 ; used.set(2)->set(1) = 3 ; // Component xz
	used.set(3)->set(0) = 2 ; used.set(3)->set(1) = 2 ; // Component yy
	used.set(4)->set(0) = 2 ; used.set(4)->set(1) = 3 ; // Component yz
	used.set(5)->set(0) = 3 ; used.set(5)->set(1) = 3 ; // Component zz
	
	// New system
	// General metric
	Metric_general mettwo (gmet) ;

	// The system
	System_of_eqs systtwo (space, 1) ;
	mettwo.set_system (systtwo, "g") ;
	
	// Pass the constants
	systtwo.add_cst ("gin", ginner) ;
	systtwo.add_cst ("gout", gouter) ;

	// Now the equations
	systtwo.add_eq_bc (1, INNER_BC, "g_ij = gin_ij") ;
	systtwo.add_eq_inside (1, "lap(g_ij) = g_ik * R_j^k + g_jk * R_i^k", used) ;
	systtwo.add_eq_bc (1, OUTER_BC, "g_ij = gout_ij") ;

	// Compute the initial errors
	// This ensure that the code properly determines the number of conditions
	Array<double> errorstwo (systtwo.sec_member()) ;
	
	// One can print the number of conditions arising from the equations
	cout << "Right number of conditions : " << systtwo.get_nbr_conditions() << endl ;
	
	return EXIT_SUCCESS ;
}

