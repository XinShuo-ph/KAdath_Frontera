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
    
	// A Cartesian tensorial basis
	Base_tensor basis_cart (space, CARTESIAN_BASIS) ;
	Metric_flat fmet_cart (space, basis_cart) ; // The associated flat metric

	// Same thing in orthonormal spherical basis
	Base_tensor basis_spher (space, SPHERICAL_BASIS) ;
	Metric_flat fmet_spher (space, basis_spher) ; // The associated flat metric

	// Construction of a vector
	Vector vecU (space, CON, basis_cart) ;
	
	// Sets (x,y,z) in the nucleus :
	for (int cmp=1 ; cmp<=3 ; cmp++)
    		vecU.set(cmp).set_domain(0) = space.get_domain(0)->get_cart(cmp) ;

	// Sets (x/r, y/r, z/r) for the other domains
	for (int d=1 ; d<ndom ; d++)
    		for (int cmp=1 ; cmp<=3 ; cmp++)
        		vecU.set(cmp).set_domain(d) = space.get_domain(d)->get_cart_surr(cmp) ;

	// Sets the appropriate spectral basis (here the standard one)
	vecU.std_base() ;
	// Pass in spherical basis
	vecU.change_basis_cart_to_spher() ;

	// The system only in the first shell
	System_of_eqs syst (space, 1) ;

	// Use the orthonormal spherical basis
	fmet_spher.set_system (syst, "f") ;
	// The string denotes the names by which the metric will be recognized by the system

	// Recover and print the inverse of the metric
	syst.add_def ("inv^ij = f^ij") ;
	//cout << syst.give_val_def("inv") << endl ;

	// Assume a contravariant vector vecU has been defined
	// Use it as a constant
	syst.add_cst ("U", vecU) ;
	// One can manipulate its indices :
	syst.add_def ("Ucov_i = U_i") ;
	// One can compute the scalar product (should be one)
	syst.add_def ("product = U_i * U^i") ;
	//cout << syst.give_val_def("product") << endl ;
	
	// The covariant derivative
	syst.add_def ("der_i^j = D_i U^j") ;
	cout << syst.give_val_def("der") << endl ;

	
	return EXIT_SUCCESS ;
}

