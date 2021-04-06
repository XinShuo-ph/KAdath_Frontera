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
	res.set(1) = 5 ; // Points in theta
	res.set(2) = 4 ; // Points in phi

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
    
	// Assume a spherical space has been previously constructed
	// Recover the number of domains if need be
	//int ndom = space.get_nbr_domains() ;

	Scalar func (space) ;

	// Puts 1-r^2 in domain 0
	func.set_domain(0) = 1 - pow(space.get_domain(0)->get_radius(),2) ;
	// Put 1/r in the other ones
	for (int d=1 ; d<ndom ; d++)
     		func.set_domain(d)  = 1/space.get_domain(d)->get_radius() ;

	// Assume a scalar field func has been previously computed (see tutorial 3)
	// One sets the basis to the standard one
	func.std_base() ;

	// Print the basis in domain 1 (the first shell)
	cout << "Basis in domain 1" << endl ;
	cout << func(1).get_base() << endl ;

	// The result is read as follows
	// Variable 2 is phi, only one base, labelled 4, which is a full sequence of sines and cosines 

	// Variable 1 is theta, an array of basis, with values 6 and 10.
	// It corresponds to either even cosines or odd sines, depending on the index of phi.

	// Variable 0 is r, a two dimensional array with values 1
	// It corresponds to Chebyshev polynomials, not matter what the indices of the angles are.

	// Print the basis in domain 0 (the nucleus)
	cout << "Basis in domain 0" << endl ;
	cout << func(0).get_base() << endl ;

	// Combination of fields
	Scalar f1 (2*func + func*func) ;
	//Base is known and is the standard one ;
	cout << "Base of 2*f + f^2 in the nucleus" << endl ;
	cout << f1(0).get_base() << endl ;
	// It is the standard base

	// Taking the derivative
	Scalar der (func.der_r()) ;
	cout << "Base of of f' in the nucleus" << endl ;
	cout << der(0).get_base() << endl ;
	// It is not the standard base ; change in the r parity

	// Operator cosine
	Scalar test (cos(func)) ;
	cout << "Base for cos(f) is not defined" << endl ;
	cout << test(0).get_base() << endl ;

	// Define a scalar field being z in the nucleus and z/r^2 everywhere else
	Scalar odd (space) ;
	odd.set_domain(0) = space.get_domain(0)->get_cart(3) ; // z coordinate
	for (int d=1 ; d<ndom ; d++)
    		 odd.set_domain(d) = space.get_domain(d)->get_cart(3) 
		/ space.get_domain(d)->get_radius() / space.get_domain(d)->get_radius() ;
	// Sets the right value at infinity
	odd.set_val_inf(0.) ;
	// Try to use std_base() 
	odd.std_base() ;
	// Compare the numerical value to the analytical one using val_point
	Point M(3) ;
	M.set(1) = 2.32 ;
	M.set(2) = 1.32 ;
	M.set(3) = 1.98 ;
	double rr = sqrt(M(1)*M(1)+M(2)*M(2)+M(3)*M(3)) ; // Get the radius ; check it not in the nucleus
	cout << "Analytical value  = " << M(3) / rr / rr << endl ;
	cout << "Wrong numerical value = " << odd.val_point(M) << endl ;

	// Put the right spectral basis :
	odd.std_anti_base() ;
	odd.set_in_conf() ; // One needs to kill the previously computed (wrong) coefficients.
	cout << "True numerical value = " << odd.val_point(M) << endl ;

	return EXIT_SUCCESS ;
}
