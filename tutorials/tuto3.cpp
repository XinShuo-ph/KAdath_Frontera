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
	// Puts 1/r in the other ones
	for (int d=1 ; d<ndom ; d++)
     		func.set_domain(d)  = 1/space.get_domain(d)->get_radius() ;

	// Prints the values at each collocation point
	cout << "Scalar function func " << endl ;
	cout << func << endl ;

	// Use the same scalar as before
	// Put r exp(-r*r/10.) in the outer domains
	for (int d=1 ; d<ndom ; d++)
    	   func.set_domain(d)  = space.get_domain(d)->get_radius()* exp(- space.get_domain(d)->get_radius()*space.get_domain(d)->get_radius()/10.) ;

	// Nans appear at infinty
	cout << "Nans at infinity" << endl ;
	cout << func(ndom-1) << endl ;

	// Sets the right value at infinity by hand
	func.set_val_inf(0.) ;

	// No more Nans
	cout << "Nans are removed"  << endl ;
	cout << func(ndom-1) << endl ;

	// Sets the standard base for scalar fields
	func.std_base() ;

	// Define a Point in a three dimensional space
	Point M (3) ;

	// Sets its coordinates
	M.set(1) = 2.93 ; // X coordinate
	M.set(2) = -3.19 ; // Y coordinate
	M.set(3) = 1.76 ; // Z coordinate

	// Prints the values of the scalar field at M
	cout << "The point M is " << M << endl ;
	cout << "Field value at M " << func.val_point(M) << endl ;


	// Compute the radial derivative
	Scalar derr (func.der_r()) ;
	cout << "Value of the radial derivative at M " << derr.val_point(M) << endl ;

	// Use the multiplication by r :
	Scalar multr (func.mult_r()) ;
	//Force the code to compute the values at the collocation points
	multr.coef_i() ;
	
	// No nans appear at infinity
	cout << "Value in the last domain" << endl ;
	cout << multr(ndom-1) << endl ;

	return EXIT_SUCCESS ;
}
