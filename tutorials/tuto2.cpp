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
    cout << space << endl ;

	// Prints the number of points in the domain 0
cout << space.get_domain(0)->get_nbr_points() << endl ;

// Prints the nulber of coefficients in the domain 1 
cout << space.get_domain(1)->get_nbr_coefs() << endl ;

// Returns the radius in the domain 0
cout << space.get_domain(0)->get_radius() << endl ;

// Returns the X coordinate in the domain 1
cout << space.get_domain(1)->get_cart(1) << endl ;


	return EXIT_SUCCESS ;
}
