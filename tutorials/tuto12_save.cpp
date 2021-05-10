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
    
	Scalar field (space) ;
	field.set_domain(0) = 1. ;
	for (int d=1 ; d<ndom ; d++)
		field.set_domain(d) = 1./space.get_domain(d)->get_radius() ;
	field.std_base() ;
	
	// Construction of a vector
	Base_tensor basis (space, CARTESIAN_BASIS) ;
	Vector vecU (space, CON, basis) ;
	
	// Sets (x,y,z) in the nucleus :
	for (int cmp=1 ; cmp<=3 ; cmp++)
    		vecU.set(cmp).set_domain(0) = space.get_domain(0)->get_cart(cmp) ;

	// Sets (x/r, y/r, z/r) for the other domains
	for (int d=1 ; d<ndom ; d++)
    		for (int cmp=1 ; cmp<=3 ; cmp++)
        		vecU.set(cmp).set_domain(d) = space.get_domain(d)->get_cart_surr(cmp) ;

	// Sets the appropriate spectral basis (here the standard one)
	vecU.std_base() ;

	// Opening the file in write mode.
	FILE* fout = fopen ("file.dat", "w") ;

	// Two values to be saved
	int valint = 2 ;
	double valdouble = 1.31732 ;

	// Write using fwrite
	fwrite_be (&valint, sizeof(int), 1, fout) ;
	fwrite_be (&valdouble, sizeof(double), 1, fout) ;

	// Assume a space, a scalar and a vector have been previously defined
	// They are saved by invoking save
	space.save(fout) ;
	field.save(fout) ;
	vecU.save(fout) ;

	// Dont forget to close the file in the end
	fclose(fout) ;

	return EXIT_SUCCESS ;
}

