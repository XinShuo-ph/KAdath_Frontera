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
    
	// Assume a spherical space is known and ndom its number of domains
	Base_tensor basis_cart (space, CARTESIAN_BASIS) ;
	// Constructor with the same tensorial basis everywhere

	// Sets an orthonormal spherical basis in the first shell
	Base_tensor basis_mixed (space, CARTESIAN_BASIS) ;
	basis_mixed.set_basis(1) = SPHERICAL_BASIS ;
	cout << "The mixed tensorial basis is " << endl ;
	cout << basis_mixed << endl ;

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

	// Prints the spectral basis in the nucleus
	for (int cmp=1 ; cmp<=3 ; cmp++) {
    		cout << "Basis of cmp " << cmp << " in the nucleus: " << endl ;
    		cout << vecU(cmp)(0).get_base() << endl ;
    	}
    	
    	// Change the tensorial basis
	vecU.change_basis_cart_to_spher() ;

	// Prints the new tensorial basis
	cout << vecU.get_basis() << endl ;


	// Check that the spectral basis has changed
	cout << "Spectral basis of component 3 in the nucleus: " << endl ;
	cout << vecU(3)(0).get_base() << endl ;
	
	// Construct a rank-two tensor, both indices being contravariant
	Tensor Tone (space, 2, CON, basis_cart) ;
	cout << Tone << endl ; /// The values of the components are not defined at this point

	// A rank 3 tensor with COV, CON, CON indices
	int valence = 3 ;
	Array<int> type_indices (valence) ;
	type_indices.set(0) = COV ; 
	type_indices.set(1) = CON ;
	type_indices.set(2) = CON ;
	Tensor Ttwo (space, valence, type_indices, basis_cart) ;
	
	// Direct accessor
	Ttwo.set(1,2,2) = 1. ; // Sets the component (x,y,y) to one

	// Accessor using an Index
	Index indices (Ttwo) ;
	indices.set(0) = 2 ; indices.set(1) = 2 ; indices.set(2) = 0 ; // Corresponds to (z,z,x)
	Ttwo.set(indices) = 2. ;

	cout << Ttwo << endl ; // Only two components are defined
	return EXIT_SUCCESS ;
}

