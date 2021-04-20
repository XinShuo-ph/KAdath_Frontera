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
    
	// Assume a spherical space is known and ndom its number of domains
	Base_tensor basis_cart (space, CARTESIAN_BASIS) ;
	
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

	// Construct a rank-two tensor, both indices being covariant
	Tensor tenT (space, 2, COV, basis_cart) ;
	// Sets some values (antisymmetric)
	for (int i=1 ; i<=3 ; i++) {
		tenT.set(i,i) = 0. ;
		for (int j=i+1 ; j<=3 ; j++) {
			tenT.set(i,j).set_domain(0) = space.get_domain(0)->get_cart(i)*space.get_domain(0)->get_cart(j) ;
			for (int d=1 ; d<ndom ; d++)
				tenT.set(i,j).set_domain(d) = space.get_domain(d)->get_cart_surr(i)*space.get_domain(d)->get_cart_surr(j) ;
			tenT.set(j,i) = -tenT(i,j) ;
		}
	}
	
	// Sets the appropriate spectral basis (here the standard one)
	tenT.std_base() ;
	
	// The system of equation in domain 1 only
	System_of_eqs syst (space, 1) ;
	// Pass V as an unknown field
	syst.add_var ("U", vecU) ;
	// Pass T as a constant
	syst.add_cst ("T", tenT) ;

	// Tensorial product 
	syst.add_def ("TV_ij^k = T_ij * U^k") ;

	// Tensorial product with contraction
	syst.add_def ("Contract_i = T_ij * U^j") ;

	// Transposition 
	syst.add_def ("Transpose_ij = T_ji") ;

	// Index manipulation produces an arror as no metric has been defined
	//syst.add_def ("Ucov_i = U_i") ;  

	// The various definitions can be accessed by give_val_def()
	// For instance for printing
	// cout << syst.give_val_def("Contract") << endl ;

	// One will solve only in domain 1 using syst
	// Define a tensor for the inner BC (using twice the initial value of vecU)

	Vector vecInner (2*vecU) ;
	syst.add_cst ("I", vecInner) ;

	//Random inner BC
	syst.add_eq_bc (1, INNER_BC, "U^k= I^k") ; 
	// Bulk equation
	syst.add_eq_inside (1, "lap(U^k)= 0") ; 
	// Outer BC
	syst.add_eq_bc (1, OUTER_BC, "U^k= 0") ; 

	// Newton-Raphson solver
	double prec = 1e-8 ;
	double error ;
	bool endloop = false ;
	while (!endloop) {
           endloop = syst.do_newton(prec, error) ;
	}
	// The problem is linear and converges in one step.
	// One can print the solution
	for (int cmp=1 ; cmp<=3 ; cmp++) {
     		cout << "Component " << cmp << endl ;
    		cout << vecU(cmp)(1) << endl ;
	}
	// At this point vecU is no longer vecInner/2 (the first one has been modified by the solver).

	return EXIT_SUCCESS ;
}

