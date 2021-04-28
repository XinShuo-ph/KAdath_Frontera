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
    
	// Assume a space is known (see tutorial 1)
	// Get the number of domains, if needed
	// int ndom = space.get_nbr_domains() ;

	// Constructor 
	System_of_eqs syst (space, 1, ndom-1) ;
	// Here one will solve the equations in the domains 1 to ndom-1 (so everywhere but in the domain 0)

	// Construct a scalar field
	Scalar field (space) ;
	field = 1 ;
	field.std_base() ;

	// The field is an unknown
	syst.add_var ("F", field) ;
	// The string is the name by which the field will be refereed to by the system.

	// Define a constant, begin the radius of the first domain
	Index pos (space.get_domain(1)->get_nbr_points()) ; // Index on the domain 1
	double rad = space.get_domain(1)->get_radius()(pos) ;
	cout << "The inner radius of the first shell is " << rad << endl ;

	// This number is a constant for the system
	syst.add_cst ("a", rad) ;
	// The string is the name by which the constant will be refereed to by the system.
	
	// The equation is the same in each domain
	for (int d=1 ; d<ndom ; d++) // Loop on the domains
        	syst.add_eq_inside (d, "Lap(F)=0") ; 
	// Lap is a reserved word of Kadath that means the flat Laplacian and we have defined F before.

	// Inner boundary condition
	syst.add_eq_bc (1, INNER_BC, "dn(F) + 0.5*F/a = 0") ;
	// It is a Robin BC
	// dn stands for the normal derivative with respect to the boundary (d/dr in this case)

	// Outer boundary condition
	syst.add_eq_bc (ndom-1, OUTER_BC, "F=1") ;
	// The field is one at spatial infinity

	// Loop on all the boundaries between the domains
	for (int d=1 ; d<ndom-1 ; d++) {
      	 	// Matching of the field between domain d and d+1
       	 syst.add_eq_matching (d, OUTER_BC, "F") ;
        	// Matching of the normal derivative of the field between domain d and d+1
        	syst.add_eq_matching (d, OUTER_BC, "dn(F)") ;
	}

	Array<double> errors (syst.check_equations()) ;
	// Verify the number of equations, being the size of errors
	int neq = errors.get_size(0) ;
	cout << "The system has " << neq << " equations." << endl ;
	// Print the errors, when they are "big"
	for (int n=0 ; n<neq ; n++)
       	if (fabs(errors(n)) > 1e-10)
               	 cout << "eq " << n << " ; error " << errors(n) << endl ;

	// Define the threshold
	double prec = 1e-8 ;
	// The current error
	double error ;
	// End of the iteration ?
	bool endloop = false ;

	while (!endloop) {
        	endloop = syst.do_newton(prec, error) ;
	}
	// At the first call the error is big and so the solver does one iteration and endloop is still false
	// At the second one, the error is smaller than the threshold and endloop is true thus ending the loop.
	// In that example the problem is linear and so converges in one iteration only.


	// Random point
	Point MM(3) ;
	MM.set(1) = 1.3873 ;
	MM.set(2) = -0.827 ;
	MM.set(3) = 0.982 ;
	// The numerical solution
	double numerical = field.val_point(MM) ;
	// The analytic one
	// Get the radius
	double rr = sqrt(MM(1)*MM(1)+MM(2)*MM(2)+MM(3)*MM(3)) ;
	double analytic = 1 + rad  / rr ;
	cout << "Numerical = " << numerical << endl ;
	cout << "Analytic  = " << analytic << endl ;
	cout << "Relative difference = " << fabs(numerical-analytic) / numerical << endl ;
	
	return EXIT_SUCCESS ;
}
