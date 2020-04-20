/*
    Copyright 2017 Philippe Grandclement

    This file is part of Kadath.

    Kadath is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Kadath is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Kadath.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "headcpp.hpp"

#include "utilities.hpp"
#include "spheric.hpp"
#include "point.hpp"
#include "val_domain.hpp"
namespace Kadath {
// Standard constructor
Domain_shell_surr::Domain_shell_surr (int num, int ttype, double rint, double rext, const Point& cr, const Dim_array& nbr) : Domain_shell(num, ttype, rint, rext, cr, nbr)  {
     alpha = (1./rext-1./rint)/2. ;
     beta = (1./rext+1./rint)/2. ; 
     assert (nbr.get_ndim()==3) ;
     assert (cr.get_ndim()==3) ;    // Affectation de type_point :
     do_coloc() ;

}


// Constructor by copy
Domain_shell_surr::Domain_shell_surr (const Domain_shell_surr& so) : Domain_shell(so) {

    // Update the alpha and beta
    double rint = (beta - alpha) ;
    double rext = (beta+alpha) ;
    alpha = (1./rext-1./rint)/2. ;
    beta = (1./rext+1./rint)/2. ; 
}

Domain_shell_surr::Domain_shell_surr (int num, FILE* fd) : Domain_shell(num, fd) {
	do_coloc() ;
}

// Destructor
Domain_shell_surr::~Domain_shell_surr() {}

void Domain_shell_surr::save (FILE* fd) const {
	nbr_points.save(fd) ;
	nbr_coefs.save(fd) ;
	fwrite_be (&ndim, sizeof(int), 1, fd) ;
	fwrite_be (&type_base, sizeof(int), 1, fd) ;	
	center.save(fd) ;
	fwrite_be (&alpha, sizeof(double), 1, fd) ;	
	fwrite_be (&beta, sizeof(double), 1, fd) ;
}

ostream& operator<< (ostream& o, const Domain_shell_surr& so) {
  o << "Shell surr" << endl ;
  o << 1./(so.beta - so.alpha)  << " < r < " << 1./(so.beta+so.alpha) << endl ;
  o << "Center  = " << so.center << endl ;
  o << "Nbr pts = " << so.nbr_points << endl ;
  o << endl ;
  return o ;
}



Val_domain Domain_shell_surr::der_normal (const Val_domain& so, int bound) const {

	if ((bound!=OUTER_BC) && (bound!=INNER_BC)) {
	    cerr << "Unknown boundary case in Domain_shell_surr::der_normal" << endl ;
	    abort() ;
	}
	return der_r (so) ;
}

// Computes the Cartesian coordinates
void Domain_shell_surr::do_absol () const  {
	for (int i=0 ; i<3 ; i++)
	   assert (coloc[i] != 0x0) ;
	for (int i=0 ; i<3 ; i++)
	   assert (absol[i] == 0x0) ;
	for (int i=0 ; i<3 ; i++) {
	   absol[i] = new Val_domain(this) ;
	   absol[i]->allocate_conf() ;
	   }
	Index index (nbr_points) ;

	do  {
		absol[0]->set(index) = 1./(alpha* ((*coloc[0])(index(0))) +beta)*
		 sin((*coloc[1])(index(1)))*cos((*coloc[2])(index(2))) + center(1) ;
		absol[1]->set(index) = 1./(alpha* ((*coloc[0])(index(0))) +beta) *
		 sin((*coloc[1])(index(1)))*sin((*coloc[2])(index(2))) + center(2) ;
		absol[2]->set(index) = 1./(alpha* ((*coloc[0])(index(0))) + beta) * cos((*coloc[1])(index(1))) + center(3) ;
			}
	while (index.inc())  ;
}

// Computes the radius
void Domain_shell_surr::do_radius () const  {
	for (int i=0 ; i<3 ; i++)
	   assert (coloc[i] != 0x0) ;
	assert (radius == 0x0) ;
	radius = new Val_domain(this) ;
	radius->allocate_conf() ;
	Index index (nbr_points) ;
	do 
		radius->set(index) = 1./(alpha* ((*coloc[0])(index(0))) + beta) ;
	while (index.inc())  ;
	radius->std_base() ;
}

// Computes the Cartesian coordinates
void Domain_shell_surr::do_cart () const  {
	for (int i=0 ; i<3 ; i++)
	   assert (coloc[i] != 0x0) ;
	for (int i=0 ; i<3 ; i++)
	   assert (cart[i] == 0x0) ;
	for (int i=0 ; i<3 ; i++) {
	   cart[i] = new Val_domain(this) ;
	   cart[i]->allocate_conf() ;
	   }
	Index index (nbr_points) ;

	do  {
		cart[0]->set(index) = 1./(alpha* ((*coloc[0])(index(0))) +beta)*
		 sin((*coloc[1])(index(1)))*cos((*coloc[2])(index(2))) + center(1) ;
		cart[1]->set(index) = 1./(alpha* ((*coloc[0])(index(0))) +beta) *
		 sin((*coloc[1])(index(1)))*sin((*coloc[2])(index(2))) + center(2) ;
		cart[2]->set(index) = 1./(alpha* ((*coloc[0])(index(0))) + beta) * cos((*coloc[1])(index(1))) + center(3) ;
			}
	while (index.inc())  ;
}

// Check if a point is inside this domain
bool Domain_shell_surr::is_in (const Point& xx, double prec) const {

	assert (xx.get_ndim()==3) ;
	
	double x_loc = xx(1) - center(1) ;
	double y_loc = xx(2) - center(2) ;
	double z_loc = xx(3) - center(3) ;
	double air_loc = sqrt (x_loc*x_loc + y_loc*y_loc + z_loc*z_loc) ;

	bool res = ((air_loc*(alpha+beta) -1  <= +prec) && (air_loc * (beta-alpha) -1 >= -prec)) ? true : false ;
	return res ;
}

//Computes the numerical coordinates, as a function of the absolute ones.
const Point Domain_shell_surr::absol_to_num(const Point& abs) const {

	assert (is_in(abs)) ;
	Point num(3) ;
	
	double x_loc = abs(1) - center(1) ;
	double y_loc = abs(2) - center(2) ;
	double z_loc = abs(3) - center(3) ;
	double air = sqrt(x_loc*x_loc+y_loc*y_loc+z_loc*z_loc) ;
	num.set(1) = (1./air-beta)/alpha ;
	double rho = sqrt(x_loc*x_loc+y_loc*y_loc) ;
	
	if (rho==0) {
	    // On the axis ?
	    num.set(2) = (z_loc>=0) ? 0 : M_PI ;
	    num.set(3) = 0 ;
	}
 	else {
	    num.set(2) = atan(rho/z_loc) ;
	    num.set(3) = atan2 (y_loc, x_loc) ;
        }	
	
	if (num(2) <0)
	    num.set(2) = M_PI + num(2) ;

	return num ;
}
 

// Convert absolute coordinates to numerical ones
const Point Domain_shell_surr::absol_to_num_bound(const Point& abs, int bound) const {

	assert ((bound==OUTER_BC) || (bound==INNER_BC)) ;
	assert (is_in(abs, 1e-3)) ;
	Point num(3) ;
	
	double x_loc = abs(1) - center(1) ;
	double y_loc = abs(2) - center(2) ;
	double z_loc = abs(3) - center(3) ;
	
	switch (bound) {
	  case INNER_BC:
	    num.set(1) = -1 ;
	    break ;
	  case OUTER_BC:
	    num.set(1) = 1 ;
	    break ;
	  default:
	      cerr << "unknown boundary in Domain_shell::absol_to_num" << endl ;
	      abort() ;
	}
	
	double rho = sqrt(x_loc*x_loc+y_loc*y_loc) ;
	
	if (rho==0) {
	    // Sur l'axe
	    num.set(2) = (z_loc>=0) ? 0 : M_PI ;
	    num.set(3) = 0 ;
	}
 	else {
	    num.set(2) = atan(rho/z_loc) ;
	    num.set(3) = atan2 (y_loc, x_loc) ;
        }	
	
	if (num(2) <0)
	    num.set(2) = M_PI + num(2) ;
	
	return num ;
}

// Computes the derivatives with respect to XYZ as a function of the numerical ones.
void Domain_shell_surr::do_der_abs_from_der_var(const Val_domain_ptr_array &der_var, Val_domain_ptr_array &der_abs) const {

	// d/dx :
	Val_domain sintdr (-der_var[0]->mult_sin_theta()/alpha/get_radius()/get_radius()) ;	
	Val_domain dtsr (*der_var[1]/get_radius()) ;
	dtsr.set_base() = der_var[1]->get_base() ;
	Val_domain dpsr (*der_var[2]/get_radius()) ;	
	dpsr.set_base() = der_var[2]->get_base() ;
	Val_domain costdtsr (dtsr.mult_cos_theta()) ;

	Val_domain dpsrssint (dpsr.div_sin_theta()) ;
	der_abs[0] = new Val_domain ((sintdr+costdtsr).mult_cos_phi() - dpsrssint.mult_sin_phi()) ;

	// d/dy :
	der_abs[1] = new Val_domain ((sintdr+costdtsr).mult_sin_phi() + dpsrssint.mult_cos_phi()) ;
	// d/dz :
	der_abs[2] = new Val_domain (-der_var[0]->mult_cos_theta()/alpha/get_radius()/get_radius() - dtsr.mult_sin_theta()) ;
}

}
