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
#include "polar.hpp"
#include "point.hpp"
#include "array.hpp"
#include "val_domain.hpp"
namespace Kadath {
// Standard constructor
Domain_polar_shell::Domain_polar_shell (int num, int ttype, double rint, double rext, const Point& cr, const Dim_array& nbr) : Domain(num, ttype, nbr),
							 alpha((rext-rint)/2.), beta((rext+rint)/2.), center(cr) {
    
     assert (nbr.get_ndim()==2) ;
     assert (cr.get_ndim()==2) ;    // Affectation de type_point :
     do_coloc() ;

}


// Constructor by copy
Domain_polar_shell::Domain_polar_shell (const Domain_polar_shell& so) : Domain(so), alpha(so.alpha), 
						beta(so.beta), center(so.center){}



Domain_polar_shell::Domain_polar_shell (int num, FILE* fd) : Domain(num, fd), center(fd) {
	fread_be (&alpha, sizeof(double), 1, fd) ;
	fread_be (&beta, sizeof(double), 1, fd) ;
	do_coloc() ;
}

// Destructor
Domain_polar_shell::~Domain_polar_shell() {}

void Domain_polar_shell::save (FILE* fd) const {
	nbr_points.save(fd) ;
	nbr_coefs.save(fd) ;
	fwrite_be (&ndim, sizeof(int), 1, fd) ;
	fwrite_be (&type_base, sizeof(int), 1, fd) ;	
	center.save(fd) ;
	fwrite_be (&alpha, sizeof(double), 1, fd) ;	
	fwrite_be (&beta, sizeof(double), 1, fd) ;
}

ostream& Domain_polar_shell::print (ostream& o) const {
  o << "Polar shell" << endl ;
  o <<  beta-alpha << " < r < " << beta+alpha << endl ;
  o << "Center  = " << center << endl ;
  o << "Nbr pts = " << nbr_points << endl ;
  o << endl ;
  return o ;
}



Val_domain Domain_polar_shell::der_normal (const Val_domain& so, int bound) const {

	Val_domain res (so.der_var(1)) ;
	switch (bound) {
		case OUTER_BC :
			res /= alpha ;
			break ;
		case INNER_BC :
			res /= alpha ;
			break ;
		
		default:
			cerr << "Unknown boundary case in Domain_polar_shell::der_normal" << endl ;
			abort() ;
		}
return res ;
}


// Computes the Cartesian coordinates
void Domain_polar_shell::do_absol () const  {
	for (int i=0 ; i<2 ; i++)
	   assert (coloc[i] != 0x0) ;
	for (int i=0 ; i<2 ; i++)
	   assert (absol[i] == 0x0) ;
	for (int i=0 ; i<2 ; i++) {
	   absol[i] = new Val_domain(this) ;
	   absol[i]->allocate_conf() ;
	   }
	Index index (nbr_points) ;

	do  {
		absol[0]->set(index) = (alpha* ((*coloc[0])(index(0))) +beta)*
		 sin((*coloc[1])(index(1))) + center(1) ;
		absol[1]->set(index) = (alpha* ((*coloc[0])(index(0))) + beta) * cos((*coloc[1])(index(1))) + center(2) ;
			}
	while (index.inc())  ;
}

// Computes the radius
void Domain_polar_shell::do_radius () const  {
	for (int i=0 ; i<2 ; i++)
	   assert (coloc[i] != 0x0) ;
	assert (radius == 0x0) ;
	radius = new Val_domain(this) ;
	radius->allocate_conf() ;
	Index index (nbr_points) ;
	do 
		radius->set(index) = alpha* ((*coloc[0])(index(0))) + beta ;
	while (index.inc())  ;
}

// Computes the Cartesian coordinates
void Domain_polar_shell::do_cart () const  {
	for (int i=0 ; i<2 ; i++)
	   assert (coloc[i] != 0x0) ;
	for (int i=0 ; i<2 ; i++)
	   assert (cart[i] == 0x0) ;
	for (int i=0 ; i<2 ; i++) {
	   cart[i] = new Val_domain(this) ;
	   cart[i]->allocate_conf() ;
	   }
	Index index (nbr_points) ;

	do  {
		cart[0]->set(index) = (alpha* ((*coloc[0])(index(0))) +beta)*
		 sin((*coloc[1])(index(1))) + center(1) ;
		cart[1]->set(index) = (alpha* ((*coloc[0])(index(0))) + beta) * cos((*coloc[1])(index(1))) + center(2) ;
			}
	while (index.inc())  ;
}

// Computes the Cartesian coordinates
void Domain_polar_shell::do_cart_surr () const  {
	for (int i=0 ; i<2 ; i++)
	   assert (coloc[i] != 0x0) ;
	for (int i=0 ; i<2 ; i++)
	   assert (cart_surr[i] == 0x0) ;
	for (int i=0 ; i<2 ; i++) {
	   cart_surr[i] = new Val_domain(this) ;
	   cart_surr[i]->allocate_conf() ;
	   }
	Index index (nbr_points) ;

	do  {
		cart_surr[0]->set(index) = sin((*coloc[1])(index(1))) + center(1) ;
		cart_surr[1]->set(index) =  cos((*coloc[1])(index(1))) + center(2) ;
			}
	while (index.inc())  ;
}

// Check if a point is inside this domain
bool Domain_polar_shell::is_in (const Point& xx, double prec) const {

	assert (xx.get_ndim()==2) ;
		
	double rho_loc = xx(1) - center(1) ;
	double z_loc = xx(2) - center(2) ;
	double air_loc = sqrt (rho_loc*rho_loc + z_loc*z_loc) ;
	
	bool res = ((air_loc <= alpha+beta+prec) && (air_loc >= beta-alpha-prec)) ? true : false ;
	return res ;
}

//Computes the numerical coordinates, as a function of the absolute ones.
const Point Domain_polar_shell::absol_to_num(const Point& abs) const {

	assert (is_in(abs)) ;
	Point num(2) ;
	
	double rho_loc = fabs(abs(1) - center(1)) ;
	double z_loc = abs(2) - center(2) ;
	double air = sqrt(rho_loc*rho_loc+z_loc*z_loc) ;
	num.set(1) = (air-beta)/alpha ;
	
	if (rho_loc==0) {
	    // On the axis ?
	    num.set(2) = (z_loc>=0) ? 0 : M_PI ;
	}
 	else {
	    num.set(2) = atan(rho_loc/z_loc) ;
        }	
	
	if (num(2) <0)
	    num.set(2) = M_PI + num(2) ;

	return num ;
}

double coloc_leg (int, int) ;
void Domain_polar_shell::do_coloc () {

	switch (type_base) {
		case CHEB_TYPE:
			nbr_coefs = nbr_points ;
			del_deriv() ;
			for (int i=0 ; i<ndim ; i++)
				coloc[i] = new Array<double> (nbr_points(i)) ;
			for (int i=0 ; i<nbr_points(0) ; i++)
				coloc[0]->set(i) = -cos(M_PI*i/(nbr_points(0)-1)) ;
			for (int j=0 ; j<nbr_points(1) ; j++)
				coloc[1]->set(j) = M_PI/2.*j/(nbr_points(1)-1) ;
			break ;
		case LEG_TYPE:
			nbr_coefs = nbr_points ;
			del_deriv() ;
			for (int i=0 ; i<ndim ; i++) 
				coloc[i] = new Array<double> (nbr_points(i)) ;
			for (int i=0 ; i<nbr_points(0) ; i++)
				coloc[0]->set(i) = coloc_leg(i, nbr_points(0)) ;
			for (int j=0 ; j<nbr_points(1) ; j++)
				coloc[1]->set(j) = M_PI/2.*j/(nbr_points(1)-1) ;
			break ;
		default :
			cerr << "Unknown type of basis in Domain_polar_shell::do_coloc" << endl ;
			abort() ;
	}
}
// standard base for a symetric function in z, using Chebyshev
void Domain_polar_shell::set_cheb_base(Base_spectral& base) const  {
  set_cheb_base_with_m(base, 0) ;
}

// standard base for a anti-symetric function in z, using Chebyshev
void Domain_polar_shell::set_anti_cheb_base(Base_spectral& base) const  {
  set_anti_cheb_base_with_m(base, 0) ;
}

// standard base for a symetric function in z, using Legendre
void Domain_polar_shell::set_legendre_base(Base_spectral& base) const {
  set_legendre_base_with_m(base, 0) ;
}

// standard base for a anti-symetric function in z, using Legendre
void Domain_polar_shell::set_anti_legendre_base(Base_spectral& base) const {
  set_anti_legendre_base_with_m(base, 0) ;
}


// standard base for a symetric function in z, using Chebyshev
void Domain_polar_shell::set_cheb_base_with_m(Base_spectral& base, int m) const  {

	assert (type_base == CHEB_TYPE) ;

	base.allocate(nbr_coefs) ;
	    
	base.def =true ;
	
	base.bases_1d[1]->set(0) = (m%2==0) ? COS_EVEN : SIN_ODD ;
	for (int j=0 ; j<nbr_coefs(1) ; j++)
	  base.bases_1d[0]->set(j) = CHEB ;
}

// standard base for a anti-symetric function in z, using Chebyshev
void Domain_polar_shell::set_anti_cheb_base_with_m(Base_spectral& base, int m) const  {

	assert (type_base == CHEB_TYPE) ;

	base.allocate(nbr_coefs) ;
	    
	base.def =true ;
	
	base.bases_1d[1]->set(0) = (m%2==0) ? COS_ODD : SIN_EVEN ;
	for (int j=0 ; j<nbr_coefs(1) ; j++)
	  base.bases_1d[0]->set(j) = CHEB ;
}

// standard base for a symetric function in z, using Legendre
void Domain_polar_shell::set_legendre_base_with_m(Base_spectral& base, int m) const {

	assert (type_base == CHEB_TYPE) ;

	base.allocate(nbr_coefs) ;
	    
	base.def =true ;
	
	base.bases_1d[1]->set(0) = (m%2==0) ? COS_EVEN : SIN_ODD ;
	for (int j=0 ; j<nbr_coefs(1) ; j++)
	  base.bases_1d[0]->set(j) = LEG ;
}

// standard base for a anti-symetric function in z, using Legendre
void Domain_polar_shell::set_anti_legendre_base_with_m(Base_spectral& base, int m) const {

	assert (type_base == CHEB_TYPE) ;

	base.allocate(nbr_coefs) ;
	    
	base.def =true ;
	
	base.bases_1d[1]->set(0) = (m%2==0) ? COS_ODD : SIN_EVEN ;
	for (int j=0 ; j<nbr_coefs(1) ; j++)
	  base.bases_1d[0]->set(j) = LEG ;
}

// Computes the derivatives with respect to XYZ as a function of the numerical ones.
void Domain_polar_shell::do_der_abs_from_der_var(Val_domain** der_var, Val_domain** der_abs) const {
	
	Val_domain dr (*der_var[0]/alpha) ;
	Val_domain dtsr (*der_var[1]/get_radius()) ;
	dtsr.set_base() = der_var[1]->get_base() ;

	// D / drho
	der_abs[0] = new Val_domain (dr.mult_sin_theta() + dtsr.mult_cos_theta()) ;
	// d/dz :
	der_abs[1] = new Val_domain (dr.mult_cos_theta() - dtsr.mult_sin_theta()) ;
}

// Rules for the multiplication of two basis.
Base_spectral Domain_polar_shell::mult (const Base_spectral& a, const Base_spectral& b) const {

	assert (a.ndim==2) ;
	assert (b.ndim==2) ;
	
	Base_spectral res(2) ;
	bool res_def = true ;

	if (!a.def)
		res_def=false ;
	if (!b.def)
		res_def=false ;
		
	if (res_def) {

	// Bases in theta :
	res.bases_1d[1] = new Array<int> (a.bases_1d[1]->get_dimensions()) ;
	switch ((*a.bases_1d[1])(0)) {
		case COS_EVEN:
			switch ((*b.bases_1d[1])(0)) {
				case COS_EVEN:
					res.bases_1d[1]->set(0) = COS_EVEN ;
					break ;
				case COS_ODD:
					res.bases_1d[1]->set(0) = COS_ODD ;
					break ;
				case SIN_EVEN:
					res.bases_1d[1]->set(0) = SIN_EVEN ;
					break ;
				case SIN_ODD:
					res.bases_1d[1]->set(0) = SIN_ODD ;
					break ;
				default:
					res_def = false ;
					break ;
				}
			break ;
		case COS_ODD:
			switch ((*b.bases_1d[1])(0)) {
				case COS_EVEN:
					res.bases_1d[1]->set(0) = COS_ODD ;
					break ;
				case COS_ODD:
					res.bases_1d[1]->set(0) = COS_EVEN ;
					break ;
				case SIN_EVEN:
					res.bases_1d[1]->set(0) =  SIN_ODD  ;
					break ;
				case SIN_ODD:
					res.bases_1d[1]->set(0) = SIN_EVEN  ;
					break ;
				default:
					res_def = false ;
					break ;
				}
			break ;
		case SIN_EVEN:
			switch ((*b.bases_1d[1])(0)) {
				case COS_EVEN:
					res.bases_1d[1]->set(0) = SIN_EVEN ;
					break ;
				case COS_ODD:
					res.bases_1d[1]->set(0) = SIN_ODD ;
					break ;
				case SIN_EVEN:
					res.bases_1d[1]->set(0) = COS_EVEN ;
					break ;
				case SIN_ODD:
					res.bases_1d[1]->set(0) = COS_ODD ;
					break ;
				default:
					res_def = false ;
					break ;
				}
			break ;
		case SIN_ODD:
			switch ((*b.bases_1d[1])(0)) {
				case COS_EVEN:
					res.bases_1d[1]->set(0) = SIN_ODD ;
					break ;
				case COS_ODD:
					res.bases_1d[1]->set(0) = SIN_EVEN ;
					break ;
				case SIN_EVEN:
					res.bases_1d[1]->set(0) = COS_ODD ;
					break ;
				case SIN_ODD:
					res.bases_1d[1]->set(0) = COS_EVEN ;
					break ;
				default:
					res_def = false ;
					break ;
				}
			break ;
		default:
			res_def = false ;
			break ;
	}



	// Base in r :
	Index index_0 (a.bases_1d[0]->get_dimensions()) ;
	res.bases_1d[0] = new Array<int> (a.bases_1d[0]->get_dimensions()) ;
	do {
	switch ((*a.bases_1d[0])(index_0)) {
		case CHEB:
			switch ((*b.bases_1d[0])(index_0)) {
				case CHEB:
					res.bases_1d[0]->set(index_0) = CHEB ;
					break ;
				default:
					res_def = false ;
					break ;
				}
			break ;
		case LEG:
			switch ((*b.bases_1d[0])(index_0)) {
				case LEG:
					res.bases_1d[0]->set(index_0) = LEG ;
					break ;
				default:
					res_def = false ;
					break ;
				}
			break ;
		default:
			res_def = false ;
			break ;
		}
	}
	while (index_0.inc()) ;
	}
		
	if (!res_def) 
		for (int dim=0 ; dim<a.ndim ; dim++)
			if (res.bases_1d[dim]!= 0x0) {
				delete res.bases_1d[dim] ;
				res.bases_1d[dim] = 0x0 ;
				}
	res.def = res_def ;
	return res ;
}

int Domain_polar_shell::give_place_var (char* p) const {
    int res = -1 ;
    if (strcmp(p,"R ")==0)
	res = 0 ;
    if (strcmp(p,"T ")==0)
	res = 1 ;
    return res ;
}

}
