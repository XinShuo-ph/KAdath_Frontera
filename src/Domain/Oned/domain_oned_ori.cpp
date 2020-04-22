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
#include "oned.hpp"
#include "point.hpp"
#include "array_math.hpp"
#include "val_domain.hpp"

namespace Kadath {
void coef_1d (int, Array<double>&) ;
void coef_i_1d (int, Array<double>&) ;
int der_1d (int, Array<double>&) ;

// Standard constructor
Domain_oned_ori::Domain_oned_ori (int num, int ttype, double lim, const Dim_array& nbr) :  
		Domain(num, ttype, nbr), alpha(lim) {
     assert (nbr.get_ndim()==1) ;
     do_coloc() ;
}

// Constructor by copy
Domain_oned_ori::Domain_oned_ori (const Domain_oned_ori& so) : Domain(so), alpha(so.alpha) {
}

Domain_oned_ori::Domain_oned_ori (int num, FILE* fd) : Domain(num, fd) {
	fread_be (&alpha, sizeof(double), 1, fd) ;
	do_coloc() ;
}

// Destructor
Domain_oned_ori::~Domain_oned_ori() {}

void Domain_oned_ori::save (FILE* fd) const {
	nbr_points.save(fd) ;
	nbr_coefs.save(fd) ;
	fwrite_be (&ndim, sizeof(int), 1, fd) ;
	fwrite_be (&type_base, sizeof(int), 1, fd) ;
	fwrite_be (&alpha, sizeof(double), 1, fd) ;
}

ostream& operator<< (ostream& o, const Domain_oned_ori& so) {
  o << "One dimensional domain containing 0" << endl ;
  o << "Rmax    = " << so.alpha << endl ;
  o << "Nbr pts = " << so.nbr_points << endl ;
  o << endl ;
  return o ;
}


Val_domain Domain_oned_ori::der_normal (const Val_domain& so, int bound) const {

	Val_domain res (so.der_var(1)) ;
	switch (bound) {
		case OUTER_BC :
			res /= alpha ;
			break ;
		default:
			cerr << "Unknown boundary case in Domain_oned_ori::der_normal" << endl ;
			abort() ;
		}
return res ;
}

// Computes the cartesian coordinates
void Domain_oned_ori::do_absol () const {
	assert (coloc[0] != 0x0) ;
	assert (absol[0] == 0x0) ;
	absol[0] = new Val_domain(this) ;
	absol[0]->allocate_conf() ;

	Index index (nbr_points) ;
	do  {
		absol[0]->set(index) = alpha* ((*coloc[0])(index(0))) ;
	}
	while (index.inc())  ;
	
}

void Domain_oned_ori::do_radius ()  const {

	assert (coloc[0] != 0x0) ;
	assert (radius == 0x0) ;
	radius = new Val_domain(this) ;
	radius->allocate_conf() ;
	Index index (nbr_points) ;
	do
		radius->set(index) = alpha* ((*coloc[0])(index(0))) ;
 	while (index.inc())  ;
}

// Is a point inside this domain ?
bool Domain_oned_ori::is_in (const Point& xx, double prec) const {

	assert (xx.get_ndim()==1) ;	
	bool res = ((xx(1)>=-prec) && (xx(1) <= alpha+prec)) ? true : false ;
	return res ;
}

// Convert absolute coordinates to numerical ones
const Point Domain_oned_ori::absol_to_num(const Point& abs) const {

	assert (is_in(abs)) ;
	Point num(1) ;
	num.set(1) = abs(1)/alpha ;
	return num ;
}

double coloc_leg_parity(int, int) ;
void Domain_oned_ori::do_coloc () {

	switch (type_base) {
		case CHEB_TYPE:
			nbr_coefs = nbr_points ;
			del_deriv() ;
			coloc[0] = new Array<double> (nbr_points(0)) ;
			for (int i=0 ; i<nbr_points(0) ; i++)
			  coloc[0]->set(i) = sin(M_PI/2.*i/(nbr_points(0)-1)) ;
			break ;
		case LEG_TYPE:
			nbr_coefs = nbr_points ;
			del_deriv() ;
			coloc[0] = new Array<double> (nbr_points(0)) ;
			for (int i=0 ; i<nbr_points(0) ; i++)
			  coloc[0]->set(i) = coloc_leg_parity(i, nbr_points(0)) ;
			break ;
		default :
			cerr << "Unknown type of basis in Domain_oned_ori::do_coloc" << endl ;
			abort() ;
	}
}

// Base for a function symetric in z, using Chebyshev
void Domain_oned_ori::set_cheb_base(Base_spectral& base) const {
	assert (type_base == CHEB_TYPE) ;
	base.allocate(nbr_coefs) ;
	
	base.def=true ;
	base.bases_1d[0]->set(0) = CHEB_EVEN ;
}

// Base for a function anti-symetric in z, using Chebyshev
void Domain_oned_ori::set_cheb_base_odd(Base_spectral& base) const {
	assert (type_base == CHEB_TYPE) ;
	base.allocate(nbr_coefs) ;
	
	base.def=true ;
	base.bases_1d[0]->set(0) = CHEB_ODD ;
}

// Base for a function symetric in z, using Legendre
void Domain_oned_ori::set_legendre_base(Base_spectral& base) const  {
	assert (type_base == LEG_TYPE) ;
	base.allocate(nbr_coefs) ;
	
	base.def=true ;
	base.bases_1d[0]->set(0) = LEG_EVEN ;
 }

// Base for a function anti-symetric in z, using Legendre
void Domain_oned_ori::set_legendre_base_odd(Base_spectral& base) const  {
	assert (type_base == LEG_TYPE) ;
	base.allocate(nbr_coefs) ;
	
	base.def=true ;
	base.bases_1d[0]->set(0) = LEG_ODD ;
 }

// Computes the derivatives with respect to rho,Z as a function of the numerical ones.
void Domain_oned_ori::do_der_abs_from_der_var(const Val_domain *const *const der_var, Val_domain **const der_abs) const {
	der_abs[0] = new Val_domain (*der_var[0]/alpha) ;
}

// Rules for the multiplication of two basis.
Base_spectral Domain_oned_ori::mult (const Base_spectral& a, const Base_spectral& b) const {

	assert (a.ndim==1) ;
	assert (b.ndim==1) ;
	
	Base_spectral res(1) ;
	bool res_def = true ;

	if (!a.def)
		res_def=false ;
	if (!b.def)
		res_def=false ;
		
	if (res_def) {


	// Bases in theta :
	res.bases_1d[0] = new Array<int> (a.bases_1d[0]->get_dimensions()) ;
	switch ((*a.bases_1d[0])(0)) {
		case CHEB_EVEN:
			switch ((*b.bases_1d[0])(0)) {
				case CHEB_EVEN:
					res.bases_1d[0]->set(0) = CHEB_EVEN ;
					break ;
				case CHEB_ODD:
					res.bases_1d[0]->set(0) = CHEB_ODD ;
					break ;
				default:
					res_def = false ;
					break ;
				}
			break ;
			
		case CHEB_ODD:
			switch ((*b.bases_1d[0])(0)) {
				case CHEB_EVEN:
					res.bases_1d[0]->set(0) = CHEB_ODD ;
					break ;
				case CHEB_ODD:
					res.bases_1d[0]->set(0) = CHEB_EVEN ;
					break ;
				default:
					res_def = false ;
					break ;
				}
			break ;
		case LEG_EVEN:
			switch ((*b.bases_1d[0])(0)) {
				case LEG_EVEN:
					res.bases_1d[0]->set(0) = LEG_EVEN ;
					break ;
				case LEG_ODD:
					res.bases_1d[0]->set(0) = LEG_ODD ;
					break ;
				default:
					res_def = false ;
					break ;
				}
			break ;
			
		case LEG_ODD:
			switch ((*b.bases_1d[0])(0)) {
				case LEG_EVEN:
					res.bases_1d[0]->set(0) = LEG_ODD ;
					break ;
				case LEG_ODD:
					res.bases_1d[0]->set(0) = LEG_EVEN ;
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
	if (!res_def) 
		for (int dim=0 ; dim<a.ndim ; dim++)
			if (res.bases_1d[dim]!= 0x0) {
				delete res.bases_1d[dim] ;
				res.bases_1d[dim] = 0x0 ;
				}
	res.def = res_def ;
	return res ;
}

int Domain_oned_ori::give_place_var (char* p) const {
    int res = -1 ;
    if (strcmp(p,"X ")==0)
	res = 0 ;
    return res ;
}}
