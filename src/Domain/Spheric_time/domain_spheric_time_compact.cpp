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
#include "spheric_time.hpp"
#include "array.hpp"
#include "val_domain.hpp"

namespace Kadath {
void coef_1d (int, Array<double>&) ;
void coef_i_1d (int, Array<double>&) ;
int der_1d (int, Array<double>&) ;

// Standard constructor
Domain_spheric_time_compact::Domain_spheric_time_compact (int numdom, int ttype, double tmmin, double tmmax, double r, const Dim_array& nbr) :  
		Domain(numdom, ttype, nbr), alpha(-0.5/r), tmin(tmmin), tmax(tmmax) {
     assert (nbr.get_ndim()==2) ;
     do_coloc() ;
}

// Constructor by copy
Domain_spheric_time_compact::Domain_spheric_time_compact (const Domain_spheric_time_compact& so) : Domain(so), alpha(so.alpha), tmin(so.tmin), tmax(so.tmax) {
}

Domain_spheric_time_compact::Domain_spheric_time_compact (int num, FILE* fd) : Domain(num, fd) {
	fread_be (&alpha, sizeof(double), 1, fd) ;
	fread_be (&tmin, sizeof(double), 1, fd) ;
	fread_be (&tmax, sizeof(double), 1, fd) ;   
	do_coloc() ;
}

// Destructor
Domain_spheric_time_compact::~Domain_spheric_time_compact() {}

void Domain_spheric_time_compact::save (FILE* fd) const {
	nbr_points.save(fd) ;
	nbr_coefs.save(fd) ;
	fwrite_be (&ndim, sizeof(int), 1, fd) ;
	fwrite_be (&type_base, sizeof(int), 1, fd) ;
	fwrite_be (&alpha, sizeof(double), 1, fd) ;
	fwrite_be (&tmin, sizeof(double), 1, fd) ;
	fwrite_be (&tmax, sizeof(double), 1, fd) ;
}

ostream& Domain_spheric_time_compact::print (ostream& o) const {
  o << "Spheric-time compact" << endl ;
  o << "time goes from " << tmin << " to " << tmax << endl ;
  o << "Rmin    = " << -0.5/alpha << endl ;
  o << "Nbr pts = " << nbr_points << endl ;
  o << endl ;
  return o ;
}


Val_domain Domain_spheric_time_compact::der_normal (const Val_domain& so, int bound) const {

	
	Val_domain res (so.der_var(1)) ;
	switch (bound) {
		case OUTER_BC :
			res = res.mult_xm1() ;
			res = res.mult_xm1() ;
			res *= -alpha ;
			break ;
		case INNER_BC :
			res = res.mult_xm1() ;
			res = res.mult_xm1() ;
			res *= -alpha ;
			break ;
		default:
			cerr << "Unknown boundary case in Domain_spheric_periodic_compact::der_normal" << endl ;
			abort() ;
		}
return res ;
}

// Computes the cartesian coordinates
void Domain_spheric_time_compact::do_absol () const {
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
		absol[0]->set(index) = 1./alpha / ((*coloc[0])(index(0))-1)  ;
		absol[1]->set(index) = (*coloc[1])(index(1)) * (tmax-tmin)/2. + (tmax+tmin)/2. ;
	}
	while (index.inc())  ;
	
}

// Computes the radius
void Domain_spheric_time_compact::do_radius ()  const {

	for (int i=0 ; i<2 ; i++)
	   assert (coloc[i] != 0x0) ;
	assert (radius == 0x0) ;
	radius = new Val_domain(this) ;
	radius->allocate_conf() ;
	Index index (nbr_points) ;
	do
		radius->set(index) =   1./alpha/ ((*coloc[0])(index(0))-1) ;
 	while (index.inc())  ;
}


// Is a point inside this domain ?
bool Domain_spheric_time_compact::is_in (const Point& xx, double prec) const {

	assert (xx.get_ndim()==2) ;
	
	bool res = (xx(1) >= -0.5/alpha+prec) ? true : false ;
	if ((xx(2)<tmin-prec) || (xx(2)>tmax + prec))
	      res = false ;
	return res ;
}

// Convert absolute coordinates to numerical ones
const Point Domain_spheric_time_compact::absol_to_num(const Point& abs) const {

	assert (is_in(abs)) ;
	Point num(2) ;
	
	num.set(1) = 1 + 1./abs(1)/alpha ;
	num.set(2) = 2./(tmax-tmin)*(abs(2)- (tmax+tmin)/2.) ;
	
	return num ;
}

double coloc_leg(int, int) ;
void Domain_spheric_time_compact::do_coloc () {

	switch (type_base) {
		case CHEB_TYPE:
			nbr_coefs = nbr_points ;
			del_deriv() ;
			for (int i=0 ; i<ndim ; i++)
				coloc[i] = new Array<double> (nbr_points(i)) ;
			for (int i=0 ; i<nbr_points(0) ; i++)
				coloc[0]->set(i) = -cos(M_PI*i/(nbr_points(0)-1)) ;
			for (int j=0 ; j<nbr_points(1) ; j++)
				coloc[1]->set(j) =  -cos(M_PI*j/(nbr_points(1)-1))  ;
			break ;
		case LEG_TYPE:
			nbr_coefs = nbr_points ;
			del_deriv() ;
			for (int i=0 ; i<ndim ; i++) 
				coloc[i] = new Array<double> (nbr_points(i)) ;
			for (int i=0 ; i<nbr_points(0) ; i++)
				coloc[0]->set(i) = coloc_leg(i, nbr_points(0)) ;
			for (int j=0 ; j<nbr_points(1) ; j++)
				coloc[1]->set(j) = coloc_leg(j, nbr_points(1)) ;
			break ;
		default :
			cerr << "Unknown type of basis in Domain_spheric_time_compact::do_coloc" << endl ;
			abort() ;
	}
}

// Base for a function symetric in z, using Chebyshev
void Domain_spheric_time_compact::set_cheb_base(Base_spectral& base) const {

	assert (type_base == CHEB_TYPE) ;
	base.allocate(nbr_coefs) ;
	
	Index index(base.bases_1d[0]->get_dimensions()) ;
	
	base.def=true ;
	base.bases_1d[1]->set(0) = CHEB ;
	for (int j=0 ; j<nbr_coefs(1) ; j++)
	    base.bases_1d[0]->set(j) = CHEB ;
}

// Base for a function symetric in z, using Legendre
void Domain_spheric_time_compact::set_legendre_base(Base_spectral& base) const  {

	assert (type_base == LEG_TYPE) ;
	base.allocate(nbr_coefs) ;
	
	Index index(base.bases_1d[0]->get_dimensions()) ;
	base.def=true ;
	base.bases_1d[1]->set(0) = LEG ;
	for (int j=0 ; j<nbr_coefs(1) ; j++)
	    base.bases_1d[0]->set(j) = LEG ;
 }
 
// Computes the derivatives with respect to rho,Z as a function of the numerical ones.
void Domain_spheric_time_compact::do_der_abs_from_der_var(Val_domain** der_var, Val_domain** der_abs) const {
	// d/dr
	der_abs[0] = new Val_domain (-der_var[0]->mult_xm1()) ;
	// d/dt :
	der_abs[1] = new Val_domain (*der_var[1]*2./(tmax-tmin)) ;
}

// Rules for the multiplication of two basis.
Base_spectral Domain_spheric_time_compact::mult (const Base_spectral& a, const Base_spectral& b) const {

	assert (a.ndim==2) ;
	assert (b.ndim==2) ;
	
	Base_spectral res(2) ;
	bool res_def = true ;

	if (!a.def)
		res_def=false ;
	if (!b.def)
		res_def=false ;
		
	if (res_def) {


	// Bases in time :
	res.bases_1d[1] = new Array<int> (a.bases_1d[1]->get_dimensions()) ;
	switch ((*a.bases_1d[1])(0)) {
		case CHEB:
			switch ((*b.bases_1d[1])(0)) {
				case CHEB:
					res.bases_1d[1]->set(0) = CHEB ;
					break ;
				default:
					res_def = false ;
					break ;
				}
			break ;	
		case LEG:
			switch ((*b.bases_1d[1])(0)) {
				case LEG:
					res.bases_1d[1]->set(0) = LEG ;
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
					res.bases_1d[0]->set(index_0) = LEG  ;
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
}
