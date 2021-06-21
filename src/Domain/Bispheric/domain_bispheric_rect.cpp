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

#include "bispheric.hpp"
#include "array.hpp"
#include "param.hpp"
#include "val_domain.hpp"
namespace Kadath {
// standard constructor
Domain_bispheric_rect::Domain_bispheric_rect (int num, int ttype, double a, double rr, double etamin, double etapl, double chimin, const Dim_array& nbr) :  Domain(num, ttype, nbr), aa(a), r_ext(rr), eta_minus(etamin), eta_plus(etapl), chi_min(chimin), p_eta(0x0), p_chi(0x0), p_phi(0x0),
p_detadx(0x0), p_detady(0x0), p_detadz(0x0), p_dchidx(0x0), p_dchidy(0x0), p_dchidz(0x0), p_dphidy(0x0), p_dphidz(0x0), p_dsint(0x0) {

     assert (nbr.get_ndim()==3) ;
     do_coloc() ;
}

// Copy constructor
Domain_bispheric_rect::Domain_bispheric_rect (const Domain_bispheric_rect& so) : Domain(so), aa(so.aa), 
		eta_minus(so.eta_minus), eta_plus(so.eta_plus), chi_min(so.chi_min) {
	p_eta = (so.p_eta!=0x0) ? new Val_domain(*so.p_eta) : 0x0 ;
	p_chi = (so.p_chi!=0x0) ? new Val_domain(*so.p_chi) : 0x0 ;
	p_phi = (so.p_phi!=0x0) ? new Val_domain(*so.p_phi) : 0x0 ;
	p_detadx = (so.p_detadx!=0x0) ? new Val_domain(*so.p_detadx) : 0x0 ;
	p_detady = (so.p_detady!=0x0) ? new Val_domain(*so.p_detady) : 0x0 ;
	p_detadz = (so.p_detadz!=0x0) ? new Val_domain(*so.p_detadz) : 0x0 ;
	p_dchidx = (so.p_dchidx!=0x0) ? new Val_domain(*so.p_dchidx) : 0x0 ;
	p_dchidy = (so.p_dchidy!=0x0) ? new Val_domain(*so.p_dchidy) : 0x0 ;
	p_dchidz = (so.p_dchidz!=0x0) ? new Val_domain(*so.p_dchidz) : 0x0 ;
	p_dphidy = (so.p_dphidy!=0x0) ? new Val_domain(*so.p_dphidy) : 0x0 ;
        p_dphidz = (so.p_dphidz!=0x0) ? new Val_domain(*so.p_dphidz) : 0x0 ;  
	p_dsint = (so.p_dsint!=0x0) ? new Val_domain(*so.p_dsint) : 0x0 ;
}

Domain_bispheric_rect::Domain_bispheric_rect (int num, FILE* fd) : Domain(num, fd) {
	fread_be (&aa, sizeof(double), 1, fd) ;
	fread_be (&eta_minus, sizeof(double), 1, fd) ;
	fread_be (&eta_plus, sizeof(double), 1, fd) ;
	fread_be (&chi_min, sizeof(double), 1, fd) ;

	p_eta = 0x0 ;
	p_chi = 0x0 ;
	p_phi = 0x0 ;
	p_detadx = 0x0 ;
	p_detady = 0x0 ;
	p_detadz = 0x0 ;
	p_dchidx = 0x0 ;
	p_dchidy = 0x0 ;
	p_dchidz = 0x0 ;
	p_dphidy = 0x0 ;
	p_dphidz = 0x0 ;
	p_dsint = 0x0 ;
	do_coloc() ;
}

// Destructor
Domain_bispheric_rect::~Domain_bispheric_rect() {
	del_deriv() ;
}

void Domain_bispheric_rect::save (FILE* fd) const {
	nbr_points.save(fd) ;
	nbr_coefs.save(fd) ;
	fwrite_be (&ndim, sizeof(int), 1, fd) ;
	fwrite_be (&type_base, sizeof(int), 1, fd) ;
	fwrite_be (&aa, sizeof(double), 1, fd) ;
	fwrite_be (&eta_minus, sizeof(double), 1, fd) ;
	fwrite_be (&eta_plus, sizeof(double), 1, fd) ;
	fwrite_be (&chi_min, sizeof(double), 1, fd) ;
}

// Deletes the derived members
void Domain_bispheric_rect::del_deriv() const  {
	for (int l=0 ; l<ndim ; l++) {
		if (coloc[l] !=0x0) delete coloc[l] ;
		if (cart[l] !=0x0) delete cart[l] ;
		coloc[l] = 0x0 ;
		cart[l] = 0x0 ;
	}
	
	if (radius !=0x0)
	    delete radius ;
	radius = 0x0 ;
	if (p_eta !=0x0)
	    delete p_eta ;
	p_eta = 0x0 ;
	if (p_chi !=0x0)
	    delete p_chi ;
	p_chi = 0x0 ;
	if (p_phi !=0x0)
	    delete p_phi ;
	p_phi = 0x0 ;
	if (p_detadx !=0x0)
	    delete p_detadx ;
	p_detadx = 0x0 ;
	if (p_detady !=0x0)
	    delete p_detady ;
	p_detady = 0x0 ;
	if (p_detadz !=0x0)
	    delete p_detadz ;
	p_detadz = 0x0 ;
	if (p_dchidx !=0x0)
	    delete p_dchidx ;
	p_dchidx = 0x0 ;
	if (p_dchidy !=0x0)
	    delete p_dchidy ;
	p_dchidy = 0x0 ;
        if (p_dchidz !=0x0)
	    delete p_dchidz ;
	p_dchidz = 0x0 ;
	if (p_dphidy !=0x0)
	    delete p_dphidy ;
	p_dphidy = 0x0 ;
	if (p_dphidz !=0x0)
	    delete p_dphidz ;
	p_dphidz = 0x0 ;	
	if (p_dsint !=0x0)
	    delete p_dsint ;
	p_dsint = 0x0 ;
}

// Display
ostream& Domain_bispheric_rect::print (ostream& o) const {

  o << "Bispherical domain, rectangular part" << endl ;
  o << "aa      = " << aa << endl ;
  o << eta_minus << " < eta < " << eta_plus << endl ;
  o << chi_min << " < chi < pi " << endl ;
  o << "Nbr pts = " << nbr_points << endl ;
  o << endl ;
  return o ;
}

// Computes eta from eta star
void Domain_bispheric_rect::do_eta() const {
	for (int i=0 ; i<3 ; i++)
	   assert (coloc[i] != 0x0) ;
	assert (p_eta==0x0) ;
	p_eta= new Val_domain(this) ;
	p_eta->allocate_conf() ;
	Index index (nbr_points) ;
	do 
	   p_eta->set(index) = (-eta_minus+eta_plus)/2.*((*coloc[0])(index(0))) + (eta_minus+eta_plus)/2. ;
	while (index.inc()) ;
}

// Computes chi from chi star
void Domain_bispheric_rect::do_chi() const {
	for (int i=0 ; i<3 ; i++)
	   assert (coloc[i] != 0x0) ;
	assert (p_chi==0x0) ;
	p_chi= new Val_domain(this) ;
	p_chi->allocate_conf() ;
	Index index (nbr_points) ;
	do 
	   p_chi->set(index) =  (chi_min-M_PI)*(*coloc[1])(index(1))  + M_PI ;
	while (index.inc()) ;
}

// Computes phi from phi star
void Domain_bispheric_rect::do_phi() const {
	for (int i=0 ; i<3 ; i++)
	   assert (coloc[i] != 0x0) ;
	assert (p_phi==0x0) ;
	p_phi= new Val_domain(this) ;
	p_phi->allocate_conf() ;
	Index index (nbr_points) ;
	do 
	   p_phi->set(index) = ((*coloc[2])(index(2))) ;
	while (index.inc()) ;
}

Val_domain Domain_bispheric_rect::get_chi() const {
	if (p_chi==0x0)
		do_chi() ;
	return *p_chi ;
}

Val_domain Domain_bispheric_rect::get_eta() const {
	if (p_eta==0x0)
		do_eta() ;
	return *p_eta ;
}

void Domain_bispheric_rect::do_absol () const  {
	for (int i=0 ; i<3 ; i++)
	   assert (coloc[i] != 0x0) ;
	for (int i=0 ; i<3 ; i++)
	   assert (absol[i] == 0x0) ;
	for (int i=0 ; i<3 ; i++) {
	   absol[i] = new Val_domain(this) ;
	   absol[i]->allocate_conf() ;
	   }

	if (p_eta==0x0)
		do_eta() ;
	if (p_chi==0x0)
		do_chi() ;
	if (p_phi==0x0)
		do_phi() ;
	
	Index index (nbr_points) ;
	do  {

		absol[0]->set(index) = aa*sinh((*p_eta)(index))/(cosh((*p_eta)(index))-cos((*p_chi)(index))) ;
		absol[1]->set(index) =  
   aa*sin((*p_chi)(index))*cos((*p_phi)(index))/(cosh((*p_eta)(index))-cos((*p_chi)(index))) ;
		absol[2]->set(index) = 
   aa*sin((*p_chi)(index))*sin((*p_phi)(index))/(cosh((*p_eta)(index))-cos((*p_chi)(index))) ;
		}
	while (index.inc())  ;	

	// Basis
	absol[0]->std_base() ;
	absol[1]->std_base() ;
	absol[2]->std_anti_base() ;
}


// Computes the generalized radius
void Domain_bispheric_rect::do_radius () const  {

	for (int i=0 ; i<3 ; i++)
	   assert (coloc[i] != 0x0) ;
	assert (radius == 0x0) ;
	radius = new Val_domain(sqrt (get_cart(1)*get_cart(1)+get_cart(2)*get_cart(2)+get_cart(3)*get_cart(3))) ;
}

// Computes the Cartesian coordinates
void Domain_bispheric_rect::do_cart () const  {
	for (int i=0 ; i<3 ; i++)
	   assert (coloc[i] != 0x0) ;
	for (int i=0 ; i<3 ; i++)
	   assert (cart[i] == 0x0) ;
	for (int i=0 ; i<3 ; i++) {
	   cart[i] = new Val_domain(this) ;
	   cart[i]->allocate_conf() ;
	   }

	if (p_eta==0x0)
		do_eta() ;
	if (p_chi==0x0)
		do_chi() ;
	if (p_phi==0x0)
		do_phi() ;
	
	Index index (nbr_points) ;
	do  {

		cart[0]->set(index) = aa*sinh((*p_eta)(index))/(cosh((*p_eta)(index))-cos((*p_chi)(index))) ;
		cart[1]->set(index) =  
   aa*sin((*p_chi)(index))*cos((*p_phi)(index))/(cosh((*p_eta)(index))-cos((*p_chi)(index))) ;
		cart[2]->set(index) = 
   aa*sin((*p_chi)(index))*sin((*p_phi)(index))/(cosh((*p_eta)(index))-cos((*p_chi)(index))) ;
		}
	while (index.inc())  ;	

	// Basis
	cart[0]->std_base() ;
	cart[1]->std_base() ;
	cart[2]->std_anti_base() ;
}


void Domain_bispheric_rect::do_dsint () const {
	double rr = aa/sinh(fabs(eta_minus)) ;
	double xc = aa*cosh(eta_minus)/sinh(eta_minus) ;

	Val_domain rho2 (get_cart(2)*get_cart(2)+get_cart(3)*get_cart(3)) ;
	Val_domain xx (this) ;
	xx = (xc<0) ? get_cart(1)-xc : xc-get_cart(1) ;
	
	assert (p_dsint==0x0) ;
	p_dsint = new Val_domain ((0.5*rho2.der_var(2)*xx - xx.der_var(2)*rho2)/rr) ;
}

// Check if a point is inside ?
bool Domain_bispheric_rect::is_in (const Point& abs, double prec) const {
	assert (abs.get_ndim()==3) ;
		
	double xx = abs(1) ;
	double yy = abs(2) ;
	double zz = abs(3) ;
	
	double air = sqrt (xx*xx+yy*yy+zz*zz) ;

	double chi ;
	double rho = sqrt(yy*yy+zz*zz) ;
	
	double x_out = aa*cosh(eta_minus)/sinh(eta_minus) ;
     
	if (rho<prec)
		chi = (fabs(xx)>fabs(x_out)) ? 0 : M_PI  ; 
	else    
		chi = atan (2*aa*rho/(air*air-aa*aa)) ;

	if (chi<0)
	    chi += M_PI ;
	    
	double eta = 0.5*log((1+(2*aa*xx)/(aa*aa+air*air))/(1-(2*aa*xx)/(aa*aa+air*air))) ;
	double eta_num = 2./(eta_plus-eta_minus)*(eta-(eta_plus+eta_minus)/2.) ;

	bool res = true ;
	
	if (chi<chi_min-prec)
		res = false ;
	if (chi>M_PI+prec)
		res = false ;

	if (eta_num<-1-prec)
		res = false ;
	if (eta_num>1+prec)
		res = false ;

	return res ;
}
// Convert absolute coordinates to numerical ones
const Point Domain_bispheric_rect::absol_to_num(const Point& abs) const {

	assert (is_in(abs, 1e-12)) ;
	Point num(3) ;
		
	double xx = abs(1) ;
	double yy = abs(2) ;
	double zz = abs(3) ;
	double air = sqrt (xx*xx+yy*yy+zz*zz) ;

	num.set(3) = atan2 (zz, yy) ;
	
	if (num(3)<0)
	    num.set(3) += 2*M_PI ;
	
	double chi ;
	double rho = sqrt(yy*yy+zz*zz) ;
	
	// On the axis ?
	if (rho<PRECISION)
		chi =  M_PI  ; 
	else    
		chi = atan (2*aa*rho/(air*air-aa*aa)) ;

	if (chi<0)
	    chi += M_PI ;
	    
	double eta = 0.5*log((1+(2*aa*xx)/(aa*aa+air*air))/(1-(2*aa*xx)/(aa*aa+air*air))) ;
	
	num.set(1) = 2./(eta_plus-eta_minus)*(eta-(eta_plus+eta_minus)/2.) ;
	num.set(2) = (chi-M_PI)/(chi_min-M_PI) ;

	return num ;
}

const Point Domain_bispheric_rect::absol_to_num_bound(const Point& abs, int bound) const {

	assert (is_in(abs, 1e-3)) ;
	Point num(3) ;
		
	double xx = abs(1) ;
	double yy = abs(2) ;
	double zz = abs(3) ;
	double air = sqrt (xx*xx+yy*yy+zz*zz) ;

	num.set(3) = atan2 (zz, yy) ;
	
	if (num(3)<0)
	    num.set(3) += 2*M_PI ;
	
	double chi ;
	double rho = sqrt(yy*yy+zz*zz) ;
	
	// On the axis ?
	if (rho<PRECISION)
		chi =  M_PI  ; 
	else    
		chi = atan (2*aa*rho/(air*air-aa*aa)) ;

	if (chi<0)
	    chi += M_PI ;
	    
	double eta = 0.5*log((1+(2*aa*xx)/(aa*aa+air*air))/(1-(2*aa*xx)/(aa*aa+air*air))) ;
	
	num.set(1) = 2./(eta_plus-eta_minus)*(eta-(eta_plus+eta_minus)/2.) ;
	num.set(2) = (chi-M_PI)/(chi_min-M_PI) ;

	switch (bound) {
	  case INNER_BC :
	    num.set(1) = -1 ;
	    break ;
	  case ETA_PLUS_BC :
	    num.set(1) = 1 ;
	    break ;
	  case CHI_ONE_BC :
	    num.set(2) = 1  ;
	    break ;
	  case OUTER_BC :
	    num.set(1) = 1 ;
	    num.set(2) = 1 ;
	    break ;
	  default :
	    cerr << "Unknown case in Domain_bispheric_rect::absol_to_num_bound" << endl ;
	    abort() ;
	}
	return num ;
}

double coloc_leg(int,int) ;
double coloc_leg_parity(int,int) ;
void Domain_bispheric_rect::do_coloc () {

	switch (type_base) {
		case CHEB_TYPE:
			nbr_coefs = nbr_points ;
			del_deriv() ;
			for (int i=0 ; i<ndim ; i++)
				coloc[i] = new Array<double> (nbr_points(i)) ;
			for (int i=0 ; i<nbr_points(0) ; i++)
				coloc[0]->set(i) = -cos(M_PI*i/(nbr_points(0)-1)) ;
			for (int j=0 ; j<nbr_points(1) ; j++)
				coloc[1]->set(j) = sin(M_PI*j/2./(nbr_points(1)-1)) ;
			for (int k=0 ; k<nbr_points(2) ; k++)
				coloc[2]->set(k) = M_PI*k/(nbr_points(2)-1) ;
			do_eta() ;
			do_chi() ;
			do_phi() ;
			break ;
		case LEG_TYPE:
			nbr_coefs = nbr_points ;
			del_deriv() ;
			for (int i=0 ; i<ndim ; i++)
				coloc[i] = new Array<double> (nbr_points(i)) ;
			for (int i=0 ; i<nbr_points(0) ; i++)
				coloc[0]->set(i) = coloc_leg(i,nbr_points(0)) ;
			for (int j=0 ; j<nbr_points(1) ; j++)
				coloc[1]->set(j) = coloc_leg_parity(j,nbr_points(1)) ;
			for (int k=0 ; k<nbr_points(2) ; k++)
				coloc[2]->set(k) = M_PI*k/(nbr_points(2)-1) ;
			do_eta() ;
			do_chi() ;
			do_phi() ;
			break ;
		default : 
			cerr << "Unknown type of basis in Domain_bispheric_rect::do_coloc" << endl ;
			abort() ;
	}
}

// Standard base fo Chebyshev
void Domain_bispheric_rect::set_cheb_base(Base_spectral& base) const {

	assert (type_base == CHEB_TYPE) ;
	base.allocate(nbr_coefs) ;
	
	Index index(base.bases_1d[0]->get_dimensions()) ;
	
	base.def=true ;
	base.bases_1d[2]->set(0) = COS ;
	for (int k=0 ; k<nbr_coefs(2) ; k++) {
		base.bases_1d[1]->set(k) = (k%2==0) ? CHEB_EVEN : CHEB_ODD ;
		for (int j=0 ; j<nbr_coefs(1) ; j++) {
		    index.set(0) = j ; index.set(1) = k ;
		    base.bases_1d[0]->set(index) = CHEB ;
		 }
	}
}

// Antisymetric base for Chebyshev
void Domain_bispheric_rect::set_anti_cheb_base(Base_spectral& base) const {

	assert (type_base == CHEB_TYPE) ;
	base.allocate(nbr_coefs) ;
	
	Index index(base.bases_1d[0]->get_dimensions()) ;
	
	base.def=true ;
	base.bases_1d[2]->set(0) = SIN ;
	for (int k=0 ; k<nbr_coefs(2) ; k++) {
		base.bases_1d[1]->set(k) = (k%2==0) ? CHEB_EVEN : CHEB_ODD ;
		for (int j=0 ; j<nbr_coefs(1) ; j++) {
		    index.set(0) = j ; index.set(1) = k ;
		    base.bases_1d[0]->set(index) = CHEB ;
		 }
	}
}

double coloc_leg(int,int) ;
double coloc_leg_parity(int,int) ;
// Standard base for Legendre
void Domain_bispheric_rect::set_legendre_base(Base_spectral& base) const {
	
	assert (type_base == LEG_TYPE) ;
	base.allocate(nbr_coefs) ;
	
	Index index(base.bases_1d[0]->get_dimensions()) ;

	base.def=true ;    
	base.bases_1d[2]->set(0) = COS ;
	for (int k=0 ; k<nbr_coefs(2) ; k++) {
		base.bases_1d[1]->set(k) = (k%2==0) ? LEG_EVEN : LEG_ODD ;
		for (int j=0 ; j<nbr_coefs(1) ; j++) {
		    index.set(0) = j ; index.set(1) = k ;
		    base.bases_1d[0]->set(index) = LEG ;
		 }
	}
}

// Antisymetric base for Legendre
void Domain_bispheric_rect::set_anti_legendre_base(Base_spectral& base) const {
	
	assert (type_base == LEG_TYPE) ;
	base.allocate(nbr_coefs) ;
	
	Index index(base.bases_1d[0]->get_dimensions()) ;

	base.def=true ;    
	base.bases_1d[2]->set(0) = SIN ;
	for (int k=0 ; k<nbr_coefs(2) ; k++) {
		base.bases_1d[1]->set(k) = (k%2==0) ? LEG_EVEN : LEG_ODD ;
		for (int j=0 ; j<nbr_coefs(1) ; j++) {
		    index.set(0) = j ; index.set(1) = k ;
		    base.bases_1d[0]->set(index) = LEG ;
		 }
	}
}

// Computes the partial derivatives of the numerical coordinates with respect to the Cartesian ones.
void Domain_bispheric_rect::do_for_der() const {

	if (cart[0]==0x0)
		do_cart() ;
	if (radius==0x0)
		do_radius() ;

	// Partial derivatives of eta 
	Val_domain denom_eta 
		((eta_plus-eta_minus)/2.*((aa*aa+(*radius)*(*radius))*((aa*aa+(*radius)*(*radius))) -
					 4.*aa*aa*(*cart[0])*(*cart[0]))) ;
	p_detadx = new Val_domain ((2.*aa*(aa*aa+(*radius)*(*radius)-2.*(*cart[0])*(*cart[0])))/denom_eta) ;
	p_detady = new Val_domain (-4.*aa*(*cart[0])*(*cart[1])/denom_eta) ;
	p_detadz = new Val_domain (-4.*aa*(*cart[0])*(*cart[2])/denom_eta) ;
 
	// Partial derivatives of chi 
	Val_domain rho (sqrt((*cart[1])*(*cart[1])+(*cart[2])*(*cart[2]))) ;
	Val_domain denom_chi ((chi_min-M_PI)*(((*radius)*(*radius)-aa*aa)*((*radius)*(*radius)-aa*aa) + 
					4.*aa*aa*rho*rho)) ;
	p_dchidx = new Val_domain (-4.*aa*(*cart[0])*rho/denom_chi) ;
	Val_domain chiyz (2.*aa*((*radius)*(*radius)-aa*aa-2*rho*rho)/denom_chi) ;
	
	// Partial derivatives of phi :
	Val_domain phiyz (((exp(*p_eta)+exp(-*p_eta))/2.-cos(*p_chi))/aa) ;

	// The basis :
	p_detadx->std_base() ;
	p_detady->std_base() ;
	p_detadz->std_anti_base() ;
	p_dchidx->base = p_detadx->der_var(2).get_base() ;
	chiyz.std_base() ;
	phiyz.std_base() ;

	p_dchidy = new Val_domain (chiyz.mult_cos_phi()) ;
	p_dchidz = new Val_domain (chiyz.mult_sin_phi()) ;
	p_dphidy = new Val_domain (-phiyz.mult_sin_phi()) ;
	p_dphidz = new Val_domain (phiyz.mult_cos_phi()) ;
}

// Computes the derivatives with respect to the absolute coordinates with respect to the numerical ones.
void Domain_bispheric_rect::do_der_abs_from_der_var(Val_domain** der_var, Val_domain** der_abs) const {
	if (p_detadx==0x0)
		do_for_der() ;

	// d/dx :
	der_abs[0] = new Val_domain((*der_var[0])*(*p_detadx) + (*der_var[1])*(*p_dchidx)) ;

	// d/dy :
	Val_domain auchi_y ((*der_var[2])*(*p_dphidy)) ;
	Val_domain part_y_phi (auchi_y.div_sin_chi()) ;
	der_abs[1] = new Val_domain((*der_var[0])*(*p_detady) + (*der_var[1])*(*p_dchidy) + part_y_phi) ;

	// d/dz :
	Val_domain auchi_z ((*der_var[2])*(*p_dphidz)) ;
	Val_domain part_z_phi (auchi_z.div_sin_chi()) ;
	der_abs[2] = new Val_domain((*der_var[0])*(*p_detadz) + (*der_var[1])*(*p_dchidz) + part_z_phi) ;
}

// Multplication rules for the basis.
Base_spectral Domain_bispheric_rect::mult (const Base_spectral& a, const Base_spectral& b) const {

	assert (a.ndim==3) ;
	assert (b.ndim==3) ;
	
	Base_spectral res(3) ;
	bool res_def = true ;

	if (!a.def)
		res_def=false ;
	if (!b.def)
		res_def=false ;
		
	if (res_def) {
	// Base in phi :
	res.bases_1d[2] = new Array<int> (a.bases_1d[2]->get_dimensions()) ;
	switch ((*a.bases_1d[2])(0)) {
		case COS:
			switch ((*b.bases_1d[2])(0)) {
				case COS:
					res.bases_1d[2]->set(0) = COS ;
					break ;
				case SIN:
					res.bases_1d[2]->set(0) = SIN ;
					break ;
				default:
					res_def = false ;
					break ;
				}
			break ;
		case SIN:
			switch ((*b.bases_1d[2])(0)) {
				case COS:
					res.bases_1d[2]->set(0) = SIN ;
					break ;
				case SIN:
					res.bases_1d[2]->set(0) = COS ;
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

	// Bases in chi :
	// On check l'alternance :
	Index index_1 (a.bases_1d[1]->get_dimensions()) ;
	res.bases_1d[1] = new Array<int> (a.bases_1d[1]->get_dimensions()) ;
	do {
	switch ((*a.bases_1d[1])(index_1)) {
		case CHEB_EVEN:
			switch ((*b.bases_1d[1])(index_1)) {
				case CHEB_EVEN:
					res.bases_1d[1]->set(index_1) = (index_1(0)%2==0) ? CHEB_EVEN : CHEB_ODD ;
					break ;
				case CHEB_ODD:
					res.bases_1d[1]->set(index_1) = (index_1(0)%2==0) ? CHEB_ODD : CHEB_EVEN ;
					break ;
				default:
					res_def = false ;
					break ;
				}
			break ;
		case CHEB_ODD:
			switch ((*b.bases_1d[1])(index_1)) {
				case CHEB_EVEN:
					res.bases_1d[1]->set(index_1) = (index_1(0)%2==0) ? CHEB_ODD : CHEB_EVEN ;
					break ;
				case CHEB_ODD:
					res.bases_1d[1]->set(index_1) = (index_1(0)%2==0) ? CHEB_EVEN : CHEB_ODD ;
					break ;
				default:
					res_def = false ;
					break ;
				}
			break ;
		case LEG_EVEN:
			switch ((*b.bases_1d[1])(index_1)) {
				case LEG_EVEN:
					res.bases_1d[1]->set(index_1) = (index_1(0)%2==0) ? LEG_EVEN : LEG_ODD ;
					break ;
				case LEG_ODD:
					res.bases_1d[1]->set(index_1) = (index_1(0)%2==0) ? LEG_ODD : LEG_EVEN ;
					break ;
				default:
					res_def = false ;
					break ;
				}
			break ;
		case LEG_ODD:
			switch ((*b.bases_1d[1])(index_1)) {
				case LEG_EVEN:
					res.bases_1d[1]->set(index_1) = (index_1(0)%2==0) ? LEG_ODD : LEG_EVEN ;
					break ;
				case LEG_ODD:
					res.bases_1d[1]->set(index_1) = (index_1(0)%2==0) ? LEG_EVEN : LEG_ODD ;
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
	while (index_1.inc()) ;

	// Base in eta :
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


Val_domain Domain_bispheric_rect::der_normal (const Val_domain& so, int bound) const { 

	double rr, xc ;
	switch (bound) {
		case INNER_BC :
			xc = aa*cosh(eta_minus)/sinh(eta_minus) ;
			rr = aa/sinh(fabs(eta_minus)) ; 

			return Val_domain(so.der_abs(1)*(get_cart(1)-xc)/rr + so.der_abs(2)*get_cart(2)/rr +
				 so.der_abs(3)*get_cart(3)/rr) ;
		case OUTER_BC :
			return Val_domain(so.der_abs(1)*get_cart(1)/r_ext + so.der_abs(2)*get_cart(2)/r_ext +
				 so.der_abs(3)*get_cart(3)/r_ext) ;
		case CHI_ONE_BC :
			return Val_domain (so.der_var(2)/ (chi_min-M_PI)) ;
		case ETA_PLUS_BC :
			return Val_domain (so.der_var(1)/(eta_plus-eta_minus)*2.) ;
		default:
			cerr << "Unknown boundary case in Domain_bispheric_rect::der_normal" << endl ;
			abort() ;
		}
}

double Domain_bispheric_rect::integ (const Val_domain& so, int bound) const {

	if (bound !=INNER_BC) {
	  cerr << "Domain_bispheric_rect::integ only defined for inner boundary yet" << endl ;
	  abort() ;
	}

	double res = 0 ;

	int basep = (*so.base.bases_1d[2]) (0) ;
	if (basep == SIN) {
		// Odd function
		return res ;
	}
	else {
		
		// Multiply by the surface element :
		if (p_dsint==0x0)
			do_dsint() ;

		Val_domain auxi (so*(*p_dsint)) ;
		auxi.coef() ;

		Index pos (get_nbr_coefs()) ;
		pos.set(2) = 0 ;

		if (type_base==CHEB_TYPE) {
			for (int j=0 ; j<nbr_coefs(1) ; j++) {
				pos.set(1) = j ;
				int base_chi = (*auxi.base.bases_1d[1]) (0) ;
				double val_cheb ;
				switch (base_chi) {
					case CHEB_EVEN :
						if (j==0)
							val_cheb = 1. ;
						else
							val_cheb = double(2*j)/double(4*j*j-1) -1./double(2*j-1) ;
						break ;
					case CHEB_ODD :
						if (j==0) 
							val_cheb = 1./2. ;
						else
							if (j%2==0)
								val_cheb = 2*double(2*j+1)/double((2*j+1)*(2*j+1)-1) - 1./double(2*j) ;
							else
								val_cheb = -1./double(2*j) ;
						break ;
					default:
						cerr << "Unknown basis in Domain_bispheric_rect::integ" << endl ;
						abort() ;
				}
				for (int i=0 ; i<nbr_coefs(0) ; i++) {
					pos.set(0) = i ;
					if (i%2==0)
						res += val_cheb * (*auxi.cf)(pos) ;
					else
						res -= val_cheb * (*auxi.cf)(pos) ;
				}
			}
		}
	
		if (type_base==LEG_TYPE) {
			int base_chi = (*auxi.base.bases_1d[1]) (0) ;

			double val_leg = 1. ;
			double val_m1 = 1. ;
			double val_m = -0.5 ;

			for (int j=0 ; j<nbr_coefs(1) ; j++) {
				pos.set(1) = j ;
			
			switch (base_chi) {
					case LEG_EVEN :
						val_leg = 0. ;
						break ;
					case LEG_ODD :
						val_leg = (val_m1 - val_m)/double(4*j+3) ;
						val_m *= -double(2*j+3)/double(2*j+4) ;
						val_m1 *= -double(2*j+1)/double(2*j+2) ;
						break ;
					default:
						cerr << "Unknown basis in Domain_bispheric_rect::integ" << endl ;
						abort() ;
				}
		
			for (int i=0 ; i<nbr_coefs(0) ; i++) {
				pos.set(0) = i ;
				if (i%2==0)
					res += val_leg*(*auxi.cf)(pos) ;
				else
					res -= val_leg*(*auxi.cf)(pos) ;
			}
		}
	}
	return res*M_PI*2 ;
	}
}

int Domain_bispheric_rect::give_place_var (char* p) const {
    int res = -1 ;
    if (strcmp(p,"ETA ")==0)
	res = 0 ;
    if (strcmp(p,"CHI ")==0)
	res = 1 ;
    if (strcmp(p,"P ")==0)
	res = 2 ;
    return res ;
}}
