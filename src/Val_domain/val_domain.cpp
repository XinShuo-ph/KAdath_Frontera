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

#include "val_domain.hpp"
#include "utilities.hpp"
namespace Kadath {
Val_domain::Val_domain(const Domain* dom) : 
    zone{dom}, base{zone->get_ndim()}, is_zero{false}, c{nullptr}, cf{nullptr}, in_conf{false},
    in_coef{false}, p_der_var{zone->get_ndim()}, p_der_abs{zone->get_ndim()}
{
	for (int i=0 ; i<zone->get_ndim() ; i++) {
	     p_der_var[i] = nullptr ;
	     p_der_abs[i] = nullptr ;
     }
}

Val_domain::Val_domain (const Val_domain& so, bool const copie) :
    zone{so.zone}, base{so.base}, is_zero{so.is_zero}, in_conf{copie && so.in_conf}, in_coef{copie && so.in_coef},
    c{nullptr}, cf{nullptr}, p_der_var{zone->get_ndim()}, p_der_abs{zone->get_ndim()}
{
    if(copie){
        c  = (so.c  ? new Array<double> {*so.c}  : nullptr);
        cf = (so.cf ? new Array<double> {*so.cf} : nullptr);
    }
	for (int i=0 ; i<zone->get_ndim() ; i++) {
	     p_der_var[i] = ((so.p_der_var[i]!=nullptr) && copie) ? new Val_domain(*so.p_der_var[i]) : nullptr ;
	     p_der_abs[i] = ((so.p_der_abs[i]!=nullptr) && copie) ? new Val_domain(*so.p_der_abs[i]) : nullptr ;
     }
}

void Val_domain::swap(Val_domain &so) noexcept {
    assert(zone == so.zone);
    base.swap(so.base);
    std::swap(is_zero,so.is_zero);
    std::swap(c,so.c);
    std::swap(cf,so.cf);
    std::swap(in_conf,so.in_conf);
    std::swap(in_coef,so.in_coef);
    p_der_var.swap(so.p_der_var);
    p_der_abs.swap(so.p_der_abs);
}

Val_domain::Val_domain(Val_domain && so) noexcept :
    zone{so.zone}, base{std::move(so.base)}, is_zero{so.is_zero}, c{so.c}, cf{so.cf}, in_conf{so.in_conf},
    in_coef{so.in_coef}, p_der_var{std::move(so.p_der_var)}, p_der_abs{std::move(so.p_der_abs)}
{
    so.c = nullptr;
    so.cf = nullptr;
}

Val_domain & Val_domain::operator=(Val_domain && so) noexcept
{
    assert(zone = so.zone);
    base = std::move(so.base);
    is_zero = so.is_zero ;
    in_conf = so.in_conf ;
    in_coef = so.in_coef ;
    std::swap(c,so.c);
    std::swap(cf,so.cf);
    p_der_var = std::move(so.p_der_var);
    p_der_abs = std::move(so.p_der_abs);
    return *this;
}

Val_domain::Val_domain (const Domain* so, FILE* fd) :
    zone (so), base(fd), is_zero{}, c{}, cf{}, in_conf{}, in_coef{}, p_der_var{}, p_der_abs{}
{
	int indic ;
	fread_be (&indic, sizeof(int), 1, fd) ;
	is_zero = (indic==0) ? false : true ;
	fread_be (&indic, sizeof(int), 1, fd) ;
	in_conf = (indic==0) ? true : false ;
	c = (in_conf) ? new Array<double>(fd) : nullptr ;
	fread_be (&indic, sizeof(int), 1, fd) ;
	in_coef = (indic==0) ? true : false ;
	cf = (in_coef) ? new Array<double>(fd) : nullptr ;

    p_der_var.resize(zone->get_ndim()) ;
    p_der_abs.resize(zone->get_ndim()) ;
	for (int i=0 ; i<zone->get_ndim() ; i++) {
	     p_der_var[i] = nullptr ;
	     p_der_abs[i] = nullptr ;
     }
}

Val_domain::~Val_domain() {
	del_deriv();
//	if (c!=nullptr) delete c ;
//	if (cf!=nullptr) delete cf ;
    safe_delete(c);
    safe_delete(cf);
}

void Val_domain::save (FILE* fd) const {
	base.save(fd) ;
	int indic = (is_zero) ? 1 : 0 ;
	fwrite_be (&indic, sizeof(int), 1, fd) ;
	indic = (in_conf) ? 0 : 1 ;
	fwrite_be (&indic, sizeof(int), 1, fd) ;
	if (in_conf)
		c->save(fd) ;
	indic = (in_coef) ? 0 : 1 ;
	fwrite_be (&indic, sizeof(int), 1, fd) ;
	if (in_coef)
		cf->save(fd) ;
}

void Val_domain::del_deriv() const {
//    auto safe_delete = [](Val_domain * & p) -> void {if(p) {delete p; p = nullptr;}};
    p_der_var.apply(Safe_deleter<Val_domain>{});
    p_der_abs.apply(Safe_deleter<Val_domain>{});
}

void Val_domain::operator=(const Val_domain& so) {
	assert (zone == so.zone) ;
	is_zero = so.is_zero ;
	in_conf = so.in_conf ;
	in_coef = so.in_coef ;
	safe_delete(c);
	c = (in_conf) ? new Array<double> {*so.c} : nullptr ;
	safe_delete(cf);
	cf = (in_coef) ? new Array<double> (*so.cf) : nullptr ;
	if (so.base.is_def()) base=so.base ;
	else base.set_non_def() ;
	del_deriv() ;
	for (int i=0 ; i<zone->get_ndim() ; i++) {
	     p_der_var[i] = (so.p_der_var[i]==nullptr) ? nullptr : new Val_domain(*so.p_der_var[i]) ;
	     p_der_abs[i] = (so.p_der_abs[i]==nullptr) ? nullptr : new Val_domain(*so.p_der_abs[i]) ;
     }
}

void Val_domain::operator=(double xx) {
	if (xx==0)
		set_zero() ;
	else {
		is_zero = false ;
		allocate_conf() ;
		del_deriv() ;
		*c = xx ;
	}
}

void Val_domain::annule_hard() {
	allocate_conf() ;
	del_deriv() ;
	*c = 0. ;
}

double& Val_domain::set(const Index& index)  {
	set_in_conf() ;
	return c->set(index) ;
}

double& Val_domain::set_coef(const Index& index)  {
	set_in_coef() ;
	return cf->set(index) ;
}

double Val_domain::get_coef(const Index& index) const {
	if (is_zero)
		return 0. ;
	else {
		assert (in_coef) ;
		return cf->set(index) ;
	}
}

double Val_domain::operator()(const Index& index) const {
	coef_i() ;
	return (*c)(index) ;
}

void Val_domain::set_in_conf()  {
	safe_delete(cf);
	in_conf = true ;
	in_coef = false ;
}

void Val_domain::set_in_coef()  {
	safe_delete(c);
	in_conf = false ;
	in_coef = true ;
}

void Val_domain::allocate_conf() {
	set_in_conf() ;
	is_zero = false ;
	safe_delete(c);
	c = new Array<double>{zone->get_nbr_points()} ;
}

void Val_domain::allocate_coef() {
	set_in_coef() ;
	is_zero = false ;
	safe_delete(cf);
	cf = new Array<double>{zone->get_nbr_coefs()} ;
}

void Val_domain::set_zero() {
	if (!is_zero) {
		del_deriv() ;

		if (in_conf) {
		    assert(c);
			delete c ;
			c = nullptr ;
			in_conf = false ;
			}

		if (in_coef) {
            assert(cf);
			delete cf ;
			cf = nullptr ;
			in_coef = false ;
		}
		base.set_non_def() ;

	}
	is_zero = true ;
}

void Val_domain::std_base() {
	// recupere le type :
	int typeb = zone->get_type_base() ;
	switch (typeb) {
		case CHEB_TYPE :
			zone->set_cheb_base(base) ;
			break ;
		case LEG_TYPE :
			zone->set_legendre_base(base) ;
			break ;
		default:
			cerr << "Unknown type of base in Val_domain::std_base" ;
			abort() ;
	}
}

void Val_domain::std_r_base() {
	// recupere le type :
	int typeb = zone->get_type_base() ;
	switch (typeb) {
		case CHEB_TYPE :
			zone->set_cheb_r_base(base) ;
			break ;
		case LEG_TYPE :
			zone->set_legendre_r_base(base) ;
			break ;
		default:
			cerr << "Unknown type of base in Val_domain::std_base" ;
			abort() ;
	}
}


void Val_domain::std_anti_base() {
	// recupere le type :
	int typeb = zone->get_type_base() ;
	switch (typeb) {
		case CHEB_TYPE :
			zone->set_anti_cheb_base(base) ;
			break ;
		case LEG_TYPE :
			zone->set_anti_legendre_base(base) ;
			break ;
		default:
			cerr << "Unknown type of base in Val_domain::std_anti_base" ;
			abort() ;
	}
}

void Val_domain::std_base(int m) {
	// recupere le type :
	int typeb = zone->get_type_base() ;
	switch (typeb) {
		case CHEB_TYPE :
			zone->set_cheb_base_with_m(base,m) ;
			break ;
		case LEG_TYPE :
			zone->set_legendre_base_with_m(base,m) ;
			break ;
		default:
			cerr << "Unknown type of base in Val_domain::std_base" ;
			abort() ;
	}
}


void Val_domain::std_anti_base(int m) {
	// recupere le type :
	int typeb = zone->get_type_base() ;
	switch (typeb) {
		case CHEB_TYPE :
			zone->set_anti_cheb_base_with_m(base, m) ;
			break ;
		case LEG_TYPE :
			zone->set_anti_legendre_base_with_m(base, m) ;
			break ;
		default:
			cerr << "Unknown type of base in Val_domain::std_anti_base" ;
			abort() ;
	}
}

void Val_domain::std_base_rt_spher()
{
   int typeb(zone->get_type_base());
   switch (typeb)
   {
      case CHEB_TYPE :
         zone->set_cheb_base_rt_spher(base);
         break;
      default:
         cerr << "Unknown type of base in Val_domain::std_base_rt_spher";
         abort();
   }
}

void Val_domain::std_base_rp_spher()
{
   int typeb(zone->get_type_base());
   switch (typeb)
   {
      case CHEB_TYPE :
         zone->set_cheb_base_rp_spher(base);
         break;
      default:
         cerr << "Unknown type of base in Val_domain::std_base_rp_spher";
         abort();
   }
}

void Val_domain::std_base_tp_spher()
{
   int typeb(zone->get_type_base());
   switch (typeb)
   {
      case CHEB_TYPE :
         zone->set_cheb_base_tp_spher(base);
         break;
      default:
         cerr << "Unknown type of base in Val_domain::std_base_tp_spher";
         abort();
   }
}

void Val_domain::std_base_xy_cart()
{
   int typeb(zone->get_type_base());
   switch (typeb)
   {
      case CHEB_TYPE :
         zone->set_cheb_base_xy_cart(base);
         break;
      default:
         cerr << "Unknown type of base in Val_domain::std_base_xy_cart";
         abort();
   }
}

void Val_domain::std_base_xz_cart()
{
   int typeb(zone->get_type_base());
   switch (typeb)
   {
      case CHEB_TYPE :
         zone->set_cheb_base_xz_cart(base);
         break;
      default:
         cerr << "Unknown type of base in Val_domain::std_base_xz_cart";
         abort();
   }
}

void Val_domain::std_base_yz_cart()
{
   int typeb(zone->get_type_base());
   switch (typeb)
   {
      case CHEB_TYPE :
         zone->set_cheb_base_yz_cart(base);
         break;
      default:
         cerr << "Unknown type of base in Val_domain::std_base_yz_cart";
         abort();
   }
}

void Val_domain::std_base_r_spher() {
	// recupere le type :
	int typeb = zone->get_type_base() ;
	switch (typeb) {
		case CHEB_TYPE :
			zone->set_cheb_base_r_spher(base) ;
			break ;
		case LEG_TYPE :
			zone->set_legendre_base_r_spher(base) ;
			break ;
		default:
			cerr << "Unknown type of base in Val_domain::std_base_r_spher" ;
			abort() ;
	}
}

void Val_domain::std_base_t_spher() {
	// recupere le type :
	int typeb = zone->get_type_base() ;
	switch (typeb) {
		case CHEB_TYPE :
			zone->set_cheb_base_t_spher(base) ;
			break ;
		case LEG_TYPE :
			zone->set_legendre_base_t_spher(base) ;
			break ;
		default:
			cerr << "Unknown type of base in Val_domain::std_base_t_spher" ;
			abort() ;
	}
}

void Val_domain::std_base_p_spher() {
	// recupere le type :
	int typeb = zone->get_type_base() ;
	switch (typeb) {
		case CHEB_TYPE :
			zone->set_cheb_base_p_spher(base) ;
			break ;
		case LEG_TYPE :
			zone->set_legendre_base_p_spher(base) ;
			break ;
		default:
			cerr << "Unknown type of base in Val_domain::std_base_p_spher" ;
			abort() ;
	}
}

void Val_domain::std_base_x_cart() {
	// recupere le type :
	int typeb = zone->get_type_base() ;
	switch (typeb) {
		case CHEB_TYPE :
			zone->set_cheb_base_x_cart(base) ;
			break ;
		case LEG_TYPE :
			zone->set_legendre_base_x_cart(base) ;
			break ;
		default:
			cerr << "Unknown type of base in Val_domain::std_base_x_cart" ;
			abort() ;
	}
}

void Val_domain::std_base_y_cart() {
	// recupere le type :
	int typeb = zone->get_type_base() ;
	switch (typeb) {
		case CHEB_TYPE :
			zone->set_cheb_base_y_cart(base) ;
			break ;
		case LEG_TYPE :
			zone->set_legendre_base_y_cart(base) ;
			break ;
		default:
			cerr << "Unknown type of base in Val_domain::std_base_y_cart" ;
			abort() ;
	}
}

void Val_domain::std_base_z_cart() {
	// recupere le type :
	int typeb = zone->get_type_base() ;
	switch (typeb) {
		case CHEB_TYPE :
			zone->set_cheb_base_z_cart(base) ;
			break ;
		case LEG_TYPE :
			zone->set_legendre_base_z_cart(base) ;
			break ;
		default:
			cerr << "Unknown type of base in Val_domain::std_base_z_cart" ;
			abort() ;
	}
}

void Val_domain::std_base_r_mtz() {
	// recupere le type :
	int typeb = zone->get_type_base() ;
	switch (typeb) {
		case CHEB_TYPE :
			zone->set_cheb_base_r_mtz(base) ;
			break ;
		case LEG_TYPE :
			zone->set_legendre_base_r_mtz(base) ;
			break ;
		default:
			cerr << "Unknown type of base in Val_domain::std_base_r_mtz" ;
			abort() ;
	}
}

void Val_domain::std_base_t_mtz() {
	// recupere le type :
	int typeb = zone->get_type_base() ;
	switch (typeb) {
		case CHEB_TYPE :
			zone->set_cheb_base_t_mtz(base) ;
			break ;
		case LEG_TYPE :
			zone->set_legendre_base_t_mtz(base) ;
			break ;
		default:
			cerr << "Unknown type of base in Val_domain::std_base_t_mtz" ;
			abort() ;
	}
}

void Val_domain::std_base_p_mtz() {
	// recupere le type :
	int typeb = zone->get_type_base() ;
	switch (typeb) {
		case CHEB_TYPE :
			zone->set_cheb_base_p_mtz(base) ;
			break ;
		case LEG_TYPE :
			zone->set_legendre_base_p_mtz(base) ;
			break ;
		default:
			cerr << "Unknown type of base in Val_domain::std_base_p_mtz" ;
			abort() ;
	}
}


void Val_domain::std_xodd_base() {
	// recupere le type :
	int typeb = zone->get_type_base() ;
	switch (typeb) {
		case CHEB_TYPE :
			zone->set_cheb_xodd_base(base) ;
			break ;
		case LEG_TYPE :
			zone->set_legendre_xodd_base(base) ;
			break ;
		default:
			cerr << "Unknown type of base in Val_domain::std_xodd_base" ;
			abort() ;
	}
}

void Val_domain::std_todd_base() {
	// recupere le type :
	int typeb = zone->get_type_base() ;
	switch (typeb) {
		case CHEB_TYPE :
			zone->set_cheb_todd_base(base) ;
			break ;
		case LEG_TYPE :
			zone->set_legendre_todd_base(base) ;
			break ;
		default:
			cerr << "Unknown type of base in Val_domain::std_todd_base" ;
			abort() ;
	}
}

void Val_domain::std_xodd_todd_base() {
	// recupere le type :
	int typeb = zone->get_type_base() ;
	switch (typeb) {
		case CHEB_TYPE :
			zone->set_cheb_xodd_todd_base(base) ;
			break ;
		case LEG_TYPE :
			zone->set_legendre_xodd_todd_base(base) ;
			break ;
		default:
			cerr << "Unknown type of base in Val_domain::std_xodd_todd_base" ;
			abort() ;
	}
}

void Val_domain::std_base_odd() {
	// recupere le type :
	int typeb = zone->get_type_base() ;
	switch (typeb) {
		case CHEB_TYPE :
			zone->set_cheb_base_odd(base) ;
			break ;
		case LEG_TYPE :
			zone->set_legendre_base_odd(base) ;
			break ;
		default:
			cerr << "Unknown type of base in Val_domain::std_base_odd" ;
			abort() ;
	}
}


void Val_domain::coef() const {
	if ((in_coef) || (is_zero)) {
		return ;
    	}
	else {
	    if (base.is_def()==false) {
		cout << "Base not defined in Val_domain::coef" << endl ;
		abort() ;
		}
	    else {
   	    assert(in_conf) ;
            if (cf !=nullptr)
	        delete cf ;
	    cf = new Array<double>(base.coef(zone->get_nbr_coefs(), *c)) ;
	    }
	    in_coef = true ;
	}
}

void Val_domain::coef_i() const {
	if ((in_conf) || (is_zero)) {
		return ;
    	}
	else {
	    if (base.is_def()==false) {
		cout << "Base not defined in Val_domain::coef" << endl ;
		abort() ;
		}
	    else {
   	    assert(in_coef) ;
            if (c !=nullptr)
	        delete c ;
	    c = new Array<double>(base.coef_i(zone->get_nbr_points(), *cf)) ;
	    }
    	in_conf = true ;
	}
}

ostream& operator<< (ostream& o, const Val_domain& so) {

    if (so.is_zero)
	o << "Null Val_domain" << endl ;

    if (so.in_conf) {
    	o << "Configuration space : " << endl ;
	o << *so.c<< endl ;
    }

    if (so.in_coef) {
    	o << "Coefficient space : " << endl ;
	o << *so.cf << endl ;
	}
    return o ;
}

Val_domain Val_domain::der_var(int var) const {
	assert ((var>0) && (var<=zone->get_ndim())) ;
	if (is_zero)
		return *this ;
	else {
		if (p_der_var[var-1] == nullptr)
			compute_der_var() ;
		return *p_der_var[var-1] ;
	}
}

Val_domain Val_domain::der_abs(int var) const
{
   assert ((var>=0) && (var<=zone->get_ndim())) ;
   if (var == 0)
   {
      Val_domain zero(get_domain());
      zero = 0.0;
      return zero;
   }
   if (is_zero)
      return *this ;
   else {
      if (p_der_abs[var-1] == nullptr)
         compute_der_abs() ;
      return *p_der_abs[var-1] ;
   }
}

Val_domain Val_domain::der_spher(int var) const
{
   assert ( (var > 0 ) and ( var <= zone->get_ndim() ) );
   if (is_zero)
      return *this;
   else
   {
      Val_domain zero(zone);
      zero = 0.0;
      if (var == 0)
         return zero;
      else if (var == 1)
         return this->der_r();
      else if (var == 2)
         return this->der_t();
      else if (var == 0)
         return this->der_p();
      else
      {
         cerr << "bad index in der_spher" << endl;
         abort();
      }
   }
}

Val_domain Val_domain::der_r() const {
	return zone->der_r (*this) ;
}

Val_domain Val_domain::der_t() const {
	return zone->der_t (*this) ;
}

Val_domain Val_domain::der_p() const {
	return zone->der_p (*this) ;
}


Val_domain Val_domain::der_r_rtwo() const {
	return zone->der_r_rtwo (*this) ;
}

Val_domain Val_domain::div_r() const {
  if (is_zero)
    return *this ;
  else
    return zone->div_r (*this) ;
}

Val_domain Val_domain::mult_r() const {
  if (is_zero)
    return *this ;
  else
    return zone->mult_r (*this) ;
}
double Val_domain::integrale() const {
	if (is_zero)
		return 0. ;
	else
	  return zone->integrale(*this) ;
}

double Val_domain::integ_volume() const {
	if (is_zero)
		return 0. ;
	else
	  return zone->integ_volume(*this) ;
}

int der_1d (int, Array<double>&) ;
void Val_domain::compute_der_var () const {
	coef() ;
	for (int var=0 ; var<zone->get_ndim() ; var++) {
	    Val_domain res(zone) ;
	    res.base = base ;
	    res.cf = new Array<double>{base.ope_1d(der_1d, var, *cf, res.base)} ;
	    res.in_coef = true ;
	    safe_delete(p_der_var[var]);
	    p_der_var[var] = new Val_domain{res} ;
	}
}

void Val_domain::compute_der_abs () const {
	for (int i=0 ; i<zone->get_ndim() ; i++)
		if (p_der_var[i]==nullptr)
			compute_der_var() ;
	zone->do_der_abs_from_der_var(p_der_var.get_data(), p_der_abs.set_data()) ;
}

}

