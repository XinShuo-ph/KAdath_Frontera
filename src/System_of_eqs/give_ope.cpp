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

#include "system_of_eqs.hpp"
#include "name_tools.hpp"
#include "term_eq.hpp" 
#include "base_tensor.hpp"
namespace Kadath {
bool System_of_eqs::isvar_double (const char* name, int& which) const {
	bool res = false ;
	char auxi[LMAX] ;
	trim_spaces (auxi, name) ;
	for (int i=0 ; i<nvar_double ; i++) {
		if (!res) 
		if (strcmp(auxi, names_var_double[i])==0) {
			res = true ;
			which = i ;
		}
	}
	return res ;
}

bool System_of_eqs::isvar (const char* name, int& which, int& valence, char*& name_ind, Array<int>*& type_ind) const {
	bool res = false ;
	for (int i=0 ; i<nvar ; i++) {
		if (name_ind !=0x0) {
		    delete [] name_ind ;
		    name_ind = 0x0 ;
		}
		if (type_ind !=0x0) {
		  delete type_ind ;
		  type_ind = 0x0 ;
		}
		res = is_tensor(name, names_var[i], valence, name_ind, type_ind) ;
		if (res) {
			which = i ;
			if (valence != var[which]->get_valence()) {
				cerr << "Bad valence for " << name << endl ;
				abort() ;
			}
			break ;
		}
	}
	return res ;
}

bool System_of_eqs::iscst (const char* name, int& which, int& valence, char*& name_ind, Array<int>*& type_ind) const {
	bool res = false ;
	for (int i=0 ; i<ncst ; i++) {

		// Type of constant
		int ind = i*ndom ;
		int type_cst = cst[ind]->get_type_data() ;
	
		switch (type_cst) {
			case TERM_T :
			  if (name_ind !=0x0) {
			    delete [] name_ind ;
			    name_ind = 0x0 ;
			}
			if (type_ind !=0x0) {
			  delete type_ind ;
			  type_ind = 0x0 ;
			}
				res = is_tensor(name, names_cst[i], valence, name_ind, type_ind) ;
				if (res) {
					which = i ;
					if (valence != cst[ind]->get_val_t().get_valence()) {
					cerr << "Bad valence for " << name << endl ;
					abort() ;
					}
				}
				break ;
			case TERM_D :
				if (strcmp(names_cst[i], name)==0) {
					res = true ;
					which = i ;
				}
				break ;
			default :
				cerr << "Unknown type of data in System_of_eqs::iscst" << endl ;
				abort() ;
			}
		if (res)
			break ;
		}
	return res ;
}

bool System_of_eqs::isdef (int dd, const char* name, int& which, int& valence, char*& name_ind, Array<int>*& type_ind) const {
	bool res = false ;
	for (int i=0 ; i<ndef ; i++) {
		if (name_ind !=0x0) {
		    delete [] name_ind ;
		    name_ind = 0x0 ;
		}
		if (type_ind !=0x0) {
		  delete type_ind ;
		  type_ind = 0x0 ;
		}
		res = is_tensor(name, names_def[i], valence, name_ind, type_ind) ;
		if ((res) && (def[i]->get_dom()==dd)) {
			which = i ;
			break ;
		}
		res = false ;
	}
	return res ;
}

bool System_of_eqs::isdef_glob (int dd, const char* name, int& which) const {
	bool res = false ;
	for (int i=0 ; i<ndef_glob ; i++) {
	    if ((strcmp(names_def_glob[i], name)==0) && (def_glob[i]->get_dom()==dd)) {
		res = true ;
		which = i ;
		break ;
	 }
	}
	return res ;
}


bool System_of_eqs::isdouble (const char*name, double& val) const {
	char* error ;
	val = strtod (name, &error);
	bool res = ((*error != ' ') || (strlen(error)>1)) ? false : true ;
	return res ;
}

bool System_of_eqs::ismet (const char* name, char*& name_ind, int& type_indice) const {
	if (met==0x0)
		return false ;
	int valence ;
	Array<int>* type_ind = 0x0 ;
	bool res = is_tensor(name, name_met, valence, name_ind, type_ind) ;
	if (res) {
		if (valence !=2) {
			cerr << "Bad valence for the metric " << name_met << " in " << name << endl ;
			abort() ;
		}

		if ((*type_ind)(0)!=(*type_ind)(1)) {
			cerr << "Indices of the metric must be of the same type" << endl ;
			abort() ;
		}
		
		type_indice = (*type_ind)(0) ;
		delete type_ind ;
	}
	return res ;
}

bool System_of_eqs::ismet (const char* name) const {

	char auxi[LMAX] ;
	trim_spaces(auxi, name) ;

	if (met==0x0)
		return false ;
	int same = strcmp(auxi, name_met) ;

	if (same!=0)
	    return false;
	else
	  return true ;
}

bool System_of_eqs::ischristo (const char* name, char*& name_ind, Array<int>*& type_ind) const {
	
	if (met==0x0)
		return false ;
	int valence ;
	bool res = is_tensor(name, "Gam ", valence, name_ind, type_ind) ;
	if (res) {
		if (valence != 3) {
			cerr << "Bad valence for Christoffel symbol in " << name << endl ;
			abort() ;
		}

		if (((*type_ind)(0)!=COV) || ((*type_ind)(1)!=COV) || ((*type_ind)(2)!=CON)) {
			cerr << "Indices of the wrong type in " << name << endl ;
			abort() ;
		}
	}
	return res ;
}

bool System_of_eqs::isriemann (const char* name, char*& name_ind, Array<int>*& type_ind) const {
	
	if (met==0x0)
		return false ;
	int valence ;
	bool res = is_tensor(name, "R ", valence, name_ind, type_ind) ;
	if (res)
	  if (valence !=4)
	    res = false ;
	if (res) {

		if (((*type_ind)(0)!=CON) || ((*type_ind)(1)!=COV) || ((*type_ind)(2)!=COV) || ((*type_ind)(3)!=COV)) {
		  delete [] name_ind ;
		  delete type_ind ;
		  res = false ;
		}
	}

	if (!res) {
		if (name_ind !=0x0)
			delete [] name_ind ;
		if (type_ind !=0x0)
			delete type_ind ;
	}
	return res ;
}


bool System_of_eqs::isricci_tensor (const char* name, char*& name_ind, Array<int>*& type_ind) const {
	
	if (met==0x0)
		return false ;
	int valence ;
	bool res = is_tensor(name, "R ", valence, name_ind, type_ind) ;
	if ((res) && (valence != 2)) {
		 delete [] name_ind ;
		delete type_ind ;
		res = false ;
	}

	if (!res) {
		if (name_ind !=0x0)
			delete [] name_ind ;
		if (type_ind !=0x0)
			delete type_ind ;
	}
	return res ;
}

bool System_of_eqs::isricci_scalar (const char* name, char*& name_ind, Array<int>*& type_ind) const {
	
	if (met==0x0)
		return false ;
	int valence ;
	bool res = is_tensor(name, "R ", valence, name_ind, type_ind) ;
	if ((res) && (valence != 0)) { 
		delete [] name_ind ;
		delete type_ind ;
		res = false ;
	}

	if (!res) {
		if (name_ind !=0x0)
			delete [] name_ind ;
		if (type_ind !=0x0)
			delete type_ind ;
	}
	return res ;
}


bool System_of_eqs::is_ope_bin (const char* name, char* part1, char* part2, char cible) const {
	assert ((cible=='+') || (cible=='-') || (cible=='*') || (cible=='/') || (cible=='=')) ;

	bool res = false ;
	int nbr_cible = nbr_char(name, cible) ;
	int conte = 0 ;

	while ((!res) && (conte<nbr_cible)) {
		get_parts(name, part1, part2, cible, conte) ;
		if ((part1[0] !='\0') && (part2[0]!='\0'))
		  if ((nbr_char(part1, '(')==nbr_char(part1, ')')) && (nbr_char(part2, '(')==nbr_char(part2, ')')))
			res = true ;
		conte ++ ;
	  }
	return res ;
}

bool System_of_eqs::is_ope_minus (const char* name, char* part) const {

  // Check the beginning
  char auxi[LMAX] ;
  trim_spaces(auxi, name) ;
  int len = static_cast<int>(strlen(auxi)) ;

  bool res = true ;
  if (auxi[0] != '-')
    res = false ;
  if (res) {
    char occi[LMAX] ;
    for (int i=1 ; i<len-1 ; i++)
      occi[i-1] = auxi[i] ;
    occi[len-2] = ' ' ;
    occi[len-1] = '\0' ;
    trim_spaces(part, occi) ;
  }
  return res ;
}

bool System_of_eqs::is_ope_uni (const char* name, char* part, const char* zeope) const {

	// Check the beginning
	char auxi[LMAX] ;
	trim_spaces(auxi, zeope) ;
	int len = static_cast<int>(strlen(auxi)) ;

	bool res = true ;
	for (int i=0 ; i<len-1 ; i++)
		if (name[i]!=auxi[i]) {
			res = false ;
			break ;
		}

	if (res)
		if (name[len-1]!='(')
			res = false ;

	// check the end
	if (res) {
		int lentot = static_cast<int>(strlen(name)) ;
		assert (name[lentot-1]==' ') ;
		if (name[lentot-2]!=')')
			res = false ;

	for (int i=len ; i<lentot-2 ; i++)
		part[i-len] = name[i] ;

	part[lentot-len-2] = ' ' ;
	part[lentot-len-1] = '\0' ;
	}

	return res ;
}

bool System_of_eqs::is_ope_uni (const char* name, char* p1, char* p2, const char* zeope) const {

	// Check the beginning
	char auxi[LMAX] ;
	trim_spaces(auxi, zeope) ;
	int len = static_cast<int>(strlen(auxi)) ;
	bool res = true ;
	for (int i=0 ; i<len-2 ; i++)
		if (name[i]!=auxi[i]) {
			res = false ;
			break ;
		}
	if (res) 
	if (name[len-1]!='(')
		res = false ;

	// check the end
	if (res) {

	char part[LMAX] ;
	int lentot = static_cast<int>(strlen(name)) ;
	assert (name[lentot-1]==' ') ;
	if (name[lentot-2]!=')')
		res = false ;

	for (int i=len ; i<lentot-2 ; i++)
		part[i-len] = name[i] ;

	part[lentot-len-2] = ' ' ;
	part[lentot-len-1] = '\0' ;

	get_parts(part, p1, p2, ',') ;
	}
	return res ;
}

bool System_of_eqs::is_ope_deriv (const char* nn, char* part, int& type_der, char& name_ind) const {
	char name[LMAX] ;
	trim_spaces(name, nn) ;
	bool res = true ;
	int len = static_cast<int>(strlen(name)) ;
	// Check first char and last
	if (name[0]!='D')
		res = false ;
	// Check type of derivative
	if (res) {
		if (name[1]=='_')
			type_der = COV ;
		else {
			if (name[1]=='^')
				type_der = CON ;
			else
				res = false ;
		      }
	}

	if (res) {
		// Get the indice :
		name_ind = name[2] ;

		// Get the rest :
		char auxi[LMAX] ;
		for (int i=3 ; i<len ; i++)
			auxi[i-3] = name[i] ;
		auxi[len-3] = '\0' ;
		trim_spaces(part, auxi) ;
	}
	return res ;
}

bool System_of_eqs::is_ope_deriv_flat (const char* nn, char* part, int& type_der, char& name_ind) const {
	char name[LMAX] ;
	trim_spaces(name, nn) ;
	bool res = true ;
	int len = static_cast<int>(strlen(name)) ;
	// Check first char and last
	if (name[0]!='D')
		res = false ;
	if (name[1]!='F')
		res = false ;
	// Check type of derivative
	if (res) {
		if (name[2]=='_')
			type_der = COV ;
		else {
			if (name[2]=='^')
				type_der = CON ;
			else
				res = false ;
		      }
	}

	if (res) {
		// Get the indice :
		name_ind = name[3] ;

		// Get the rest :
		char auxi[LMAX] ;
		for (int i=4 ; i<len ; i++)
			auxi[i-4] = name[i] ;
		auxi[len-4] = '\0' ;
		trim_spaces(part, auxi) ;
	}
	return res ;
}

bool System_of_eqs::is_ope_deriv_background (const char* nn, char* part, int& type_der, char& name_ind) const {
	char name[LMAX] ;
	trim_spaces(name, nn) ;
	bool res = true ;
	int len = static_cast<int>(strlen(name)) ;
	// Check first char and last
	if (name[0]!='D')
		res = false ;
	if (name[1]!='B')
		res = false ;
	// Check type of derivative
	if (res) {
		if (name[2]=='_')
			type_der = COV ;
		else {
			if (name[2]=='^')
				type_der = CON ;
			else
				res = false ;
		      }
	}

	if (res) {
		// Get the indice :
		name_ind = name[3] ;

		// Get the rest :
		char auxi[LMAX] ;
		for (int i=4 ; i<len ; i++)
			auxi[i-4] = name[i] ;
		auxi[len-4] = '\0' ;
		trim_spaces(part, auxi) ;
	}
	return res ;
}

bool System_of_eqs::is_ope_partial (const char* nn, char* part, char& name_ind) const {
	char name[LMAX] ;
	trim_spaces(name, nn) ;
	bool res = true ;
	int len = static_cast<int>(strlen(name)) ;

	// Check name
	if (strncmp(name,"partial", 7) != 0)
		res = false ;
	// Check type of derivative
	if (res) 
		if (name[7]!='_')
			res = false ;
	
	if (res) {
		// Get the indice :
		name_ind = name[8] ;

		// Get the rest :
		char auxi[LMAX] ;
		for (int i=9 ; i<len ; i++)
			auxi[i-9] = name[i] ;
		auxi[len-9] = '\0' ;
		trim_spaces(part, auxi) ;
	}
	return res ;
}

bool System_of_eqs::is_ope_pow (const char* name, char* part, int& expo) const {
	char auxi[LMAX] ;
	trim_spaces(auxi, name) ;

	bool res = (nbr_char(auxi, '^')>=1) ? true : false ;

	if (res) {
		
		char p2[LMAX] ;
		get_parts(auxi, part, p2, '^') ;

		// Is p2 an integer ?
		char* error ;
		expo = static_cast<int>(strtol (p2, &error, 0));
		res = ((*error != ' ') || (strlen(error)>1)) ? false : true ;
	}

	return res ;	
}

bool System_of_eqs::is_ope_der_var (int dd, const char*name, char*part, int& numvar) const {

    char auxi[LMAX] ;
    trim_spaces(auxi, name) ;

    bool res = (nbr_char(name, ',')==1) ? true : false ;
    if (res) {
      char p2[LMAX] ;
      get_parts(auxi, part, p2, ',') ;
      
      // Check if p2 is the name of a variable 
     int isitvar = espace.get_domain(dd)->give_place_var (p2) ;
    if (isitvar==-1)
      res = false ;
    else
      numvar = isitvar ;
    }
  return res ;
}

Ope_eq* System_of_eqs::give_ope (int dd, const char* name, int bound) const {
	Ope_eq* p_ope = 0x0 ;
	bool indic ;
	int which = -1 ;
	char p1[LMAX] ;
	char p2[LMAX] ;

	// Check if addition :
	indic = is_ope_bin(name, p1, p2, '+') ;
	if (indic) {
		p_ope = new Ope_add(this, give_ope(dd,p1, bound), give_ope(dd,p2, bound)) ;
		return p_ope ;
		}

	// Check if substraction :
	indic = is_ope_bin(name, p1, p2, '-') ;
	if (indic) {
		p_ope = new Ope_sub(this, give_ope(dd,p1,bound), give_ope(dd,p2, bound)) ;
		return p_ope ;
		}

	// Check if minus
	indic = is_ope_minus(name, p1) ;
	if (indic) {
		p_ope = new Ope_minus (this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	// Check if multiplication :
	indic = is_ope_bin(name, p1, p2, '*') ;
	if (indic) {
		p_ope = new Ope_mult(this, give_ope(dd,p1, bound), give_ope(dd,p2, bound)) ;
		return p_ope ;
		}

	// Check if division :
	indic = is_ope_bin(name, p1, p2, '/') ;
	if (indic) {
		p_ope = new Ope_div(this, give_ope(dd,p1, bound), give_ope(dd,p2, bound)) ;
		return p_ope ;
		}

	// Check if laplacian :
	indic = is_ope_uni(name, p1, "Lap") ;
	if (indic) {
		p_ope = new Ope_lap(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	indic = is_ope_uni(name, p1, "lap") ;
	if (indic) {
		p_ope = new Ope_lap(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	indic = is_ope_uni(name, p1, "Lap2") ;
	if (indic) {
		p_ope = new Ope_lap2(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	indic = is_ope_uni(name, p1, "lap2") ;
	if (indic) {
		p_ope = new Ope_lap2(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}
	// Check if dn :
	indic = is_ope_uni(name, p1, "dn") ;
	if (indic) {	
		p_ope = new Ope_dn(this, bound, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}
		
	  // Check if mult by x :
	indic = is_ope_uni(name, p1, "multx") ;
	if (indic) {
		p_ope = new Ope_mult_x(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	// Check if mult by r :
	indic = is_ope_uni(name, p1, "multr") ;
	if (indic) {
		p_ope = new Ope_mult_r(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	// Check if mult by 1 - r/L :
	indic = is_ope_uni(name, p1, "mult1mrsL") ;
	if (indic) {
		p_ope = new Ope_mult_1mrsL(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	// Check if div by 1 - r/L :
	indic = is_ope_uni(name, p1, "div1mrsL") ;
	if (indic) {
		p_ope = new Ope_div_1mrsL(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	indic = is_ope_uni(name, p1, "valori") ;
	if (indic) {
		p_ope = new Ope_val_ori(this, dd, give_ope(0,p1, bound)) ;
		return p_ope ;
		}

	// Check if srdr :
	indic = is_ope_uni(name, p1, "srdr") ;
	if (indic) {
		p_ope = new Ope_srdr(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	// Check if ddp :
	indic = is_ope_uni(name, p1, "ddp") ;
	if (indic) {
		p_ope = new Ope_ddp(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	// Check if ddt :
	indic = is_ope_uni(name, p1, "ddt") ;
	if (indic) {
		p_ope = new Ope_ddt(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	// Check if dt :
	indic = is_ope_uni(name, p1, "dt") ;
	if (indic) {
		p_ope = new Ope_dt(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	// Check if dt :
	indic = is_ope_uni(name, p1, "dtime") ;
	if (indic) {
		p_ope = new Ope_dtime(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}
	
	// Check if ddt
	indic = is_ope_uni(name, p1, "ddtime") ;
	if (indic) {
		p_ope = new Ope_ddtime (this, give_ope(dd, p1, bound)) ;
		return p_ope ;
		}


	// Check if ddr :
	indic = is_ope_uni(name, p1, "ddr") ;
	if (indic) {
		p_ope = new Ope_ddr(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	// Check if dr :
	indic = is_ope_uni(name, p1, "dr") ;
	if (indic) {
		p_ope = new Ope_dr(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	// Check if division by r:
	indic = is_ope_uni(name, p1, "divr") ;
	if (indic) {
		p_ope = new Ope_div_r(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	// Check if multiplication by r sint:
	indic = is_ope_uni(name, p1, "multrsint") ;
	if (indic) {
		p_ope = new Ope_mult_rsint(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	// Check if division by r sint:
	indic = is_ope_uni(name, p1, "divrsint") ;
	if (indic) {
		p_ope = new Ope_div_rsint(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	// Check if division by sint:
	indic = is_ope_uni(name, p1, "divsint") ;
	if (indic) {
		p_ope = new Ope_div_sint(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	// Check if division by cost:
	indic = is_ope_uni(name, p1, "divcost") ;
	if (indic) {
		p_ope = new Ope_div_cost(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	// Check if multiplication by sint:
	indic = is_ope_uni(name, p1, "multsint") ;
	if (indic) {
		p_ope = new Ope_mult_sint(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	// Check if division by xpone:
	indic = is_ope_uni(name, p1, "divxpone") ;
	if (indic) {
		p_ope = new Ope_div_xpone(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	// Check if division by xpone:
	indic = is_ope_uni(name, p1, "div1mx2") ;
	if (indic) {
		p_ope = new Ope_div_1mx2(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	// Check if integ :
	indic = is_ope_uni(name, p1, "integ") ;
	if (indic) {
		p_ope = new Ope_int(this, bound, give_ope(dd,p1,bound)) ;
		return p_ope ;
		}

	// Check if integ :
	indic = is_ope_uni(name, p1, "integvolume") ;
	if (indic) {
		p_ope = new Ope_int_volume(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}
	
	// Check if grad :
	indic = is_ope_uni(name, p1, "grad") ;
	if (indic) {
		p_ope = new Ope_grad(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	// Check if sqrt :
	indic = is_ope_uni(name, p1, "sqrt") ;
	if (indic) {
		p_ope = new Ope_sqrt(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	// Check if sqrt rho basis :
	indic = is_ope_uni(name, p1, "sqrtrho") ;
	if (indic) {
		p_ope = new Ope_sqrt_nonstd(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	// Check if sqrt anti basis :
	indic = is_ope_uni(name, p1, "sqrtanti") ;
	if (indic) {
		p_ope = new Ope_sqrt_anti(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}


	// Check if exp :
	indic = is_ope_uni(name, p1, "exp") ;
	if (indic) {
		p_ope = new Ope_exp(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	// Check if log :
	indic = is_ope_uni(name, p1, "log") ;
	if (indic) {
		p_ope = new Ope_log(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}
	
	// Check if atanh :
	indic = is_ope_uni(name, p1, "atanh") ;
	if (indic) {
		p_ope = new Ope_atanh(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	// Check if cos :
	indic = is_ope_uni(name, p1, "cos") ;
	if (indic) {
		p_ope = new Ope_cos(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	// Check if sin :
	indic = is_ope_uni(name, p1, "sin") ;
	if (indic) {
		p_ope = new Ope_sin(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	// Check if cosh :
	indic = is_ope_uni(name, p1, "cosh") ;
	if (indic) {
		p_ope = new Ope_cosh(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	// Check if sinh :
	indic = is_ope_uni(name, p1, "sinh") ;
	if (indic) {
		p_ope = new Ope_sinh(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	// Check if scalar product :
	indic = is_ope_uni(name, p1, p2, "scal") ;
	if (indic) {
		p_ope = new Ope_scal(this, give_ope(dd,p1, bound), give_ope(dd,p2, bound)) ;
		return p_ope ;
		}

	// Check if complex conjugate
	indic = is_ope_uni(name, p1, "conjug") ;
	if (indic) {
		p_ope = new Ope_conjug(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}
		
	// Check if determinant :
	indic = is_ope_uni(name, p1, "determinant") ;
	if (indic) {
		p_ope = new Ope_determinant(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	// Check if inverse :
	indic = is_ope_uni(name, p1, "inverse") ;
	if (indic) {
		p_ope = new Ope_inverse(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	// Check if inversenodet :
	indic = is_ope_uni(name, p1, "inversenodet") ;
	if (indic) {
		p_ope = new Ope_inverse_nodet(this, give_ope(dd,p1, bound)) ;
		return p_ope ;
		}

	// Check if fit waves
	indic = is_ope_uni(name, p1, p2, "fitwaves") ;
	if (indic) {
	    p_ope = new Ope_fit_waves(this, give_ope(dd,p1, bound), give_ope(dd,p2, bound)) ;
	    return p_ope ;
	}

	// Check if import
	indic = is_ope_uni(name, p1, "import") ;
	if (indic) {
	    p_ope = new Ope_import(this, dd, bound, p1) ;
	    return p_ope ;
	}

	// Check if to spherical tensorial basis
	indic = is_ope_uni(name, p1, "tospherical") ;
	if (indic) {
	    p_ope = new Ope_change_basis(this, SPHERICAL_BASIS , give_ope(dd, p1, bound)) ;
	    return p_ope ;
	}
	
	// Check if to cartesian tensorial basis
	indic = is_ope_uni(name, p1, "tocartesian") ;
	if (indic) {
	    p_ope = new Ope_change_basis(this, CARTESIAN_BASIS , give_ope(dd, p1, bound)) ;
	    return p_ope ;
	}
	
	// Check if derivative :
	int type_der ;
	char ind_der ;
	indic = is_ope_deriv(name, p1, type_der, ind_der) ;
	if (indic) {
		 p_ope = new Ope_der (this, type_der, ind_der, give_ope(dd,p1, bound)) ;
		return p_ope ;
	}
	
	// Check if flat derivative :
	indic = is_ope_deriv_flat(name, p1, type_der, ind_der) ;
	if (indic) {
		 p_ope = new Ope_der_flat (this, type_der, ind_der, give_ope(dd,p1, bound)) ;
		return p_ope ;
	}
		
	// Check if background derivative :
	indic = is_ope_deriv_background(name, p1, type_der, ind_der) ;
	if (indic) {
		 p_ope = new Ope_der_background (this, type_der, ind_der, give_ope(dd,p1, bound)) ;
		return p_ope ;
	}
	

	// Check if partial :
	indic = is_ope_partial(name, p1, ind_der) ;
	if (indic) {
		p_ope = new Ope_partial (this, ind_der, give_ope(dd,p1, bound)) ;
		return p_ope ;
	}

	// Check if partial derivative with respect to one particular variable :
	int ind_var ;
	indic = is_ope_der_var (dd, name, p1, ind_var) ;
	if (indic) {
	    p_ope = new Ope_partial_var (this, ind_var, give_ope(dd, p1, bound)) ;
	    return p_ope ;
	}

	// Check if exponent :
	int val_exp ;
	indic = is_ope_pow (name, p1, val_exp) ;
	if (indic) {
		p_ope = new Ope_pow(this, val_exp, give_ope(dd, p1, bound)) ;
		return p_ope ;
	}

	// Check if double variable :
	indic = isvar_double (name, which) ;
	if (indic) {
		p_ope = new Ope_id (this, give_term_double(which, dd)) ;
		return p_ope ;
		}

	// Check if variable :
	int valence ;
	char* name_ind = 0x0 ;
	Array<int>* type_ind = 0x0 ;
	indic = isvar (name, which, valence, name_ind, type_ind) ;
	if (indic) {
		p_ope = new Ope_id (this, give_term(which, dd), valence, name_ind, type_ind) ;
		return p_ope ;
		}


	// Check if metric 
	name_ind = 0x0 ;
	int type_indice ;
	indic = ismet (name, name_ind, type_indice) ;
	if (indic) {
		type_ind = new Array<int> (2) ;
		type_ind->set(0) = type_indice ; type_ind->set(1) = type_indice ;
		p_ope = new Ope_id (this, met->give_term (dd, type_indice), 2, name_ind, type_ind) ;
		return p_ope ;
	}

	// Check if Christoffel :
	name_ind = 0x0 ;
	type_ind = 0x0 ;
	indic = ischristo (name, name_ind, type_ind) ;
	if (indic) {
		p_ope = new Ope_id (this, met->give_christo (dd), 3, name_ind, type_ind) ;
		return p_ope ;
	}

	// Check if Riemann :
	name_ind = 0x0 ;
	type_ind = 0x0 ;
	indic = isriemann (name, name_ind, type_ind) ;
	if (indic) {
		p_ope = new Ope_id (this, met->give_riemann (dd), 4, name_ind, type_ind) ;
		return p_ope ;
	}

	// Check if Ricci tensor :
	name_ind = 0x0 ;
	type_ind = 0x0 ;
	indic = isricci_tensor (name, name_ind, type_ind) ;
	if (indic) {
		p_ope = new Ope_id (this, met->give_ricci_tensor (dd), 2, name_ind, type_ind) ;
		return p_ope ;
	}

	// Check if Ricci scalar :
	name_ind = 0x0 ;
	type_ind = 0x0 ;
	indic = isricci_scalar (name, name_ind, type_ind) ;
	if (indic) {
		p_ope = new Ope_id (this, met->give_ricci_scalar (dd)) ;
		return p_ope ;
	}

	// Check if Dirac :	
	name_ind = 0x0 ;
	type_ind = 0x0 ;
	indic = is_tensor(name, "dirac ", valence, name_ind, type_ind) ;
	if (indic) {
	    if (valence != 1) {
	      cerr << "Dirac operator must be of valence one" << endl ;
	      abort() ;
	    }
	    p_ope = new Ope_id (this, met->give_dirac(dd), valence, name_ind, type_ind) ;
	    return p_ope ;
	}
	
	// Check if normal vector :	
	name_ind = 0x0 ;
	type_ind = 0x0 ;
	indic = is_tensor(name, "normal ", valence, name_ind, type_ind) ;
	if (indic) {
	    if (valence != 1) {
	      cerr << "Normal vector operator must be of valence one" << endl ;
	      abort() ;
	    }
	    if (met==0x0) {
	      cerr << "Metric must be passed to call normal" << endl ;
	      abort() ;
	    }
	    int typebase = met->give_type(dd) ;
	    p_ope = new Ope_id (this, espace.get_domain(dd)->give_normal(bound, typebase), valence, name_ind, type_ind) ;
	    return p_ope ;
	}
	
	// Check if passed cst :
	name_ind = 0x0 ;
	type_ind = 0x0 ;
	indic = iscst (name, which, valence, name_ind, type_ind) ;
	if (indic) {
		if (name_ind==0x0)
			p_ope = new Ope_id (this, give_cst(which, dd)) ;
		else
			p_ope = new Ope_id (this, give_cst(which, dd), valence, name_ind, type_ind) ;
		return p_ope ;
		}

	// Check if cst passed in hard (i.e. a double)
	double val = 0. ;
	indic = isdouble (name, val) ;	
	if (indic) {
		valence = 0 ;
		p_ope = new Ope_id (this, give_cst_hard(val, dd)) ;
		return p_ope ;
		}

	// Check if definition :
	name_ind = 0x0 ;
	type_ind = 0x0 ;
	indic = isdef (dd, name, which, valence, name_ind, type_ind) ;
	if (indic) {
		p_ope = new Ope_id (this, give_def(which)->get_res(), valence, name_ind, type_ind) ;
		return p_ope ;
	}
	
	// Check if global definition :;
	indic = isdef_glob (dd, name, which) ;
	if (indic) {
		p_ope = new Ope_id (this, give_def_glob(which)->get_res()) ;
		return p_ope ;
	}

	// Check if user defined operator
	for (int i=0 ; i<nopeuser ; i++) {
	  indic = is_ope_uni (name, p1, names_opeuser[i]) ;
	  if (indic) {
	      p_ope = new Ope_user (this, opeuser[i], paruser[i], give_ope(dd,p1, bound)) ; 
	      return p_ope ;
	  }
	}

      // Check if user defined operator
	for (int i=0 ; i<nopeuser_bin ; i++) {
	  indic = is_ope_uni (name, p1, p2, names_opeuser_bin[i]) ;
	  if (indic) {
	      p_ope = new Ope_user_bin (this, opeuser_bin[i], paruser_bin[i], give_ope(dd,p1, bound), give_ope(dd, p2, bound)) ; 
	      return p_ope ;
	  }
	}
	
	assert (p_ope==0x0) ;
	cerr << "Unknown operator " << name << endl ;
	abort() ;
}
}
