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

#include "scalar.hpp"
#include "tensor_impl.hpp"
#include "tensor.hpp"
namespace Kadath {
void affecte_one_dom (int dd, Tensor* res, const Tensor* so) {
    if ((so->is_name_affected()) && (!res->is_name_affected())) {
		res->set_name_affected() ;
		for (int i=0 ; i<res->get_valence() ; i++)
			res->set_name_ind(i, so->get_name_ind()[i]) ;
	      }
	      
    if (!so->is_name_affected()) {
      for (int i=0 ; i<res->get_n_comp() ; i++) {
	Array<int> ind (res->indices(i)) ;
	res->set(ind).set_domain(dd) = (*so)(ind)(dd) ;
      }
    }
    else {
 	Array<int> perm (so->valence) ;
	bool same_ind = res->find_indices(*so, perm) ;
	if (!same_ind) {
		cerr << "Indices do not match in affecte_one_dom" << endl ;
		abort() ;
		}
	Array<int> ind_so (so->valence) ;
	
      for (int i=0 ; i<res->get_n_comp() ; i++) {
	  Array<int> ind (res->indices(i)) ;
		for (int j=0 ; j<so->valence ; j++)
			ind_so.set(perm(j)) = ind(j) ;
	res->set(ind).set_domain(dd) = (*so)(ind_so)(dd) ;
      }
    }  
}

Tensor add_one_dom (int dd, const Tensor & t1, const Tensor & t2) {
    assert (&t1.espace==&t2.espace) ;

    assert (t1.valence == t2.valence) ;
    if (t1.valence != 0) {
	assert (t1.basis.get_basis(dd) == t2.basis.get_basis(dd)) ;
    }    	
  // Check same type :
   bool same = ((t1.give_place_array==t2.give_place_array) && (t1.give_place_index==t2.give_place_index) &&  (t1.give_indices==t2.give_indices)) ? 
			true : false ;

     Tensor aux (t1.espace, t1.valence, t1.type_indice, t1.basis) ;
     int ncmp = (same) ? t1.n_comp : aux.n_comp ;
     Tensor res (t1.espace, t1.valence, t1.type_indice, ncmp, t1.basis) ;
 
     if (same) {
	res.give_place_array = t1.give_place_array ;
	res.give_place_index = t1.give_place_index ;
	res.give_indices = t1.give_indices ;
	}
     else {
	res.give_place_array = aux.give_place_array ;
	res.give_place_index = aux.give_place_index ;
	res.give_indices = aux.give_indices ;
	}


  // Case when one does not care about the indices :
    if ((!t1.name_affected) || (!t2.name_affected)) {  

	    for (int i=0 ; i<t1.valence ; i++)
		assert(t1.type_indice(i) == t2.type_indice(i)) ;

    	for (int i=0 ; i<res.n_comp ; i++) {
      	Array<int> ind (res.indices(i)) ;
      	res.set(ind).set_domain(dd) = t1(ind)(dd) + t2(ind)(dd) ;
    	}
    }
    else {

	// Res with the same indices as t1 
	res.name_affected= true ;
	for (int i=0 ; i<t1.valence ; i++)
		res.name_indice[i] = t1.name_indice[i] ;

	// Find the permutation :
	Array<int> perm (t1.valence) ;
	bool same_ind = t1.find_indices(t2, perm) ;
	if (!same_ind) {
		cerr << "Indices do not match in add_one_dom" << endl ;
		abort() ;
		}
	else 
		for (int i=0 ; i<res.n_comp ; i++) {
			Array<int> ind (res.indices(i)) ;
			Array<int> ind_t2 (t1.valence) ;
			for (int j=0 ; j<t1.valence ; j++)
				ind_t2.set(perm(j)) = ind(j) ;
			res.set(ind).set_domain(dd) = t1(ind)(dd) + t2(ind_t2)(dd) ;
		}
	}

    // Put parameters :
    int m_res = add_m_quant (t1.get_parameters(), t2.get_parameters()) ;
    if (m_res!=0) {
      res.set_parameters().set_m_quant() = m_res ;
    }

    return res ;
}

Tensor add_one_dom (int dd, const Tensor& t, double xx) {
	Tensor res(t.espace, t.valence, t.type_indice, t.n_comp, t.basis) ;
	res.give_place_array = t.give_place_array ;
	res.give_place_index = t.give_place_index ;
	res.give_indices = t.give_indices ;
	for (int i=0 ; i<res.n_comp ; i++)
		res.cmp[i]->set_domain(dd) = t(res.indices(i))(dd)+xx ;

	// Copy parameters in this case...
	if (t.parameters)
	  res.parameters = t.get_parameters() ;
	return res ;
}

Tensor add_one_dom (int dd, double xx, const Tensor& t) {
	Tensor res(t.espace, t.valence, t.type_indice, t.n_comp, t.basis) ;
	res.give_place_array = t.give_place_array ;
	res.give_place_index = t.give_place_index ;
	res.give_indices = t.give_indices ;
	for (int i=0 ; i<res.n_comp ; i++)
		res.cmp[i]->set_domain(dd) = t(res.indices(i))(dd)+xx ;

	// Copy parameters in this case...
	if (t.parameters)
	  res.parameters = t.get_parameters() ;
	return res ;
}

Tensor sub_one_dom (int dd, const Tensor & t1, const Tensor & t2) {
    assert (&t1.espace==&t2.espace) ;

    assert (t1.valence == t2.valence) ;
    if (t1.valence != 0) {
	assert (t1.basis.get_basis(dd) == t2.basis.get_basis(dd)) ;
    }

   // Check same type :
   bool same = ((t1.give_place_array==t2.give_place_array) && (t1.give_place_index==t2.give_place_index) &&  (t1.give_indices==t2.give_indices)) ? 
			true : false ;
   Tensor aux (t1.espace, t1.valence, t1.type_indice, t1.basis) ;
     int ncmp = (same) ? t1.n_comp : aux.n_comp ;
     Tensor res (t1.espace, t1.valence, t1.type_indice, ncmp, t1.basis) ;
 
     if (same) {
	res.give_place_array = t1.give_place_array ;
	res.give_place_index = t1.give_place_index ;
	res.give_indices = t1.give_indices ;
	}
     else {
	res.give_place_array = aux.give_place_array ;
	res.give_place_index = aux.give_place_index ;
	res.give_indices = aux.give_indices ;
	}

    // Case when one does not care about the indices :
    if ((!t1.name_affected) || (!t2.name_affected)) {  

	    for (int i=0 ; i<t1.valence ; i++)
		assert(t1.type_indice(i) == t2.type_indice(i)) ;

    	for (int i=0 ; i<res.n_comp ; i++) {
      	Array<int> ind (res.indices(i)) ;
      	res.set(ind).set_domain(dd) = t1(ind)(dd) - t2(ind)(dd) ;
    	}
    }
    else {

	// Res with the same indices as t1 
	res.name_affected= true ;
	for (int i=0 ; i<t1.valence ; i++)
		res.name_indice[i] = t1.name_indice[i] ;

	// Find the permutation :
	Array<int> perm (t1.valence) ;
	bool same_ind = t1.find_indices(t2, perm) ;
	if (!same_ind) {
		cerr << "Indices do not match in sub_one_dom" << endl ;
		abort() ;
		}
	else 
		for (int i=0 ; i<res.n_comp ; i++) {
			Array<int> ind (res.indices(i)) ;
			Array<int> ind_t2 (t1.valence) ;
			for (int j=0 ; j<t1.valence ; j++)
				ind_t2.set(perm(j)) = ind(j) ;
			res.set(ind).set_domain(dd) = t1(ind)(dd) - t2(ind_t2)(dd) ;
		}
	}

  // Put parameters :
    int m_res = add_m_quant (t1.get_parameters(), t2.get_parameters()) ;
      if (m_res!=0) {
      res.set_parameters().set_m_quant() = m_res ;
    }

    return res ;
}

Tensor sub_one_dom (int dd, const Tensor& t, double xx) {
	Tensor res(t.espace, t.valence, t.type_indice, t.n_comp, t.basis) ;
	res.give_place_array = t.give_place_array ;
	res.give_place_index = t.give_place_index ;
	res.give_indices = t.give_indices ;
	for (int i=0 ; i<res.n_comp ; i++)
		res.cmp[i]->set_domain(dd) = t(res.indices(i))(dd)-xx ;
	// Copy parameters in this case...
	if (t.parameters)
	    res.parameters = t.get_parameters() ;
	return res ;
}

Tensor sub_one_dom (int dd, double xx, const Tensor& t) {
	Tensor res(t.espace, t.valence, t.type_indice, t.n_comp, t.basis) ;
	res.give_place_array = t.give_place_array ;
	res.give_place_index = t.give_place_index ;
	res.give_indices = t.give_indices ;
	for (int i=0 ; i<res.n_comp ; i++)
		res.cmp[i]->set_domain(dd) = -t(res.indices(i))(dd)+xx ;
	// Copy parameters in this case...
	if (t.parameters) res.parameters = t.get_parameters() ;
	return res ;
}

Tensor mult_one_dom (int dd, const Tensor& t1, const Tensor& t2) {
	assert (&t1.espace==&t2.espace) ;

	if ((t1.valence!=0) && (t2.valence!=0))
	    assert (t1.basis.get_basis(dd)==t2.basis.get_basis(dd)) ;

	bool do_name = true ;
	if ((!t1.name_affected) && (t1.valence!=0))
		do_name = false ;
	if ((!t2.name_affected) && (t2.valence!=0))
		do_name = false ;
	
	if (!do_name) {
	// Standard product
	int val_res = t1.valence+t2.valence ;
	Array<int> ind_res (val_res) ;
	for (int i=0 ; i<t1.valence ; i++)
		ind_res.set(i) = t1.type_indice(i) ;
	for (int i=0 ; i<t2.valence ; i++)
		ind_res.set(i+t1.valence) = t2.type_indice(i) ;

	const Base_tensor& base = (t1.valence!=0) ? t1.basis : t2.basis ;
	Tensor res (t1.espace, val_res, ind_res, base) ;

	Array<int> ind1 (t1.valence) ;
	Array<int> ind2 (t2.valence) ;
	for (int i=0 ; i<res.n_comp ; i++) {
		ind_res = res.indices(i) ;
		for (int j=0 ; j<t1.valence ; j++)
			ind1.set(j) = ind_res(j) ;
		for (int j=0 ; j<t2.valence ; j++)
			ind2.set(j) = ind_res(j+t1.valence) ;
		res.set(ind_res).set_domain(dd) = (t1(ind1))(dd) * (t2(ind2))(dd) ;
	}
	  
	// Put parameters :
      int m_res = mult_m_quant (t1.get_parameters(), t2.get_parameters()) ;
        if (m_res!=0) {
	res.set_parameters().set_m_quant() = m_res ;
      }

	return res ;
	}
	else {

		// Check if one needs to sum on some indices :
		Array<int> sum_1 (t1.valence) ;
		Array<int> sum_2 (t2.valence) ;
		for (int i=0 ; i<t2.valence ; i++)
			sum_2.set(i) = -1 ;

		for (int i=0 ; i<t1.valence ; i++) {
			sum_1.set(i) = -1 ;
			for (int j=0 ; j<t2.valence ; j++)
				if (t1.name_indice[i]==t2.name_indice[j]){
					sum_1.set(i) = j ;
					sum_2.set(j) = i ;
				}
			if ((sum_1(i)!=-1) && (t1.type_indice(i)==t2.type_indice(sum_1(i)))) {
				cerr << "Can not sum on indices of the same type in operator*" << endl ;
				abort() ;
			}
		}

		// Valence of the result
		int val_res = 0 ;
		for (int i=0 ; i<t1.valence ; i++)
			if (sum_1(i)==-1)
				val_res ++ ;
		for (int i=0 ; i<t2.valence ; i++)
			if (sum_2(i)==-1)
				val_res ++ ;

		if (val_res>0) {
		  
		// Type on indices
		Array<int> ind_res (val_res) ;
		int conte = 0 ;
		for (int i=0 ; i<t1.valence ; i++)
			if (sum_1(i)==-1) {
				ind_res.set(conte)=t1.type_indice(i) ;
				conte ++ ;
			}
		for (int i=0 ; i<t2.valence ; i++)
			if (sum_2(i)==-1) {
				ind_res.set(conte) = t2.type_indice(i) ;
				conte ++ ;
			}

		// Le tenseur
		const Base_tensor& base = (t1.valence!=0) ? t1.basis : t2.basis ;
		Tensor res (t1.espace, val_res, ind_res, base) ;

		// Name of the remaining variables :
		res.name_affected = true ;
		conte = 0 ;
		for (int i=0 ; i<t1.valence ; i++)
			if (sum_1(i)==-1) {
				res.name_indice[conte] = t1.name_indice[i] ;
				conte ++ ;
			}
		for (int i=0 ; i<t2.valence ; i++)
			if (sum_2(i)==-1) {
				res.name_indice[conte] = t2.name_indice[i] ;
				conte ++ ;
			}
		
		Array<bool> first (res.get_n_comp()) ;
		first = true ;

		Index pos_t1 (t1) ;
		Index pos_t2 (t2) ;
		Index pos_res (res) ;

		do {
			pos_t2.set_start() ;
			do {
				// Check if the summation indices are the same 
				bool same = true ;
				for (int i=0 ; i<t1.valence ; i++)
					if (sum_1(i)!=-1)
						if (pos_t1(i)!=pos_t2(sum_1(i)))
							same = false ;
				// if the same :
				if (same) {
					conte = 0 ;
					for (int i=0 ; i<t1.valence ; i++)
						if (sum_1(i)==-1) {
							pos_res.set(conte) = pos_t1(i) ;
							conte ++ ;
						}
					for (int i=0 ; i<t2.valence ; i++)
						if (sum_2(i)==-1) {
							pos_res.set(conte) = pos_t2(i) ;
							conte ++ ;
						}
					int ind (res.position(pos_res)) ;
					if (first(ind)) {
						res.set(pos_res).set_domain(dd) = t1(pos_t1)(dd)*t2(pos_t2)(dd) ;
						first.set(ind) = false ;
					}
					else
						res.set(pos_res).set_domain(dd) += t1(pos_t1)(dd)*t2(pos_t2)(dd) ;
				}
			}
			while (pos_t2.inc()) ;
		}
		while (pos_t1.inc()) ;

		// Put parameters :
      int m_res = mult_m_quant (t1.get_parameters(), t2.get_parameters()) ;
        if (m_res!=0) {
	res.set_parameters().set_m_quant() = m_res ;
      }
		return res ;
		}
	else {
			// Result is a scalar :
			Scalar res (t1.espace) ;
			Index pos_t1 (t1) ;
			Index pos_t2 (t2) ;
			int first = true ;

			do {
				pos_t2.set_start() ;
				do {
					// Check if the summation indices are the same 
					bool same = true ;
					for (int i=0 ; i<t1.valence ; i++)
						if (sum_1(i)!=-1)
							if (pos_t1(i)!=pos_t2(sum_1(i)))
								same = false ;
				// if the same :
				if (same) {
					if (first) {
						res.set_domain(dd) = t1(pos_t1)(dd)*t2(pos_t2)(dd) ;
						first = false ;
					}
					else
						res.set_domain(dd) += t1(pos_t1)(dd)*t2(pos_t2)(dd) ;
				}
			}
			while (pos_t2.inc()) ;
		}
		while (pos_t1.inc()) ;
      // Put parameters :
      int m_res = mult_m_quant (t1.get_parameters(), t2.get_parameters()) ;
        if (m_res!=0) {
	res.set_parameters().set_m_quant() = m_res ;
      }
		return res ;
		}
	}
}

Tensor mult_one_dom (int dd, const Tensor& t, double xx) {
	Tensor res(t.espace, t.valence, t.type_indice, t.basis) ;
	for (int i=0 ; i<res.n_comp ; i++)
		res.cmp[i]->set_domain(dd) = t(res.indices(i))(dd)*xx ;

	  // Les indices if needed :
	   if (t.name_affected) {
		res.name_affected= true ;
		for (int i=0 ; i<res.valence ; i++)
			res.set_name_ind(i, t.name_indice[i]) ;
   	}
// Copy parameters in this case...
	if (t.parameters)
	  res.parameters = t.get_parameters() ;
	return res ;
}

Tensor mult_one_dom (int dd, double xx, const Tensor& t) {
	Tensor res(t.espace, t.valence, t.type_indice, t.basis) ;
	for (int i=0 ; i<res.n_comp ; i++)
		res.cmp[i]->set_domain(dd) =  t(res.indices(i))(dd)*xx ;  

	// Les indices if needed :
	   if (t.name_affected) {
		res.name_affected= true ;
		for (int i=0 ; i<res.valence ; i++)
			res.set_name_ind(i, t.name_indice[i]) ;
   	}
      // Copy parameters in this case...
	if (t.parameters)
	  res.parameters = t.get_parameters() ;
	return res ;
}

Tensor mult_one_dom (int dd, const Tensor& t, int mm) {
	Tensor res(t.espace, t.valence, t.type_indice, t.basis) ;
	for (int i=0 ; i<res.n_comp ; i++)
		res.cmp[i]->set_domain(dd) = t(res.indices(i))(dd)*mm ;

	  // Les indices if needed :
	   if (t.name_affected) {
		res.name_affected= true ;
		for (int i=0 ; i<res.valence ; i++)
			res.set_name_ind(i, t.name_indice[i]) ;
   	}
// Copy parameters in this case...
	if (t.parameters)
	  res.parameters = t.get_parameters() ;
	return res ;
}

Tensor mult_one_dom (int dd, int mm, const Tensor& t) {
	Tensor res(t.espace, t.valence, t.type_indice, t.basis) ;
	for (int i=0 ; i<res.n_comp ; i++)
		res.cmp[i]->set_domain(dd) = t(res.indices(i))(dd)*mm ;  

	// Les indices if needed :
	   if (t.name_affected) {
		res.name_affected= true ;
		for (int i=0 ; i<res.valence ; i++)
			res.set_name_ind(i, t.name_indice[i]) ;
   	}
	// Copy parameters in this case...
	if (t.parameters)
	  res.parameters = t.get_parameters() ;
	return res ;
}

Tensor div_one_dom (int dd, const Tensor& t1, const Tensor& t2) {
   if (t2.valence !=0) {
	cerr << "Division only defined for a scalar" << endl ;
	abort() ;
	}
   Tensor res(t1.espace, t1.valence, t1.type_indice, t1.basis) ;

   for (int i=0 ; i<res.n_comp ; i++)
	   res.cmp[i]->set_domain(dd) = t1(res.indices(i))(dd)/(*t2.cmp[0])(dd) ;

   // Les indices if needed :
   if (t1.name_affected) {
	res.name_affected= true ;
	for (int i=0 ; i<res.valence ; i++)
		res.set_name_ind(i, t1.name_indice[i]) ;
   }
 
   // Put parameters :
      int m_res = div_m_quant (t1.get_parameters(), t2.get_parameters()) ;
        if (m_res!=0) {
	res.set_parameters().set_m_quant() = m_res ;
      }
   return res ;
}


Tensor div_one_dom (int dd, double x, const Tensor& t) {
   if (t.valence !=0) {
	cerr << "Division only defined for a scalar" << endl ;
	abort() ;
	}
   Tensor res(t.espace, t.valence, t.type_indice, t.basis) ;

   res.cmp[0]->set_domain(dd) = x/(*t.cmp[0])(dd) ;
   
   // Copy parameters in this case...
      int m_res = inv_m_quant (t.get_parameters()) ;
       if (t.is_m_quant_affected()) {
	res.set_parameters().set_m_quant() = m_res ;
      }
   return res ;
}

Tensor div_one_dom (int dd, const Tensor& t, double xx) {
	Tensor res(t.espace, t.valence, t.type_indice, t.basis) ;
	for (int i=0 ; i<res.n_comp ; i++)
		res.cmp[i]->set_domain(dd) = t(res.indices(i))(dd)/xx ;  

	// Les indices if needed :
       if (t.name_affected) {
		res.name_affected= true ;
		for (int i=0 ; i<res.valence ; i++)
			res.set_name_ind(i, t.name_indice[i]) ;
   }

      // Copy parameters in this case...
	if (t.parameters)
	  res.parameters = t.get_parameters() ;
	return res ;
}

Tensor scal_one_dom (int dd, const Tensor & t1, const Tensor & t2) {
    assert (&t1.espace==&t2.espace) ;
    assert (t1.valence == t2.valence) ;
    assert (t1.valence!=0) ;
    assert (t1.basis.get_basis(dd) == t2.basis.get_basis(dd)) ;

    Array<int> ind (t1.indices(0)) ;
    Tensor auxi (mult_one_dom(dd, t1(ind), t2(ind))) ;
    for (int i=1 ; i<t1.n_comp ; i++) {
      ind = t1.indices(i) ;
      auxi = add_one_dom(dd, auxi, mult_one_dom(dd, t1(ind), t2(ind))) ;
    }

    return auxi ;
}


Tensor Tensor::do_summation_one_dom(int dd) const {

	if (!name_affected) {
		cerr << "Names of indices must be affected in Tensor::do_summation" << endl ;
		abort() ;
	}

	Array<int> sum_ind (valence) ;
	Array<int> nsame (valence) ;
	sum_ind = -1 ;
	nsame = 0 ;
	for (int i=0 ; i<valence ; i++)
		for (int j=i+1 ; j<valence ; j++) {
			if (name_indice[i]==name_indice[j]) {
					nsame.set(i) ++ ;
					nsame.set(j) ++ ;
					if ((nsame(i)>1) || (nsame(j)>1)) {
						cerr << "Too many identical indices in Tensor::do_summation_one_dom" << endl ;
						abort() ;
					}
					sum_ind.set(i) = j ;
					sum_ind.set(j) = i ;
					if (type_indice(i)==type_indice(j)) {
						cerr << "Can not sum on indices of the same type in Tensor::do_summation_one_dom" << endl ;
						abort() ;
					}
				}
		}

	int valence_res = 0 ;
	for (int i=0 ; i<valence ; i++)
		if (nsame(i)==0)
			valence_res ++ ;
	
	// difference Scalar or not :
	if (valence_res > 0) {
		Array<int> type_ind_res (valence_res) ;
		char* name_ind_res = new char[valence_res] ;
		int conte = 0 ;
		for (int i=0 ; i<valence ; i++)
			if (nsame(i)==0) {
				type_ind_res.set(conte) = type_indice(i) ;
				name_ind_res[conte] = name_indice[i] ;
				conte ++ ;
		}

		Tensor res (espace, valence_res, type_ind_res, basis) ;
		res.name_affected = true ;
		for (int i=0 ; i<valence_res ; i++)
			res.name_indice[i] = name_ind_res[i] ;
		delete [] name_ind_res ;


		Index pos_res (res) ;
		do {
			res.set(pos_res).set_domain(dd) = 0 ;
		}
		while (pos_res.inc()) ;

		Array<bool> first (res.get_n_comp()) ;
		first = true ;

		// Loop on the various components :
		Index pos (*this) ;
		do {
			bool take = true ;
			for (int i=0 ; i<valence ; i++)
				if ((nsame(i)!=0) && (pos(sum_ind(i))!=pos(i)))
					take = false ;

			if (take) {
				conte = 0 ;
				for (int i=0 ; i<valence ; i++) 
					if (nsame(i)==0) {
						pos_res.set(conte)=pos(i) ;
						conte ++ ;
					}
				int ind (res.position(pos_res)) ;
				if (first(ind)) {
					res.set(pos_res).set_domain(dd).set_base() = (*this)(pos)(dd).get_base() ;
					first.set(ind) = false ;
					}
				res.set(pos_res).set_domain(dd) += (*this)(pos)(dd) ;
			}
		}
		while (pos.inc()) ;
		return res ;
	}
	else {
		// Scalar case :
		Scalar res(espace) ;
		res.set_domain(dd) = 0 ;
		Index pos (*this) ;
		bool first = true ;
		do {
			bool take = true ;
			for (int i=0 ; i<valence ; i++)
				if (pos(sum_ind(i))!=pos(i))
					take = false ;

			if (take) {
				if (first) {
					res.set_domain(dd).set_base() = (*this)(pos)(dd).get_base() ;
					first = false ;
				}
				res.set_domain(dd) += (*this)(pos)(dd) ;
			}
		}
		while (pos.inc()) ;
		return res ;
	}
}

Tensor partial_one_dom (int dd, char ind_partial, const Tensor& so) {
	//Check triad :
	if ((so.valence>0) && (so.basis.get_basis(dd) != CARTESIAN_BASIS)) {
		cerr << "Partial_one_dom only defined with respect to cartesian basis so far" << endl  ;
		abort() ;
	}	

	bool dosum = false ;
	if (so.name_affected)
		for (int i=0 ; i<so.valence ; i++)
			if (ind_partial==so.name_indice[i])
				dosum = true ;

	Array<int> type_res (so.valence+1) ;
	type_res.set(0) = COV ;
	for (int i=0 ; i<so.valence ; i++)
		type_res.set(i+1) = so.type_indice(i) ;
	
	Base_tensor basis (so.get_space(), CARTESIAN_BASIS) ;
	
	Tensor res(so.espace, so.valence+1, type_res, basis) ;

	Index pos (res) ;
	Index pos_so(so) ;
	do {
		for (int i=0 ; i<so.valence  ; i++)
			pos_so.set(i) = pos(i+1) ;
		res.set(pos).set_domain(dd) = so(pos_so)(dd).der_abs(pos(0)+1) ;
	}
	while (pos.inc()) ;

	if ((so.name_affected) || (so.valence==0)) {
		res.name_affected = true ;
		res.name_indice[0] = ind_partial ;
		for (int i=0 ; i<so.valence ; i++)
			res.name_indice[i+1] = so.name_indice[i] ;
	}

	if (dosum)
		return res.do_summation_one_dom(dd) ;
	else
		return res ;
}

Tensor sqrt_one_dom (int dd, const Tensor& t) {
	if (t.get_valence() !=0) {
		cerr << "sqrt_one_dom only defined for scalars" << endl ;
		abort() ;
	}
	Scalar res (t.get_space()) ;
	res.set_domain(dd) = sqrt(t()(dd)) ;
	return res ;
}

}
