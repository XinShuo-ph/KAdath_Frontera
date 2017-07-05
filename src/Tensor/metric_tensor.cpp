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
#include "metric_tensor.hpp"
namespace Kadath {
double sign(gsl_permutation* p)     // compute permutation signature
{
   int card(static_cast<int>(gsl_permutation_size(p)));
   double signature(1.0);
   for (int i(0 ); i < card - 1 ; ++i)
      for (int j(i+1) ; j < card ; ++j)
      {
         double pj(static_cast<double>(gsl_permutation_get(p, j)));        // convert to int just because I don't like compiler warnings...
         double pi(static_cast<double>(gsl_permutation_get(p, i)));
         signature *= (pj - pi) / (j - i);   // compute signature
      }
   return signature;
}

vector<int> ind_com(int i, int j, int k, int l)  // introduce matrix indices into the transposed comatrix components
{
   vector<int> res(2);
   if (k < j and l < i)
   {
      res[0] = k;
      res[1] = l;
   }
   if (k >= j and l < i)
   {
      res[0] = k + 1;
      res[1] = l;
   }
   if (k < j and l >= i)
   {
      res[0] = k;
      res[1] = l + 1;
   }
   if (k >= j and l >= i)
   {
      res[0] = k + 1;
      res[1] = l + 1;
   }
   return res;
}
// Storage functions :
int metric_tensor_position_array (const Array<int>& idx, int ndim) {

    assert (idx.get_ndim() == 1) ;
    assert (idx.get_size(0) == 2) ;
    
    for (int i=0 ; i<2 ; i++)
	assert ((idx(i)>=1) && (idx(i)<=ndim)) ;

    int res = 0 ;
    bool found = false ;
   
    if (idx(1)>idx(0)) {
	    for (int i=0 ; i<ndim ; i++)
		for (int j=i ; j<ndim ; j++) {
			if ((i+1==idx(0)) && (j+1==idx(1)))
				found = true ;
			if (!found)
				res ++ ;
    		}
    }
    else {
	    for (int i=0 ; i<ndim ; i++)
		for (int j=i ; j<ndim ; j++) {
			if ((i+1==idx(1)) && (j+1==idx(0)))
				found = true ;
			if (!found)
				res ++ ;
    		}
    }
    return res;
}

int metric_tensor_position_index (const Index& idx, int ndim) {
    assert (idx.get_ndim() == 2) ;
    for (int i=0 ; i<2 ; i++)
	assert ((idx(i)>=0) && (idx(i)<ndim)) ;

    int res = 0 ;
    bool found = false ;
   
    if (idx(1)>idx(0)) {
	    for (int i=0 ; i<ndim ; i++)
		for (int j=i ; j<ndim ; j++) {
			if ((i==idx(0)) && (j==idx(1)))
				found = true ;
			if (!found)
				res ++ ;
    		}
    }
    else {
	    for (int i=0 ; i<ndim ; i++)
		for (int j=i ; j<ndim ; j++) {
			if ((i==idx(1)) && (j==idx(0)))
				found = true ;
			if (!found)
				res ++ ;
    		}
    }
    return res ;
}

Array<int> metric_tensor_indices (int pos, int valence, int ndim) {

	assert (valence==2) ;
        Array<int> res (valence) ;
	int current = 0 ;
	for (int i=1 ; i<=ndim ; i++)
		for (int j=i ; j<=ndim ; j++) {
			if (current==pos) {
				res.set(0) = i ;
				res.set(1) = j ;
			}
			current ++ ;
	}
	return res ;
}

		//////////////
		// Members
		//////////////

Metric_tensor::Metric_tensor (const Space& sp, int type_descr, const Base_tensor& bb) : 
	Tensor (sp, 2, type_descr, 6, bb, 3) {

	give_place_array = metric_tensor_position_array ;
	give_place_index = metric_tensor_position_index ;
	give_indices = metric_tensor_indices ;
}

Metric_tensor::Metric_tensor (const Metric_tensor& so, bool copie) :
	Tensor (so, copie) {
}

Metric_tensor::Metric_tensor (const Space& sp, FILE* ff) :
	Tensor (sp, 3, ff) {

	assert (valence==2) ;
	assert (type_indice(0)==type_indice(1)) ;
	assert (n_comp==6) ;

	// Overwrite the storage functions :
	give_place_array = metric_tensor_position_array ;
	give_place_index = metric_tensor_position_index ;
	give_indices = metric_tensor_indices ;
}

Metric_tensor::~Metric_tensor () {
}

void Metric_tensor::operator= (const Metric_tensor& so) {
    assert (&espace==&so.espace) ;
    basis = so.basis ;
    assert(so.get_type()==get_type()) ;

    for (int i=0 ; i<n_comp ; i++)
      *cmp[i] = *so.cmp[i] ;
}


void Metric_tensor::operator= (const Tensor& so) {
	assert (so.valence==2) ; 
	assert (&espace==&so.espace) ;
	basis = so.basis ;
	assert (so.type_indice(0)==get_type()) ;
	assert (so.type_indice(1)==get_type()) ;
	
	for (int i=1 ; i<=ndim ; i++)
		for (int j=i ; j<=ndim ; j++)
			set(i,j) = (so(i,j) + so(j,i))/2. ; // On symetrise in case non sym
}

void Metric_tensor::operator= (double xx) {
	for (int i=0 ; i<n_comp ; i++)
		*cmp[i] = xx ;
}

Metric_tensor Metric_tensor::inverse(int ndom_max)
{
   
   assert(ndim == 3);
   Metric_tensor res(*this);
   res.set_index_type(0) = -this->get_index_type(0);
   res.set_index_type(1) = -this->get_index_type(1);
   for (int d(0) ; d <= ndom_max ; ++d)      // sum over components, only the upper triangular part
   {
      Val_domain detval(espace.get_domain(d));     // let's compute the determinant...
      detval = 0.0;
      gsl_permutation* p(gsl_permutation_calloc(ndim));  // initialise to identity 0, 1, 2
      do
      {
         int p1(static_cast<int>(gsl_permutation_get(p, 0)) + 1);   // for tensor indices, you need to add one to get a permutation of 1, 2, 3
         int p2(static_cast<int>(gsl_permutation_get(p, 1)) + 1);
         int p3(static_cast<int>(gsl_permutation_get(p, 2)) + 1);
         detval += sign(p)*res(p1, 1)(d)*res(p2, 2)(d)*res(p3, 3)(d);   // sum on permutations
      }while(gsl_permutation_next(p) == GSL_SUCCESS);
      gsl_permutation_free(p);

      // compute the transposed comatrix
      Val_domain cmpval(espace.get_domain(d));
      for (int i(1) ; i <= ndim ;  ++i)      // sum over components, only the upper triangular part
      {
         for (int j(i) ; j <= ndim ; ++j)
         {
            cmpval = 0.0;                   // initialise to 0 before each computation
            p = gsl_permutation_calloc(ndim - 1);              // reinitialise to identity 0, 1, ready to be used for next component
            do
            {
               int p1(static_cast<int>(gsl_permutation_get(p, 0)) + 1);  // always add one for tensor indices, now we have a permutation of 1, 2
               int p2(static_cast<int>(gsl_permutation_get(p, 1)) + 1);
               vector<int> index1, index2;                      // manage the indices of the transposed comatrix
               index1 = ind_com(i, j, p1, 1);                   // we need to play a bit with indices to get transpose(COM(metric))
               index2 = ind_com(i, j, p2, 2);
               cmpval += pow(-1.0, i + j)*sign(p)*(*this)(index1[0], index1[1])(d)*(*this)(index2[0],index2[1])(d);      // sum on permutations
            }while(gsl_permutation_next(p) == GSL_SUCCESS);
            gsl_permutation_free(p);
            res.set(i,j).set_domain(d) = cmpval / detval;              // now we have the component of the inverse matrix
         }
      }
   }
   return res;
}

void Metric_tensor::std_base2()
{
   for (int d(0) ; d < ndom ; ++d)
   {
      cmp[0]->std_base_domain(d);    // xx ou rr
      cmp[3]->std_base_domain(d);    // yy ou tt
      cmp[5]->std_base_domain(d);    // zz ou pp
      if (basis.get_basis(d) == CARTESIAN_BASIS)
      {
         cmp[1]->std_base_xy_cart_domain(d);
         cmp[2]->std_base_xz_cart_domain(d);
         cmp[4]->std_base_yz_cart_domain(d);
      }
      if (basis.get_basis(d) == SPHERICAL_BASIS)
      {
         cmp[1]->std_base_rt_spher_domain(d);
         cmp[2]->std_base_rp_spher_domain(d);
         cmp[4]->std_base_tp_spher_domain(d);
      }
   }
}}
