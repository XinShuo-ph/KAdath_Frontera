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
#include "base_spectral.hpp"
#include "array.hpp"

namespace Kadath{
Base_spectral::Base_spectral (int dimensions) :  def(false), ndim(dimensions) {
	bases_1d = new Array<int>* [ndim] ;
	for (int i=0 ; i<ndim ; i++)
	    bases_1d[i] = 0x0 ;
}

Base_spectral::Base_spectral(const Base_spectral& so) : def(so.def), ndim(so.ndim) {
	bases_1d = new Array<int>* [ndim] ;
	for (int l=0 ; l<ndim ; l++)
		bases_1d[l] = (so.bases_1d[l] == 0x0) ? 0x0 : new Array<int>(*so.bases_1d[l]) ;
}

Base_spectral::Base_spectral (FILE* fd) {
	int indic ;
	fread_be (&indic, sizeof(int), 1, fd) ;
	def = (indic==0) ? true : false ;
	fread_be (&ndim, sizeof(int), 1, fd) ;	
	bases_1d = new Array<int>* [ndim] ;
	if (def)
		for (int i=0 ; i<ndim ; i++)
	    		bases_1d[i] = new Array<int>(fd) ;
	else
		for (int i=0 ; i<ndim ; i++)
	    	bases_1d[i] = 0x0 ;
}

#ifdef ENABLE_MOVE_SEMANTIC
    Base_spectral::Base_spectral(Base_spectral && so) : def{so.def}, ndim{0}, bases_1d{nullptr}
    {
        std::swap(ndim,so.ndim);
        std::swap(bases_1d,so.bases_1d);
        so.def = false;
    }
    Base_spectral& Base_spectral::operator=(Base_spectral &&so)
    {
        std::swap(def,so.def);
        ndim = so.ndim;
        std::swap(bases_1d,so.bases_1d);
        return *this;
    }
#endif

Base_spectral::~Base_spectral() {
    if(bases_1d)
    {
        for (int l = 0; l < ndim; l++)
            if (bases_1d[l] != 0x0) delete bases_1d[l];
        delete[] bases_1d;
    }
}

void Base_spectral::set_non_def() {
	for (int l=0 ; l<ndim ; l++)
	     if (bases_1d[l] !=0x0) {
		delete bases_1d[l] ;
		bases_1d[l] = 0x0 ;
	}
	def = false ;
}

void Base_spectral::save (FILE* fd) const {
	int indic = (def) ? 0 : 1 ;
	fwrite_be (&indic, sizeof(int), 1, fd) ;
	fwrite_be (&ndim, sizeof(int), 1, fd) ;
	if (def)
		for (int i=0 ; i<ndim ; i++)
			bases_1d[i]->save(fd) ;
}

void Base_spectral::operator= (const Base_spectral& so) {
	assert (ndim==so.ndim) ;
	def = so.def ;
	for (int i=0 ; i<ndim ; i++) {
	    if (bases_1d[i] !=0x0)
	    	delete bases_1d[i] ;
	    if ((so.def) && (so.bases_1d[i] !=0x0))
		bases_1d[i] = new Array<int>(*so.bases_1d[i]) ;
	    else 
	    	bases_1d[i] = 0x0 ;
	    }
}

void Base_spectral::allocate(const Dim_array& nbr_coef) {
	
	for (int i=0 ; i<ndim ; i++)
	    if (bases_1d[i] != 0x0) delete bases_1d[i] ;
	    
	bases_1d[ndim-1] = new Array<int>(1) ;
	for (int i=ndim-2 ; i>=0 ; i--) {
		Dim_array taille (ndim-1-i) ;
		for (int k=0 ; k<ndim-1-i ; k++)
		    taille.set(k) = nbr_coef(i+k+1) ;
		bases_1d[i] = new Array<int>(taille) ;
	}
}	

bool operator== (const Base_spectral& a, const Base_spectral& b) {
	bool res {(a.def) && (b.def) && (a.ndim == b.ndim)} ;
	for (int i=0 ; i<a.ndim ; i++) {
		Array_iterator index (a.bases_1d[i]->get_dimensions()) ;
		do
			res = ((*a.bases_1d[i])(index) == (*b.bases_1d[i])(index)) ;
		while (index.inc() && res) ;
	}
	return res ;
} 

ostream& operator<< (ostream& o, const Base_spectral& so) {


	o << so.ndim << "-dimensional spectral base" << endl ;
	if (so.def) {
	for (int l=0 ; l<so.ndim ; l++) {
	    o << "Variable " << l << endl ;
	    o << *so.bases_1d[l] << endl ;
	}
	}
	else
		o << "Base not defined" << endl ;
return o ;
}

void Base_spectral::set(Dim_array const& nbr_coefs, int BASEPHI, int BASETHETA, int BASER)
{
 
   assert (nbr_coefs.get_ndim()==3) ;

   allocate(nbr_coefs);
   def = true;
   bases_1d[2]->set(0) = BASEPHI;
//   Array_iterator index(bases_1d[0]->get_dimensions());
    Index index(bases_1d[0]->get_dimensions());
   for (int k(0) ; k < nbr_coefs(2) ; ++k) 
   {
      bases_1d[1]->set(k) = BASETHETA;
      for (int j(0) ; j < nbr_coefs(1) ; ++j) 
      {
          index.set(0) = j ; index.set(1) = k;
//          index.inc();
          bases_1d[0]->set(index) = BASER;
       }
   }
}
}
