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
#include "matrice.hpp"
#include "array.hpp"
#include "array_math.hpp"
namespace Kadath {
//Destructeur des quantites derivees

void Matrice::del_deriv() {
    if (band != 0x0) delete band ;
    if (lu != 0x0) delete lu ;
    if (permute != 0x0) delete permute ;
    band = 0x0 ;
    lu = 0x0 ;
    permute = 0x0 ;
}

// sets to zero in std
void Matrice::annule() {
    del_deriv() ;
    Index index(sizes) ;
    do std->set(index) = 0 ;
	while (index.inc()) ;
}
	 
// Standard Constructor
Matrice::Matrice (int i, int j) : sizes(2) {
    sizes.set(0) = i ; sizes.set(1) = j ;
    std = new Array<double>(sizes) ;
    kl = 0 ;
    ku = 0 ;
    band = 0x0 ;
    lu = 0x0 ;
    permute = 0x0 ;
}


// Copy constructor
Matrice::Matrice (const Matrice & source) : sizes(source.sizes) {
    kl = source.kl ;
    ku = source.ku ;
    std = new Array<double>(*source.std) ;
    if (source.band != 0x0) band = new Array<double>(*source.band) ;
    else band = 0x0 ;
    if (source.lu != 0x0) lu = new Array<double>(*source.lu) ;
    else lu = 0x0 ;
    if (source.permute != 0x0) permute = new Array<int>(*source.permute) ;
    else permute = 0x0 ;
}


// Copy from a 2 dimensional array
Matrice::Matrice (const Array<double> & source) : sizes(source.get_dimensions()) {
    std = new Array<double> (source) ;
    assert (sizes.get_ndim()==2) ;
    kl = 0 ;
    ku = 0 ;
    band = 0x0 ;
    lu = 0x0 ;
    permute = 0x0 ;
}

// destructor
Matrice::~Matrice() {
    del_deriv() ;
    delete std ;
}

// Assignement to a double in the std representation
void Matrice::operator= (double x) {
    del_deriv() ;
    *std = x ;
}

// Assignement to another matrix
void Matrice::operator= (const Matrice &source) {
    
    assert (sizes==source.sizes) ;
    del_deriv() ;
    delete std ;
    std = new Array<double>(*source.std) ;
    kl = source.kl ;
    ku = source.ku ;
    if (source.band !=0x0)
        band = new Array<double>(*source.band) ;
    if (source.lu !=0x0)
        lu = new Array<double>(*source.lu) ;
    if (source.permute !=0x0)
        permute = new Array<int>(*source.permute) ;
}

// Assignement to a 2dimensional array
void Matrice::operator= (const Array<double> &source) {
    
    assert (sizes == source.get_dimensions()) ;
    del_deriv() ;
    delete std ;
    std = new Array<double>(source) ;
    kl = 0 ;
    ku = 0 ;
}

// Read/write of an element
double& Matrice::set (int i, int j) {
	del_deriv() ;
	Index index (sizes) ;
	index.set(0) = i ;
	index.set(1) = j ;
	return std->set(index) ;
}

void Matrice::copy_inside (int i, int j, const Matrice& so) {
	del_deriv() ;
	Index index_so (so.sizes) ;
	Index index (sizes) ;
	do {
		index.set(0) = index_so(0) + i ;
		index.set(1) = index_so(1) + j ;
		std->set(index) = (*so.std)(index_so) ;
	}
	while (index_so.inc()) ;
}


// Read only of an element
double Matrice::operator() (int i, int j) const {
	Index index (sizes) ;
	index.set(0) = i ;
	index.set(1) = j ;
	return (*std)(index) ;
}

//Display
ostream& operator<< (ostream& flux, const Matrice & source) {
    
	flux << "Matrix " << source.sizes(0) << " * " << source.sizes(1) << endl ;
	    for (int i=0 ; i<source.sizes(0) ; i++) {
		for (int j=0 ; j<source.sizes(1) ; j++)
		    flux << source(i, j) << "  " ;
		flux << endl ;	
	    }
    flux << endl ;
    return flux ;
}

// Computes the banded representation : LAPACK storage
void Matrice::set_band (int u, int l) const {
    if (band != 0x0) return ;
    else {
	int n = sizes(0) ;
	assert (n == sizes(1)) ;
	
	ku = u ; kl = l ;
	int ldab = 2*l+u+1 ;
	band = new Array<double>(ldab*n) ;
	*band = 0 ; 
	
	for (int i=0 ; i<u ; i++)
	    for (int j=u-i ; j<n ; j++)
		band->set(j*ldab+i+l) = (*this)(j-u+i, j) ;
	
	for (int j=0 ; j<n ; j++)
	    band->set(j*ldab+u+l) = (*this)(j, j) ;
	
	for (int i=u+1 ; i<u+l+1 ; i++)
	    for (int j=0 ; j<n-i+u ; j++)
		band->set(j*ldab+i+l) = (*this) (i+j-u, j) ;

    }
    return ; 
}

//LU decomposition : LAPACK storage
void Matrice::set_lu() const {
    if (lu != 0x0) {
        assert (permute != 0x0) ;
        return ;
    }
    else {
        // LU decomposition
        int n = sizes(0) ;
        int ldab, info ;
        permute = new Array<int>(n) ;

        // Case of a banded matrix
        if (band != 0x0) {
            ldab = 2*kl+ku+1 ;
            lu = new Array<double>(*band) ;

            F77_dgbtrf(&n, &n, &kl, &ku, lu->set_data(), &ldab, permute->set_data(), &info) ;
        }
        else { // General matrix
            ldab = n ;
            lu = new Array<double>(*std) ;

            F77_dgetrf(&n, &n, lu->set_data(), &ldab, permute->set_data(), &info) ;
        }
    }
    return ;
}

// Solution of Ax = B : use LAPACK et the LU decomposition.
Array<double> Matrice::solve (const Array<double>& source) const {
    
    assert(lu != 0x0) ;
    assert(permute != 0x0) ;
       
    int n = source.get_size(0) ;
    assert (sizes(1) == n) ;
    int ldab, info ;
    char trans ;
    int nrhs = 1 ;
    int ldb = n ;
        
    Array<double> res(source) ;
    
    if (band != 0x0) { //Banded matrix
        ldab = 2*kl+ku+1 ;
        trans = 'N' ;
        F77_dgbtrs(&trans, &n, &kl, &ku, &nrhs, lu->set_data(),
               &ldab, permute->set_data(), res.set_data(), &ldb, &info);
    }
    else { // General case
        ldab = n ;
        trans = 'T' ;
        F77_dgetrs(&trans, &n, &nrhs, lu->set_data(), &ldab, permute->set_data(),
               res.set_data(), &ldb, &info) ;
    }
    
    return res ;
}


// Eighenvalues of the matrix, use of LAPACK :
Array<double> Matrice::val_propre() const {
    
    char jobvl = 'N' ;
    char jobvr = 'N' ;
    
    int n = sizes(0) ;
    assert (n == sizes(1)) ;
    
    double* a = new double [n*n] ;
    Index index(sizes) ;
    for (int i=0 ; i<n*n ; i++) {
	a[i] = (*std)(index) ;
	index.inc() ;
	}
	
    int lda = n ;
    double* wr = new double[n] ;
    double* wi = new double[n] ;
    
    int ldvl = 1 ;
    double* vl = 0x0 ;
    int ldvr = 1 ;
    double* vr = 0x0 ;
    
    int ldwork = 3*n ;
    double* work = new double[ldwork] ;
    
    int info ;
    
    F77_dgeev(&jobvl, &jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
	    work, &ldwork, &info) ;
    
   
    Dim_array res_out (2) ;
    res_out.set(0) = 2 ;
    res_out.set(1) = n ;
    Array<double> result(res_out) ;
    Index index_out(res_out) ;
    for (int i=0 ; i<n ; i++) {
	result.set(index_out) = wr[n-i-1] ;
	index_out.inc() ;
	result.set(index_out) = wi[n-i-1] ;
	index_out.inc() ;
	}
    
    delete [] wr ;
    delete [] wi ;
    delete [] a ;
    delete [] work ;
    
    return result ; 
    
}

// Eighenvectors of the matrix (use of LAPACK) :
Matrice Matrice::vect_propre() const {
    
    assert (std != 0x0) ;
    
    char jobvl = 'V' ;
    char jobvr = 'V' ;

    int n = sizes(0) ;
    assert (n == sizes(1)) ;
    
    double* a = new double [n*n] ;
    Index index(sizes) ;
    for (int i=0 ; i<n*n ; i++) {
	a[i] = (*std)(index) ;
	index.inc() ;
	}

    int lda = n ;
    double* wr = new double[n] ;
    double* wi = new double[n] ;
    
    int ldvl = n ;
    double* vl = new double[ldvl*ldvl] ;
    int ldvr = n ;
    double* vr = new double[ldvl*ldvl] ;
    
    int ldwork = 4*n ;
    double* work = new double[ldwork] ;
    
    int info ;
    
    F77_dgeev(&jobvl, &jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
        work, &ldwork, &info) ;
    
   
    Matrice res (n,n) ;

    int conte = 0 ;
    for (int i=0 ; i<n ; i++)
      for (int j=0 ; j<n ; j++) {
	res.set(j,n-i-1) = vr[conte] ;
	conte ++ ;
      }
    
    delete [] wr ;
    delete [] wi ;
    delete [] a ;
    delete [] work ;
    delete [] vl ;

    return res ;
}

// Transposed matrix
Matrice Matrice::transpose() const {

	int nbl = sizes(1) ;
	int nbc = sizes(0) ;

	Matrice resu(nbc, nbl) ;
	
	Index index (std->get_dimensions()) ;
	for (int i=0; i<nbc; i++)
		for (int j=0; j<nbl; j++) {
			index.set(0) = j ;
			index.set(1) = i ;
			resu.set(i,j) = (*std)(index) ;
			}
	return resu ;
}


void Matrice::operator+=(const Matrice& a) {
    assert((std != 0x0)&&(a.std != 0x0)) ;
    std->operator+=(*a.std) ;
}

void Matrice::operator-=(const Matrice& a) {
    assert((std != 0x0)&&(a.std != 0x0)) ;
    std->operator-=(*a.std) ;
}

void Matrice::operator+=(double x) {
    assert(std != 0x0);
    std->operator+=(x) ;
}

void Matrice::operator-=(double x) {
    assert(std != 0x0);
    std->operator-=(x) ;
}

void Matrice::operator*=(double x) {
    assert(std != 0x0);
    std->operator*=(x) ;
}

void Matrice::operator/=(double x) {
    assert(std != 0x0);
    assert(x != 0) ;
    std->operator/=(x) ;
}

Matrice operator+ (const Matrice& a, const Matrice& b) {
    assert((a.std != 0x0) && (b.std != 0x0)) ;
    Matrice res(*a.std+*b.std) ;
    return res ;
}

Matrice operator- (const Matrice& a, const Matrice& b) {
    assert((a.std != 0x0) && (b.std != 0x0)) ;
    Matrice res(*a.std-*b.std) ;
    return res ;
}

Matrice operator- (const Matrice& a) {
    assert(a.std != 0x0) ;
    Matrice res(-*a.std) ;
    return res ;
}


Matrice operator* (const Matrice& a, double x) {
    assert(a.std != 0x0) ;
    Matrice res(*a.std*x);
    return res ;
}

Matrice operator* (double x, const Matrice& a) {
    assert(a.std != 0x0) ;
    Matrice res(*a.std*x);
    return res ;
}

Matrice operator* (const Matrice& aa, const Matrice& bb) {

	int nbla = aa.sizes(1) ;
	int nbca = aa.sizes(0) ;
	int nbcb = bb.sizes(0) ;
	
	assert( nbca ==  bb.sizes(1)) ;
	
	Matrice resu(nbla, nbcb) ;

	for (int i=0; i<nbla; i++)
		for (int j=0; j<nbcb; j++) {
			double sum = 0 ;
			for (int k=0; k<nbca; k++) {
				sum += aa(i,k) * bb(k, j) ;
			}
			resu.set(i,j) = sum ;
		}
    return resu ;
}

Matrice operator/ (const Matrice& a, double x) {
    assert (x != 0) ;
    assert(a.std != 0x0) ;
    Matrice res(*a.std/x);
    return res ;
}}
