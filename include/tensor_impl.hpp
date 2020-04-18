//
// Created by sauliac on 18/04/2020.
//

#include "headcpp.hpp"
#include "tensor.hpp"
#include "scalar.hpp"
#include "vector.hpp"
#include "metric_tensor.hpp"

#ifndef __TENSOR_IMPL_HPP_
#define __TENSOR_IMPL_HPP_
namespace Kadath {
//    int std_position_array (const Array<int>& , int );
//    int std_position_index (const Index& , int );
//    Array<int> std_indices (int , int , int );

    inline void Tensor::annule_hard() {
        for (int i=0 ; i<n_comp ; i++)
            cmp[i]->annule_hard() ;
    }

    //-------------
    // Accessors
    //-------------

    // Affectation d'un tenseur d'ordre 1 :
    inline Scalar& Tensor::set(int i) {

        assert (valence == 1) ;

        Array<int> ind (valence) ;
        ind.set(0) = i ;

        int place = position(ind) ;

        return *cmp[place] ;
    }

    // Affectation d'un tenseur d'ordre 2 :
    inline Scalar& Tensor::set(int ind1, int ind2) {

        assert (valence == 2) ;

        Array<int> ind (valence) ;
        ind.set(0) = ind1 ;
        ind.set(1) = ind2 ;

        int place = position(ind) ;

        return *cmp[place] ;
    }

    // Affectation d'un tenseur d'ordre 3 :
    inline Scalar& Tensor::set(int ind1, int ind2, int ind3) {

        assert (valence == 3) ;

        Array<int> idx(valence) ;
        idx.set(0) = ind1 ;
        idx.set(1) = ind2 ;
        idx.set(2) = ind3 ;
        int place = position(idx) ;

        return *cmp[place] ;
    }


    // Affectation d'un tenseur d'ordre 4 :
    inline Scalar& Tensor::set(int ind1, int ind2, int ind3, int ind4) {

        assert (valence == 4) ;

        Array<int> idx(valence) ;
        idx.set(0) = ind1 ;
        idx.set(1) = ind2 ;
        idx.set(2) = ind3 ;
        idx.set(3) = ind4 ;
        int place = position(idx) ;

        return *cmp[place] ;
    }


    // Affectation cas general
    inline Scalar& Tensor::set(const Array<int>& idx) {

        assert (idx.get_ndim() == 1) ;
        assert (idx.get_size(0) == valence) ;

        int place = position(idx) ;
        return *cmp[place] ;
    }

    // Affectation cas general from an Index
    inline Scalar& Tensor::set(const Index& idx) {

        Array<int> ind (valence) ;
        for (int i=0 ; i<valence ; i++)
            ind.set(i) = idx(i)+1 ;

        return set(ind) ;
    }

    inline const Scalar& Tensor::operator()() const {

        assert(valence == 0) ;

        return *cmp[0] ;

    }

    inline const Scalar& Tensor::operator()(int indice) const {

        assert(valence == 1) ;

        Array<int> idx(1) ;
        idx.set(0) = indice ;
        return *cmp[position(idx)] ;

    }

    inline const Scalar& Tensor::operator()(int indice1, int indice2) const {

        assert(valence == 2) ;

        Array<int> idx(2) ;
        idx.set(0) = indice1 ;
        idx.set(1) = indice2 ;
        return *cmp[position(idx)] ;

    }

    inline const Scalar& Tensor::operator()(int indice1, int indice2, int indice3) const {

        assert(valence == 3) ;

        Array<int> idx(3) ;
        idx.set(0) = indice1 ;
        idx.set(1) = indice2 ;
        idx.set(2) = indice3 ;
        return *cmp[position(idx)] ;
    }


    inline const Scalar& Tensor::operator()(int indice1, int indice2, int indice3,
                                     int indice4) const {

        assert(valence == 4) ;

        Array<int> idx(4) ;
        idx.set(0) = indice1 ;
        idx.set(1) = indice2 ;
        idx.set(2) = indice3 ;
        idx.set(3) = indice4 ;
        return *cmp[position(idx)] ;
    }

    inline const Scalar& Tensor::at(int indice1, int indice2) const {
        return operator()(indice1,indice2);

    }

    inline const Scalar& Tensor::operator()(const Array<int>& ind) const {

        assert (ind.get_ndim() == 1) ;
        assert (ind.get_size(0) == valence) ;
        return *cmp[position(ind)] ;

    }

    inline const Scalar& Tensor::operator()(const Index& idx) const {
        Array<int> ind (valence) ;
        for (int i=0 ; i<valence ; i++)
            ind.set(i) = idx(i)+1 ;

        return operator()(ind) ;

    }

    inline void Tensor::set_name_ind (int pos, char name) {
        assert ((pos>=0) && (pos<valence)) ;
        if (!name_affected)
            name_affected = true ;
        name_indice[pos] = name ;
    }

    inline void Tensor::coef() const
    {
        Index pos(*this);
        do
        {
            (*this)(pos).coef();
        }while(pos.inc());
    }

    inline void Tensor::coef_i() const
    {
        Index pos(*this);
        do
        {
            (*this)(pos).coef_i();
        }while(pos.inc());
    }

    inline void Tensor::filter_phi(int dom, int ncf)
    {
        Index pos(*this);
        do
        {
            set(pos).filter_phi(dom, ncf);
        }while(pos.inc());
    }


    inline void Tensor::filter (double threshold)  {

        for (int d=0 ; d<espace.get_nbr_domains() ; d++) {
            espace.get_domain(d)->filter(*this, d, threshold) ;
        }

    }
}
#endif //__TENSOR_IMPL_HPP_
