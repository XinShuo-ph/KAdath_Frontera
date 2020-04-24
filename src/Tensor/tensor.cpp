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

/*
 *  Methods of class Tensor
 *

 */


#include "headcpp.hpp"
#include "tensor.hpp"
#include "scalar.hpp"
#include "tensor_impl.hpp"
#include "vector.hpp"
#include "metric_tensor.hpp"

namespace Kadath {
// Functions standard for storage :
    int std_position_array (const Array<int>& idx, int ndim) {
        assert (idx.get_ndim() == 1) ;
        int valence =idx.get_size(0) ;

        for (int i=0 ; i<valence ; i++)
        assert ((idx(i)>=1) && (idx(i)<=ndim)) ;
        int res = 0 ;
        for (int i=0 ; i<valence ; i++)
            res = ndim*res+(idx(i)-1) ;

        return res;
    }

    int std_position_index (const Index& idx, int ndim) {
       int valence =idx.get_ndim() ;
       for (int i=0 ; i<valence ; i++)
        assert ((idx(i)>=0) && (idx(i)<ndim)) ;
        int res = 0 ;
        for (int i=0 ; i<valence ; i++)
            res = ndim*res+(idx(i)) ;

        return res;
    }

    Array<int> std_indices (int place, int valence, int ndim) {
        Array<int> res(valence) ;
        for (int i=valence-1 ; i>=0 ; i--) {
            res.set(i) = div(place, ndim).rem ;
            place = int((place-res(i))/ndim) ;
            res.set(i)++ ;
        }
        return res ;
    }


    //--------------//
    // Constructors //
    //--------------//

// Standard constructor
// --------------------
    Tensor::Tensor(const Space& sp, int val, const Array<int>& tipe, const Base_tensor& bb)
            : espace(sp), ndom(espace.get_nbr_domains()), ndim(espace.get_ndim()),
              valence(val), basis(bb), type_indice(tipe), name_affected{false}, name_indice{valence},
              n_comp(int(pow(espace.get_ndim(), val))), cmp{n_comp},  parameters{}
    {
        assert (valence >= 0) ;
        assert (tipe.get_ndim() == 1) ;
        assert (valence == tipe.get_size(0)) ;
#ifdef DEBUG_MODE
        for (int i=0 ; i<valence ; i++)
            assert ((tipe(i) == COV) || (tipe(i) == CON)) ;
#endif
        for (int i=0 ; i<n_comp ; i++) cmp[i] = new Scalar{espace};

        // Storage methods :
        give_place_array = std_position_array ;
        give_place_index = std_position_index ;
        give_indices = std_indices ;
    }

    Tensor::Tensor(const Space& sp, int val, int tipe, const Base_tensor& bb)
            : espace(sp), ndom(espace.get_nbr_domains()), ndim(espace.get_ndim()),
              valence(val), basis(bb), type_indice(val), name_affected{false}, name_indice{valence},
              n_comp(int(pow(espace.get_ndim(), val))), cmp{n_comp}, parameters{}
    {
        assert (valence >= 0) ;
        assert ((tipe == COV) || (tipe == CON)) ;
        type_indice = tipe;
        for (int i=0 ; i<n_comp ; i++) cmp[i] = new Scalar{espace} ;

        // Storage methods :
        give_place_array = std_position_array ;
        give_place_index = std_position_index ;
        give_indices = std_indices ;
    }

    Tensor::Tensor(const Space& sp, int val, const Array<int>& tipe, const Base_tensor& bb, int dim)
            : espace(sp), ndom(espace.get_nbr_domains()), ndim(dim),
              valence(val), basis(bb), type_indice(tipe), name_affected{false}, name_indice{valence},
              n_comp(int(pow(ndim, val))), cmp{n_comp}, parameters{}
    {
        assert (valence >= 0) ;
        assert (tipe.get_ndim() == 1) ;
        assert (valence == tipe.get_size(0)) ;
#ifdef DEBUG_MODE
        for (int i=0 ; i<valence ; i++)
            assert ((tipe(i) == COV) || (tipe(i) == CON)) ;
#endif
        for (int i=0 ; i<n_comp ; i++) cmp[i] = new Scalar{espace} ;

        // Storage methods :
        give_place_array = std_position_array ;
        give_place_index = std_position_index ;
        give_indices = std_indices ;
    }

    Tensor::Tensor(const Space& sp, int val, int tipe, const Base_tensor& bb, int dim)
            : espace(sp), ndom(espace.get_nbr_domains()), ndim(dim),
              valence(val), basis(bb), type_indice(val), name_affected{false}, name_indice{valence},
              n_comp(int(pow(ndim, val))), cmp{n_comp}, parameters{}
    {
        assert (valence >= 0) ;
        assert ((tipe == COV) || (tipe == CON)) ;
        type_indice = tipe;
        for (int i=0 ; i<n_comp ; i++) cmp[i] = new Scalar{espace} ;

        // Storage methods :
        give_place_array = std_position_array ;
        give_place_index = std_position_index ;
        give_indices = std_indices ;
    }



// Copy constructor
// ----------------
    Tensor::Tensor (const Tensor& source, bool copy) :
            espace(source.espace), ndom(source.ndom), ndim(source.ndim), valence(source.valence), basis(source.basis),
            type_indice(source.type_indice), name_affected{copy && source.name_affected}, name_indice{valence},
            n_comp (source.n_comp), cmp{n_comp}, parameters{source.parameters}
    {
        for (int i=0 ; i<n_comp ; i++) cmp[i] = new Scalar{*source.cmp[i], copy} ;
        if (name_affected) {
            for (int i=0 ; i<valence ; i++)
                name_indice[i] = source.name_indice[i] ;
        }

        // Storage methods :
        give_place_array = source.give_place_array ;
        give_place_index = source.give_place_index ;
        give_indices = source.give_indices ;
    }

//  Constructor for a scalar field: to be used by the derived
//  class {\tt Scalar}
//-----------------------------------------------------------
    Tensor::Tensor(const Space& sp) : espace(sp), ndom(espace.get_nbr_domains()), ndim(espace.get_ndim()),
                                      valence(0), basis(sp), type_indice(0), name_affected{false}, name_indice{},
                                      n_comp{1}, cmp{1}, parameters{}
    {
        cmp[0] = nullptr ;

        // Storage methods :
        give_place_array = std_position_array ;
        give_place_index = std_position_index ;
        give_indices = std_indices ;
    }

    Tensor::Tensor(const Space& sp, int val, int tipe, int ncompi, const Base_tensor& bb)
            : espace(sp), ndom(espace.get_nbr_domains()), ndim(espace.get_ndim()),
              valence(val), basis(bb), type_indice(val), name_affected{false}, name_indice{valence},
              n_comp(ncompi), cmp{n_comp}, parameters{}
    {
        assert (valence >= 0) ;
        assert ((tipe == COV) || (tipe == CON)) ;
        type_indice = tipe;
        for (int i=0 ; i<n_comp ; i++) cmp[i] = new Scalar{espace} ;
    }

    Tensor::Tensor(const Space& sp, int val, const Array<int>& tipe, int ncompi, const Base_tensor& bb)
            : espace(sp), ndom(espace.get_nbr_domains()), ndim(espace.get_ndim()),
              valence(val), basis(bb), type_indice(tipe), name_affected{false}, name_indice{valence},
              n_comp(ncompi), cmp{n_comp}, parameters{}
    {
        // Des verifs :
        assert (valence >= 0) ;
        assert (tipe.get_ndim() == 1) ;
        assert (valence == tipe.get_size(0)) ;
#ifdef DEBUG_MODE
        for (int i=0 ; i<valence ; i++)
            assert ((tipe(i) == COV) || (tipe(i) == CON)) ;
#endif
        for (int i=0 ; i<n_comp ; i++) cmp[i] = new Scalar{espace} ;
    }

    Tensor::Tensor(const Space& sp, int val, int tipe, int ncompi, const Base_tensor& bb, int dim)
            : espace(sp), ndom(espace.get_nbr_domains()), ndim(dim),
              valence(val), basis(bb), type_indice(val), name_affected{false}, name_indice{valence},
              n_comp(ncompi), cmp{n_comp}, parameters{}
    {
        assert (valence >= 0) ;
        assert ((tipe == COV) || (tipe == CON)) ;
        type_indice = tipe;
        for (int i=0 ; i<n_comp ; i++) cmp[i] = new Scalar{espace} ;
    }

    Tensor::Tensor(const Space& sp, int val, const Array<int>& tipe, int ncompi, const Base_tensor& bb, int dim)
            : espace(sp), ndom(espace.get_nbr_domains()), ndim(dim),
              valence(val), basis(bb), type_indice(tipe), name_affected{false}, name_indice{valence},
              n_comp(ncompi), cmp{n_comp}, parameters{}
    {
        // Des verifs :
        assert (valence >= 0) ;
        assert (tipe.get_ndim() == 1) ;
        assert (valence == tipe.get_size(0)) ;
#ifdef DEBUG_MODE
        for (int i=0 ; i<valence ; i++)
            assert ((tipe(i) == COV) || (tipe(i) == CON)) ;
#endif
        for (int i=0 ; i<n_comp ; i++) cmp[i] = new Scalar{espace} ;
    }


    Tensor::Tensor (const Space& sp, FILE* fd) :
        espace(sp), ndom(espace.get_nbr_domains()), ndim(espace.get_ndim()),
        basis(sp, fd), type_indice(fd), name_affected{false}, name_indice{valence}, n_comp{}, cmp{}, parameters{}
    {
        fread_be (&valence, sizeof(int), 1, fd) ;
        assert (type_indice.get_size(0) == valence) ;
        fread_be (&n_comp, sizeof(int), 1, fd) ;
        cmp.resize(n_comp);
        for (int i=0 ; i<n_comp ; i++) cmp[i] = new Scalar{espace, fd} ;
        // Storage methods (overwritten in case of non std things)
        give_place_array = std_position_array ;
        give_place_index = std_position_index ;
        give_indices = std_indices ;
    }

    Tensor::Tensor (const Space& sp, int dim, FILE* fd) :
        espace(sp), ndom(espace.get_nbr_domains()), ndim(dim),
        basis(sp, fd), type_indice(fd), name_affected{false}, name_indice{valence}, n_comp{}, cmp{}, parameters{}
    {
        fread_be (&valence, sizeof(int), 1, fd) ;
        assert (type_indice.get_size(0) == valence) ;
        fread_be (&n_comp, sizeof(int), 1, fd) ;
        cmp.resize(n_comp);
        for (int i=0 ; i<n_comp ; i++) cmp[i] = new Scalar{espace, fd} ;
        // Storage methods (overwritten in case of non std things)
        give_place_array = std_position_array ;
        give_place_index = std_position_index ;
        give_indices = std_indices ;
    }

    void Tensor::swap(Tensor & so) noexcept {
        assert(&espace == &so.espace);
        std::swap(ndom,so.ndom);
        std::swap(ndim,so.ndim);
        std::swap(valence,so.valence);
        basis.swap(so.basis);
        type_indice.swap(so.type_indice);
        std::swap(name_affected,so.name_affected);
        name_indice.swap(so.name_indice);
        std::swap(n_comp,so.n_comp);
        cmp.swap(so.cmp);
        parameters.swap(so.parameters);
        std::swap(give_place_array,so.give_place_array);
        std::swap(give_place_index,so.give_place_index);
        std::swap(give_indices,so.give_indices);
    }

#ifdef TENSOR_MOVE_SEMANTIC
    Tensor::Tensor(Tensor&& so) noexcept: espace{so.espace}, ndom{so.ndom}, ndim{so.ndim}, valence{so.valence},
    basis{std::move(so.basis)}, type_indice{std::move(so.type_indice)}, name_affected{so.name_affected},
    name_indice{std::move(so.name_indice)}, n_comp{so.n_comp}, cmp{std::move(so.cmp)},
    parameters{std::move(so.parameters)}, give_place_array{so.give_place_array},
    give_place_index{so.give_place_index}, give_indices{so.give_indices}
    {}

    void Tensor::do_move(Tensor &&so, bool move_cmp) noexcept
    {
        std::swap(ndom,so.ndom);
        ndim = so.ndim;
        valence = so.valence;
        basis = std::move(so.basis);
        type_indice = std::move(so.type_indice);
        name_affected = so.name_affected;
        std::swap(name_indice,so.name_indice);
        std::swap(n_comp,so.n_comp);
        if(move_cmp) cmp = std::move(so.cmp);
        std::swap(parameters,so.parameters);
        give_place_array = so.give_place_array;
        give_place_index = so.give_place_index;
        give_indices = so.give_indices;
    }

    Tensor & Tensor::operator=(Tensor && so) noexcept
    {
        this->do_move(std::move(so),true);
        return *this;
    }
#endif //#ifdef TENSOR_MOVE_SEMANTIC
                //--------------//
                //  Destructor  //
                //--------------//


    Tensor::~Tensor () {
        for(auto & v : cmp) if(v != nullptr) delete v;
    }

    void Tensor::save (FILE* fd) const {
        basis.save(fd) ;
        type_indice.save(fd) ;
        fwrite_be (&valence, sizeof(int), 1, fd) ;
        assert (type_indice.get_size(0) == valence) ;
        fwrite_be (&n_comp, sizeof(int), 1, fd) ;
        for (int i=0 ; i<n_comp ; i++)
            cmp[i]->save(fd) ;
    }


    // Le cout :
    ostream& operator<<(ostream& flux, const Tensor &source ) {

        flux << '\n' ;
        flux << "Class : " << typeid(source).name()
            << "           Valence : " << source.valence << '\n' ;

        if (source.valence!=0) {
            flux << source.get_basis() << endl ;
        }

        if (source.valence != 0) {
            flux <<    "Type of the indices : " ;
            for (int i=0 ; i<source.valence ; i++) {
                flux << "index " << i << " : " ;
                if (source.type_indice(i) == CON)
                    flux << " contravariant." << '\n' ;
                else
                    flux << " covariant." << '\n' ;
                if ( i < source.valence-1 ) flux << "                      " ;
            }
            flux << '\n' ;
        }

        for (int i=0 ; i<source.n_comp ; i++) {

            if (source.valence == 0) {
                flux <<
                "===================== Scalar field ========================= \n" ;
            }
            else {
                flux << "================ Component " ;
                Array<int> num_indices (source.indices(i)) ;
                for (int j=0 ; j<source.valence ; j++) {
                    flux << " " << num_indices(j) ;
                }
                flux << " ================ \n" ;
            }
            flux << '\n' ;

            flux << *source.cmp[i] << '\n' ;
        }

        return flux ;
    }

    // Sets the standard spectal bases of decomposition for each component
    void Tensor::std_base() {
        // Loop on the domains
        for (int d=0 ; d<ndom ; d++) {

        switch (valence) {

            case 0 : {
                if (!is_m_quant_affected())
                  cmp[0]->std_base_domain(d) ;
                else
                  cmp[0]->std_base_domain(d, parameters.get_m_quant()) ;
                break ;
            }
            case 1 : {
                bool done = false ;
                if (basis.get_basis(d) ==CARTESIAN_BASIS) {
                  cmp[0]->std_base_x_cart_domain(d) ;
                  cmp[1]->std_base_y_cart_domain(d) ;
                  cmp[2]->std_base_z_cart_domain(d) ;
                  done = true ;
                }
                if (basis.get_basis(d) == SPHERICAL_BASIS) {
                  cmp[0]->std_base_r_spher_domain(d) ;
                  cmp[1]->std_base_t_spher_domain(d) ;
                  cmp[2]->std_base_p_spher_domain(d) ;
                  done = true ;
                  }
                if (basis.get_basis(d) == MTZ_BASIS) {
                  cmp[0]->std_base_r_mtz_domain(d) ;
                  cmp[1]->std_base_t_mtz_domain(d) ;
                  cmp[2]->std_base_p_mtz_domain(d) ;
                  done = true ;
                  }
#ifndef REMOVE_ALL_CHECKS
                if (!done) {
                    cerr << "Tensor::std_base not yet implemented for " << basis << endl ;
                    abort() ;
                }
#endif
                break ;
                }
            default : {
                Vector auxi (espace, CON, basis) ;
                auxi.std_base() ;

                for (int cc=0 ; cc<n_comp ; cc++) {
                      Array<int> ind (indices(cc)) ;
                      Base_spectral result (auxi(ind(0))(d).get_base()) ;
                      for (int i=1 ; i<valence ; i++)
                        result =  espace.get_domain(d)->mult(result, auxi(ind(i))(d).get_base()) ;
                      cmp[cc]->set_domain(d).set_base() = result ;

                  }
                }
                break ;
            }
        }
    }


    bool Tensor::find_indices (const Tensor& tt, Array<int>& perm) const {

        assert (name_affected) ;
        bool res = true ;
        Array<bool> found (valence) ;
        found = false ;
        for (int ncmp=0 ; ncmp<valence ; ncmp++) {
            int pos=0 ;
            bool finloop = false ;
            do {
                if ((!found(pos)) && (tt.name_indice[pos]==name_indice[ncmp])) {
                    found.set(pos) = true ;
                    perm.set(ncmp) = pos ;
                    finloop = true ;
                    // Check valence :
                    if (type_indice(ncmp) != tt.type_indice(pos))
                        res = false ;
                }
                pos ++ ;
                if ((pos == valence) && (!finloop)) {
                    finloop = true ;
                    res = false ;
                }
            }
            while ((!finloop) && (res));
        }
        return res ;
    }

    void Tensor::change_basis_spher_to_cart()
    {
       for (int d(0) ; d <= espace.get_nbr_domains() - 1 ; ++d)
       {
          Tensor auxi(espace, get_valence(), get_index_type(), get_basis());
          auxi = espace.get_domain(d)->change_basis_spher_to_cart(d, *this);
          Index pos(auxi);
          do
          {
             set(pos).set_domain(d) = auxi(pos)(d);
          }while(pos.inc());
          set_basis(d) = CARTESIAN_BASIS;
       }
    }

    void Tensor::change_basis_cart_to_spher()
    {
       for (int d(0) ; d <= espace.get_nbr_domains() - 1 ; ++d)
       {
          Tensor auxi(espace, get_valence(), get_index_type(), get_basis());
          auxi = espace.get_domain(d)->change_basis_cart_to_spher(d, *this);
          Index pos(auxi);
          do
          {
             set(pos).set_domain(d) = auxi(pos)(d);
          }while(pos.inc());
          set_basis(d) = SPHERICAL_BASIS;
       }
    }



}
