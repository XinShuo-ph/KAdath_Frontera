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

#ifndef __BASE_SPECTRAL_HPP_
#define __BASE_SPECTRAL_HPP_

#include "dim_array.hpp"
#include "array.hpp"

enum : int {
    NBR_MAX_BASE    = 30,
    CHEB         = 1,
    CHEB_EVEN    = 2,
    CHEB_ODD     = 3,
    COSSIN       = 4,
    COS          = 5,
    COS_EVEN     = 6,
    COS_ODD      = 7,
    SIN          = 8,
    SIN_EVEN     = 9,
    SIN_ODD     = 10,
    LEG         = 11,
    LEG_EVEN    = 12,
    LEG_ODD     = 13,
    COSSIN_EVEN = 14,
    COSSIN_ODD  = 15
};

namespace Kadath {

    void coef_1d (int, Array<double>&) ;

    class Point ;

    /**
    * Class for storing the basis of decompositions of a field.
    *
    * It mainly consists in a list of integers, encoding the various basis.
    * \li \c CHEB : Chebyshev polynomials.
    * \li \c CHEB_EVEN : even Chebyshev polynomials.
    * \li \c CHEB_ODD : odd Chebyshev polynomials.
    * \li \c COSSIN : series of sines and cosines.
    * \li \c COS : series of cosines.
    * \li \c COS_EVEN : series of even cosines.
    * \li \c COS_ODD : series of odd cosines.
    * \li \c SIN : series of sines.
    * \li \c SIN_EVEN : series of even sines.
    * \li \c SIN_ODD : series of odd sines.
    * \li \c LEG : Legendre polynomials.
    * \li \c LEG_EVEN : even Legendre polynomials.
    * \li \c LEG_ODD : odd Legendre polynomials.
    * \li \c COSSIN_EVEN : series of sines and cosines with even harmonics only.
    * \li \c COSSIN_ODD : series of sines and cosines with odd harmonics only.
    * \ingroup spectral
    **/
    class Base_spectral : public Memory_mapped {
    public:
	//! Sylvain'stuff
        using Bases_container = Memory_mapped_array<Array<int>*,int>;

    protected:
        bool def; ///< \c true if the \c Base_spectral is defined and \c false otherwise.
        int ndim; ///< Number of dimensions.
        /**
        * Arrays containing the various basis of decomposition.
        * The size of each array depends on the order of the various numerical coordinates.
        */
        //	Array<int>** bases_1d ;
        Bases_container bases_1d;

    public:
        /**
        * Standard constructor, the \c Base_spectral is not defined.
        * @param dimensions [input] : number of dimensions.
        */
        explicit Base_spectral(int dimensions) : def{false}, ndim{dimensions}, bases_1d{ndim} {
            for (int i = 0; i < ndim; i++) bases_1d[i] = nullptr;
        }

        Base_spectral(const Base_spectral &so); ///< Copy constructor.
        Base_spectral(FILE *);    ///< Constructor from a file
        Base_spectral(Base_spectral &&so); ///< Move constructor.
        Base_spectral& operator=(Base_spectral && so) ; ///< Move assignment operator.

    public:
        ~Base_spectral() {for (int l = 0; l < ndim; l++) if (bases_1d[l] != nullptr) delete bases_1d[l];} ///< Destructor.
        void save(FILE*) const ; ///< Saving function
    public:
        /**
         * Swaps data between the source and this.
         */
        void swap(Base_spectral & so) ;
        /**
         * @returns : \c true if the basis is defined, \c false otherwise.
         */
        bool is_def () const {return def;} ;

        /**
        * Sets all the basis to the undefined state.
        */
        void set_non_def() ;

        Base_spectral & operator= (const Base_spectral&) ; ///< Assignement operator.

	/// Returns one of the 1d base array.
        const Array<int>* get_base_1d (int i) const {return bases_1d[i] ;} ; 
        /**
        * Allocates the various arrays, for a given number of coefficients.
        * @param nbr_coefs [input] : a \c Dim_array storing the number of coefficients in each dimenions.
        */
        void allocate (const Dim_array& nbr_coefs) ;
        /**
        * Allocates the various arrays, for a given number of coefficients and sets basis to some values (same for all harmonics).
        * It assumes one is working in 3D.
        * @param nbr_coefs [input] : a \c Dim_array storing the number of coefficients in each dimenions.
        * @param basephi : basis for \f$\varphi\f$
        * @param basetheta : basis for \f$\varphi\f$
        * @param baser : basis for \f$\varphi\f$
        */
        void set(Dim_array const& nbr_coefs, int basephi, int basetheta, int baser);

        /**
            * performs the coefficient transformation for one particular variable.
        * @param var [input] : the variable for which the coefficients are to be computed.
        * @param nbr [input] : number of coefficients for \c var.
        * @param tab [input/output] : values of the field, before and after the computation of the coefficients.
        */
        void coef_dim (int var, int nbr, Array<double> *& tab) const ;
        /**
            * Computes the coefficients for all the variables.
        * @param nbr_coefs [input] : number of coefficients.
            * @param so [input] : values of the field in the configuration space.
        * @returns the coefficients of the field.
        */
        Array<double> coef (const Dim_array& nbr_coefs, const Array<double> & so) const ;
       /**
        * Performs the inverse coefficient transformation for one particular variable.
        * @param var [input] : the variable for which the coefficients are to be computed.
        * @param nbr [input] : number of points in the configuration space for \c var.
        * @param tab [input/output] : values of the field, before and after the computation of the coefficients.
        */
        void coef_i_dim (int var, int nbr, Array<double> *& tab) const ;
        /**
            * Computes the values in the configuration space.
        * @param nbr_points [input] : number of points.
            * @param so [input] : values of the field in the coefficients space.
        * @returns the coefficients of the field.
        */
        Array<double> coef_i (const Dim_array& nbr_points, const Array<double> & so) const ;
	/**
	* Computes the spectral summation.
	* @param num [input] : numerical coordinates used in the summation.
        * @param tab [input] : spectral coefficients of the field.
        * @returns the summation.
        */
        double summation (const Point& num, const Array<double>& tab) const ;
        /**
        * One-dimensional operator acting in the coefficient space.
        * @param function [input] : function pointing to the particular operation to be performed (i.e. like
        * \c mult_sin_x_1d )
        * @param var [input] : variable on which the operatio is to be performed.
        * @param so [input] : coefficients of the field before the operation.
        * @param base [input/output] : basis of the field. After the operation, the basis associated with \c var
        * may is changed but not the others basis which should be changed by hand, if needed.
        * @returns the coefficients of the result.
        */
        Array<double> ope_1d (int (*function) (int, Array<double>&),
                    int var , const Array<double>& so, Base_spectral& base) const ;

         friend bool operator== (const Base_spectral&, const Base_spectral&) ; ///< Comparison operator
         friend ostream& operator<< (ostream& o, const Base_spectral&) ; ///< Display
         friend class Space ;
         friend class Space_spheric ;
         friend class Scalar ;
         friend class Domain_nucleus ;
         friend class Domain_shell ;
         friend class Domain_shell_log ;
         friend class Domain_compact ;
         friend class Domain_bispheric_rect ;
         friend class Domain_bispheric_chi_first ;
         friend class Domain_bispheric_eta_first ;
         friend class Domain_critic_inner ;
         friend class Domain_critic_outer ;
         friend class Domain_polar_nucleus ;
         friend class Domain_polar_shell ;
         friend class Domain_polar_shell_inner_adapted ;
         friend class Domain_polar_shell_outer_adapted ;
         friend class Domain_polar_compact ;
         friend class Domain_oned_ori ;
         friend class Domain_oned_qcq ;
         friend class Domain_oned_inf ;
         friend class Domain_spheric_periodic_nucleus ;
         friend class Domain_spheric_periodic_shell ;
         friend class Domain_spheric_periodic_compact ;
         friend class Domain_spheric_time_nucleus ;
         friend class Domain_spheric_time_shell ;
         friend class Domain_spheric_time_compact ;
         friend class Domain_shell_inner_adapted ;
         friend class Domain_shell_outer_adapted ;
         friend class Domain_shell_inner_homothetic ;
         friend class Domain_shell_outer_homothetic ;
         friend class Domain_nucleus_symphi ;
         friend class Domain_shell_symphi ;
         friend class Domain_compact_symphi ;
         friend class Domain_polar_periodic_nucleus ;
         friend class Domain_polar_periodic_shell ;
         friend class Domain_fourD_periodic_nucleus ;
         friend class Domain_fourD_periodic_shell ;
    } ;

    inline Base_spectral::Base_spectral(const Base_spectral &so) : def{so.def}, ndim{so.ndim}, bases_1d{ndim} {
        for (int l = 0; l < ndim; l++)
            bases_1d[l] = (so.bases_1d[l] == nullptr) ? nullptr : new Array<int>(*so.bases_1d[l]);
    }

    inline Base_spectral::Base_spectral(Base_spectral &&so) : def{so.def}, ndim{0}, bases_1d{} {
        std::swap(ndim, so.ndim);
        std::swap(bases_1d, so.bases_1d);
        so.def = false;
    }

    inline Base_spectral& Base_spectral::operator=(Base_spectral && so) {
        std::swap(def,so.def);
        ndim = so.ndim;
        std::swap(bases_1d,so.bases_1d);
        return *this;
    }

    inline void Base_spectral::set_non_def() {
        for (int l=0 ; l<ndim ; l++)
            if (bases_1d[l] !=0x0) {
                delete bases_1d[l] ;
                bases_1d[l] = 0x0 ;
            }
        def = false ;
    }

    inline Base_spectral & Base_spectral::operator= (const Base_spectral& so) {
        assert (ndim==so.ndim) ;
        def = so.def ;
        for (int i=0 ; i<ndim ; i++) {
            if (bases_1d[i] != nullptr) delete bases_1d[i] ;
            if ((so.def) && (so.bases_1d[i] != nullptr))
                bases_1d[i] = new Array<int>(*so.bases_1d[i]) ;
            else
                bases_1d[i] = nullptr ;
        }
        return *this;
    }

    inline void Base_spectral::swap(Base_spectral & so) {
        std::swap(def,so.def);
        std::swap(ndim,so.ndim);
        bases_1d.swap(so.bases_1d);
    }

    inline bool operator== (const Base_spectral& a, const Base_spectral& b) {
        bool res {(a.def) && (b.def) && (a.ndim == b.ndim)} ;
        for (int i=0 ; i<a.ndim ; i++) {
            Array_iterator index (a.bases_1d[i]->get_dimensions()) ;
            do
                res = ((*a.bases_1d[i])(index) == (*b.bases_1d[i])(index)) ;
            while (index.inc() && res) ;
            if (!res)
            	break ;
        }
        return res ;
    }
}
#endif
