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

#ifndef __METRIC_symphi_HPP_
#define __METRIC_symphi_HPP_
#include "metric.hpp"

namespace Kadath {
class Space ;
class Tensor ;
class Term_eq ;
class Metric_tensor ;
class System_of_eqs ;
class Base_tensor ;

/**
* Class that deals with flat metric assuming a axisymmetric setting (nothing depends on \f$varphi\f$).
* \ingroup metric
*/
class Metric_flat_symphi : public Metric {

  protected:
    const Base_tensor& basis ; ///< The tensorial basis used.
  
	public:
		/**
		* Standard constructor
		* @param sp : the associated \c Space.
		*/
		Metric_flat_symphi (const Space& sp, const Base_tensor& bb) ;
		Metric_flat_symphi (const Metric_flat_symphi& ) ; ///< Copy constructor

	protected:
		virtual void compute_con (int) const ;
		virtual void compute_cov (int) const ;		
		virtual void compute_christo (int) const ;
		
  private:
	/**
	* Computes the partial derivative part of the covariant derivative, in orthonormal spherical coordinates.
	*  The \f$varphi\f$ component is set to zero.
	* The result is \f$ \partial_r, \frac{1}{r}\partial_\theta, 0\f$.
	* If need be inner summation of the result is performed and index manipulated.
	* @param tder : type of the result (COV or CON)
	* @param indder : name of the index corresponding to the derivative.
	* @param so : the field to be derived.
	* @return : the result as a \c Term_eq.
	*/
	Term_eq derive_partial_spher (int tder, char indder, const Term_eq& so) const ;
	
	/**
	* Computes the flat covariant derivative, in orthonormal spherical coordinates.
	* If need be inner summation of the result is performed and index manipulated.
	* @param tder : type of the result (COV or CON)
	* @param indder : name of the index corresponding to the derivative.
	* @param so : the field to be derived.
	* @return : the result as a \c Term_eq.
	*/
	Term_eq derive_spher (int tder, char indder, const Term_eq& so) const ;

	/**
	* Computes the flat covariant derivative, in orthonormal spherical coordinates.
	* If need be inner summation of the result is performed.
	* If the CON derivative is asked for, the index is NOT manipulated by the flat metric but by an arbitrary one.
	* @param tder : type of the result (COV or CON)
	* @param indder : name of the index corresponding to the derivative.
	* @param so : the field to be derived.
	* @param othermet : pointer on the \c Metric used for index manipulation.
	* @return : the result as a \c Term_eq.
	*/
	Term_eq derive_with_other_spher (int tder, char indder, const Term_eq& so, const Metric* othermet) const ;

	public:
		virtual ~Metric_flat_symphi() ;	
		virtual void update() ;	
		virtual void update(int) ;
		virtual void manipulate_ind (Term_eq&, int) const ;
		virtual Term_eq derive_partial (int, char, const Term_eq&) const ;
		virtual Term_eq derive (int, char, const Term_eq&) const ;

		/**
		* Computes the flat covariant derivative.
		* If need be inner summation of the result is performed.
		* If the CON derivative is asked for, the index is NOT manipulated by the flat metric but by an arbitrary one.
		* @param tder : type of the result (COV or CON)
		* @param indder : name of the index corresponding to the derivative.
		* @param so : the field to be derived.
		* @param othermet : pointer on the \c Metric used for index manipulation.
		* @return : the result as a \c Term_eq.
		*/
		Term_eq derive_with_other (int tder, char indder, const Term_eq& so, const Metric* othermet) const ;

		/**
		* Associate the metric to a given system of equations.
		* @param syst : the \c System_of_eqs.
		* @param name : name by which the metric will be known in the system (like "g", "f"...)
		*/
		virtual void set_system (System_of_eqs& syst, const char* name) ;	
		
		virtual int give_type (int) const ;
} ;

/**
 * Class to deal with a metric independant of \f$\varphi\f$.
 * \ingroup metric
 */
class Metric_symphi : public Metric {

	protected:
	  Metric_tensor* p_met ; ///< Pointer on the \c Metric_tensor describing the metric.
	  const Base_tensor& basis ; ///< The tensorial basis used.
	  Metric_flat_symphi fmet ; ///< Associated flat metric.
	  int place_syst ; ///< Gives the location of the metric amongst the various unknowns of the associated \c System_of_eqs.

	public:
		Metric_symphi (Metric_tensor&) ; ///< Constructor from a \c Metric_tensor.
		Metric_symphi (const Metric_symphi& ) ; ///< Copy constructor.

	protected:
		virtual void compute_con (int) const ;
		virtual void compute_cov (int) const ;		
		virtual void compute_christo (int) const ;		
		virtual void compute_riemann (int) const ;		
		virtual void compute_ricci_tensor (int) const ;
		
	public:
		virtual ~Metric_symphi() ;
		virtual Term_eq derive (int, char, const Term_eq&) const ;
		/**
		* Associate the metric to a given system of equations.
		* @param syst : the \c System_of_eqs.
		* @param name : name by which the metric will be known in the system (like "g", "f"...)
		*/
		virtual void set_system (System_of_eqs& syst, const char* name) ;
		virtual Term_eq derive_flat (int, char, const Term_eq&) const ;
		
		virtual int give_type (int) const ;
} ;

/**
 * Class to deal with a metric independant of \f$\varphi\f$ and constant. By this one means that it is a given and is not modified by the \c System_of_eqs.
 * \ingroup metric
 */
class Metric_symphi_const : public Metric_symphi {

	public:
		Metric_symphi_const (Metric_tensor&) ; ///< Constructor from a \c Metric_tensor.
		Metric_symphi_const (const Metric_symphi_const& ) ; ///< Copy constructor
	
	protected:
		virtual void compute_con (int) const ;
		virtual void compute_cov (int) const ;


	public:
		virtual ~Metric_symphi_const() ;

		/**
		* Associate the metric to a given system of equations.
		* @param syst : the \c System_of_eqs.
		* @param name : name by which the metric will be known in the system (like "g", "f"...)
		*/
		virtual void set_system (System_of_eqs& syst, const char* name) ;
} ;
}
#endif
