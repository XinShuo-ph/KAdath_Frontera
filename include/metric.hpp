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

#ifndef __METRIC_HPP_
#define __METRIC_HPP_
#include "metric_tensor.hpp"
#include "vector.hpp"
#include "scalar.hpp"
#include "memory.hpp"

namespace Kadath {
class Space ;
class Tensor ;
class Term_eq ;
class Metric_tensor ;
class System_of_eqs ;
class Base_tensor ;


/**
* Purely abstract class for metric handling. Can not be instanciated directly. One must use the various derived classes.
* \ingroup metric
*/
class Metric : public MemoryMappable {

	protected:
		const Space& espace ; ///< The associated \c Space

		const System_of_eqs* syst ; ///< Pointer of the system of equations where the metric is used (only one for now).

		/**
		* Array of pointers on various \c Term_eq.
		* Each one points onto one component of the covariant representation of the \c Metric, in a given \c Domain.
		*/
		mutable Term_eq** p_met_cov ;
		/**
		* Array of pointers on various \c Term_eq.
		* Each one points onto one component of the contravariant representation of the \c Metric, in a given \c Domain.
		*/
		mutable Term_eq** p_met_con ;
		/**
		* Array of pointers on various \c Term_eq.
		* Each one points onto one component of the Christoffel symbols, in a given \c Domain.
		*/
		mutable Term_eq** p_christo ;
		/**
		* Array of pointers on various \c Term_eq.
		* Each one points onto one component of the Riemann tensor, in a given \c Domain.
		*/
		mutable Term_eq** p_riemann ;
		/**
		* Array of pointers on various \c Term_eq.
		* Each one points onto one component of the Ricci tensor, in a given \c Domain.
		*/
		mutable Term_eq** p_ricci_tensor ;
		/**
		* Array of pointers on various \c Term_eq.
		* Each one points onto the Ricci scalar, in a given \c Domain.
		*/
		mutable Term_eq** p_ricci_scalar ;
		/**
		* Array of pointers on various \c Term_eq.
		* Each one points onto the potential of the Dirac gauge, in a given \c Domain.
		*/
		mutable Term_eq** p_dirac ;
		/**
		* Array of pointers on various \c Term_eq.
		* Each one points onto determinant of the covariant representation, in a given \c Domain.
		*/
		mutable Term_eq** p_det_cov ;
		int type_tensor ; ///< States if one works in the CON or COV representation.

	protected:
		Metric (const Space&) ; ///< Constructor from a \c Space only (internal use only).
		virtual void compute_con (int) const ; ///< Computes the contravariant representation, in a given \c Domain.
		virtual void compute_cov (int) const ; ///< Computes the covariant representation, in a given \c Domain.
		virtual void compute_christo (int) const ; ///< Computes the Christoffel symbols, in a given \c Domain.
		virtual void compute_riemann (int) const ; ///< Computes the Riemann tensor, in a given \c Domain.
		virtual void compute_ricci_tensor (int) const ; ///< Computes the Ricci tensor, in a given \c Domain.
		virtual void compute_ricci_scalar (int) const ; ///< Computes the Ricci scalar, in a given \c Domain.
		virtual void compute_dirac (int) const ; ///< Computes the Dirac gauge term, in a given \c Domain.
		virtual void compute_det_cov (int) const ; ///< Computes the determinant of the covariant representation, in a given \c Domain.

	public:
		Metric (const Metric&) ; ///< Copy constructor

		virtual ~Metric() ; ///< Destructor
		 /** Updates the derived quantities (Christoffels etc..)
		* This is done only for the ones that are needed, i.e. for the ones that have already been computed.
		*/
		virtual void update() ;
		/** Updates the derived quantities (Christoffels etc..) in a given \c Domain
		* This is done only for the ones that are needed, i.e. for the ones that have already been computed.
		* @param dd : the index of the \c Domain.
		*/
		virtual void update(int dd) ;

		/**
		* Uses the \c Metric to manipulate one of the index of a \c Term_eq (i.e. moves the index up or down)
		* @param so : the \c Term_eq to be manipulated.
		* @param ind : index number.
		*/
		virtual void manipulate_ind (Term_eq& so, int ind) const ;

		/**
		* Computes the partial derivative of a \c Term_eq (assumes Cartesian basis of decomposition).
		* The index corresponding to the derivation is given a name.
		* If need be inner summation is performed.
		* @param typeder : the result is either \f$\partial_i\f$ or \f$\partial^i\f$
		* @param nameder : the name given as the derivation index.
		* @param so : the \c Term_eq to be derived.
		* @returns : the derivative.
		*/
		virtual Term_eq derive_partial (int typeder, char nameder, const Term_eq& so) const ;
		/**
		* Computes the covariant derivative of a \c Term_eq (assumes Cartesian basis of decomposition).
		* The index corresponding to the derivation is given a name.
		* If need be inner summation is performed.
		* @param typeder : the result is either \f$D_i\f$ or \f$D^i\f$
		* @param nameder : the name given as the derivation index.
		* @param so : the \c Term_eq to be derived.
		* @returns : the covariant derivative.
		*/
		virtual Term_eq derive (int typeder, char nameder, const Term_eq& so) const ;

		/**
		* Computes the covariant flat derivative of a \c Term_eq.
		* The index corresponding to the derivation is given a name.
		* If need be inner summation is performed.
		* @param typeder : the result is either \f$\bar{D}_i\f$ or \f$\bar{D}^i\f$
		* @param nameder : the name given as the derivation index.
		* @param so : the \c Term_eq to be derived.
		* @returns : the covariant flat derivative.
		*/
		virtual Term_eq derive_flat (int typeder, char nameder, const Term_eq& so) const ;

		virtual const Metric* get_background() const ; ///< Returns a pointer on the background metric, if it exists.

		virtual int give_type (int) const ; ///< Returns the type of tensorial basis of the covariant representation, in a given \c Domain.

		// Extraction :
		/**
		* Gives one representation of the metric, in a \c Domain. This is the \c Term_eq version.
		* @param typemet : CON or COV.
		* @param dd : in which \c Domain.
		* @return : pointer on the result (being a \c Term_eq).
		*/
 		const Term_eq* give_term (int typemet, int dd) const ;
		/**
		* Gives the Christoffel symbols, in a \c Domain. This is the \c Term_eq version.
		* @param dd : in which \c Domain.
		* @return : pointer on the result (being a \c Term_eq).
		*/
		const Term_eq* give_christo (int dd) const ;
		/**
		* Gives the Riemann tensor, in a \c Domain. This is the \c Term_eq version.
		* @param dd : in which \c Domain.
		* @return : pointer on the result (being a \c Term_eq).
		*/
		const Term_eq* give_riemann (int dd) const ;
		/**
		* Gives the Ricci tensor, in a \c Domain. This is the \c Term_eq version.
		* @param dd : in which \c Domain.
		* @return : pointer on the result (being a \c Term_eq).
		*/
		const Term_eq* give_ricci_tensor (int dd) const ;
		/**
		* Gives the Ricci scalar, in a \c Domain. This is the \c Term_eq version.
		* @param dd : in which \c Domain.
		* @return : pointer on the result (being a \c Term_eq).
		*/
		const Term_eq* give_ricci_scalar (int dd) const ;
		/**
		* Gives the potential of the Dirac gauge, in a \c Domain. This is the \c Term_eq version.
		* @param dd : in which \c Domain.
		* @return : pointer on the result (being a \c Term_eq).
		*/
		const Term_eq* give_dirac (int dd) const ;
		/**
		* Gives the determinant of the covariant metric, in a \c Domain. This is the \c Term_eq version.
		* @param dd : in which \c Domain.
		* @return : pointer on the result (being a \c Term_eq).
		*/
		const Term_eq* give_det_cov (int dd) const ;

		friend class System_of_eqs ;
		friend class Metric_flat ;
		friend class Metric_dirac ;
		friend class Metric_conf ;
		friend class Domain_nucleus ;
		friend class Domain_shell ;
		friend class Domain_compact ;
		friend class Domain_bispheric_rect ;
		friend class Domain_bispheric_chi_first ;
		friend class Domain_bispheric_eta_first ;
		friend class Domain_shell_outer_adapted ;
		friend class Domain_shell_inner_adapted ;
		friend class Domain_nucleus_symphi ;
		friend class Domain_shell_symphi ;
		friend class Domain_compact_symphi ;
} ;

/**
* Class that deals with flat metric.
* \ingroup metric
*/
class Metric_flat : public Metric {

  protected:
    const Base_tensor& basis ; ///< The tensorial basis used.

	public:
		/**
		* Standard constructor
		* @param sp : the associated \c Space.
		* @param bb : the tensorial basis.
		*/
		Metric_flat (const Space& sp, const Base_tensor& bb) ;
		Metric_flat (const Metric_flat& ) ; ///< Copye constructor

	protected:
		virtual void compute_con (int) const ;
		virtual void compute_cov (int) const ;
		virtual void compute_christo (int) const ;
		virtual void compute_ricci_tensor (int) const ;
		virtual void compute_ricci_scalar (int) const ;

  private:
	/**
	* Computes the partial derivative part of the covariant derivative, in orthonormal spherical coordinates.
	* The result is \f$ \partial_r, \frac{1}{r}\partial_\theta, \frac{1}{r\sin\theta}\partial_\varphi\f$.
	* If need be inner summation of the result is performed and index manipulated.
	* @param tder : type of the result (COV or CON)
	* @param indder : name of the index corresponding to the derivative.
	* @param so : the field to be derived.
	* @return : the result as a \c Term_eq.
	*/
	Term_eq derive_partial_spher (int tder, char indder, const Term_eq& so) const ;
	/**
	* Computes the partial derivative part of the covariant derivative, in Cartesian coordinates.
	* If need be inner summation of the result is performed and index manipulated.
	* @param tder : type of the result (COV or CON)
	* @param indder : name of the index corresponding to the derivative.
	* @param so : the field to be derived.
	* @return : the result as a \c Term_eq.
	*/
	Term_eq derive_partial_cart (int tder, char indder, const Term_eq& so) const ;

	/**
	* Computes the partial derivative part of the covariant derivative, in orthonormal coordinates
	* where the constant radii sections have a negative curvature.
	* The result is \f$ \partial_r, \frac{\cos \theta}{r}\partial_\theta, \frac{\cos \theta}{r\sin\theta}\partial_\varphi\f$.
	* If need be inner summation of the result is performed and index manipulated.
	* @param tder : type of the result (COV or CON)
	* @param indder : name of the index corresponding to the derivative.
	* @param so : the field to be derived.
	* @return : the result as a \c Term_eq.
	*/
	Term_eq derive_partial_mtz (int tder, char indder, const Term_eq& so) const ;
	/**

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
	* Computes the flat covariant derivative, in Cartesian coordinates.
	* If need be inner summation of the result is performed and index manipulated.
	* @param tder : type of the result (COV or CON)
	* @param indder : name of the index corresponding to the derivative.
	* @param so : the field to be derived.
	* @return : the result as a \c Term_eq.
	*/
	Term_eq derive_cart (int, char, const Term_eq&) const ;

        /**
	* Computes the flat covariant derivative, in orthonormal coordinates
        * where the constant radii sections have a negative curvature.
	* If need be inner summation of the result is performed and index manipulated.
	* @param tder : type of the result (COV or CON)
	* @param indder : name of the index corresponding to the derivative.
	* @param so : the field to be derived.
	* @return : the result as a \c Term_eq.
	*/
	Term_eq derive_mtz (int tder, char indder, const Term_eq& so) const ;


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

	/**
	* Computes the flat covariant derivative, in orthonormal coordinates
        * where the constant radii sections have a negative curvature.
	* If need be inner summation of the result is performed.
	* If the CON derivative is asked for, the index is NOT manipulated by the flat metric but by an arbitrary one.
	* @param tder : type of the result (COV or CON)
	* @param indder : name of the index corresponding to the derivative.
	* @param so : the field to be derived.
	* @param othermet : pointer on the \c Metric used for index manipulation.
	* @return : the result as a \c Term_eq.
	*/
	Term_eq derive_with_other_mtz (int tder, char indder, const Term_eq& so, const Metric* othermet) const ;

	/**
	* Computes the flat covariant derivative, in Cartesian coordinates.
	* If need be inner summation of the result is performed.
	* If the CON derivative is asked for, the index is NOT manipulated by the flat metric but by an arbitrary one.
	* @param tder : type of the result (COV or CON)
	* @param indder : name of the index corresponding to the derivative.
	* @param so : the field to be derived.
	* @param othermet : pointer on the \c Metric used for index manipulation.
	* @return : the result as a \c Term_eq.
	*/
	Term_eq derive_with_other_cart (int tder, char indder, const Term_eq& so, const Metric* othermet) const ;

	public:
		virtual ~Metric_flat() ;
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
 * Class to deal with arbitrary type of metric.
 * \ingroup metric
 */
class Metric_general : public Metric {

	protected:
	  Metric_tensor* p_met ; ///< Pointer on the \c Metric_tensor describing the metric.
	  const Base_tensor& basis ; ///< The tensorial basis used.
	  Metric_flat fmet ; ///< Associated flat metric.
	  int place_syst ; ///< Gives the location of the metric amongst the various unknowns of the associated \c System_of_eqs.

	public:
		Metric_general (Metric_tensor&) ; ///< Constructor from a \c Metric_tensor.
		Metric_general (const Metric_general& ) ; ///< Copy constructor.

	protected:
		virtual void compute_con (int) const ;
		virtual void compute_cov (int) const ;
		virtual void compute_christo (int) const ;
		virtual void compute_riemann (int) const ;
		virtual void compute_ricci_tensor (int) const ;
		virtual void compute_dirac (int) const ;

	public:
		virtual ~Metric_general() ;
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
 * Class to deal with arbitrary type of metric but that is constant. By this one means that it is a given and is not modified by the \c System_of_eqs.
 * \ingroup metric
 */
class Metric_const : public Metric_general {

	public:
		Metric_const (Metric_tensor&) ; ///< Constructor from a \c Metric_tensor.
		Metric_const (const Metric_const& ) ; ///< Copy constructor

	protected:
		virtual void compute_con (int) const ;
		virtual void compute_cov (int) const ;


	public:
		virtual ~Metric_const() ;

		/**
		* Associate the metric to a given system of equations.
		* @param syst : the \c System_of_eqs.
		* @param name : name by which the metric will be known in the system (like "g", "f"...)
		*/
		virtual void set_system (System_of_eqs& syst, const char* name) ;
} ;

/**
 * Class to deal with a metric which determinant is 1.
 * It is however not enforced but only assumed is some computation (like the inverse).
 * \ingroup metric
 */
class Metric_conf : public Metric {

	protected:
	  Metric_tensor* p_met ; ///< Pointer on the \c Metric_tensor describing the metric.
	  const Base_tensor& basis ;  ///< The tensorial basis used.
	  Metric_flat fmet ;  ///< Associated flat metric.
	  int place_syst ; ///< Gives the location of the metric amongst the various unknowns of the associated \c System_of_eqs.



	public:
		Metric_conf (Metric_tensor&) ; ///< Constructor from a \c Metric_tensor.
		Metric_conf (const Metric_conf& ) ; ///< Copy constructor

	protected:
		virtual void compute_con (int) const ;
		virtual void compute_cov (int) const ;
		virtual void compute_christo (int) const ;
		virtual void compute_riemann (int) const ;
		virtual void compute_ricci_tensor (int) const ;
		virtual void compute_dirac (int) const ;

	public:
		virtual ~Metric_conf() ;
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
 * Class to deal with a conformal metric in the Dirac gauge.
 * The Dirac gauge is not enforced at this level but only assumed to be verified.
 * It has to be part of the equations solved for things to be consistent.
 * \ingroup metric
 */
 class Metric_dirac : public Metric_conf {

	public:
		Metric_dirac (Metric_tensor&) ; ///< Constructor from a \c Metric_tensor.
		Metric_dirac (const Metric_dirac& ) ;///< Copy constructor
		/**
		* @returns : the associated flat metric.
		*/
		const Metric_flat& get_fmet () const {return fmet;} ;

	protected:
		virtual void compute_ricci_tensor (int) const ;
		virtual void compute_ricci_scalar (int) const ;


	public:
		virtual ~Metric_dirac() ;
} ;

/**
 * Class to deal with a conformal metric in the Dirac gauge.
 * The metric is a constant, which means that it is a given and is not modified by the \c System_of_eqs.
 * The Dirac gauge is not enforced at this level but only assumed to be verified.
 * \ingroup metric
 */
class Metric_dirac_const : public Metric_dirac {

	public:
		Metric_dirac_const (Metric_tensor&) ; ///< Constructor from a \c Metric_tensor.
		Metric_dirac_const (const Metric_dirac_const& ) ; ///< Copy constructor

	protected:
		virtual void compute_con (int) const ;
		virtual void compute_cov (int) const ;

	public:
		virtual ~Metric_dirac_const() ;
		/**
		* Associate the metric to a given system of equations.
		* @param syst : the \c System_of_eqs.
		* @param name : name by which the metric will be known in the system (like "g", "f"...)
		*/
		virtual void set_system (System_of_eqs& syst, const char* name) ;
} ;
/**
 * Class to deal with a metric in the spatial harmonic gauge.
 * The gauge not enforced at this level but only assumed to be verified.
 * It has to be part of the equations solved for things to be consistent.
 * \ingroup metric
*/
 class Metric_harmonic : public Metric_general {

	public:
		Metric_harmonic (Metric_tensor&) ; ///< Constructor from a \c Metric_tensor.
		Metric_harmonic (const Metric_harmonic& ) ;///< Copy constructor
		/**
		* @returns : the associated flat metric.
		*/
		const Metric_flat& get_fmet () const {return fmet;} ;

	protected:
		virtual void compute_ricci_tensor (int) const ;
		virtual void compute_ricci_scalar (int) const ;
		virtual void compute_det_cov (int) const ;

	public:
		virtual ~Metric_harmonic() ;
} ;

/**
 * Class to deal with a metric with a conformal decomposition (mainly used for AADS spacetimes)
 * The true metric and the conformal are related via
 * \f$\gamma_{ij} = \frac{1}{\Omega^2}\tilde{\gamma}_{ij}\f$.
 * The conformal factor vanishes at some boundary so that the various quantities (Christoffels) are multiplied by appropriate factors
 * of $\Omega$ to ensure regularity.
 * \ingroup metric
 */
class Metric_conf_factor : public Metric {

	protected:
	  Metric_tensor* p_met ; ///< Pointer on the \c Metric_tensor describing the coformal metric.
	  const Base_tensor& basis ; ///< The tensorial basis used.
	  Metric_flat fmet ; ///< Associated flat metric.
          Scalar conformal ; ///< The conformal factor $\Omega$ (must be a purely radial function)
	  Vector grad_conf ; ///< flat gradient of the conformal factor
	  int place_syst ; ///< Gives the location of the metric amongst the various unknowns of the associated \c System_of_eqs.

	public:
		Metric_conf_factor (Metric_tensor&, const Scalar& conf) ; ///< Constructor from a \c Metric_tensor and a conformal factor.
		Metric_conf_factor (const Metric_conf_factor& ) ; ///< Copy constructor.

	protected:
		virtual void compute_con (int) const ;
		virtual void compute_cov (int) const ;
		virtual void compute_christo (int) const ;
		virtual void compute_riemann (int) const ;
		virtual void compute_ricci_tensor (int) const ;

	public:
		virtual ~Metric_conf_factor() ;
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
 * * Class to deal with a metric with a conformal decomposition (mainly used for AADS spacetimes)
 * The true metric and the conformal are related via
 * \f$\gamma_{ij} = \frac{1}{\Omega^2}\tilde{\gamma}_{ij}\f$.
 * The conformal factor vanishes at some boundary so that the various quantities (Christoffels) are multiplied by appropriate factors
 * of $\Omega$ to ensure regularity.
 * The metric is assumed to be constant.
 * \ingroup metric
 */
class Metric_conf_factor_const : public Metric_conf_factor {

	public:
		Metric_conf_factor_const (Metric_tensor&, const Scalar& conf) ; ///< Constructor from a \c Metric_tensor and a conformal factor.
		Metric_conf_factor_const (const Metric_conf_factor_const& ) ; ///< Copy constructor.

	protected:
		virtual void compute_con (int) const ;
		virtual void compute_cov (int) const ;


	public:
		virtual ~Metric_conf_factor_const() ;

		/**
		* Associate the metric to a given system of equations.
		* @param syst : the \c System_of_eqs.
		* @param name : name by which the metric will be known in the system (like "g", "f"...)
		*/
		virtual void set_system (System_of_eqs& syst, const char* name) ;
} ;

/**
 * Class to deal with a metric with a conformaly flat metric.
 * \ingroup metric
 */
class Metric_cfc : public Metric {

	protected:
	  Scalar* p_conf ; ///< Pointer on the \c Scalar being the conformal factor.
	  const Base_tensor& basis ;  ///< The tensorial basis used.
	  Metric_flat fmet ;  ///< Associated flat metric.
	  int place_syst ; ///< Gives the location of the metric amongst the various unknowns of the associated \c System_of_eqs.

	public:
		Metric_cfc (Scalar&, const Base_tensor&) ; ///< Constructor from a \c Metric_tensor.
		Metric_cfc (const Metric_cfc& ) ; ///< Copy constructor

	protected:
		virtual void compute_con (int) const ;
		virtual void compute_cov (int) const ;
		virtual void compute_christo (int) const ;
		virtual void compute_riemann (int) const ;
		virtual void compute_ricci_tensor (int) const ;

	public:
		virtual ~Metric_cfc() ;

		/**
		* Associate the metric to a given system of equations.
		* @param syst : the \c System_of_eqs.
		* @param name : name by which the metric will be known in the system (like "g", "f"...)
		*/
		virtual void set_system (System_of_eqs& syst, const char* name) ;

		virtual Term_eq derive (int, char, const Term_eq&) const ;
		virtual Term_eq derive_flat (int, char, const Term_eq&) const ;

		virtual int give_type (int) const ;
} ;


}
#endif
