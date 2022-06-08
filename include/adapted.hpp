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

#ifndef __ADAPTED_HPP_
#define __ADAPTED_HPP_

#include "space.hpp"
#include "term_eq.hpp"
#include "spheric.hpp"
#include "metric.hpp"
#include <vector>

namespace Kadath {

/**
* Class for a spherical-like domain, having a symmetry with respect to the plane \f$ z=0 \f$.
* \li 3 dimensions.
* \li centered on the point \c center \f$(X_c, Y_c, Z_c)\f$
*
* The outer surface is spherical whereas the inner boundary is deformed (typically to accomodate things like surfaces of stars).
*
* The link between numerical and physical coordinates is as follows :
*
* \f$ r = (R - r_i(\theta, \varphi))/2. x + (R + r_i(\theta, \varphi))/2. \f$
*
* \f$ \theta = \theta^\star \f$
*
* \f$ \varphi = \varphi^\star \f$
*
* \f$ R \f$ is the fixed outer boundary and \f$ r_i \f$ is an angular function giving the value of the inner variable boundary.
*
* \ingroup domain
*/


class Domain_shell_inner_adapted : public Domain {

 protected:
  /**
   * The corresponding \c Space ; required for updating fields whene the mapping changes.
   */
  const Space& sp ;
  Val_domain* inner_radius ; ///< Pointer on the inner boundary \f$r_i\f$, as a \c Val_domain
  mutable Term_eq* inner_radius_term_eq ; ///< Pointer on the inner boundary \f$r_i\f$, as a \c Term_eq
  double outer_radius ; ///< The outer radius \f$ R\f$.
  /**
   * Pointer on the \c Term_eq containing the normal vector to the inner boundary, in orthonormal spherical coordinates.
   */
  mutable Term_eq* normal_spher ;
 /**
   * Pointer on the \c Term_eq containing the normal vector to the inner boundary, in Cartesian coordinates.
   */
  mutable Term_eq* normal_cart ;
 /**
   * Pointer on the \c Term_eq containing the radius.
   */
  mutable Term_eq* rad_term_eq ;
  /**
   * Pointer on the \c Term_eq containing the \f$ {\rm d} r / {\rm d} x\f$.
   */
  mutable Term_eq* der_rad_term_eq ;
   /**
   * Pointer on the \c Term_eq containing the \f$ {\rm d} r / {\rm d} \theta\f$.
   */
  mutable Term_eq* dt_rad_term_eq ;
  /**
   * Pointer on the \c Term_eq containing the \f$ {\rm d} r / {\rm d} \varphi\f$.
   */
  mutable Term_eq* dp_rad_term_eq ;

  Point center ; ///< Absolute coordinates of the center.

 public:
  /**
  * Constructor :
  * @param sp [input] : the associated \c Space.
  * @param num : number of the domain (used by the \c Space).
  * @param ttype [input] : Chebyshev or Legendre type of spectral expansion.
  * @param rin [input] : inner radius (constant with this constructor).
  * @param rout [input] : outer radius.
  * @param cr [input] : center of the spherical coordinates.
  * @param nbr [nbr] : number of points in each dimension.
  */
  Domain_shell_inner_adapted (const Space& sp, int num, int ttype, double rin, double rout, const Point& cr, const Dim_array& nbr) ;
 /**
  * Constructor :
  * @param sp [input] : the associated \c Space.
  * @param num : number of the domain (used by the \c Space).
  * @param ttype [input] : Chebyshev or Legendre type of spectral expansion.
  * @param rin [input] : inner radius.
  * @param rout [input] : outer radius.
  * @param cr [input] : center of the spherical coordinates.
  * @param nbr [nbr] : number of points in each dimension.
  */
  Domain_shell_inner_adapted (const Space& sp, int num, int ttype, const Val_domain& rin, double rout, const Point& cr, const Dim_array& nbr) ;
  Domain_shell_inner_adapted (const Domain_shell_inner_adapted & so) ; ///< Copy constructor.
 /**
  * Constructor from a file
  * @param sp [input] : the associated \c Space.
  * @param num : number of the domain (used by the \c Space).
  * @param fd : the fiel (generated by the save function.
  */
  Domain_shell_inner_adapted (const Space& sp, int num, FILE* fd) ;

  virtual ~Domain_shell_inner_adapted () ; ///< Destructor
  virtual void del_deriv() const ;
  virtual void save (FILE*) const ;
  virtual double integ (const Val_domain&, int) const ;
   /**
    * Returns the inner variable boundary.
    */
  Val_domain get_inner_radius() const {return *inner_radius;} ;

  private:
    virtual void do_absol ()  const ;
    virtual void do_radius () const ;
    virtual void do_cart ()  const ;
    virtual void do_cart_surr () const ;


  protected:

     virtual void set_cheb_base(Base_spectral&) const ;
     virtual void set_legendre_base(Base_spectral&) const ;

     virtual void set_anti_cheb_base(Base_spectral&) const ;
     virtual void set_anti_legendre_base(Base_spectral&) const ;

     virtual Tensor change_basis_cart_to_spher (int, const Tensor&) const ;
     virtual Tensor change_basis_spher_to_cart (int, const Tensor&) const ;

     virtual void set_cheb_base_r_spher(Base_spectral&) const ;
     virtual void set_cheb_base_t_spher(Base_spectral&) const ;
     virtual void set_cheb_base_p_spher(Base_spectral&) const ;
     virtual void set_legendre_base_r_spher(Base_spectral&) const ;
     virtual void set_legendre_base_t_spher(Base_spectral&) const ;
     virtual void set_legendre_base_p_spher(Base_spectral&) const ;

     virtual void set_cheb_r_base(Base_spectral&) const ;
     virtual void set_legendre_r_base(Base_spectral&) const ;

     virtual void do_coloc () ;
     virtual int give_place_var (char*) const ;

     virtual int nbr_unknowns_from_adapted() const ;
     virtual void vars_to_terms() const ;
     virtual void affecte_coef (int&, int, bool&) const ;
     virtual void xx_to_vars_from_adapted (Val_domain&, const Array<double>&, int&) const ;
     virtual void xx_to_ders_from_adapted (const Array<double>&, int&) const ;
     virtual void update_term_eq (Term_eq*) const ;
     virtual void update_variable (const Val_domain&, const Scalar&, Scalar&) const ;
     virtual void update_constante (const Val_domain&, const Scalar&, Scalar&) const ;
     virtual void update_mapping(const Val_domain&) ;

   public:
    /** Affects the inner radius.
     * @param so : the value to be used.
     */
     void set_mapping(const Val_domain& so) const ;
    /**
     * Updates all the quantities that depend on the inner radius (like the normal vectors).
     */
     void update() const ;

  public:
     virtual Point get_center () const {return center ;} ;
     virtual bool is_in(const Point&xx, double prec=1e-13) const ;
     virtual const Point absol_to_num(const Point&) const;

     virtual const Point absol_to_num_bound(const Point&, int) const;


     virtual void do_der_abs_from_der_var(Val_domain** der_var, Val_domain** der_abs) const ;
     virtual Base_spectral mult (const Base_spectral&, const Base_spectral&) const ;

  public:

     virtual Val_domain mult_cos_phi (const Val_domain&) const ;
     virtual Val_domain mult_sin_phi (const Val_domain&) const ;
     virtual Val_domain mult_cos_theta (const Val_domain&) const ;
     virtual Val_domain mult_sin_theta (const Val_domain&) const ;
     virtual Val_domain div_sin_theta (const Val_domain&) const ;
     virtual Val_domain div_cos_theta (const Val_domain&) const ;
     virtual Val_domain ddp (const Val_domain&) const ;
     virtual Val_domain der_r (const Val_domain&) const ;
     virtual Val_domain div_r (const Val_domain&) const ;
     virtual Val_domain laplacian (const Val_domain&, int) const ;
     virtual Val_domain laplacian2 (const Val_domain&, int) const ;

     virtual double val_boundary (int, const Val_domain&, const Index&) const ;
     virtual void find_other_dom (int, int, int&, int&) const ;
     virtual Val_domain der_normal (const Val_domain&, int) const ;

     virtual int nbr_unknowns (const Tensor&, int) const ;
	/**
	* Computes the number of true unknowns of a \c Val_domain.
	* It takes into account the various symmetries and regularity conditions to determine the precise number of degrees of freedom.
	* @param so : the field.
	* @param mlim: limit for the regularity (quantum number wrt \f$\varphi\f$).
	* @returns the number of true unknowns.
	*/
     int nbr_unknowns_val_domain (const Val_domain& so, int mlim) const ;
     virtual Array<int> nbr_conditions (const Tensor&, int, int, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Computes number of discretized equations associated with a given tensorial equation in the bulk.
	* It takes into account the various Galerkin basis used.
	* @param so : the residual of the equation.
	* @param mlim : limit for the regularity (quantum number wrt \f$\varphi\f$).
	* @param order : order of the equation (i.e. 2 for a Laplacian for instance)
	* @returns the number of true unknowns.
	*/
     int nbr_conditions_val_domain (const Val_domain& so, int mlim, int order) const ;
     virtual Array<int> nbr_conditions_boundary (const Tensor&, int, int, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Computes number of discretized equations associated with a given equation on a boundary.
	* It takes into account the various Galerkin basis used.
	* It is used for implementing boundary conditions and matching ones.
	* @param eq : the residual of the equation.
	* @param mlim : mimit for the regularity (quantum number wrt \f$\varphi\f$).
	* @returns the number of true conditions.
	*/
     int nbr_conditions_val_domain_boundary (const Val_domain& eq, int mlim) const ;
     virtual void export_tau (const Tensor&, int, int, Array<double>&, int&, const Array<int>&,  int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Exports a residual equation in the bulk.
	* It makes use of the various Galerkin basis used.
	* @param eq : the residual of the equation.
	* @param mlim : limit for the regularity (quantum number wrt \f$\varphi\f$).
	* @param order : describes the order of the equation (2 for a Laplacian for instance).
	* @param res : The \c Array where the discretized equations are stored.
	* @param pos_res : current position in res.
	* @param ncond :  the corresponding number of equations. It is used when the equation is null.
	*/
     void export_tau_val_domain (const Val_domain& eq, int mlim, int order, Array<double>& res, int& pos_res, int ncond) const ;
     virtual void export_tau_boundary (const Tensor&, int, int, Array<double>&, int&, const Array<int>&, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Exports all the residual equations corresponding to a tensorial one on a given boundary
	* It makes use of the various Galerkin basis used.
	* @param eq : the residual of the equation.
	* @param mlim : limit for the regularity (quantum number wrt \f$\varphi\f$).
	* @param bound : the boundary at which the equation is enforced.
	* @param res : The \c Array where the discretized equations are stored.
	* @param pos_res : current position in res.
	* @param ncond : the corresponding number of equations. It is used when the residual is null.
	*/
     void export_tau_val_domain_boundary (const Val_domain& eq, int mlim, int bound, Array<double>& res, int& pos_res, int ncond) const ;
     virtual void affecte_tau (Tensor&, int, const Array<double>&, int&) const ;
	/**
	* Affects some coefficients to a \c Val_domain.
	* It takes into account the various symmetries and regularity conditions (by means of Garlekin basis).
	* @param so : the field to be affected.
	* @param mlim : limit for the regularity (quantum number wrt \f$\varphi\f$).
	* @param cf : \c Array of the coefficients used.
	* @param pos_cf : current position in the array of coefficients.
	*/
     void affecte_tau_val_domain (Val_domain& so, int mlim, const Array<double>& cf, int& pos_cf) const ;
     virtual void affecte_tau_one_coef (Tensor&, int, int, int&) const ;
	/**
	* Sets at most one coefficient of a \c Val_domain to 1.
	* It takes into account the various symmetries and regularity conditions (by means of Garlekin basis).
	* @param so : the \c Val_domain to be affected. It is set to zero if cc does not corresponds to another field.
	* @param mlim : limit for the regularity (quantum number wrt \f$\varphi\f$).
	* @param cc : location, in the overall system, of the coefficient to be set to 1.
	* @param pos_cf : current position.
	*/
     void affecte_tau_one_coef_val_domain (Val_domain& so, int mlim, int cc, int& pos_cf) const ;

     virtual int nbr_points_boundary (int, const Base_spectral&) const ;
     virtual void do_which_points_boundary (int, const Base_spectral&, Index**, int) const ;

     /**
      * Computes the flat gradient of a field, in orthonormal spherical coordinates.
      * @param so : the input field.
      * @returns : the gradient.
      */
     Term_eq flat_grad_spher (const Term_eq&) const ;

     virtual Term_eq partial_spher (const Term_eq&) const ;
     virtual Term_eq partial_cart (const Term_eq&) const ;
     virtual Term_eq connection_spher (const Term_eq&) const ;
     virtual const Term_eq* give_normal(int, int) const ;

	/**
	* Computes \f$ \partial_r \f$.
	* @param so : the source \c Term_eq.
	* @returns : the derivative.
	*/
     Term_eq derive_r (const Term_eq& so) const ;
	/**
	* Computes \f$ \partial_\theta \f$.
	* @param so : the source \c Term_eq.
	* @returns : the derivative.
	*/
     Term_eq derive_t (const Term_eq& so) const ;
	/**
	* Computes \f$ \partial_\varphi \f$.
	* @param so : the source \c Term_eq.
	* @returns : the derivative.
	*/
     Term_eq derive_p (const Term_eq& so) const ;
	/** Computes the normal wrt the inner boundary, in orthonormal spherical coordinates.
	* The result is stored in \c normal_spher
	*/
     void do_normal_spher () const ;
	/** Computes the normal wrt the inner boundary, in Cartesian coordinates.
	* The result is stored in \c normal_cart
	*/
     void do_normal_cart () const ;

     virtual Term_eq der_normal_term_eq (const Term_eq&, int) const ;
     virtual Term_eq dr_term_eq (const Term_eq&) const ;
     virtual Term_eq lap_term_eq (const Term_eq&, int) const ;
     virtual Term_eq mult_r_term_eq (const Term_eq&) const ;
     virtual Term_eq integ_volume_term_eq (const Term_eq&) const ;

     virtual Term_eq derive_flat_spher (int, char, const Term_eq&, const Metric*) const ;
     virtual Term_eq derive_flat_cart (int, char, const Term_eq&, const Metric*) const ;

     virtual Tensor import (int, int, int, const Array<int>&,  Tensor**) const ;

     virtual double integ_volume (const Val_domain&) const ;

public:
     virtual ostream& print (ostream& o) const ;

     friend class Space_spheric_adapted ;
     friend class Space_bin_ns ;
     friend class Space_bhns ;
     friend class Space_bin_bh ;
     friend class Space_adapted_bh ;
     friend class Space_KerrSchild_bh ;
     friend class Space_bbh ;
     friend class Space_Kerr_bbh;
} ;

/**
* Class for a spherical-like domain, having a symmetry with respect to the plane \f$ z=0 \f$.
* \li 3 dimensions.
* \li centered on the point \c center \f$(X_c, Y_c, Z_c)\f$
*
* The inner surface is spherical whereas the outer boundary is deformed (typically to accomodate things like surfaces of stars).
*
* The link between numerical and physical coordinates is as follows :
*
* \f$ r = (r_o(\theta, \varphi) - R)/2. x + (r_o(\theta, \varphi) + R)/2. \f$
*
* \f$ \theta = \theta^\star \f$
*
* \f$ \varphi = \varphi^\star \f$
*
* \f$ R \f$ is the fixed inner boundary and \f$ r_o \f$ is an angular function giving the value of the outer variable boundary.
*
* \ingroup domain
*/
class Domain_shell_outer_adapted : public Domain {

 protected:
 /**
   * The corresponding \c Space ; required for updating fields whene the mapping changes.
   */
  const Space& sp ;
  Val_domain* outer_radius ; ///< Pointer on the outer boundary \f$r_o\f$, as a \c Val_domain
  mutable Term_eq* outer_radius_term_eq ; ///< Pointer on the outer boundary \f$r_i\f$, as a \c Term_eq
  double inner_radius ; ///< The inner radius \f$ R\f$.
 /**
   * Pointer on the \c Term_eq containing the normal vector to the outer boundary, in orthonormal spherical coordinates.
   */
  mutable Term_eq* normal_spher ;
 /**
   * Pointer on the \c Term_eq containing the normal vector to the outer boundary, in Cartesian coordinates.
   */
  mutable Term_eq* normal_cart ;
/**
   * Pointer on the \c Term_eq containing the radius.
   */
  mutable Term_eq* rad_term_eq ;
 /**
   * Pointer on the \c Term_eq containing the \f$ {\rm d} r / {\rm d} x\f$.
   */
  mutable Term_eq* der_rad_term_eq ;
  /**
   * Pointer on the \c Term_eq containing the \f$ {\rm d} r / {\rm d} \theta\f$.
   */
  mutable Term_eq* dt_rad_term_eq ;
  /**
   * Pointer on the \c Term_eq containing the \f$ {\rm d} r / {\rm d} \varphi\f$.
   */
  mutable Term_eq* dp_rad_term_eq ;

  Point center ; ///< Absolute coordinates of the center.

 public:
  /**
  * Constructor :
  * @param sp [input] : the associated \c Space.
  * @param num : number of the domain (used by the \c Space).
  * @param ttype [input] : Chebyshev or Legendre type of spectral expansion.
  * @param rin [input] : inner radius.
  * @param rout [input] : outer radius (constant with this constructor).
  * @param cr [input] : center of the spherical coordinates.
  * @param nbr [nbr] : number of points in each dimension.
  */
  Domain_shell_outer_adapted (const Space& sp, int num, int ttype, double rin, double rout, const Point& cr, const Dim_array& nbr) ;
  /**
  * Constructor :
  * @param sp [input] : the associated \c Space.
  * @param num : number of the domain (used by the \c Space).
  * @param ttype [input] : Chebyshev or Legendre type of spectral expansion.
  * @param rin [input] : inner radius.
  * @param rout [input] : outer radius.
  * @param cr [input] : center of the spherical coordinates.
  * @param nbr [nbr] : number of points in each dimension.
  */
  Domain_shell_outer_adapted (const Space& sp, int num, int ttype, double rin, const Val_domain& rout, const Point& cr, const Dim_array& nbr) ;
  Domain_shell_outer_adapted (const Domain_shell_outer_adapted & so) ; ///< Copy constructor.
 /**
  * Constructor from a file
  * @param sp [input] : the associated \c Space.
  * @param num : number of the domain (used by the \c Space).
  * @param fd : the file (generated by the save function.
  */
  Domain_shell_outer_adapted (const Space& sp, int num, FILE* fr) ;

  virtual ~Domain_shell_outer_adapted () ;
  virtual void del_deriv() const ;
  virtual void save (FILE*) const ;
   /**
    * Returns the outer variable boundary.
    */
  Val_domain get_outer_radius() const {return *outer_radius;} ;

  private:
    virtual void do_absol ()  const ;
    virtual void do_radius () const ;
    virtual void do_cart ()  const ;
    virtual void do_cart_surr () const ;


  protected:

     virtual void set_cheb_base(Base_spectral&) const ;
     virtual void set_legendre_base(Base_spectral&) const ;

     virtual void set_anti_cheb_base(Base_spectral&) const ;
     virtual void set_anti_legendre_base(Base_spectral&) const ;

     virtual void set_cheb_base_r_spher(Base_spectral&) const ;
     virtual void set_cheb_base_t_spher(Base_spectral&) const ;
     virtual void set_cheb_base_p_spher(Base_spectral&) const ;
     virtual void set_legendre_base_r_spher(Base_spectral&) const ;
     virtual void set_legendre_base_t_spher(Base_spectral&) const ;
     virtual void set_legendre_base_p_spher(Base_spectral&) const ;

     virtual void set_cheb_r_base(Base_spectral&) const ;
     virtual void set_legendre_r_base(Base_spectral&) const ;

     virtual void do_coloc () ;
     virtual int give_place_var (char*) const ;

     virtual int nbr_unknowns_from_adapted() const ;
     virtual void vars_to_terms() const ;
     virtual void affecte_coef (int&, int, bool&) const ;
     virtual void xx_to_vars_from_adapted (Val_domain&, const Array<double>&, int&) const ;
     virtual void xx_to_ders_from_adapted (const Array<double>&, int&) const ;
     virtual void update_term_eq (Term_eq*) const ;
     virtual void update_variable (const Val_domain&, const Scalar&, Scalar&) const ;
     virtual void update_constante (const Val_domain&, const Scalar&, Scalar&) const ;
     virtual void update_mapping(const Val_domain&) ;
  public:
    /** Affects the outer radius.
     * @param so : the value to be used.
     */
     void set_mapping(const Val_domain& so) const ;
    /**
     * Updates all the quantities that depend on the inner radius (like the normal vectors).
     */
     void update() const ;

  public:
     virtual Point get_center () const {return center ;} ;

     virtual bool is_in(const Point&xx, double prec=1e-13) const ;
     virtual const Point absol_to_num(const Point&) const;
     virtual const Point absol_to_num_bound(const Point&, int) const;
     virtual void do_der_abs_from_der_var(Val_domain** der_var, Val_domain** der_abs) const ;
     virtual Base_spectral mult (const Base_spectral&, const Base_spectral&) const ;

  public:
     virtual Val_domain mult_cos_phi (const Val_domain&) const ;
     virtual Val_domain mult_sin_phi (const Val_domain&) const ;
     virtual Val_domain mult_cos_theta (const Val_domain&) const ;
     virtual Val_domain mult_sin_theta (const Val_domain&) const ;
     virtual Val_domain div_sin_theta (const Val_domain&) const ;
     virtual Val_domain div_cos_theta (const Val_domain&) const ;
     virtual Val_domain laplacian (const Val_domain&, int) const ;
     virtual Val_domain laplacian2 (const Val_domain&, int) const ;

     virtual Tensor change_basis_cart_to_spher (int dd, const Tensor&) const ;
     virtual Tensor change_basis_spher_to_cart (int dd, const Tensor&) const ;

     virtual Val_domain ddp (const Val_domain&) const ;
     virtual Val_domain der_r (const Val_domain&) const ;
     virtual Val_domain div_r (const Val_domain&) const ;

     virtual double val_boundary (int, const Val_domain&, const Index&) const ;
     virtual void find_other_dom (int, int, int&, int&) const ;
     virtual Val_domain der_normal (const Val_domain&, int) const ;

     virtual int nbr_unknowns (const Tensor&, int) const ;
	/**
	* Computes the number of true unknowns of a \c Val_domain.
	* It takes into account the various symmetries and regularity conditions to determine the precise number of degrees of freedom.
	* @param so : the field.
	* @param mlim: limit for the regularity (quantum number wrt \f$\varphi\f$).
	* @returns the number of true unknowns.
	*/
     int nbr_unknowns_val_domain (const Val_domain& so, int mlim) const ;
     virtual Array<int> nbr_conditions (const Tensor&, int, int, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Computes number of discretized equations associated with a given tensorial equation in the bulk.
	* It takes into account the various Galerkin basis used.
	* @param so : the residual of the equation.
	* @param mlim : limit for the regularity (quantum number wrt \f$\varphi\f$).
	* @param order : order of the equation (i.e. 2 for a Laplacian for instance)
	* @returns the number of true unknowns.
	*/
     int nbr_conditions_val_domain (const Val_domain& so, int mlim, int order) const ;
     virtual Array<int> nbr_conditions_boundary (const Tensor&, int, int, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Computes number of discretized equations associated with a given equation on a boundary.
	* It takes into account the various Galerkin basis used.
	* It is used for implementing boundary conditions and matching ones.
	* @param eq : the residual of the equation.
	* @param mlim : mimit for the regularity (quantum number wrt \f$\varphi\f$).
	* @returns the number of true conditions.
	*/
     int nbr_conditions_val_domain_boundary (const Val_domain& eq, int mlim) const ;
     virtual void export_tau (const Tensor&, int, int, Array<double>&, int&, const Array<int>&,  int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Exports a residual equation in the bulk.
	* It makes use of the various Galerkin basis used.
	* @param eq : the residual of the equation.
	* @param mlim : limit for the regularity (quantum number wrt \f$\varphi\f$).
	* @param order : describes the order of the equation (2 for a Laplacian for instance).
	* @param res : The \c Array where the discretized equations are stored.
	* @param pos_res : current position in res.
	* @param ncond :  the corresponding number of equations. It is used when the equation is null.
	*/
     void export_tau_val_domain (const Val_domain& eq, int mlim, int order, Array<double>& res, int& pos_res, int ncond) const ;
     virtual void export_tau_boundary (const Tensor&, int, int, Array<double>&, int&, const Array<int>&, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Exports all the residual equations corresponding to a tensorial one on a given boundary
	* It makes use of the various Galerkin basis used.
	* @param eq : the residual of the equation.
	* @param mlim : limit for the regularity (quantum number wrt \f$\varphi\f$).
	* @param bound : the boundary at which the equation is enforced.
	* @param res : The \c Array where the discretized equations are stored.
	* @param pos_res : current position in res.
	* @param ncond : the corresponding number of equations. It is used when the residual is null.
	*/
     void export_tau_val_domain_boundary (const Val_domain& eq, int mlim, int bound, Array<double>& res, int& pos_res, int ncond) const ;
     virtual void affecte_tau (Tensor&, int, const Array<double>&, int&) const ;
	/**
	* Affects some coefficients to a \c Val_domain.
	* It takes into account the various symmetries and regularity conditions (by means of Garlekin basis).
	* @param so : the field to be affected.
	* @param mlim : limit for the regularity (quantum number wrt \f$\varphi\f$).
	* @param cf : \c Array of the coefficients used.
	* @param pos_cf : current position in the array of coefficients.
	*/
     void affecte_tau_val_domain (Val_domain& so, int mlim, const Array<double>& cf, int& pos_cf) const ;
     virtual void affecte_tau_one_coef (Tensor&, int, int, int&) const ;
	/**
	* Affects some coefficients to a \c Val_domain.
	* It takes into account the various symmetries and regularity conditions (by means of Garlekin basis).
	* @param so : the field to be affected.
	* @param mlim : limit for the regularity (quantum number wrt \f$\varphi\f$).
	* @param cf : \c Array of the coefficients used.
	* @param pos_cf : current position in the array of coefficients.
	*/
     void affecte_tau_one_coef_val_domain (Val_domain& so, int mlim, int cf, int& pos_cf) const ;

     virtual int nbr_points_boundary (int, const Base_spectral&) const ;
     virtual void do_which_points_boundary (int, const Base_spectral&, Index**, int) const ;

     /**
      * Computes the flat gradient of a field, in orthonormal spherical coordinates.
      * @param so : the input field.
      * @returns : the gradient.
      */
     Term_eq flat_grad_spher (const Term_eq&) const ;
     virtual Term_eq partial_spher (const Term_eq&) const ;
     virtual Term_eq partial_cart (const Term_eq&) const ;
     virtual Term_eq connection_spher (const Term_eq&) const ;
     virtual const Term_eq* give_normal(int, int) const ;

	/**
	* Computes \f$ \partial_r \f$.
	* @param so : the source \c Term_eq.
	* @returns : the derivative.
	*/
     Term_eq  derive_r (const Term_eq& so) const ;
	/**
	* Computes \f$ \partial_\theta \f$.
	* @param so : the source \c Term_eq.
	* @returns : the derivative.
	*/
     Term_eq derive_t (const Term_eq& so) const ;
	/**
	* Computes \f$ \partial_\varphi \f$.
	* @param so : the source \c Term_eq.
	* @returns : the derivative.
	*/
     Term_eq derive_p (const Term_eq& so) const ;
	/** Computes the normal wrt the inner boundary, in orthonormal spherical coordinates.
	* The result is stored in \c normal_spher
	*/
    void do_normal_spher () const ;
	/** Computes the normal wrt the inner boundary, in Cartesian coordinates.
	* The result is stored in \c normal_cart
	*/
       void do_normal_cart () const ;

     virtual Term_eq der_normal_term_eq (const Term_eq&, int) const ;
     virtual Term_eq dr_term_eq (const Term_eq&) const ;
     virtual Term_eq lap_term_eq (const Term_eq&, int) const ;
     virtual Term_eq mult_r_term_eq (const Term_eq&) const ;
     virtual Term_eq integ_volume_term_eq (const Term_eq&) const ;
     virtual Term_eq integ_term_eq (const Term_eq&, int) const ;
     virtual Term_eq derive_flat_spher (int, char, const Term_eq&, const Metric*) const ;
     virtual Term_eq derive_flat_cart (int, char, const Term_eq&, const Metric*) const ;

     virtual double integ_volume (const Val_domain&) const ;
     virtual Tensor import (int, int, int, const Array<int>&,  Tensor**) const ;
 
public:
     virtual ostream& print (ostream& o) const ;
     
     friend class Space_spheric_adapted ;
     friend class Space_bin_ns ;
     friend class Space_bhns ;
     friend class Space_bin_bh ;
     friend class Space_adapted_bh ;     
     friend class Space_bbh ;
     friend class Space_KerrSchild_bh ;
     friend class Space_Kerr_bbh ;
} ;

/**
 * The \c Space_spheric_adapted class fills the space with one nucleus, one shell adapted on the outside, one shell adapted on the inside, several standard shells and a compactified domain, all centered on the same point.
 * \ingroup domain
 */
class Space_spheric_adapted : public Space {
     public:
	/**
     	* Standard constructor ; all the shells are initially not deformed.
     	* @param ttype [input] : the type of basis.
	* @param cr [input] : absolute coordinates of the center.
	* @param nbr [input] : number of points in each domain.
	* @param bounds [input] : radii of the various shells (and also determines the total number of domains).
	*/
	Space_spheric_adapted (int ttype, const Point& cr, const Dim_array& nbr, const Array<double>& bounds) ;
     Space_spheric_adapted (int ttype, const Point& cr, const Dim_array& nbr, const std::vector<double>& bounds) ;
	Space_spheric_adapted (FILE*) ; ///< Constructor from a file
	virtual ~Space_spheric_adapted() ; ///< Destructor
	virtual void save(FILE*) const ;

	virtual int nbr_unknowns_from_variable_domains() const ;
	virtual void affecte_coef_to_variable_domains(int& , int, Array<int>&) const ;
	virtual void xx_to_ders_variable_domains(const Array<double>&, int&) const ;
	virtual void xx_to_vars_variable_domains(System_of_eqs*, const Array<double>&, int&) const ;

	/**
	* Adds an equation being the value of some field at the origin.
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the quantity that must be zero at the origin
	*/
	void add_eq_ori (System_of_eqs& syst, const char* eq) ;
	/**
	* Adds a bulk equation and two matching conditions.
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the bulk equation.
	* @param rac : the string describing the first matching condition.
	* @param rac_der : the string describing the second matching condition.
	* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
	* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
	*/
	void add_eq (System_of_eqs& syst, const char* eq, const char* rac, const char* rac_der, int nused=-1, Array<int>** pused=0x0)  ;
	/**
	* Adds an equation being a surface integral at infinity.
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the equation (should contain something like integ(f)=b)
	*/
	void add_eq_int_inf (System_of_eqs& syst, const char* eq) ;
  
	/**
	* Adds an equation being a surface integral at an arbitrary domain and boundary so long as it is defined.
	* @param syst : the \c System_of_eqs.
  * @param dom : domain to apply equation
  * @param bc : boundary to integrate at
	* @param eq : the string describing the equation (should contain something like integ(f)=b)
	*/
  void add_eq_int (System_of_eqs& sys, const int dom, const int bc, const char* eq);
  
  /**
	* Adds an equation being a volume integral in the domains below a given number.
	* @param syst : the \c System_of_eqs.
	* @param nz : the integral is taken for all the \c Domains which number is \f$ <= nz \f$.
	* @param eq : the string describing the equation (should contain something like integvolume(f)=b)
	*/
	void add_eq_int_volume (System_of_eqs& syst, int nz, const char* eq) ;

	virtual Array<int> get_indices_matching_non_std(int, int) const ;
} ;
}
#endif
