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

#ifndef __ONED_HPP_
#define __ONED_HPP_

#include "space.hpp"
#include "val_domain.hpp" 

namespace Kadath {
/**
* Class for a 1-dimensional spherical domain containing the origin.
* It is intended to deal with systems with spherical symmetry.
* \li The numerical coordinate is :
*
* \f$ 0 \leq x \leq 1 \f$
*
* \li Standard spherical coordinates :
*
* \f$ r = \alpha x \f$
*
* \ingroup domain
*/
class Domain_oned_ori : public Domain {

 private:
  double alpha ; ///< Relates the numerical radius to the physical one.
  
 public: 
   /**
  * Standard constructor :
  * @param num : number of the domain (used by the \c Space).
  * @param ttype : Chebyshev or Legendre type of spectral expansion.
  * @param radius : radius of the nucleus.
  * @param nbr : number of points in each dimension.
  */
  Domain_oned_ori (int num, int ttype, double radius, const Dim_array& nbr) ;
  Domain_oned_ori (const Domain_oned_ori& so) ; ///< Copy constructor.
 /**
  * Constructor from a file
  * @param num : number of the domain (used by the \c Space).
  * @param ff: file containd the domain, generated by the save function.
  */
  Domain_oned_ori (int num, FILE* ff) ;

  virtual ~Domain_oned_ori() ; ///< Destructor
  virtual void save (FILE*) const ;

  private:    
    virtual void do_absol ()  const ;

  private:
     virtual void set_cheb_base(Base_spectral&) const ;       
     virtual void set_legendre_base(Base_spectral&) const ;
     virtual void set_cheb_base_odd(Base_spectral&) const ;       
     virtual void set_legendre_base_odd(Base_spectral&) const ;
     virtual void do_coloc () ;
     virtual int give_place_var (char*) const ;
     virtual void do_radius () const ;

  public:
     virtual bool is_in(const Point&xx, double prec=1e-13) const ;
     virtual const Point absol_to_num(const Point&) const;
     virtual void do_der_abs_from_der_var(Val_domain** der_var, Val_domain** der_abs) const ;
     virtual Base_spectral mult (const Base_spectral&, const Base_spectral&) const ;
   
     virtual double val_boundary (int, const Val_domain&, const Index&) const ;
     virtual void find_other_dom (int, int, int&, int&) const ;
     virtual Val_domain der_normal (const Val_domain&, int) const ;
     virtual Val_domain der_partial_var (const Val_domain&, int) const ;
   virtual Val_domain div_r (const Val_domain&) const ;
     virtual int nbr_unknowns (const Tensor&, int) const ;
	/**
	* Computes the number of true unknowns of a \c Val_domain.
	* It takes into account the various symmetries and regularity conditions to determine the precise number of degrees of freedom.
	* @param so : the field.
	* @returns the number of true unknowns.
	*/
     int nbr_unknowns_val_domain (const Val_domain& so) const ;
     virtual Array<int> nbr_conditions (const Tensor&, int, int, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;

	/**
	* Computes number of discretized equations associated with a given tensorial equation in the bulk.
	* It takes into account the various Galerkin basis used.
	* @param so : the residual of the equation.
	* @param order : order of the equation (i.e. 2 for a Laplacian for instance)
	* @returns the number of true unknowns.
	*/
     int nbr_conditions_val_domain (const Val_domain& so, int order) const ;
     virtual Array<int> nbr_conditions_boundary (const Tensor&, int, int, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;

	/**
	* Computes number of discretized equations associated with a given equation on a boundary.
	* It takes into account the various Galerkin basis used.
	* It is used for implementing boundary conditions and matching ones.
	* @param eq : the residual of the equation.
	* @returns the number of true conditions.
	*/
     int nbr_conditions_val_domain_boundary (const Val_domain& eq) const ;
     virtual void export_tau (const Tensor&, int, int, Array<double>&, int&, const Array<int>&,int n_cmp=-1, Array<int>** p_cmp=0x0) const ;

	/**
	* Exports a residual equation in the bulk.
	* It makes use of the various Galerkin basis used.
	* @param eq : the residual of the equation.
	* @param order : describes the order of the equation (2 for a Laplacian for instance).
	* @param res : The \c Array where the discretized equations are stored.
	* @param pos_res : current position in res.
	* @param ncond :  the corresponding number of equations. It is used when the equation is null.
	*/
     void export_tau_val_domain (const Val_domain& eq, int order, Array<double>& res, int& pos_res, int ncond) const ;
     virtual void export_tau_boundary (const Tensor&, int, int, Array<double>&, int&, const Array<int>&, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;

	/**
	* Exports all the residual equations corresponding to a tensorial one on a given boundary
	* It makes use of the various Galerkin basis used.
	* @param eq : the residual of the equation.
	* @param bound : the boundary at which the equation is enforced.
	* @param res : The \c Array where the discretized equations are stored.
	* @param pos_res : current position in res.
	* @param ncond :  the corresponding number of equations. It is used when the residual is null.
	*/
     void export_tau_val_domain_boundary (const Val_domain& eq, int bound, Array<double>& res, int& pos_res, int ncond) const ;
     virtual void affecte_tau (Tensor&, int, const Array<double>&, int&) const ;

	/**
	* Affects some coefficients to a \c Val_domain.
	* It takes into account the various symmetries and regularity conditions (by means of Garlekin basis).
	* @param so : the field to be affected.
	* @param cf : \c Array of the coefficients used.
	* @param pos_cf : current position in the array of coefficients.
	*/
     void affecte_tau_val_domain (Val_domain& so, const Array<double>& cf, int& pos_cf) const ;
     virtual void affecte_tau_one_coef (Tensor&, int, int, int&) const ;
	/**
	* Sets at most one coefficient of a \c Val_domain to 1.
	* It takes into account the various symmetries and regularity conditions (by means of Garlekin basis).
	* @param so : the \c Val_domain to be affected. It is set to zero if cc does not corresponds to another field.
	* @param cc : location, in the overall system, of the coefficient to be set to 1.
	* @param pos_cf : current position.
	*/
      void affecte_tau_one_coef_val_domain (Val_domain& so, int cc, int& pos_cf) const ;

     virtual double integrale(const Val_domain&) const ; 
    
     friend ostream& operator<< (ostream& o, const Domain_oned_ori& so) ; ///< Display
} ;

/**
* Class for a 1-dimensional spherical domain bounded between two raii.
* It is intended to deal with systems with spherical symmetry.
* \li The numerical coordinate is :
*
* \f$ -1 \leq x \leq 1 \f$
*
* \li Standard spherical coordinates :
*
* \f$ r = \alpha x + \beta \f$
*
* \ingroup domain
*/
class Domain_oned_qcq : public Domain {

 private:
  double alpha ; ///< Relates the numerical radius to the physical one.
  double beta ; ///< Relates the numerical radius the physical one.
  
 public:  
   /**
  * Standard constructor :
  * @param num : number of the domain (used by the \c Space).
  * @param ttype : Chebyshev or Legendre type of spectral expansion.
  * @param x_int : inner radius.
  * @param x_exy : outer radius.
  * @param nbr : number of points in each dimension.
  */
  Domain_oned_qcq (int num, int ttype, double x_int, double x_ext, const Dim_array& nbr) ;
  Domain_oned_qcq (const Domain_oned_qcq& so) ; ///< Copy constructor.
 /**
  * Constructor from a file
  * @param num : number of the domain (used by the \c Space).
  * @param ff: file containd the domain, generated by the save function.
  */
  Domain_oned_qcq (int num, FILE* ff) ;

  virtual ~Domain_oned_qcq() ;
  virtual void save (FILE*) const ;

  private:    
    virtual void do_absol ()  const ;
   
  private :
     virtual void set_cheb_base(Base_spectral&) const ;  
     virtual void set_legendre_base(Base_spectral&) const ;
     virtual void set_cheb_base_odd(Base_spectral&) const ;       
     virtual void set_legendre_base_odd(Base_spectral&) const ;
     
     virtual void do_coloc () ;
     virtual int give_place_var (char*) const ;
     virtual void do_radius () const ;

  public:     
     virtual bool is_in(const Point& xx, double prec=1e-13) const ;
     virtual const Point absol_to_num(const Point&) const;
     virtual void do_der_abs_from_der_var(Val_domain** der_var, Val_domain** der_abs) const ;
   
     virtual Base_spectral mult (const Base_spectral&, const Base_spectral&) const ;

  public:
     virtual void find_other_dom (int, int, int&, int&) const ;     
     virtual Val_domain der_normal (const Val_domain&, int) const ;
     virtual double val_boundary (int, const Val_domain&, const Index&) const ;
     virtual Val_domain der_partial_var (const Val_domain&, int) const ;
   
     virtual int nbr_unknowns (const Tensor&, int) const ;
	/**
	* Computes the number of true unknowns of a \c Val_domain.
	* It takes into account the various symmetries and regularity conditions to determine the precise number of degrees of freedom.
	* @param so : the field.
	* @returns the number of true unknowns.
	*/
     int nbr_unknowns_val_domain (const Val_domain& so) const ;
     virtual Array<int> nbr_conditions (const Tensor&, int, int, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Computes number of discretized equations associated with a given tensorial equation in the bulk.
	* It takes into account the various Galerkin basis used.
	* @param so : the residual of the equation.
	* @param order : order of the equation (i.e. 2 for a Laplacian for instance)
	* @returns the number of true unknowns.
	*/
     int nbr_conditions_val_domain (const Val_domain& so, int order) const ;
     virtual Array<int> nbr_conditions_boundary (const Tensor&, int, int, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Computes number of discretized equations associated with a given equation on a boundary.
	* It takes into account the various Galerkin basis used.
	* It is used for implementing boundary conditions and matching ones.
	* @param eq : the residual of the equation.
	* @returns the number of true conditions.
	*/
     int nbr_conditions_val_domain_boundary (const Val_domain& eq) const ;
     virtual void export_tau (const Tensor&, int, int, Array<double>&, int&, const Array<int>&,int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Exports a residual equation in the bulk.
	* It makes use of the various Galerkin basis used.
	* @param eq : the residual of the equation.
	* @param order : describes the order of the equation (2 for a Laplacian for instance).
	* @param res : The \c Array where the discretized equations are stored.
	* @param pos_res : current position in res.
	* @param ncond :  the corresponding number of equations. It is used when the equation is null.
	*/
     void export_tau_val_domain (const Val_domain& eq, int order, Array<double>& res, int& pos_res, int ncond) const ;
     virtual void export_tau_boundary (const Tensor&, int, int, Array<double>&, int&, const Array<int>&, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Exports all the residual equations corresponding to a tensorial one on a given boundary
	* It makes use of the various Galerkin basis used.
	* @param eq : the residual of the equation.
	* @param bound : the boundary at which the equation is enforced.
	* @param res : The \c Array where the discretized equations are stored.
	* @param pos_res : current position in res.
	* @param ncond :  the corresponding number of equations. It is used when the residual is null.
	*/
     void export_tau_val_domain_boundary (const Val_domain& eq, int bound, Array<double>& res, int& pos_res, int ncond) const ;
     virtual void affecte_tau (Tensor&, int, const Array<double>&, int&) const ;
	/**
	* Affects some coefficients to a \c Val_domain.
	* It takes into account the various symmetries and regularity conditions (by means of Garlekin basis).
	* @param so : the field to be affected.
	* @param cf : \c Array of the coefficients used.
	* @param pos_cf : current position in the array of coefficients.
	*/
     void affecte_tau_val_domain (Val_domain& so, const Array<double>& cf, int& pos_cf) const ;
     virtual void affecte_tau_one_coef (Tensor&, int, int, int&) const ;
	/**
	* Sets at most one coefficient of a \c Val_domain to 1.
	* It takes into account the various symmetries and regularity conditions (by means of Garlekin basis).
	* @param so : the \c Val_domain to be affected. It is set to zero if cc does not corresponds to another field.
	* @param cc : location, in the overall system, of the coefficient to be set to 1.
	* @param pos_cf : current position.
	*/
     void affecte_tau_one_coef_val_domain (Val_domain& so, int cc, int& pos_cf) const ;

     virtual Val_domain div_xp1 (const Val_domain&) const ;
     virtual double integrale(const Val_domain&) const ; 
  
     friend ostream& operator<< (ostream& o, const Domain_oned_qcq& so) ; ///< Display
} ;

/**
* Class for a 1-dimensional compactified spherical domain.
* It is intended to deal with systems with spherical symmetry.
* \li The numerical coordinate is :
*
* \f$ 0 \leq x \leq 1 \f$
*
* \li Standard spherical coordinates :
*
* \f$ r = \frac{1}{\left(\alpha x -1\right)} \f$
*
* \ingroup domain
*/
class Domain_oned_inf : public Domain {

 private:
  double alpha ; ///< Relates the numerical radius to the physical one.

 public:
   /**
  * Standard constructor :
  * @param num : number of the domain (used by the \c Space).
  * @param ttype : Chebyshev or Legendre type of spectral expansion.
  * @param radius : inner radius.
  * @param nbr : number of points in each dimension.
  */
  Domain_oned_inf (int num, int ttype, double radius, const Dim_array& nbr) ;
  Domain_oned_inf (const Domain_oned_inf& so) ; ///< Copy constructor.
 /**
  * Constructor from a file
  * @param num : number of the domain (used by the \c Space).
  * @param ff: file containd the domain, generated by the save function.
  */
  Domain_oned_inf (int num, FILE* ff) ;

  virtual ~Domain_oned_inf() ;
  virtual void save(FILE*) const ;
 
  private:
     virtual void set_cheb_base(Base_spectral&) const ;     
     virtual void set_legendre_base(Base_spectral&) const ;     
     virtual void set_cheb_base_odd(Base_spectral&) const ;       
     virtual void set_legendre_base_odd(Base_spectral&) const ;
     
     virtual void do_coloc()  ;     

     virtual void do_absol ()  const ; 
     virtual int give_place_var (char*) const ;
     virtual void do_radius () const ;

 public:    
	/**
	* Returns the \f$ \alpha \f$ of the mapping.
	*/   
     double get_alpha() const {return alpha ;} ;
   
     virtual bool is_in(const Point& xx, double prec=1e-13) const ;    
     virtual const Point absol_to_num(const Point&) const;
     virtual void do_der_abs_from_der_var(Val_domain** der_var, Val_domain** der_abs) const ; 
     virtual Base_spectral mult (const Base_spectral&, const Base_spectral&) const ;

  public:      
    
     virtual void set_val_inf (Val_domain& so, double xx) const ;    

     virtual void find_other_dom (int, int, int&, int&) const ;     
     virtual Val_domain der_normal (const Val_domain&, int) const ;
     virtual double val_boundary (int, const Val_domain&, const Index&) const ;
     virtual Val_domain der_partial_var (const Val_domain&, int) const ;
     virtual Val_domain mult_xm1 (const Val_domain&) const ;
     virtual Val_domain div_xm1 (const Val_domain&) const ;
     virtual Val_domain mult_x (const Val_domain&) const ;
     virtual Val_domain mult_r (const Val_domain& so) const {return mult_x(so) ;} ;

     virtual int nbr_unknowns (const Tensor&, int) const ;
	/**
	* Computes the number of true unknowns of a \c Val_domain.
	* It takes into account the various symmetries and regularity conditions to determine the precise number of degrees of freedom.
	* @param so : the field.
	* @returns the number of true unknowns.
	*/
     int nbr_unknowns_val_domain (const Val_domain& so) const ;
     virtual Array<int> nbr_conditions (const Tensor&, int, int, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Computes number of discretized equations associated with a given tensorial equation in the bulk.
	* It takes into account the various Galerkin basis used.
	* @param so : the residual of the equation.
	* @param order : order of the equation (i.e. 2 for a Laplacian for instance)
	* @returns the number of true unknowns.
	*/
     int nbr_conditions_val_domain (const Val_domain& so, int order) const ;
     virtual Array<int> nbr_conditions_boundary (const Tensor&, int, int, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Computes number of discretized equations associated with a given equation on a boundary.
	* It takes into account the various Galerkin basis used.
	* It is used for implementing boundary conditions and matching ones.
	* @param eq : the residual of the equation.
	* @returns the number of true conditions.
	*/
     int nbr_conditions_val_domain_boundary (const Val_domain& eq) const ;
     virtual void export_tau (const Tensor&, int, int, Array<double>&, int&, const Array<int>&,int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Exports a residual equation in the bulk.
	* It makes use of the various Galerkin basis used.
	* @param eq : the residual of the equation.
	* @param order : describes the order of the equation (2 for a Laplacian for instance).
	* @param res : The \c Array where the discretized equations are stored.
	* @param pos_res : current position in res.
	* @param ncond :  the corresponding number of equations. It is used when the equation is null.
	*/
     void export_tau_val_domain (const Val_domain& eq, int order, Array<double>& res, int& pos_res, int ncond) const ;
     virtual void export_tau_boundary (const Tensor&, int, int, Array<double>&, int&, const Array<int>&, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Exports all the residual equations corresponding to a tensorial one on a given boundary
	* It makes use of the various Galerkin basis used.
	* @param eq : the residual of the equation.
	* @param bound : the boundary at which the equation is enforced.
	* @param res : The \c Array where the discretized equations are stored.
	* @param pos_res : current position in res.
	* @param ncond :  the corresponding number of equations. It is used when the residual is null.
	*/
     void export_tau_val_domain_boundary (const Val_domain& eq, int bound, Array<double>& res, int& pos_res, int ncond) const ;
     virtual void affecte_tau (Tensor&, int, const Array<double>&, int&) const ;
	/**
	* Affects some coefficients to a \c Val_domain.
	* It takes into account the various symmetries and regularity conditions (by means of Garlekin basis).
	* @param so : the field to be affected.
	* @param cf : \c Array of the coefficients used.
	* @param pos_cf : current position in the array of coefficients.
	*/
     void affecte_tau_val_domain (Val_domain& so, const Array<double>& cf, int& pos_cf) const ;
     virtual void affecte_tau_one_coef (Tensor&, int, int, int&) const ;
	/**
	* Sets at most one coefficient of a \c Val_domain to 1.
	* It takes into account the various symmetries and regularity conditions (by means of Garlekin basis).
	* @param so : the \c Val_domain to be affected. It is set to zero if cc does not corresponds to another field.
	* @param cc : location, in the overall system, of the coefficient to be set to 1.
	* @param pos_cf : current position.
	*/
     void affecte_tau_one_coef_val_domain (Val_domain& so, int cc, int& pos_cf) const ;

     virtual double integrale(const Val_domain&) const ; 
  
   friend ostream& operator<< (ostream& o, const Domain_oned_inf& so) ; ///< Display
} ;


/**
 * The \c Space_oned class fills the space with 1-dimensional domains (spherical ones), intended for spherically symmetric problems.
 * \ingroup domain
 */
class Space_oned : public Space {
     public:
    	/**
     	* Standard constructor 
     	* @param ttype [input] : the type of basis.
	* @param nbr [input] : number of points in each domain.
	* @param bounds [input] : radii of the various shells (and also determines the total number of domains).
	*/
	Space_oned (int ttype, const Dim_array& nbr, const Array<double>& bounds) ;
	Space_oned (FILE*) ; ///< Constructor from a file.
	virtual ~Space_oned() ;        
	virtual void save(FILE*) const ;

	/**
	* Adds an equation being the value of some field at the origin.
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the quantity that must be zero at the origin
	*/
	void add_eq_ori (System_of_eqs& syst, const char* eq) ;
} ;
}
#endif
