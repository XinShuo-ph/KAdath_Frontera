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

#ifndef __SCALAR_HPP_
#define __SCALAR_HPP_

#include "space.hpp"
#include "val_domain.hpp"
#include <memory>
#include "tensor.hpp"

using std::unique_ptr;

namespace Kadath {
//Versions in all space
Scalar operator+ (const Scalar&) ;
Scalar operator- (const Scalar&) ;
Scalar operator+ (const Scalar&, const Scalar&) ;
Scalar operator+ (const Scalar&, double) ;
Scalar operator+ (double, const Scalar&) ;
Scalar operator- (const Scalar&, const Scalar&) ;
Scalar operator- (const Scalar&, double) ;
Scalar operator- (double, const Scalar&) ;
Scalar operator* (const Scalar&, const Scalar&) ;
Scalar operator* (const Scalar&, double) ;
Scalar operator* (double, const Scalar&) ;
Scalar operator* (const Scalar&, int) ;
Scalar operator* (int, const Scalar&) ;
Scalar operator* (const Scalar&, long int) ;
Scalar operator* (long int, const Scalar&) ;
Scalar operator/ (const Scalar&, const Scalar&) ;
Scalar operator/ (const Scalar&, double) ;
Scalar operator/ (double, const Scalar&) ;
Scalar pow (const Scalar&, int) ;
Scalar pow (const Scalar&, double) ;
Scalar sqrt (const Scalar&) ;
Scalar exp (const Scalar&) ;
Scalar sin (const Scalar&) ;
Scalar cos (const Scalar&) ;
Scalar atan (const Scalar&) ;
double diffmax (const Scalar&, const Scalar&) ;

class Vector ;

/**
* The class \c Scalar does not really implements scalars in the mathematical sense but rather tensorial coordinates of tensors.
* This class is mainly an array of various \c Val_domain. It also stores some quantities like the derivatives of 
* the field with respect to the absolute Cartesian coordinates.
* \ingroup fields
*/
class Scalar : public Tensor {

    protected:
      Val_domain** val_zones; ///< Pointers on the various \c Val_domain describing the field in each \c Domain.

    public:
      /**
       * Standard constructor. The value of the field is not initialized.
       * @param sp [input] : the space on which the field is defined.
       */
      Scalar (const Space&) ;
       /**
        * Copy constructor.
        * @param so [input] : \c Scalar to be copied.
	* @param copy [input] : the values of \c so are only copied if \c copy is \c true. Otherwise, the values of 
	* the field are left uninitialized.
	*/
      Scalar (const Scalar& so, bool copy=true) ;

	/**
        * Constructor from a \c Tensor.
        * @param so [input] : \c Tensor to be copied : it must be a valence zero tensor.
	* @param copy [input] : the values of \c so are only copied if \c copy is \c true. Otherwise, the values of 
	* the field are left uninitialized.
	*/
      Scalar (const Tensor& a, bool copy=true) ;
	/**
	* Constructor from a file.
	* @param sp : the \c Space.
	* @param : fd file (generated by the saving function).
	*/
      Scalar (const Space& sp, FILE* fd) ;

      virtual ~Scalar () ; ///< Destructor.
      
      virtual void save (FILE*) const ; ///< Saving function
    public:
	/**
	* @returns the number of dimensions.
	*/
      int get_ndim() const {return ndim ;} ;
	/**
	* @returns the number of domains.
	*/
      int get_nbr_domains() const {return ndom ;} ; 
	/**
	* @param : index of the desired \c Domain
	* @returns Pointer on a particular domain.
	*/
      const Domain* get_domain(int i) const {return val_zones[i]->get_domain() ;} ;
	/**
	* @returns the \c Space.
	*/
      const Space& get_space () const {return espace ;} ;
      
    public:
      Val_domain& set_domain (int) ; ///< Read/write of a particular \c Val_domain.
      const Val_domain& operator() (int) const; ///< Read only of a particular \c Val_domain.
      void operator= (const Scalar&) ; ///< Assignement to another \c Scalar.
      virtual void operator= (const Tensor&) ; ///< Assignement to a \c Tensor (must be scalar)
      void operator= (double) ; ///< Assignment to a double (sets all the values in the configuration space to that value.
      virtual void annule_hard() ; ///< Sets the value to zero evetywhere (the logical state of the \c Val_domain is NOT zero).

      Scalar der_var (int) const ; ///< Returns the derivative with respect to one particular numerical coordinate.
      Scalar der_abs (int) const ;///< Returns the derivative with respect to one particular absolute Cartesian coordinate.
      Scalar der_spher (int) const ;///< Returns the derivative with respect to one particular absolute Cartesian coordinate.
      Scalar der_r () const ; ///< Returns the radial derivative.
      Scalar div_r () const ; ///< Returns the division by \f$ r\f$.
      Scalar div_rsint () const ; ///< Returns the division by \f$ r\sin \theta\f$.
      Scalar div_1mx2 () const ; ///< Returns the division by \f$ 1-x^2\f$.
      Scalar mult_cos_theta () const ; ///< Returns the multiplication by \f$\cos \theta\f$.
      Scalar mult_sin_theta () const ; ///< Returns the multiplication by \f$\sin \theta\f$.
      Scalar mult_cos_phi () const ;///< Returns the multiplication by \f$\cos \varphi\f$.
      Scalar mult_sin_phi () const ;///< Returns the multiplication by \f$\sin \varphi\f$.
      double integrale() const ;///< Returns the integral in the whole space.

      unique_ptr<Scalar> clone() const; ///< Copy using unique_ptr
      Vector grad() const ; ///< Computes the gradient (in Cartesian coordinates).

    public:
	/**
	* Destroys the values in the coefficient space.
	*/
      void set_in_conf() ; 
	/**
	* Destroys the values in the configuration space.
	*/
      void set_in_coef() ;
  	/**
	* Allocates the values in the configuration space and destroys the values in the coefficients space.
	*/
      void allocate_conf() ;
	/**
	* Allocates the values in the coefficient space and destroys the values in the configuration space.
	*/
      void allocate_coef() ;
      void std_base() ;  ///< Sets the standard basis of decomposition
      void std_anti_base() ;  ///< Sets the standard, anti-symetric, basis of decomposition      
      void std_base(int m) ;  ///< Sets the standard basis of decomposition assuming  a given harmonic wrt \f$\varphi\f$.
      void std_base(int l, int m);///< Sets the standard basis of decomposition assuming  a given harmonic wrt \f$\varphi\f$ and \f$\theta\f$.
      void std_anti_base(int m) ;  ///< Sets the standard, anti-symetric, basis of decomposition assuming  a given harmonic wrt \f$\varphi\f$.
      void std_base_r_spher() ; ///< Sets the basis for the radial component of a vector in orthonormal spherical coordinates.
      void std_base_t_spher() ; ///< Sets the basis for the \f$\theta\f$ component of a vector in orthonormal spherical coordinates.
      void std_base_p_spher() ; ///< Sets the basis for the \f$\varphi\f$ component of a vector in orthonormal spherical coordinates.
      
      void std_base_domain(int) ; ///< Sets the standard basis of decomposition, in a given \c Domain
      void std_anti_base_domain(int) ;  ///< Sets the standard, anti-symetric, basis of decomposition, in a given \c Domain     
      void std_base_domain(int, int m) ;  ///< Sets the standard basis of decomposition assuming  a given harmonic wrt \f$\varphi\f$, in a given \c Domain.
      void std_base_domain(int d, int l, int m) ;///< Sets the standard basis of decomposition assuming  a given harmonic wrt \f$\varphi\f$ and \f$\theta\f$, in a given \c Domain.
      void std_base_r_spher_domain(int) ; ///< Sets the basis for the radial component of a vector in orthonormal spherical coordinates, in a given \c Domain.
      void std_base_t_spher_domain(int) ; ///< Sets the basis for the \f$\theta\f$ component of a vector in orthonormal spherical coordinates, in a given \c Domain.
      void std_base_p_spher_domain(int) ; ///< Sets the basis for the \f$\varphi\f$ component of a vector in orthonormal spherical coordinates, in a given \c Domain.
      void std_base_x_cart_domain(int) ; ///< Sets the basis for the X-component of a vector in Cartesian coordinates, in a given \c Domain.
      void std_base_y_cart_domain(int) ; ///< Sets the basis for the Y-component of a vector in Cartesian coordinates, in a given \c Domain.
      void std_base_z_cart_domain(int) ; ///< Sets the basis for the Z-component of a vector in Cartesian coordinates, in a given \c Domain.

      void std_base_rt_spher_domain(int d) ;///< Sets the basis for the \f$(r,\theta)\f$ component of a 2-tensor in orthonormal spherical coordinates, in a given \c Domain.
      void std_base_rp_spher_domain(int d) ; ///< Sets the basis for the \f$(r,\varphi)\f$ component of a 2-tensor in orthonormal spherical coordinates, in a given \c Domain.
      void std_base_tp_spher_domain(int d) ;///< Sets the basis for the \f$(\theta, \varphi)\f$ component of a 2-tensor in orthonormal spherical coordinates, in a given \c Domain.
      void std_base_xy_cart_domain(int d) ;///< Sets the basis for the XY component of a 2-tensor in Cartesian coordinates, in a given \c Domain.
      void std_base_xz_cart_domain(int d) ;///< Sets the basis for the XZ component of a 2-tensor in Cartesian coordinates, in a given \c Domain.
      void std_base_yz_cart_domain(int d) ;///< Sets the basis for the YZ component of a 2-tensor in Cartesian coordinates, in a given \c Domain.
   
      // For critic solution
      void std_xodd_base() ;  ///< Sets the basis for an odd function in \f$X\f$ (Critic case).
      void std_todd_base() ; ///< Sets the basis for an odd function in \f$T\f$ (Critic case).
      void std_xodd_todd_base() ; ///< Sets the basis for an odd function in \f$X\f$ and \f$T\f$ (Critic case).
    
      void std_base_odd() ; ///< Sets the basis in odd polynomials.

      void set_val_inf (double xx) ; ///< Sets the value at infinity (in the last domain) to \c xx.
      void set_val_inf (double xx, int l) ; ///< Sets the value at infinity (in the domain \c l) to \c xx.

      double integ_volume () const ; ///< @returns integral in the whole space.
    /**
      * Computes the value of the field at a given point, by doing the spectral summation.
      * @param xxx [input] : absolute Cartesian coordinates of the point.
      * @param sens : looks for the point starting from the origin (+1) or infinity (-1).
      * @returns the value of the field.
      */
      double val_point (const Point& xxx, int sens=-1) const ;
      /**
      * Computes the value of the field at a given point, by doing the spectral summation. Returns zero if the \c Point is not found in the computational domain.
      * @param xxx [input] : absolute Cartesian coordinates of the point.
      * @param sens : looks for the point starting from the origin (+1) or infinity (-1).
      * @returns the value of the field.
      */
      double val_point_zeronotdef (const Point& xxx, int sens=-1) const ;
   
      /**
	* Affects all the values to the one of another scalar.
	* This is done by using spectral summation and so does not requires the two fields to be have the same collocation points.
	* @param so : the source \c Scalar.
	*/
      void import (const Scalar& so) ;

    public:
      void coef() const ; ///< Computes the coefficients.
      void coef_i() const ; ///< Computes the values in the configuration space.
 	/**
      * Sets to zero all the coefficients above a given order, for the \f$ \varphi\f$ coefficients, in a gicen \c Domain.
      * Takes into account the various Galerkin basis to maintain regularity.
      * @param dom : the \c Domain where the filter is applied.
      * @param ncf : the coefficients which index is above this are set to zero.
      */
      void filter_phi (int dom, int ncf) ;
      /**
	* @returns : a \c Scalar containing zero. 
	* The result is logically zero.
	*/
      static Scalar zero(Space const& espace);

      void operator+= (const Scalar&) ; ///< Operator +=
      void operator-= (const Scalar&) ; ///< Operator -=
      void operator*= (const Scalar&) ; ///< Operator *=
      void operator/= (const Scalar&) ; ///< Operator /=
      void operator+= (double) ; ///< Operator +=
      void operator-= (double) ; ///< Operator -=
      void operator*= (double) ; ///< Operator *=
      void operator/= (double) ; ///< Operator /=

      friend ostream& operator<< (ostream& o, const Scalar&) ; ///< Display
      friend class Space ;  
      friend class Space_spheric ;
    
	friend Scalar operator+ (const Scalar&) ;
	friend Scalar operator- (const Scalar&) ;
	friend Scalar operator+ (const Scalar&, const Scalar&) ;
	friend Scalar operator+ (const Scalar&, double) ;
	friend Scalar operator+ (double, const Scalar&) ;
	friend Scalar operator- (const Scalar&, const Scalar&) ;
	friend Scalar operator- (const Scalar&, double) ;
	friend Scalar operator- (double, const Scalar&) ;
	friend Scalar operator* (const Scalar&, const Scalar&) ;
	friend Scalar operator* (const Scalar&, double) ;
	friend Scalar operator* (double, const Scalar&) ;
	friend Scalar operator/ (const Scalar&, const Scalar&) ;
	friend Scalar operator/ (const Scalar&, double) ;
	friend Scalar operator/ (double, const Scalar&) ;
	friend Scalar pow (const Scalar&, int) ;
	friend Scalar pow (const Scalar&, double) ;
	friend Scalar sqrt (const Scalar&) ;
	friend Scalar exp (const Scalar&) ;
	friend Scalar sin (const Scalar&) ;
	friend Scalar cos (const Scalar&) ;
	friend double diffmax (const Scalar&, const Scalar&) ;

    friend Scalar operator+(const Tensor&, const Scalar&) ;
    friend Scalar operator+(const Scalar&, const Tensor&) ;
    friend Scalar operator-(const Tensor&, const Scalar&) ;
    friend Scalar operator-(const Scalar&, const Tensor&) ;
    friend Tensor operator*(const Scalar&, const Tensor&) ;
    friend Tensor operator*(const Tensor&, const Scalar&) ;
    friend Tensor operator/(const Tensor&, const Scalar&) ; 
} ;

void des_coupe (const Scalar& uu, const Point& x0, 
		int num_un, double var_un_min, double var_un_max, 
		int num_deux, double var_deux_min, double var_deux_max, 
		const char* title = 0x0, const char* axis_one=0x0, const char* axis_two = 0x0, int ncour=15, int n_un=100, int n_deux=100) ;

void des_coupe_zeronotdef (const Scalar& uu, const Point& x0, 
		int num_un, double var_un_min, double var_un_max, 
		int num_deux, double var_deux_min, double var_deux_max, 
		const char* title = 0x0, const char* axis_one=0x0, const char* axis_two = 0x0, int ncour=15, int n_un=100, int n_deux=100) ;


void des_sphere (const Scalar& uu, const Point& x0, double rad, const char* title = 0x0, int ncour=15, int n_un=100, int n_deux=100) ;
}
#endif
