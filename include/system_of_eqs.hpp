/*
    Copyright 2017 Philippe Grandclement & Gregoire Martinon

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

#ifndef __SYSTEM_OF_EQS_HPP_
#define __SYSTEM_OF_EQS_HPP_

#include "headcpp.hpp"
#ifdef PAR_VERSION
#include "mpi.h"
#endif
#include "profiled_object.hpp"
#include "tensor.hpp"
#include "ope_eq.hpp"
#include "term_eq.hpp"
#include "metric.hpp"
#include "param.hpp"
#include "matrice.hpp"
#include "list_comp.hpp"

#include <vector>
using std::vector;

constexpr int VARMAX{1000};

#define GMRES_SPARSE 0
#define BICGSTAB_SPARSE 1

#define SAVED_DER  1

namespace Kadath {
	class Equation ;
	class Eq_int ;

	/**
	 * Class used to describe and solve a system of equations. It is the central object of Kadath.
	 * The equations are solved between the domains \c dom_min and \c dom_max .
	 * The various quantities are given names (char*) that are used when passing the equations to the \c System_of_eqs.
	 * \ingroup systems
	 */

	class System_of_eqs : public Profiled_object<System_of_eqs> {
	public:
		static constexpr std::size_t default_block_size {64};
		static constexpr std::size_t nb_core_per_node{24};
		//! Dummy variable for the purpose of better readability.
		static constexpr int ALL_COLUMNS {-1};
		//! Dummy names for the purpose of better readability.
		enum : bool { DO_NOT_TRANSPOSE = false, TRANSPOSE = true};
		enum Matrix_computation_parallel_paradigm: unsigned short {sequential, multi_thread, mpi, hybrid};

	protected:
	    std::ostream * output_stream; ///< Default output stream for log messages.
		const Space& espace ; ///< Associated \c Space
		int dom_min ; ///< Smallest domain number
		int dom_max ; ///< Highest domain number
		int ndom ; ///< Number of domains used.

		int nvar_double ; ///< Number of unknowns that are numbers (i.e. not fields)
		MMPtr_array<double> var_double ; ///< Pointer on the unknowns that are numbers  (i.e. not fields)
		MMPtr_array<char> names_var_double ; ///< Names of the unknowns that are numbers (i.e. not fields)

		int nvar ; ///< Number of unknown fields.
		MMPtr_array<Tensor> var ; ///< Pointer on the unknown fields.
		MMPtr_array<char> names_var ; ///< Names of the unknown fields.

		int nterm_double ; ///< Number of \c Term_eq corresponding to the unknowns that are numbers.
		MMPtr_array<Term_eq> term_double ; ///< Pointers on the \c Term_eq corresponding to the unknowns that are numbers.
		Array<int> assoc_var_double ; ///< Array giving the correspondance with the \c var_double pointers.

		int nterm ; ///< Number of \c Term_eq corresponding to the unknown fields.
		MMPtr_array<Term_eq> term ; ///< Pointers on the \c Term_eq corresponding to the unknown fields.
		Array<int> assoc_var ;///< Array giving the correspondance with the \c var pointers.

		int ncst ; ///< Number of constants passed by the user.
		int nterm_cst ; ///< Number of \c Term_eq coming from the constants passed by the user.
		MMPtr_array<Term_eq> cst ; ///< Pointers on the \c Term_eq coming from the constants passed by the user.
		MMPtr_array<char> names_cst ; ///< Names of the constants passed by the user.

		mutable int ncst_hard ; ///< Number of constants generated on the fly (when encoutering things like "2.2" etc...)
		mutable MMPtr_array<Term_eq> cst_hard ;///< Pointers on the \c Term_eq coming from the constants generated on the fly (when encoutering things like "2.2" etc...)
		mutable Array<double> val_cst_hard ; ///< Values of the constants generated on the fly (when encoutering things like "2.2" etc...)

		// Definitions (liste of operators I guess)
		int ndef ; ///< Number of definitions.
		MMPtr_array<Ope_def> def ; ///< Pointers on the definition (i.e. on the \c Ope_def that is needed to compute the result).
		MMPtr_array<char> names_def ; ///< Names of the definitions.

		// Definitions globale (liste of operators I guess)
		int ndef_glob ; ///< Number of global definitions (the one that require the knowledge of the whole space to give the result, like integrals).
		MMPtr_array<Ope_def_global> def_glob ; ///< Pointers on the global definitions.
		MMPtr_array<char> names_def_glob ;///< Names of the global definitions.

		// ser defined operators
		int nopeuser ; ///< Number of operators defined by the user (single argument).
		Term_eq (*opeuser[VARMAX]) (const Term_eq&, Param*) ; ///< Pointers on the functions used by the user defined operators (single argument).
		MMPtr_array<Param> paruser ; ///< Parameters used by the user defined operators (single argument).
		MMPtr_array<char> names_opeuser ; ///< Names of the user defined operators (single argument).

		int nopeuser_bin ;///< Number of operators defined by the user (with two arguments).
		Term_eq (*opeuser_bin[VARMAX]) (const Term_eq&, const Term_eq&, Param*) ;///< Pointers on the functions used by the user defined operators (with two arguments).
		MMPtr_array<Param> paruser_bin ; ///< Parameters used by the user defined operators (with two arguments).
		MMPtr_array<char> names_opeuser_bin ;///< Names of the user defined operators (with two arguments).

		// The metric (only one for now)
		Metric* met ; ///< Pointer on the associated \c Metric, if defined.
		char* name_met ; ///< Name by which the metric is recognized.

		int neq_int ; ///< Number of integral equations (i.e. which are doubles)
		MMPtr_array<Eq_int> eq_int ; ///< Pointers onto the integral equations.

		int neq ; ///< Number of field equations.
		MMPtr_array<Equation> eq ; ///< Pointers onto the field equations.

		MMPtr_array<Term_eq> results ; ///< Pointers on the residual of the various equations.

		int nbr_unknowns ; ///< Number of unknowns (basically the number of coefficients of all the unknown fields, once regularities are taken into account).
		int nbr_conditions ; ///< Total number of conditions (the number of coefficients of all the equations, once regularities are taken into account).

		MMPtr_array<Index> which_coef ; ///< Stores the "true" coefficients on some boundaries (probably deprecated).

		unsigned niter{0u}; ///< Counter toward the number of times the \c do_newton method has been called.
		bool display_newton_data{true};///< Boolean used to enable or disable the display of data during the Newton iterations.
		int mpi_world_size{1};
		int mpi_proc_rank{0};
		void init_proc_data() {
#ifdef PAR_VERSION
            MPI_Comm_rank(MPI_COMM_WORLD,&mpi_proc_rank);
            MPI_Comm_size(MPI_COMM_WORLD,&mpi_world_size);
#endif
		}

	public:
		/**
		* Standard constructor nothing is done. The space is affected and the equations are to be solved in all space.
		* @param so [input] : associated space.
		*/
		System_of_eqs (const Space& so) ;
		/**
		* Constructor, nothing is done. The space is affected and the equations are solved only between two domains.
		* @param so [input] : associated space.
		* @param i [input] : smallest domain number.
		* @param j [input] : highest domain number.
		*/
		System_of_eqs (const Space& so, int i, int j) ;
		/**
		* Constructor, nothing is done. The space is affected and the equations are solved only in one domain.
		* @param so [input] : associated space.
		* @param i [input] : the domain number.
		**/
		System_of_eqs (const Space& so, int i) : System_of_eqs{so,i,i} {}
		System_of_eqs (const System_of_eqs&) = delete ; ///< Constructor by copy.
		~System_of_eqs() ; ///< Destructor.

		const Metric* get_met() const ; ///< Returns a pointer on the \c Metric.

		// Accessors
		/**
		 * Returns the default output stream reference.
		 */
		std::ostream & get_output_stream() const {return *output_stream;}
		/**
		 * Sets a new output stream reference.
		 */
		void set_output_stream(std::ostream & new_output_stream) {output_stream = &new_output_stream;}
		/**
		* Returns the space.
		*/
		const Space& get_space() const {return espace ;} ;
		/**
		* Returns the smallest index of the domains.
		*/
		int get_dom_min() const {return dom_min ;} ;
		/**
		* Returns the highest index of the domains.
		*/
		int get_dom_max() const {return dom_max ;} ;
		/**
		* Returns the number of conditions.
		*/
		int get_nbr_conditions() const {return nbr_conditions ;} ;
		/**
		* Returns the number of unknowns.
		*/
		int get_nbr_unknowns() const {return nbr_unknowns ;} ;
		/**
		 * Returns the current iteration number.
		 */
		unsigned get_niter() const {return niter;}
		System_of_eqs & enable_data_display() {display_newton_data = true; return *this;}
		System_of_eqs & disable_data_display() {display_newton_data = false; return *this;}
        int get_mpi_proc_rank() const {return mpi_proc_rank;}
        int get_mpi_world_size() const {return mpi_world_size;}

		/**
		* Returns a pointer on a \c Term_eq corresponding to an unknown number.
		* @param which : index of the unknown.
		* @param dd : the index of the \c Domain.
		*/
		Term_eq* give_term_double (int which, int dd) const ;
		/**
		* Returns a pointer on a \c Term_eq corresponding to an unknown field.
		* @param which : index of the unknown.
		* @param dd : the index of the \c Domain.
		*/
		Term_eq* give_term (int which, int dd) const ;
		/**
		* Returns a pointer on a \c Term_eq corresponding to a constant.
		* @param which : index of the unknown.
		* @param dd : the index of the \c Domain.
		*/
		Term_eq* give_cst (int which, int dd) const ;
		/**
		* Returns a pointer on a \c Term_eq corresponding to a constant generated on the fly.
		* @param xx : the value of the constant
		* @param dd : the index of the \c Domain.
		*/
		Term_eq* give_cst_hard (double xx, int dd) const ;
		/**
		* Returns a pointer on a definition (better to use \c give_val_def if one wants to access the result of some definition).
		* @param  i : the index of the definition (one has to manage the different domains properly...)
		*/
		Ope_def* give_def (int i) const ;
		/**
		* Returns a pointer on a global definition.
		* @param  i : the index of the definition (one has to manage the different domains properly...)
		*/
		Ope_def_global* give_def_glob (int i) const ;

		/**
		* Gives the result of a definition. It is set to zero on domains where the definition is undefined.
		* @param name : name of the definition.
		*/
		Tensor give_val_def (const char* name) const ;


		/**
		 * Addition of a variable (number case)
		 * @param name : name of the variable (used afterwards by \c System_of_eqs)
		 * @param var : variable.
		 */
		virtual void add_var (const char* name, double& var) ;
		/**
		 * Addition of a variable (field case)
		 * @param name : name of the variable (used afterwards by \c System_of_eqs)
		 * @param var : variable.
		 */
		virtual void add_var (const char* name, Tensor& var) ;
		/**
		 * Addition of a constant (number case)
		 * @param name : name of the constant (used afterwards by \c System_of_eqs)
		 * @param cst : variable.
		 */
		virtual void add_cst (const char* name, double cst) ;
		/**
		 * Addition of a constant (field case)
		 * @param name : name of the constant (used afterwards by \c System_of_eqs)
		 * @param cst : constant.
		 */
		virtual void add_cst (const char* name, const Tensor& cst) ;

		/**
		 * Addition of a definition
		 * @param name : string describing the definition (like "A=...")
		 */
		virtual void add_def (const char* name) ;
		/**
		 * Addition of a definition in a single domain.
		 * @param dd : number of the \c Domain.
		 * @param name : string describing the definition (like "A=...")
		 */
		virtual void add_def (int dd, const char* name) ;
		/**
		 * Addition of a global definition
		 * @param name : string describing the definition (like "A=...")
		 */
		virtual void add_def_global (const char* name) ;
		/**
		 * Addition of a global definition in a single domain.
		 * @param dd : number of the \c Domain.
		 * @param name : string describing the definition (like "A=...")
		 */
		virtual void add_def_global (int dd, const char* name) ;
		/**
		* Addition of a user defined operator (one argument version)
		* @param name : name of the operator (used afterwards by \c System_of_eqs)
		* @param pope : pointer on the function describing the action of the  operator.
		* @param par : parameters of the operator.
		*/
		virtual void add_ope (const char* name, Term_eq (*pope) (const Term_eq&, Param*), Param* par) ;
		/**
		* Addition of a user defined operator (two arguments version)
		* @param name : name of the operator (used afterwards by \c System_of_eqs)
		* @param pope : pointer on the function describing the action of the  operator.
		* @param par : parameters of the operator.
		*/
		virtual void add_ope (const char* name, Term_eq (*pope) (const Term_eq&, const Term_eq&, Param*), Param* par) ;

		/**
		* Check if a string is an unknown (number).
		* @param target : the string to be tested.
		* @param which : the index of the found variable (if found).
		*/
		bool isvar_double (const char* target, int& which) const ;
		/**
		* Check if a string is an unknown field.
		* @param target : the string to be tested.
		* @param which : the index of the found variable (if found).
		*/
		bool isvar (const char*, int&, int&, char*&, Array<int>*&) const ;
		/**
		* Check if a string is a constant (can required indices manipulation and/or inner contraction).
		* @param target : the string to be tested.
		* @param which : the index of the found constant (if found).
		* @param valence : valence of the result.
		* @param name_ind : name of the indices of the result.
		* @param type_ind : type of the indices of the result (COV or CON).
		*/
		bool iscst (const char* target, int& which, int& valence, char*& name_ind, Array<int>*& type_ind) const ;
		/**
		* Check if a string is a definition (can required indices manipulation and/or inner contraction).
		* @param dd : index of the \c Domain.
		* @param target : the string to be tested.
		* @param which : the index of the found definition (if found).
		* @param valence : valence of the result.
		* @param name_ind : name of the indices of the result.
		* @param type_ind : type of the indices of the result (COV or CON).
		*/
		bool isdef (int dd, const char* target, int& which, int& valence, char*& name_ind, Array<int>*& type_ind) const ;
		/**
		* Check if a string is a global definition.
		* @param dd : index of the \c Domain.
		* @param target : the string to be tested.
		* @param which : the index of the found definition (if found).
		*/
		bool isdef_glob (int dd, const char* target, int& which) const ;


		/**
		* Checks if a string contains the operator minus
		* @param input : the string to be tested.
		* @param output : returns the argument of the minus operator ("A" if one started with "-A")
		*/
		bool is_ope_minus (const char* input, char* output) const ;
		/**
		* Checks if a string is a double
		* @param input : the string to be tested.
		* @param output : returns the double
		*/
		bool isdouble (const char* input, double& output) const ;
		/**
		* Checks if a string is a metric
		* @param input : the string to be tested.
		* @param name_ind : name of the indices of the result.
		* @param type_ind : type of the indices of the result (COV or CON).
		*/
		bool ismet (const char* input, char*& name_ind, int& type_ind) const ;
		/**
		* Checks if a string is a metric (without arguments, probably deprecated)
		* @param input : the string to be tested.
		*/
		bool ismet (const char* input) const ;
		/**
		* Checks if a string represents the Christoffel symbols.
		* The reserved word is "Gam"
		* @param input : the string to be tested.
		* @param name_ind : name of the indices of the result.
		* @param type_ind : type of the indices of the result (COV or CON).
		*/
		bool ischristo (const char* input, char*& name_ind, Array<int>*& type_ind) const ;
		/**
		* Checks if a string represents the Riemann tensor.
		* The reserved word is "R"
		* @param input : the string to be tested.
		* @param name_ind : name of the indices of the result.
		* @param type_ind : type of the indices of the result (COV or CON).
		*/
		bool isriemann (const char* input, char*& name_ind, Array<int>*&type_ind) const ;
		/**
		* Checks if a string represents the Ricci tensor.
		* The reserved word is "R"
		* @param input : the string to be tested.
		* @param name_ind : name of the indices of the result.
		* @param type_ind : type of the indices of the result (COV or CON).
		*/
		bool isricci_tensor (const char* input, char*& name_ind, Array<int>*& type_ind) const ;
		/**
		* Checks if a string represents the Ricci scalar.
		* The reserved word is "R"
		* @param input : the string to be tested.
		* @param name_ind : name of the indices of the result (not used in this case).
		* @param type_ind : type of the indices of the result (COV or CON) (not used in this case).
		*/
		bool isricci_scalar (const char* input, char*& name_ind, Array<int>*& type_ind) const ;
		/**
		* Checks if a string represents an operator of the type "a + b".
		* @param input : the string to be tested.
		* @param p1 : the returned first argument.
		* @param p2 : the returned second argument.
		* @param symb : the symbol representing the operator (can be +, -; * and so on...)
		*/
		bool is_ope_bin (const char* input, char* p1, char* p2, char symb) const ;
		/**
		* Checks if a string represents an operator of the type "ope(a)".
		* @param input : the string to be tested.
		* @param p1 : the returned  argument.
		* @param nameope : name of the operator
		*/
		bool is_ope_uni (const char* input, char* p1, const char* nameope) const ;
		/**
		* Checks if a string represents an operator of the type "ope(a,b)".
		* @param input : the string to be tested.
		* @param p1 : the returned  first argument.
		* @param p2 : the returned  first argument.
		* @param nameope : name of the operator
		*/
		bool is_ope_uni (const char* input, char* p1, char* p2, const char* nameope) const ;
		/**
		* Checks if a string represents the covariant derivative.
		* @param input : the string to be tested.
		* @param p1 : the returned argument.
		* @param typeder : the type of derivative (COV or CON).
		* @param nameind : name of the index corresponding to the derivative.
		*/
		bool is_ope_deriv (const char* input, char* p1, int& typeder, char& nameind) const ;
		/**
		* Checks if a string represents the flat covariant derivative.
		* @param input : the string to be tested.
		* @param p1 : the returned argument.
		* @param typeder : the type of derivative (COV or CON).
		* @param nameind : name of the index corresponding to the derivative.
		*/
		bool is_ope_deriv_flat (const char* input, char* p1, int& typeder, char& nameind) const ;
		/**
		* Checks if a string represents the covariant derivative wrt a background metric.
		* @param input : the string to be tested.
		* @param p1 : the returned argument.
		* @param typeder : the type of derivative (COV or CON).
		* @param nameind : name of the index corresponding to the derivative.
		*/
		bool is_ope_deriv_background (const char* input, char* p1, int& typeder, char& nameind) const ;
		/**
		* Checks if a string represents the power of something (like "a^2").
		* @param input : the string to be tested.
		* @param p1 : the returned argument.
		* @param expo : the returned power.
		*/
		bool is_ope_pow (const char* input, char* p1, int& expo) const ;
		/**
		* Checks if a string represents the partial derivative (like "partial_i a")
		* @param input : the string to be tested.
		* @param p1 : the returned argument.
		* @param nameind : name of the index corresponding to the derivative.
		*/
		bool is_ope_partial (const char* input, char* p1, char& nameind) const ;
		/**
		* Checks if a string represents the derivative wrt a numerical coordinate of a given \c Domain (like "a,T")
		* @param dd : number of the \c Domain
		* @param input : the string to be tested.
		* @param p1 : the returned argument.
		* @param which : the variable concerned by the derivative.
		*/
		bool is_ope_der_var (int dd, const char* input, char* p1, int& which) const ;

		/**
		* Function that reads a string and returns a pointer on the generated \c Ope_eq.
		* It is higly recursive, calling itslef unless the full operator is generated or an error encoutered.
		* @param dom : number of the \c Domain
		* @param name : the sting to be read.
		* @param bb : boundary, if a boundary is needed (depending on the type of equation considered : bulk vs matching for instance).
		*/
		Ope_eq* give_ope (int dom, const char* name, int bb=0) const ;

		/**
		* Addition of an equation to be solved inside a domain (assumed to be second order).
		* @param dom : number of the \c Domain.
		* @param eq : string defining the equation.
		* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
		* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
		*/
		virtual void add_eq_inside (int dom, const char* eq, int n_cmp = -1, Array<int>** p_cmp=nullptr) ;

		/**
		* Addition of an equation to be solved inside a domain (assumed to be second order).
		* Version with a list of components
		* @param dom : number of the \c Domain.
		* @param eq : string defining the equation.
		* @param list : list of the components to be considered
		*/
		virtual void add_eq_inside (int dom, const char* eq, const List_comp& list) ;

		/**
		* Addition of an equation to be solved inside a domain (of arbitrary order).
		* @param dom : number of the \c Domain.
		* @param order : order of the equation.
		* @param eq : string defining the equation.
		* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
		* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
		*/
		virtual void add_eq_order (int dom, int order, const char* eq, int n_cmp = -1, Array<int>** p_cmp=nullptr) ;

		/**
		* Addition of an equation to be solved inside a domain (of arbitrary order).
		* @param dom : number of the \c Domain.
		* @param order : order of the equation.
		* @param eq : string defining the equation.
		* @param list : list of the components to be considered.
		*/
		virtual void add_eq_order (int dom, int order, const char* eq,  const List_comp& list) ;

		/**
		* Addition of an equation describing a boundary condition.
		* @param dom : number of the \c Domain.
		* @param bb : the boundary.
		* @param eq : string defining the equation.
		* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
		* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
		*/
		virtual void add_eq_bc (int dom, int bb, const char* eq, int n_cmp = -1, Array<int>** p_cmp=nullptr) ;

		/**
		* Addition of an equation describing a boundary condition.
		* @param dom : number of the \c Domain.
		* @param bb : the boundary.
		* @param eq : string defining the equation.
		* @param list : list of the components to be considered.
		*/
		virtual void add_eq_bc (int dom, int bb, const char* eq, const List_comp& list) ;

		/**
		* Addition of an equation describing a matching condition between two domains (standard setting)
		* @param dom : number of the \c Domain.
		* @param bb : the boundary.
		* @param eq : string defining the equation.
		* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
		* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
		*/
		virtual void add_eq_matching (int dom, int bb, const char* eq, int n_cmp = -1, Array<int>** p_cmp=nullptr) ;

		/**
		* Addition of an equation describing a matching condition between two domains (standard setting)
		* @param dom : number of the \c Domain.
		* @param bb : the boundary.
		* @param eq : string defining the equation.
		* @param list : list of the components to be considered.
		*/
		virtual void add_eq_matching (int dom, int bb, const char* eq, const List_comp& list) ;

		/**
		* Addition of an equation describing a matching condition between two domains (specialized function for time evolution).
		* @param dom : number of the \c Domain.
		* @param bb : the boundary.
		* @param eq : string defining the equation.
		* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
		* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
		*/
		virtual void add_eq_matching_one_side (int dom, int bb, const char* eq, int n_cmp = -1, Array<int>** p_cmp=nullptr) ;

		/**
		* Addition of an equation describing a matching condition between two domains (specialized function for time evolution).
		* @param dom : number of the \c Domain.
		* @param bb : the boundary.
		* @param eq : string defining the equation.
		* @param list : list of the components to be considered.
		*/
		virtual void add_eq_matching_one_side (int dom, int bb, const char* eq, const List_comp& list) ;

		/**
		* Addition of an equation describing a matching condition between domains.
		* The matching is performed in the configuration space.
		* It is intended where the collocations points are different at each side of the boundary.
		* It can happen when there are more than one touching domain (bispheric vs spheric) and when the number of points is different.
		* @param dom : number of the \c Domain.
		* @param bb : the boundary.
		* @param eq : string defining the equation.
		* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
		* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
		*/
		virtual void add_eq_matching_non_std (int dom, int bb, const char* eq, int n_cmp = -1, Array<int>** p_cmp=nullptr) ;

		/**
		* Addition of an equation describing a matching condition between domains.
		* The matching is performed in the configuration space.
		* It is intended where the collocations points are different at each side of the boundary.
		* It can happen when there are more than one touching domain (bispheric vs spheric) and when the number of points is different.
		* @param dom : number of the \c Domain.
		* @param bb : the boundary.
		* @param eq : string defining the equation.
		* @param list : list of the components to be considered.
		*/
		virtual void add_eq_matching_non_std (int dom, int bb, const char* eq, const List_comp& list) ;

		/**
		* Addition of an equation describing a matching condition between domains using the ("import" setting)
		* The matching is performed in the configuration space.
		* It is intended where the collocations points are different at each side of the boundary.
		* It can happen when there are more than one touching domain (bispheric vs spheric) and when the number of points is different.
		* @param dom : number of the \c Domain.
		* @param bb : the boundary.
		* @param eq : string defining the equation.
		* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
		* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
		*/
		virtual void add_eq_matching_import (int dom, int bb, const char* eq, int n_cmp = -1, Array<int>** p_cmp=nullptr) ;

		/**
		* Addition of an equation describing a matching condition between domains using the ("import" setting)
		* The matching is performed in the configuration space.
		* It is intended where the collocations points are different at each side of the boundary.
		* It can happen when there are more than one touching domain (bispheric vs spheric) and when the number of points is different.
		* @param dom : number of the \c Domain.
		* @param bb : the boundary.
		* @param eq : string defining the equation.
		* @param list : list of the components to be considered.
		*/
		virtual void add_eq_matching_import (int dom, int bb, const char* eq, const List_comp& list) ;

		/**
		* Addition of an equation to be solved inside a domain (assumed to be zeroth order i.e. with no derivatives).
		* @param dom : number of the \c Domain.
		* @param eq : string defining the equation.
		* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
		* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
		*/
		virtual void add_eq_full (int dom, const char* eq, int n_cmp = -1, Array<int>** p_cmp=nullptr) ;

		/**
		* Addition of an equation to be solved inside a domain (assumed to be zeroth order i.e. with no derivatives).
		* @param dom : number of the \c Domain.
		* @param eq : string defining the equation.
		* @param list : list of the components to be considered.
		*/
		virtual void add_eq_full (int dom, const char* eq, const List_comp& list) ;

		/**
		* Addition of an equation to be solved inside a domain (assumed to be first order).
		* @param dom : number of the \c Domain.
		* @param eq : string defining the equation.
		* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
		* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
		*/
		virtual void add_eq_one_side (int dom, const char* eq, int n_cmp = -1, Array<int>** p_cmp=nullptr) ;

		/**
		* Addition of an equation to be solved inside a domain (assumed to be first order).
		* @param dom : number of the \c Domain.
		* @param eq : string defining the equation.
		* @param list : list of the components to be considered.
		*/
		virtual void add_eq_one_side (int dom, const char* eq, const List_comp& list) ;

		/**
		* Addition of a matching condition, except for one coefficient where an alternative condition is enforced (highly specialized usage).
		* @param dom : number of the \c Domain.
		* @param eq : string defining the equation.
		* @param par : parameters for the exceptional condition (i.e. which coefficient is concerned basically).
		* @param eq_exception : the excpetionnal equation used.
		* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
		* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
		*/
		virtual void add_eq_matching_exception (int dom, int bb, const char* eq, const Param& par, const char* eq_exception, int n_cmp = -1, Array<int>** p_cmp=nullptr) ;

		/**
		* Addition of a matching condition, except for one coefficient where an alternative condition is enforced (highly specialized usage).
		* @param dom : number of the \c Domain.
		* @param eq : string defining the equation.
		* @param par : parameters for the exceptional condition (i.e. which coefficient is concerned basically).
		* @param eq_exception : the excpetionnal equation used.
		* @param list : list of the components to be considered.
		*/
		virtual void add_eq_matching_exception (int dom, int bb, const char* eq, const Param& par, const char* eq_exception, const List_comp& list) ;

		/**
		* Addition of an equation to be solved inside a domain of arbitrary order.
		* The order can be different for each variable (first order in time and second in \f$r\f$ for instance).
		* @param dom : number of the \c Domain.
		* @param orders : orders of the equation, for each variable.
		* @param eq : string defining the equation.
		* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
		* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
		*/
		virtual void add_eq_order (int dom, const Array<int>& orders, const char* eq, int n_cmp = -1, Array<int>** p_cmp=nullptr) ;

		/**
		* Addition of an equation to be solved inside a domain of arbitrary order.
		* The order can be different for each variable (first order in time and second in \f$r\f$ for instance).
		* @param dom : number of the \c Domain.
		* @param orders : orders of the equation, for each variable.
		* @param eq : string defining the equation.
		* @param list : list of the components to be considered.
		*/
		virtual void add_eq_order (int dom, const Array<int>& orders, const char* eq, const List_comp& list) ;

		/**
		* Addition of an equation for the velocity potential of irrotational binaries.
		* @param dom : number of the \c Domain.
		* @param order : order of the equation.
		* @param eq : string defining the equation.
		* @param const_part : constant par
		*/
		virtual void add_eq_vel_pot (int dom, int order, const char* eq, const char* const_part) ;

		/**
		* Addition of an boundary equation with an exception for \f$l=m=0\f$
		* @param dom : number of the \c Domain.
		* @param bound : boundary index.
		* @param eq : string defining the equation.
		* @param const_part : constant par
		*/
		void add_eq_bc_exception (int dom, int bound, const char* eq, const char* const_part) ;

		/**
		* Addition of an equation a boundary condition of arbitrary orders.
		* The order can be different for each variable. It is irrelevant for the variable corresponding to the boundary.
		* @param dom : number of the \c Domain.
		* @param bb : the boundary.
		* @param orders : orders of the equation, for each variable.
		* @param eq : string defining the equation.
		* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
		* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
		*/
		virtual void add_eq_bc (int dom, int bb, const Array<int>& orders, const char* eq, int n_cmp = -1, Array<int>** p_cmp=nullptr) ;

		/**
		* Addition of an equation a boundary condition of arbitrary orders.
		* The order can be different for each variable. It is irrelevant for the variable corresponding to the boundary.
		* @param dom : number of the \c Domain.
		* @param bb : the boundary.
		* @param orders : orders of the equation, for each variable.
		* @param eq : string defining the equation.
		* @param list : list of the components to be considered.
		*/
		virtual void add_eq_bc (int dom, int bb, const Array<int>& orders, const char* eq, const List_comp& list) ;


		/**
		* Addition of an equation a matching condition of arbitrary orders.
		* The order can be different for each variable. It is irrelevant for the variable corresponding to the boundary.
		* @param dom : number of the \c Domain.
		* @param bb : the boundary.
		* @param orders : orders of the equation, for each variable.
		* @param eq : string defining the equation.
		* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
		* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
		*/
		virtual void add_eq_matching (int dom, int bb, const Array<int>& orders, const char* eq, int n_cmp = -1, Array<int>** p_cmp=nullptr) ;

		/**
		* Addition of an equation a matching condition of arbitrary orders.
		* The order can be different for each variable. It is irrelevant for the variable corresponding to the boundary.
		* @param dom : number of the \c Domain.
		* @param bb : the boundary.
		* @param orders : orders of the equation, for each variable.
		* @param eq : string defining the equation.
		* @param list : list of the components to be considered.
		*/
		virtual void add_eq_matching (int dom, int bb, const Array<int>& orders, const char* eq, const List_comp& list) ;

		/**
		* Addition of an equation representing a first integral.
		* @param dom : number of the \c Domain.
		* @param eq : string defining the equation.
		* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
		* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
		*/
	//	virtual void add_eq_first_integral (int dom, const char* eq, int n_cmp = -1, Array<int>** p_cmp=nullptr) ;

		/**
		* Addition of an equation representing a first integral.
		* @param dommin : index of the first \d Domain
		* @param dommax : index of the last \d Domain
		* @param integ_part : name of the integral quantity
		* @param const_part : equation fixing the value of the integral.
		*/
		virtual void add_eq_first_integral (int dom_min, int dom_max, const char* integ_part, const char* const_part) ;

		/**
		* Addition of an equation prescribing the value of one coefficient of a scalar field, on a given boundary.
		* @param dom : number of the \c Domain.
		* @param bb : the boundary.
		* @param eq : string defining the scalar field.
		* @param pos_cf : which coefficient is used.
		* @param val : the value the coefficient must have.
		*/
		virtual void add_eq_mode (int dom, int bb, const char* eq, const Index& pos_cf, double val) ;
		/**
		* Addition of an equation prescribing the value of one coefficient of a scalar field.
		* @param dom : number of the \c Domain.
		* @param eq : string defining the scalar field.
		* @param pos_cf : which coefficient is used.
		* @param val : the value the coefficient must have.
		*/
		virtual void add_eq_val_mode (int dom, const char* eq, const Index& pos_cf, double val) ;
		/**
		* Addition of an equation saying that the value of a field must be zero at one collocation point.
		* @param dom : number of the \c Domain.
		* @param eq : string defining the scalar field that must vanish.
		* @param pos : which collocation point is used.
		*/
		virtual void add_eq_val (int dom, const char* eq, const Index& pos) ;
		/**
		* Addition of an equation saying that the value of a field must be zero at one point (arbitrary).
		* @param dom : number of the \c Domain.
		* @param eq : string defining the scalar field that must vanish.
		* @param MM : which point is used.
		*/
		virtual void add_eq_point (int dom, const char* eq, const Point& MM) ;

		// Various computational stuff :
		/**
		* Computes the residual of all the equations.
		* This is essentially the value of the biggest coefficient.
		* @returns an array of the error, on all the equations.
		*/
		Array<double> check_equations() ;
		/**
		* Copies the various unknowns (doubles and \c Tensors) into their \c Term_eq counterparts.
		*/
		void vars_to_terms() ;
		virtual void vars_to_terms_impl();

		virtual void compute_nbr_of_conditions();
		/**
		* Computes the second member of the Newton-Raphson equations.
		* It is essentially the coefficients of the residual of all the equations (plus some Galerking issues).
		*/
		Array<double> sec_member() ;
		/**
		* Sets the values the variation of the fields.
		* @param vder : values of all the variations (essentially the coefficients of all the variations of the unknowns).
		*/
		void xx_to_ders(const Array<double>&) ;
		/**
		* Sets the values the  of the fields.
		* @param val : values of all the fields (essentially the coefficients of all the unknowns).
		* @param conte : current position in the \c Array \c val.
		*/
		void xx_to_vars(const Array<double>& val, int& conte) ;
		/**
		* Computes the product \f$ J \times x\f$ where \f$J\f$ is the Jacobian of the system.
		* @param xx : the vector \f$ x\f$.
		* @return the product.
		*/
		Array<double> do_JX (const Array<double>& xx) ;
		/**
		* Computes one column of the Jacobian.
		* @param i : column number.
		* @return the column.
		*/
		Array<double> do_col_J (int i) ;

		/**
		 * Compute some columns of the jacobian of the system of equations. Default arguments are to be used in the case
		 * of a non-MPI version. For the MPI case, each process has to know which columns he has to compute.
		 * @param matrix 2D array storing the result (need to be allocated before being passed).
		 * @param n size of the matrix (which is square).
		 * @param first_col index of the first column to compute.
		 * @param n_col number of columns to compute.
		 * @param num_proc number of process.
		 * @param transpose does the result has to be transposed (this has to be the case for the scalapack linear solve).
		 */
		virtual void compute_matrix_cyclic(Array<double> &matrix, int n, int first_col= 0, int n_col= ALL_COLUMNS, int num_proc = 1,
										   bool transpose = DO_NOT_TRANSPOSE);


		virtual void compute_matrix_adjacent(Array<double> &matrix, int n, int first_col= 0, int n_col= ALL_COLUMNS, int num_proc = 1,
											 bool transpose = DO_NOT_TRANSPOSE, std::vector<std::vector<std::size_t>> *dm = nullptr);

		/**
		* Does one step of the Newton-Raphson iteration.
		* @tparam computational_model template parameter for the computational model (mpi/sequential/GPU) selection.
		* @param prec : required precision.
		* @param error : achieved precision.
		* @return true if the required precision is achieved, false otherwise.
		*/
		template<Computational_model computational_model = default_computational_model>
		bool do_newton (double, double&) ;

		/**
		 * Update the values of \c var and \c var_double from the solution of the linear system of the last Newton
		 * iteration.
		 * @param xx array storing the solution of the linear system from the current Newton iteration.
		 */
		void newton_update_vars(Array<double> const & xx);

		/**
		* Updates the variations of the \c Term_eq that comes from the fact that some \c Domains are variable (i.e. their shape).
		* @param zedoms : the number of all the variable domains.
		*/
		void update_terms_from_variable_domains(const Array<int>& zedoms) ;

		/**
		* Does one step of the Newton-Raphson iteration with a linesearch algorithm.
		* Only implemented in the parallel version.
		* @param prec : required precision.
		* @param error : achieved precision.
		* @param ntrymax : first linesearch parameter.
		* @param stepmax : second linesearch parameter.
		* @return true if the required precision is achieved, false otherwise.
		*/
		template<Computational_model computational_model = default_computational_model>
		bool do_newton_with_linesearch (double precision, double& error, int ntrymax = 10, double stepmax = 1.0);


		// Parts for implementing gmres
		void do_arnoldi (int n, Array<double>& qi, Matrice& Hmat) ;
		void update_gmres (const Array<double>&) ;

		private:
		/**
		* Tests the value of the number of unknowns.
		* Used by do_newton_with_linesearch ; only implemented in parallel version.
		*/
		template<Computational_model computational_model = default_computational_model>
		void check_size_VS_unknowns(int n);
		/**
		* Tests the not too many processors are used.
		* Used by do_newton_with_linesearch ; only implemented in parallel version.
		*/
		template<Computational_model computational_model = default_computational_model>
		void check_bsize(int bsize);
		/**
		* Computes the local part of the Jacobian
		* Used by do_newton_with_linesearch ; only implemented in parallel version.
		*/
		template<Computational_model computational_model = default_computational_model>
		void compute_matloc(Array<double>& matloc_in, int nn, int bsize);
		/**
		* Distributes the second member of Newton-Raphson accross the various processors.
		* Used by do_newton_with_linesearch ; only implemented in parallel version.
		*/
		template<Computational_model computational_model = default_computational_model>
		void translate_second_member(Array<double>& secloc, Array<double> const& second, int nn, int bsize, int nprow, int myrow, int mycol);
		/**
		* Solves the linear problem in Newton-Raphson.
		* Used by do_newton_with_linesearch ; only implemented in parallel version.
		*/
		template<Computational_model computational_model = default_computational_model>
		void get_global_solution(Array<double>& auxi, Array<double> const& secloc, int nn, int bsize, int nprow, int myrow, int mycol);
		/**
		* Update the fields after a Newton-Raphson iteration.
		* Used by do_newton_with_linesearch ; only implemented in parallel version.
		*/
		template<Computational_model computational_model = default_computational_model>
		void update_fields(double lambda, vector<double> const& old_var_double, vector<Tensor> const& old_var_fields, vector<double> const& p_var_double, vector<Tensor> const& p_var_fields);
		/**
		* Update the fields when some domains have been modified.
		* Used by do_newton_with_linesearch ; only implemented in parallel version.
		*/
		template<Computational_model computational_model = default_computational_model>
		void compute_old_and_var(Array<double> const& xx, vector<double>& old_var_double, vector<Tensor>& old_var_fields, vector<double>& p_var_double, vector<Tensor>& p_var_fields);
		/**
		* Inner routine for the linesearch algorithm
		* Used by do_newton_with_linesearch ; only implemented in parallel version.
		*/
		template<Computational_model computational_model = default_computational_model>
		double compute_f(Array<double> const& second);
		/**
		* Inner routine for the linesearch algorithm
		* Used by do_newton_with_linesearch ; only implemented in parallel version.
		*/
		template<Computational_model computational_model = default_computational_model>
		void compute_p(Array<double>& xx, Array<double> const& second, int nn);
		/**
		* Tests the positivity of  \f$\delta\f$
		* Used by do_newton_with_linesearch ; only implemented in parallel version.
		*/
		template<Computational_model computational_model = default_computational_model>
		void check_positive(double delta);
		/**
		* Tests the positivity of  \f$\delta\f$
		* Used by do_newton_with_linesearch ; only implemented in parallel version.
		*/
		template<Computational_model computational_model = default_computational_model>
		void check_negative(double delta);

		friend class Space_spheric ;
		friend class Space_bispheric ;
		friend class Space_critic ;
		friend class Space_polar ;
		friend class Space_spheric_adapted ;
		friend class Space_polar_adapted ;
		friend class Space_bin_ns ;
		friend class Space_bin_bh ;
		friend class Space_bin_fake ;
		friend class Space_polar_periodic ;
		friend class Space_adapted_bh ;
		friend class Space_bbh ;
		friend class Metric ;
		friend class Metric_general ;
		friend class Metric_flat ;
		friend class Metric_dirac ;
		friend class Metric_dirac_const ;
		friend class Metric_conf ;
		friend class Metric_relax ;
		friend class Metric_ADS ;
		friend class Metric_AADS ;
		friend class Metric_const ;
		friend class Metric_flat_nophi ;
		friend class Metric_nophi ;
		friend class Metric_nophi_AADS ;
		friend class Metric_nophi_const ;
		friend class Metric_nophi_AADS_const ;
		friend class Metric_conf_factor ;
		friend class Metric_conf_factor_const ;
		friend class Metric_cfc ;

	protected:
		struct Newton_iteration_data{
			unsigned n_iter; int problem_size; double current_error;
			Duration t_load_matrix; Duration t_trans_matrix; Duration t_inv_matrix; Duration t_newton_update;
		};
		/**
		 * Sends the header lines for the Newton-Raphson iterations data.
		 * @param os stream to send the header to.
		 * @param precision is the tolerance stopping criterion for the algorithm.
		 */
		static void display_do_newton_report_header(std::ostream & os,double precision);
		/**
		 * Sends the last lines for the Newton-Raphson iterations data.
		 * @param os stream to send the end lines to.
		 * @param reached_precision value of the reached precision.
		 */
		static void display_do_newton_ending_line(std::ostream &os, double precision, double reached_precision);
		/**
		 * Format end sends the Newton-Raphson iteration data.
		 * @param os stream to send the data to.
		 * @param data data to display.
		 */
		static void display_do_newton_iteration(std::ostream &os,const Newton_iteration_data & data);
	};

	/**
	 * Class implementing an equation. This is a purely abstract class that can not be instanciated.
	 * \ingroup systems
	 */
	class Equation : public Memory_mapped {
    protected:
        const Domain* dom ; ///< Pointer on the \c Domain where the equation is defined.
        int ndom ; ///< Number of the domain
        int n_ope ; ///< Number of terms involved in the equation (one for bulk, two or more fot matching...).
        MMPtr_array<Ope_eq> parts ; ///< Array of pointers on the various terms.

        /**
        * Indicator checking whther the result has been computed already once.
        * If not the quantities \c n_cond must be computed.
        */
        bool called ;
        int n_comp ; ///< Number of components of the residual (1 for a scalar, 6 for a symmetric rank-2 tensor etc).
        int n_cond_tot ; ///< Total number of discretized equations (essentially the number of all coefficients of the residual).
        Array<int>* n_cond ; ///< Number of discretized equations, component by component.

        int n_cmp_used ; ///< Number of components used (by default the same thing as \c n_comp).
        Array<int>** p_cmp_used ; ///< Array of pointer on the indices of the used components

        /**
         * Constructor
         * @param dom : Pointer on the \d Domain
         * @param nd : index of the \d Domain (consistence is not checked).
         * @param nope : number of operators.
         * @param n_cmp : number of components of \c eq to be considered. All the components are used of it is -1.
         * @param p_cmp : pointer on the indexes of the components to be considered. Not used of nused = -1 .
         */
        Equation (const Domain *dom, int nd, int nope, int n_cmp=-1, Array<int>** p_cmp=nullptr) ;


    public:
       virtual ~Equation() ; ///< Destructor

		/**
		* @return the total number of discretized conditions.
		*/
       int get_n_cond_tot() const {return n_cond_tot;} ;

		/**
		* Computes the terms involved in computing the residual of the equations.
		* @param conte : current position in the array of terms.
		* @param res : array of pointers on the various terms.
		*/
       virtual void apply(int& conte, Term_eq** res) ;
		/**
		* Generates the discretized errors, from the various \c Term_eq computed by the equation.
		*  Basically used when computing the second member of the Newton-Raphson algorithm.
		* @param conte : current position in the array of terms.
		* @param residuals : array of pointers on the various terms.
		* @param sec : array of the discretized errors.
		* @param pos_sec : current position in \c sec.
		*/
       virtual void export_val(int& conte, Term_eq** residuals, Array<double>& sec, int& pos_sec) const = 0 ;
		/**
		* Generates the discretized variations, from the various \c Term_eq computed by the equation.
		*  Basically used when computing the Jacobian of the Newton-Raphson algorithm.
		* @param conte : current position in the array of terms.
		* @param residuals : array of pointers on the various terms.
		* @param sec : array of the discretized errors.
		* @param pos_sec : current position in \c sec.
		*/
       virtual void export_der(int& conte, Term_eq** residuals, Array<double>& sec, int& pos_sec) const = 0 ;
		/**
		* Computes the number of conditions associated with the equation.
		* @param tt : the residual of the equation.
		*/
       virtual Array<int> do_nbr_conditions (const Tensor& tt) const = 0 ;
		/**
		* Check whether the variation of the residual has to be taken into account when computing a given column.
		* @param target : domain involved in the computation of the given column.
		*/
		virtual bool take_into_account (int) const = 0 ;

		friend class System_of_eqs ;
	} ;

	/**
	 * Class for bulk equations that are solved stricly inside a given domain.
	 * By this one means that they have to be given matching/boundary conditions at all their boundaries.
	 * Typically used for second order PDE.
	 * \ingroup systems
	 */
	class Eq_inside : public Equation {

    public:
		/**
		* Constructor
		* @param dom : Pointer on the \d Domain
		* @param nd : number of the \d Domain (consistence is not checked).
		* @param ope : pointer on the operator describing the equation.
		* @param n_cmp : number of components of \c eq to be considered. All the components are used of it is -1.
		* @param p_cmp : pointer on the indexes of the components to be considered. Not used of nused = -1 .
		*/
        Eq_inside(const Domain* dom, int nd, Ope_eq* op, int n_cmp = -1, Array<int>** p_cmp=nullptr) ;
        virtual ~Eq_inside() ; ///< Destructor.

        void export_val(int&, Term_eq**, Array<double>&, int&) const override;
        void export_der(int&, Term_eq**, Array<double>&, int&) const override;
        Array<int> do_nbr_conditions (const Tensor& tt) const override;
        bool take_into_account (int) const override;
	} ;

	/**
	 * Class for the velocity potential in irrotational binray neutron stars.
	 * It is solved for all harmonics byt l=m=0 for which an exceptionnal condition is enforced.
	 * Typically used for second order PDE.
	 * \ingroup systems
	 */
	class Eq_vel_pot : public Equation {

		public:
			int order ; ///< Order of the equation.

		/**
		* Constructor
		* @param dom : Pointer on the \d Domain
		* @param nd : number of the \d Domain (consistence is not checked).
		* @param ord : order of the equation (probably only safe with 0, 1 or 2).
		* @param ope : pointer on the operator describing the equation.
		* @param ope_constant : condition for the constant part
		*/
			Eq_vel_pot(const Domain* dom, int nd, int ord, Ope_eq* op, Ope_eq* op_constant) ;
			virtual ~Eq_vel_pot() ; ///< Destructor.

			void export_val(int&, Term_eq**, Array<double>&, int&) const override;
			void export_der(int&, Term_eq**, Array<double>&, int&) const override;
			Array<int> do_nbr_conditions (const Tensor& tt) const override;
			bool take_into_account (int) const override;
	} ;

	/**
	 * Class for enforcing boundary condition.
	 * It is solved for all harmonics byt l=m=0 for which an exceptionnal condition is enforced.
	 * \ingroup systems
	 */
	class Eq_bc_exception : public Equation {

		public:
			int bound ; ///< Order of the equation.

		/**
		* Constructor
		* @param dom : Pointer on the \d Domain
		* @param nd : number of the \d Domain (consistence is not checked).
		* @param bound : boundary where the condition is enforced
		* @param ope : pointer on the operator describing the equation.
		* @param ope_constant : condition for the constant part
		*/
			Eq_bc_exception (const Domain* dom, int nd, int bound, Ope_eq* op, Ope_eq* op_constant) ;
			virtual ~Eq_bc_exception() ; ///< Destructor.

			void export_val(int&, Term_eq**, Array<double>&, int&) const override;
			void export_der(int&, Term_eq**, Array<double>&, int&) const override;
			Array<int> do_nbr_conditions (const Tensor& tt) const override;
			bool take_into_account (int) const override;
	} ;

	/**
	 * Class for bulk equation which order is passed as a parameter.
	 * \ingroup systems
	 */
	class Eq_order : public Equation {

		public:
			int order ; ///< Order of the equation.

		/**
		* Constructor
		* @param dom : Pointer on the \d Domain
		* @param nd : number of the \d Domain (consistence is not checked).
		* @param ord : order of the equation (probably only safe with 0, 1 or 2).
		* @param ope : pointer on the operator describing the equation.
		* @param n_cmp : number of components of \c eq to be considered. All the components are used of it is -1.
		* @param p_cmp : pointer on the indexes of the components to be considered. Not used of nused = -1 .
		*/

			Eq_order(const Domain* dom, int nd, int ord, Ope_eq* ope, int n_cmp = -1, Array<int>** p_cmp=nullptr) ;
			virtual ~Eq_order() ; ///< Destructor.

			void export_val(int&, Term_eq**, Array<double>&, int&) const override;
			void export_der(int&, Term_eq**, Array<double>&, int&) const override;
			Array<int> do_nbr_conditions (const Tensor& tt) const override;
			bool take_into_account (int) const override;
	} ;

	/**
	 * Class for an equation representing a boundary condition on some surface.
	 * \ingroup systems
	 */
	class Eq_bc : public Equation {

		public:
			int bound ; ///< The boundary
		/**
		* Constructor
		* @param dom : Pointer on the \d Domain
		* @param nd : number of the \d Domain (consistence is not checked).
		* @param bb : boundary where the equation is enforced.
		* @param ope : pointer on the operator describing the equation.
		* @param n_cmp : number of components of \c eq to be considered. All the components are used of it is -1.
		* @param p_cmp : pointer on the indexes of the components to be considered. Not used of nused = -1 .
		*/
			Eq_bc(const Domain* dom, int nd, int bb, Ope_eq* ope, int n_cmp = -1, Array<int>** p_cmp=nullptr) ;
			virtual ~Eq_bc() ; ///< Destructor.

			void export_val(int&, Term_eq**, Array<double>&, int&) const override;
			void export_der(int&, Term_eq**, Array<double>&, int&) const override;
			Array<int> do_nbr_conditions (const Tensor& tt) const override;
			bool take_into_account (int) const override;
	} ;

	/**
	 * Class for an equation representing the matching of quantities accross a boundary.
	 * \ingroup systems
	 */
	class Eq_matching : public Equation {

		public:
			int bound ; ///< Name of the boundary in the domain of the equation.
			int other_dom ; ///< Number of the other domain.
			int other_bound ; ///< Name of the boundary in the other domain.

		/**
		* Constructor
		* @param dom : Pointer on the \d Domain
		* @param nd : number of the \d Domain (consistence is not checked).
		* @param bb : name of the boundary in the domain of the equation.
		* @param oz_nd : number of the other domain.
		* @param oz_bb : name of the boundary in the other domain.
		* @param ope : pointer on the operator describing the quantity to be matched in the current domain.
		* @param oz_ope : pointer on the operator describing the quantity to be matched in the other domain.
		* @param n_cmp : number of components of \c eq to be considered. All the components are used of it is -1.
		* @param p_cmp : pointer on the indexes of the components to be considered. Not used of nused = -1 .
		*/

			Eq_matching(const Domain* dom, int nd, int bb, int oz_nd, int oz_bb, Ope_eq* ope, Ope_eq* oz_ope, int n_cmp = -1, Array<int>** p_cmp=nullptr) ;
			virtual ~Eq_matching() ; ///< Destructor.

			void export_val(int&, Term_eq**, Array<double>&, int&) const override;
			void export_der(int&, Term_eq**, Array<double>&, int&) const override;
			Array<int> do_nbr_conditions (const Tensor& tt) const override;
			bool take_into_account (int) const override;
	} ;

	/**
	 * Class for an equation representing the matching of quantities accross a boundary.
	 * Used when the matching surface also has boundaries.
	 * \ingroup systems
	 */
	class Eq_matching_one_side : public Equation {

		public:
			int bound ; ///< Name of the boundary in the domain of the equation.
			int other_dom ; ///< Number of the other domain.
			int other_bound ; ///< Name of the boundary in the other domain.

		/**
		* Constructor
		* @param dom : Pointer on the \d Domain
		* @param nd : number of the \d Domain (consistence is not checked).
		* @param bb : name of the boundary in the domain of the equation.
		* @param oz_nd : number of the other domain.
		* @param oz_bb : name of the boundary in the other domain.
		* @param ope : pointer on the operator describing the quantity to be matched in the current domain.
		* @param oz_ope : pointer on the operator describing the quantity to be matched in the other domain.
		* @param n_cmp : number of components of \c eq to be considered. All the components are used of it is -1.
		* @param p_cmp : pointer on the indexes of the components to be considered. Not used of nused = -1 .
		*/
		Eq_matching_one_side(const Domain* dom, int nd , int bb, int oz_nd, int oz_bb, Ope_eq* ope, Ope_eq* oz_ope, int n_cmp = -1, Array<int>** p_cmp=nullptr) ;
			virtual ~Eq_matching_one_side() ; ///< Destructor.

			void export_val(int&, Term_eq**, Array<double>&, int&) const override;
			void export_der(int&, Term_eq**, Array<double>&, int&) const override;
			Array<int> do_nbr_conditions (const Tensor& tt) const override;
			bool take_into_account (int) const override;
	} ;

	/**
	 * Class for an equation representing the matching of quantities accross a boundary.
	 * The matching is performed in the configuration space.
	 * It is intended where the collocations points are different at each side of the boundary.
	 * It can happen when there are more than one touching domain (bispheric vs spheric) and when the number of points is different.
	 * \ingroup systems
	 */
	class Eq_matching_non_std : public Equation {

		public :
			int bound ; ///< Name of the boundary in the domain of the equation.
			Array<int> other_doms ; ///< Array containing the number of the domains being on the other side of the surface.
			Array<int> other_bounds ; ///< Names of the boundary, as seen in the other domains.

			Index** which_points ; ///< Lists the collocation points on the boundary (probably...)

		/**
		* Constructor
		* @param dom : Pointer on the \d Domain
		* @param nd : number of the \d Domain (consistence is not checked).
		* @param bb : name of the boundary in the domain of the equation.
		* @param ozers : numbers of the others domains in (0,*) and names of the boundary as seen on the other side in (1,*).
		* @param n_cmp : number of components of \c eq to be considered. All the components are used of it is -1.
		* @param p_cmp : pointer on the indexes of the components to be considered. Not used of nused = -1 .
		*/
			Eq_matching_non_std (const Domain* dom, int nd, int bb, const Array<int>& ozers, int n_cmp = -1, Array<int>** p_cmp=nullptr) ;
			virtual ~Eq_matching_non_std() ; ///< Destructor.

		   void apply(int&, Term_eq**) override;
		   void export_val(int&, Term_eq**, Array<double>&, int&) const override;
		   void export_der(int&, Term_eq**, Array<double>&, int&) const override;

		   /**
		   * Computes the collocation points used
		   * @param base : the spectral base of the residual (in order to get the symmetries).
		   * @param start : starting index
		   */
		   void do_which_points (const Base_spectral& base, int start) ;
		   Array<int> do_nbr_conditions (const Tensor& tt) const override;
		   bool take_into_account (int) const override;

		friend class System_of_eqs ;
	} ;

	/**
	 * Class for an equation representing the matching of quantities accross a boundary using the "import" reserved word.
	 * The matching is performed in the configuration space.
	 * It is intended where the collocations points are different at each side of the boundary.
	 * It can happen when there are more than one touching domain (bispheric vs spheric) and when the number of points is different.
	 * \ingroup systems
	 */
	class Eq_matching_import : public Equation {

		public:
			int bound ; ///< Name of the boundary in the domain of the equation.
			Array<int> other_doms ; ///< Array containing the number of the domains being on the other side of the surface.
			Array<int> other_bounds ; ///< Names of the boundary, as seen in the other domains.


		/**
		* Constructor
		* @param dom : Pointer on the \d Domain
		* @param nd : number of the \d Domain (consistence is not checked).
		* @param bb : name of the boundary in the domain of the equation.
		* @param eq : the matching condition. Should be of the form "a = import(b)".
		* @param ozers : numbers of the others domains in (0,*) and names of the boundary as seen on the other side in (1,*).
		* @param n_cmp : number of components of \c eq to be considered. All the components are used of it is -1.
		* @param p_cmp : pointer on the indexes of the components to be considered. Not used of nused = -1 .
		*/
			Eq_matching_import (const Domain* dom, int nd, int bb, Ope_eq* eq, const Array<int>& ozers, int n_cmp = -1, Array<int>** p_cmp=nullptr) ;
			virtual ~Eq_matching_import() ; ///< Destructor.

			void export_val(int&, Term_eq**, Array<double>&, int&) const override;
			void export_der(int&, Term_eq**, Array<double>&, int&) const override;
			Array<int> do_nbr_conditions (const Tensor& tt) const override;
			bool take_into_account (int) const override;
	} ;

	/**
	 * Class implementing an integral equation.
	 * \ingroup systems
	 */
	class Eq_int {

		protected:
			int n_ope ; ///< Number of terms.
			Ope_eq** parts ; ///< Array of pointers on the various terms.

		public:
			/**
			 * Constructor just sets n_ope.
			 * @param nop : the number of operators involved.
			 */
			Eq_int (int nop) ;

		public:
		   ~Eq_int() ; ///< Destructor

		   double get_val() const ; ///< Return the value of the equation.
		   double get_der() const ; ///< Return the variation of the equation.
		/**
		* Sets one of the \c Ope_eq needed by the equation.
		* @param i : the part to be set.
		* @param ope : pointer on the operator to be used.
		*/
		   void set_part (int, Ope_eq*) ;

		   friend class System_of_eqs ;
	} ;

	/**
	 * Class for a zeroth order equation in a \c Domain.
	 * Should be used for equations without derivatives.
	 * \ingroup systems
	 */
	class Eq_full : public Equation {

		public:
		/**
		* Constructor
		* @param dom : Pointer on the \d Domain
		* @param nd : number of the \d Domain (consistence is not checked).
		* @param ope : pointer on the operator describing the equation.
		* @param n_cmp : number of components of \c eq to be considered. All the components are used of it is -1.
		* @param p_cmp : pointer on the indexes of the components to be considered. Not used of nused = -1 .
		*/
			Eq_full(const Domain* dom, int nd, Ope_eq* ope, int n_cmp = -1, Array<int>** p_cmp=nullptr) ;
			virtual ~Eq_full() ; ///< Destructor.

			void export_val(int&, Term_eq**, Array<double>&, int&) const override;
			void export_der(int&, Term_eq**, Array<double>&, int&) const override;
			Array<int> do_nbr_conditions (const Tensor& tt) const override;
			bool take_into_account (int) const override;
	} ;

	/**
	 * Class for a first order equation in a \c Domain.
	 * \ingroup systems
	 */
	class Eq_one_side : public Equation {

		public:
		/**
		* Constructor
		* @param dom : Pointer on the \d Domain
		* @param nd : number of the \d Domain (consistence is not checked).
		* @param ope : pointer on the operator describing the equation.
		* @param n_cmp : number of components of \c eq to be considered. All the components are used of it is -1.
		* @param p_cmp : pointer on the indexes of the components to be considered. Not used of nused = -1 .
		*/
			Eq_one_side(const Domain* dom, int nd, Ope_eq* ope, int n_cmp = -1, Array<int>** p_cmp=nullptr) ;
			virtual ~Eq_one_side() ; ///< Destructor.

			void export_val(int&, Term_eq**, Array<double>&, int&) const override;
			void export_der(int&, Term_eq**, Array<double>&, int&) const override;
			Array<int> do_nbr_conditions (const Tensor& tt) const override;
			bool take_into_account (int) const override;
	} ;

	/**
	 * Class for an equation in a \c Domain which order is passed, for each variable.
	 * \ingroup systems
	 */
	class Eq_order_array : public Equation {

		public:
			const Array<int>& order ; ///< Orders of the equation wrt each variable.

		/**
		* Constructor
		* @param dom : Pointer on the \d Domain
		* @param nd : number of the \d Domain (consistence is not checked).
		* @param ord : orders of the equation wrt each variable.
		* @param ope : pointer on the operator describing the equation.
		* @param n_cmp : number of components of \c eq to be considered. All the components are used of it is -1.
		* @param p_cmp : pointer on the indexes of the components to be considered. Not used of nused = -1 .
		*/
			Eq_order_array(const Domain* dom, int nd, const Array<int>& ord, Ope_eq* ope, int n_cmp = -1, Array<int>** p_cmp=nullptr) ;
			virtual ~Eq_order_array() ; ///< Destructor.

			void export_val(int&, Term_eq**, Array<double>&, int&) const override;
			void export_der(int&, Term_eq**, Array<double>&, int&) const override;
			Array<int> do_nbr_conditions (const Tensor& tt) const override;
			bool take_into_account (int) const override;
	} ;

	/**
	 * Class for a boundary condition.
	 * The order is specified for each variable. The one corresponding to the boundary is irrelevant.
	 * \ingroup systems
	 */
	class Eq_bc_order_array : public Equation {

		public:
			int bound ; ///< The boundary.
			const Array<int>& order;  ///< Orders of the equation wrt each variable.

		/**
		* Constructor
		* @param dom : Pointer on the \d Domain
		* @param nd : number of the \d Domain (consistence is not checked).
		* @param bb : the boundary.
		* @param ord : orders of the equation wrt each variable.
		* @param ope : pointer on the operator describing the equation.
		* @param n_cmp : number of components of \c eq to be considered. All the components are used of it is -1.
		* @param p_cmp : pointer on the indexes of the components to be considered. Not used of nused = -1 .
		*/
			Eq_bc_order_array(const Domain* dom, int nd, int bb, const Array<int>& ord, Ope_eq* ope, int n_cmp = -1, Array<int>** p_cmp=nullptr) ;
			virtual ~Eq_bc_order_array() ; ///< Destructor.

			void export_val(int&, Term_eq**, Array<double>&, int&) const override;
			void export_der(int&, Term_eq**, Array<double>&, int&) const override;
			Array<int> do_nbr_conditions (const Tensor& tt) const override;
			bool take_into_account (int) const override;
	} ;

	/**
	 * Class for a matching condition.
	 * The order is specified for each variable. The one corresponding to the boundary is irrelevant.
	 * \ingroup systems
	 */
	class Eq_matching_order_array : public Equation {

		public:
			int bound ;  ///< The boundary.
			int other_dom ; ///< Number of the \c Domain on the other side of the boundary.
			int other_bound ; ///< Name of the boundary as seen from the other domain.
			const Array<int>& order ; ///< Orders of the equation wrt each variable.

		/**
		* Constructor
		* @param dom : Pointer on the \d Domain
		* @param nd : number of the \d Domain (consistence is not checked).
		* @param bb : the boundary.
		* @param oz_nd : number of the \c Domain on the other side of the boundary.
		* @param oz_bb : the boundary.
		* @param ord : orders of the equation wrt each variable.
		* @param ope : pointer on the operator.
		* @param oz_ope : pointer on the operator from the other domain.
		* @param n_cmp : number of components of \c eq to be considered. All the components are used of it is -1.
		* @param p_cmp : pointer on the indexes of the components to be considered. Not used of nused = -1 .
		*/
			Eq_matching_order_array(const Domain* dom, int nd, int bb, int oz_nd, int oz_bb, const Array<int>& ord, Ope_eq* ope, Ope_eq* oz_ope, int n_cmp = -1, Array<int>** p_cmp=nullptr) ;
			virtual ~Eq_matching_order_array() ; ///< Destructor.

			void export_val(int&, Term_eq**, Array<double>&, int&) const override;
			void export_der(int&, Term_eq**, Array<double>&, int&) const override;
			Array<int> do_nbr_conditions (const Tensor& tt) const override;
			bool take_into_account (int) const override;
	} ;

	/**
	* Equation for a matching condition, except for one coefficient where an alternative condition is enforced (highly specialized usage).
	* \ingroup systems
	*/
	class Eq_matching_exception : public Equation {

		public:
			int bound ;///< The boundary.
			int other_dom ;///< Number of the \c Domain on the other side of the boundary.
			int other_bound ;///< Name of the boundary as seen from the other domain.
			const Param& parameters ; ///< Parameters needed for describing the exception.

		/**
		* Constructor
		* @param dom : Pointer on the \d Domain
		* @param nd : number of the \d Domain (consistence is not checked).
		* @param bb : the boundary.
		* @param oz_nd : number of the \c Domain on the other side of the boundary.
		* @param oz_bb : the boundary.
		* @param ope : pointer on the operator.
		* @param oz_ope : pointer on the operator from the other domain.
		* @param par : parameters needed for describing the exception.
		* @param ope_exc : operator for the exceptionnal condition.
		* @param n_cmp : number of components of \c eq to be considered. All the components are used of it is -1.
		* @param p_cmp : pointer on the indexes of the components to be considered. Not used of nused = -1 .
		*/
			Eq_matching_exception(const Domain* dom, int nd, int bb, int oz_nd, int oz_bb, Ope_eq* ope, Ope_eq* oz_ope, const Param& par, Ope_eq* ope_exc, int n_cmp = -1, Array<int>** p_cmp=nullptr) ;
			virtual ~Eq_matching_exception() ; ///< Destructor.

			void export_val(int&, Term_eq**, Array<double>&, int&) const override;
			void export_der(int&, Term_eq**, Array<double>&, int&) const override;
			Array<int> do_nbr_conditions (const Tensor& tt) const override;
			bool take_into_account (int) const override;
	} ;

	/**
	 * Equation for describing a first integral equation (i.e. a constant quantity in some domains).
	 * \ingroup systems
	 */
	class Eq_first_integral : public Equation {

		public:

			int dom_min ; ///< Index of the first \c Domain
			int dom_max ; ///< Index of the last \c Domain

		/**
		* Constructor
		* @param syst : Pointer on the associated \c System_of_eqs
		* @param dom : Pointer on the first \d Domain (needed for \c Equation constructor)
		* @param dommin : index of the first \d Domain (consistence is not checked).
		* @param dommax : index of the last \d Domain
		* @param integ_part : name of the integral quantity
		* @param const_part : equation fixing the value of the integral.
		*/
			Eq_first_integral(const System_of_eqs* syst, const Domain* dom, int dommin, int dommax, const char* integ_part, const char* const_part) ;
			virtual ~Eq_first_integral() ; ///< Destructor.

			void export_val(int&, Term_eq**, Array<double>&, int&) const override;
			void export_der(int&, Term_eq**, Array<double>&, int&) const override;
			Array<int> do_nbr_conditions (const Tensor& tt) const override;
			bool take_into_account (int) const override;
	} ;

#ifdef PAR_VERSION
    extern "C" {
    void sl_init_ (int*, int*, int*);
    void blacs_gridexit_ (int*);
    void blacs_gridinfo_ (int*, int*, int*, int*, int*);
    int numroc_ (int*, int*, int*, int*, int*);
    void descinit_ (int*, int*, int*, int*, int* , int*, int*, int*, int*, int*);
    void pdgesv_ (int*, int*, double*, int*, int*, int*, int*, double*, int*, int*, int*, int*);
    int blacs_pnum_ (int*, int*, int*);
    void Cpdgemr2d (int, int, double*, int, int, int*, double*, int, int, int*, int);
    }
#endif

    inline void split_two_d (int target, int& low, int& up) {
        int start = int(sqrt(target));
        low = start;
        up = start;

        bool endloop = (low*up==target) ? true : false;

        bool first = true;
        while (!endloop) {
            if (!first)
                low--;
            else
                first= false;

            up = int(target/low);

            if (low*up==target)
                endloop = true;
        }

    }

    template<> bool System_of_eqs::do_newton<Computational_model::sequential>(double, double &);
    template<> bool System_of_eqs::do_newton<Computational_model::mpi_parallel>(double, double &);
    template<> bool System_of_eqs::do_newton<Computational_model::gpu_sequential>(double, double &);
    template<> bool System_of_eqs::do_newton<Computational_model::gpu_mpi_parallel>(double, double &);
    template<> bool System_of_eqs::do_newton_with_linesearch<Computational_model::mpi_parallel>(double , double &, int , double );

// default : not implemented - the implemention is made through the specializations (see the seq_do_newton.cpp,
// mpi_do_newton.cpp... files).
    template<Computational_model computational_model>
    bool System_of_eqs::do_newton(double, double &)
    {
        cerr << "Error: System_of_eq::do_newton is not implemented for the "
             << computational_model_name(computational_model)
             << " computational model." << endl;
        return true;
    }

    template<Computational_model computational_model>
    bool System_of_eqs::do_newton_with_linesearch(double precision, double &error, int ntrymax, double stepmax)
    {
        cerr << "Error: System_of_eq::do_newton is not implemented for the "
             << computational_model_name(computational_model)
             << " computational model." << endl;
        return true;
    }

    template<Computational_model computational_model>
    void check_size_VS_unknowns(int n)
    {
        cerr << "Error: System_of_eq::check_size_VS_unknowns  is not implemented for the "
             << computational_model_name(computational_model)
             << " computational model." << endl;
    }

    template<Computational_model computational_model>
    void check_bsize(int bsize)
    {
        cerr << "Error: System_of_eq::check_bsize  is not implemented for the "
             << computational_model_name(computational_model)
             << " computational model." << endl;
    }

    template<Computational_model computational_model>
    void compute_matloc(Array<double>& matloc_in, int nn, int bsize)
    {
        cerr << "Error: System_of_eq::compute_matloc  is not implemented for the "
             << computational_model_name(computational_model)
             << " computational model." << endl;
    }

    template<Computational_model computational_model>
    void translate_second_member(Array<double>& secloc, Array<double> const& second, int nn, int bsize, int nprow, int myrow, int mycol)
    {
        cerr << "Error: System_of_eq::translate_second_member  is not implemented for the "
             << computational_model_name(computational_model)
             << " computational model." << endl;
    }

    template<Computational_model computational_model>
    void get_global_solution(Array<double>& auxi, Array<double> const& secloc, int nn, int bsize, int nprow, int myrow, int mycol)
    {
        cerr << "Error: System_of_eq::get_global_solution  is not implemented for the "
             << computational_model_name(computational_model)
             << " computational model." << endl;
    }

    template<Computational_model computational_model>
    void update_fields(double lambda, vector<double> const& old_var_double, vector<Tensor> const& old_var_fields, vector<double> const& p_var_double, vector<Tensor> const& p_var_fields)
    {
        cerr << "Error: System_of_eq::update_fields  is not implemented for the "
             << computational_model_name(computational_model)
             << " computational model." << endl;
    }

    template<Computational_model computational_model>
    void compute_old_and_var(Array<double> const& xx, vector<double>& old_var_double, vector<Tensor>& old_var_fields, vector<double>& p_var_double, vector<Tensor>& p_var_fields)
    {
        cerr << "Error: System_of_eq::compute_old_and_var  is not implemented for the "
             << computational_model_name(computational_model)
             << " computational model." << endl;
    }

    template<Computational_model computational_model>
    double compute_f(Array<double> const& second)
    {
        cerr << "Error: System_of_eq::compute_f  is not implemented for the "
             << computational_model_name(computational_model)
             << " computational model." << endl;
        return double{};
    }

    template<Computational_model computational_model>
    void compute_p(Array<double>& xx, Array<double> const& second, int nn)
    {
        cerr << "Error: System_of_eq::compute_p  is not implemented for the "
             << computational_model_name(computational_model)
             << " computational model." << endl;
    }

    template<Computational_model computational_model>
    void check_positive(double delta)
    {
        cerr << "Error: System_of_eq::check_positive  is not implemented for the "
             << computational_model_name(computational_model)
             << " computational model." << endl;
    }

    template<Computational_model computational_model>
    void check_negative(double delta)
    {
        cerr << "Error: System_of_eq::check_negative  is not implemented for the "
             << computational_model_name(computational_model)
             << " computational model." << endl;
    }
}
#endif
