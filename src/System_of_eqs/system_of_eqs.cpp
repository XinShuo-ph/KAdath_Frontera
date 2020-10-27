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

#include "space.hpp"
#include "system_of_eqs.hpp"
#include "name_tools.hpp"
#include "scalar.hpp"
#include "tensor_impl.hpp"
namespace Kadath {

    template<> Profiled_object_base<std::chrono::duration<double>>::Stat_map
            Profiled_object_base<std::chrono::duration<double>>::statistic_map{};

    std::size_t System_of_eqs::default_block_size {64};
// Constructor

System_of_eqs::System_of_eqs (const Space& sp) : output_stream{&std::cout}, user_stream{},
    espace{sp}, dom_min{0}, dom_max{espace.get_nbr_domains()-1}, ndom{dom_max-dom_min+1},
    nvar_double{0}, var_double{VARMAX}, names_var_double{VARMAX},
    nvar{0}, var{VARMAX}, names_var{VARMAX},
    nterm_double{0}, term_double{VARMAX * ndom}, assoc_var_double{VARMAX*ndom},
    nterm{0}, term{VARMAX * ndom}, assoc_var(VARMAX*ndom),
    ncst{0}, nterm_cst{0}, cst{VARMAX*ndom}, names_cst{VARMAX * ndom},
    ncst_hard{0}, cst_hard{VARMAX}, val_cst_hard{VARMAX},
    ndef{0}, def {VARMAX*ndom} , names_def {VARMAX*ndom},
    ndef_glob{0}, def_glob {VARMAX*ndom}, names_def_glob {VARMAX*ndom},
    nopeuser{0}, opeuser{}, paruser{VARMAX}, names_opeuser{VARMAX},
    nopeuser_bin{0}, opeuser_bin{}, paruser_bin{VARMAX}, names_opeuser_bin{VARMAX},
    met{nullptr}, name_met{nullptr},
    neq_int{0}, eq_int{VARMAX}, neq{0}, eq{VARMAX}, results{VARMAX},
    nbr_unknowns{sp.nbr_unknowns_from_variable_domains()}, nbr_conditions{-1}, which_coef{nullptr}
{
    // we may use the initialize tag in constructors, but it would restart a loop for each data member, better to make
    //only one. That said, this constructor is not likely to be extensively called, so that won't matter much.
	for (int i=0 ; i<VARMAX ; i++) {
		var_double[i] = nullptr ;
		var[i]=nullptr ;
		names_var[i]=nullptr ;
		
		cst_hard[i] = nullptr ;
		names_cst[i] = nullptr ;

		opeuser[i] = nullptr ;
		paruser[i] = nullptr ;
		
		opeuser_bin[i] = nullptr ;
		paruser_bin[i] = nullptr ;
		
		eq_int[i] = nullptr ;
		
		eq[i] = nullptr ;
		results[i] = nullptr ;
    }
	for (int i=0 ; i<VARMAX*ndom ; i++) {
		term[i] = nullptr ;
		cst[i] = nullptr ;
		term_double[i] = nullptr ;
	}
	init_proc_data();
}

// same between two bounds
System_of_eqs::System_of_eqs (const Space& sp, int dmin, int dmax) : output_stream{&std::cout},
    user_stream{}, espace{sp}, dom_min{dmin}, dom_max{dmax},ndom{dom_max-dom_min+1},
    nvar_double{0}, var_double{VARMAX}, names_var_double{VARMAX},
    nvar{0}, var{VARMAX}, names_var{VARMAX},
    nterm_double{0}, term_double{VARMAX * ndom}, assoc_var_double{VARMAX*ndom},
    nterm{0}, term{VARMAX * ndom}, assoc_var(VARMAX*ndom),
    ncst{0}, nterm_cst{0}, cst{VARMAX*ndom}, names_cst{VARMAX * ndom},
    ncst_hard{0}, cst_hard{VARMAX}, val_cst_hard{VARMAX},
    ndef{0}, def {VARMAX*ndom} , names_def {VARMAX*ndom},
    ndef_glob{0}, def_glob {VARMAX*ndom}, names_def_glob {VARMAX*ndom},
    nopeuser{0}, opeuser{}, paruser{VARMAX}, names_opeuser{VARMAX},
    nopeuser_bin{0}, opeuser_bin{}, paruser_bin{VARMAX}, names_opeuser_bin{VARMAX},
    met{nullptr}, name_met{nullptr},
    neq_int{0}, eq_int{VARMAX}, neq{0}, eq{VARMAX}, results{VARMAX},
    nbr_unknowns{sp.nbr_unknowns_from_variable_domains()}, nbr_conditions{-1}, which_coef{nullptr}
{
    // we may use the initialize tag in constructors, but it would restart a loop for each data member, better to make
    //only one. That said, this constructor is not likely to be extensively called, so that won't matter much.
    for (int i=0 ; i<VARMAX ; i++) {
        var_double[i] = nullptr ;
        var[i]=nullptr ;
        names_var[i]=nullptr ;

        cst_hard[i] = nullptr ;
        names_cst[i] = nullptr ;

        opeuser[i] = nullptr ;
        paruser[i] = nullptr ;

        opeuser_bin[i] = nullptr ;
        paruser_bin[i] = nullptr ;

        eq_int[i] = nullptr ;

        eq[i] = nullptr ;
        results[i] = nullptr ;
    }
    for (int i=0 ; i<VARMAX*ndom ; i++) {
        term[i] = nullptr ;
        cst[i] = nullptr ;
        term_double[i] = nullptr ;
    }
    init_proc_data();
}


System_of_eqs::~System_of_eqs() {
//	delete [] var_double ;
	for (int i=0 ; i<nvar_double ; i++) delete [] names_var_double[i] ;
	for (int i=0 ; i<nterm_double ; i++) delete term_double[i] ;
//	delete [] var ;
	for (int i=0 ; i<nvar ; i++) delete [] names_var[i] ;
	for (int i=0 ; i<nterm ; i++) delete term[i] ;
	for (int i=0 ; i<nterm_cst ; i++) safe_delete(cst[i]);
	for (int i=0 ; i<ncst ; i++) delete [] names_cst[i] ;
	for (int i=0 ; i<ncst_hard ; i++) delete cst_hard[i] ;
//	delete [] cst_hard ;
	for (int i=0 ; i<ndef ; i++) delete def[i] ;
//	delete [] def ;
	for (int i=0 ; i<ndef ; i++) delete [] names_def[i] ;
//	delete [] names_def ;
    for (int i=0 ; i<ndef_glob ; i++) delete def_glob[i] ;
//	delete [] def_glob ;
	for (int i=0 ; i<ndef_glob ; i++) delete [] names_def_glob[i] ;
//	delete [] names_def_glob ;
	for (int i=0 ; i<nopeuser ; i++) delete names_opeuser[i] ;
//	delete [] names_opeuser ;
//	delete [] paruser ;
	for (int i=0 ; i<nopeuser_bin ; i++) delete names_opeuser_bin[i] ;
//	delete [] names_opeuser_bin ;
//	delete [] paruser_bin ;
	if (name_met !=nullptr) delete [] name_met ;
	for (int i=0 ; i<neq_int ; i++) delete eq_int[i] ;
//	delete [] eq_int ;
	for (int i=0 ; i<neq ; i++) delete eq[i] ;
//	delete [] eq ;
	for (int i=0 ; i<VARMAX ; i++) if (results[i]!=nullptr)
			delete results[i] ;
    if(!which_coef.empty()) for (int i=0 ; i<nbr_conditions ; i++) delete which_coef[i] ;
}

void System_of_eqs::display_do_newton_report_header(std::ostream &os, double precision)
{
    os <<
       "============================================================================================================================\n"
       "|      |            |       ||b||      |                              Computational Times                                  |\n"
       "| Iter | Syst. Size |   Initial Error  |-----------------------------------------------------------------------------------|\n"
       "|      |            | (tol=" << std::setw(10) << std::setprecision(9) << precision;
    os << ") | Matrix Computation | Matrix Translation |      Linear Solver |      Newton Update |\n"
          "|======|============|==================|====================|====================|====================|====================|\n";
}

void System_of_eqs::display_do_newton_ending_line(std::ostream &os, double precision, double reached_precision)
{
    os << "===================================================================================================="
          "=======================\n";
    os << "Success: tolerance reached with ||b|| = " << reached_precision << " / " << precision << "\n" ;
}

void System_of_eqs::display_do_newton_iteration(std::ostream &os,const Newton_iteration_data &data)
{
    static constexpr int dds {16};
    os << '|';
    os << ' ' << std::setw(4) << data.n_iter << " |";
    os << ' ' << std::setw(10) << data.problem_size << " |";
    os << ' ' << std::setw(dds) << data.current_error << " |";
    os << ' ' << std::setw(dds) << to_seconds(data.t_load_matrix) << " s |";
    os << ' ' << std::setw(dds) << to_milliseconds(data.t_trans_matrix) << "ms |";
    os << ' ' << std::setw(dds) << to_seconds(data.t_inv_matrix) << " s |";
    os << ' ' << std::setw(dds) << to_milliseconds(data.t_newton_update) << "ms |";
    os << "\n";

}

const Metric* System_of_eqs::get_met() const {
	if (met==nullptr) {
		cerr << "Undefined metric in System_of_eqs" << endl ;
		abort() ;
	}
	return met ;
}

Term_eq* System_of_eqs::give_term_double (int which, int dd)  const {

	assert ((which>=0) && (which<nvar_double)) ;
	assert ((dd>=dom_min) && (dd<=dom_max)) ;

	int pos = which*ndom + (dd-dom_min) ;
	return term_double[pos] ;
}

Term_eq* System_of_eqs::give_term (int which, int dd)  const {

	assert ((which>=0) && (which<nvar)) ;
	assert ((dd>=dom_min) && (dd<=dom_max)) ;

	int pos = which*ndom + (dd-dom_min) ;
	return term[pos] ;
}

Term_eq* System_of_eqs::give_cst (int which, int dd) const {

	assert ((which>=0) && (which<ncst)) ;
	assert ((dd>=dom_min) && (dd<=dom_max)) ;
	int pos = which*ndom + (dd-dom_min) ;
	return cst[pos] ;
}

Term_eq* System_of_eqs::give_cst_hard (double xx, int dd) const {
	// Check if it has already been stored :
	for (int i=0 ; i<ncst_hard ; i++)
		if ((fabs(val_cst_hard(i) - xx) <PRECISION) && (cst_hard[i]->get_dom()==dd))
			return cst_hard[i] ;
	
	// Term not found, one needs to put it 
	cst_hard[ncst_hard] = new Term_eq (dd, xx) ;
	cst_hard[ncst_hard]->set_der_zero() ;
	val_cst_hard.set(ncst_hard) = xx ;
	ncst_hard++ ;
	return cst_hard[ncst_hard-1] ;
}

Ope_def* System_of_eqs::give_def (int which) const {

	assert ((which>=0) && (which<ndef)) ;
	return def[which] ;
}

Ope_def_global* System_of_eqs::give_def_glob (int which) const {

	assert ((which>=0) && (which<ndef_glob)) ;
	return def_glob[which] ;
}

void System_of_eqs::vars_to_terms () {

	for (int d=dom_min ; d<=dom_max ; d++)
	  espace.get_domain(d)->vars_to_terms() ;

	vars_to_terms_impl();
}

void System_of_eqs::vars_to_terms_impl()
{
    for (int vv=0 ; vv<nvar_double ; vv++)
        for (int dd=dom_min ; dd<=dom_max ; dd++)
            give_term_double(vv, dd)->set_val_d (*var_double[vv]) ;

    for (int vv=0 ; vv<nvar ; vv++)
        for (int dd=dom_min ; dd<=dom_max ; dd++)
            give_term(vv, dd)->set_val_t (*var[vv]) ;
}

// Add an uknown
void System_of_eqs::add_var (const char*nom, double& vv) {

	if ((nbr_char(nom, '_')!=0) || (nbr_char(nom, '^')!=0)) {
		cerr << "No indices allowed in names of variables" << endl ;
		abort() ;
	}
	
	var_double[nvar_double] = &vv ;
	names_var_double[nvar_double] = new char[LMAX];
	trim_spaces(names_var_double[nvar_double],nom) ;

	for (int dd=dom_min ; dd<=dom_max ; dd++) {
		term_double[nterm_double] = new Term_eq(dd, vv) ;
		assoc_var_double.set(nterm_double) = nvar_double ;
		nterm_double++ ;
	}
	nvar_double++ ;

	nbr_unknowns ++ ;
}

// Add an uknown
void System_of_eqs::add_var (const char*nom, Tensor& tt) {

	if (nom!=nullptr) {
		if ((nbr_char(nom, '_')!=0) || (nbr_char(nom, '^')!=0)) {
			cerr << "No indices allowed in names of variables" << endl ;
			abort() ;
		}
	

	names_var[nvar] = new char[LMAX];
	trim_spaces(names_var[nvar],nom) ;
	}
 	else {
		names_var[nvar] = new char[LMAX] ;
		trim_spaces (names_var[nvar], "$") ;
	}

	assert (&tt.espace==&espace) ;
	var[nvar] = &tt ;
	for (int dd=dom_min ; dd<=dom_max ; dd++) {
		term[nterm] = new Term_eq(dd, tt) ;
		assoc_var.set(nterm) = nvar ;
		nterm++ ;
	}

	nvar++ ;
	
	for (int dd=dom_min ; dd<=dom_max ; dd++) 
	    nbr_unknowns += espace.get_domain(dd)->nbr_unknowns (tt, dd) ;
}

// Add a  constant
void System_of_eqs::add_cst (const char*nom, const Tensor& so) {

	if (nom!=nullptr) {
	if ((nbr_char(nom, '_')!=0) || (nbr_char(nom, '^')!=0)) {
		cerr << "No indices allowed in names of constants" << endl ;
		abort() ;
	}

	char nomcst[LMAX];
	trim_spaces(nomcst,nom) ;

	names_cst[ncst] = new char[LMAX] ;
	trim_spaces (names_cst[ncst], nomcst) ;
	}
	else {
	  names_cst[ncst] = new char[LMAX] ;
	  trim_spaces (names_cst[ncst], "$") ;
	}
	
	ncst++ ;

	for (int dd=dom_min ; dd<=dom_max ; dd++) {
		cst[nterm_cst] = new Term_eq(dd, so) ;
		cst[nterm_cst]->set_der_zero() ;
		nterm_cst ++ ;
	}
}

void System_of_eqs::add_cst (const char*nom, double xx) {
	char  nomcst[LMAX];
	trim_spaces(nomcst,nom) ;

	names_cst[ncst] = new char[LMAX] ;
	trim_spaces (names_cst[ncst], nomcst) ;
	ncst++ ;

	for (int dd=dom_min ; dd<=dom_max ; dd++) {
		cst[nterm_cst] = new Term_eq(dd, xx) ;
		cst[nterm_cst]->set_der_zero() ;
		nterm_cst++ ;
	}
}

void System_of_eqs::add_ope (const char* nom, Term_eq (*pope) (const Term_eq&, Param*), Param* par) {
   
    names_opeuser[nopeuser] = new char [LMAX] ;
    trim_spaces (names_opeuser[nopeuser], nom) ;
    opeuser[nopeuser] = pope ;
    paruser[nopeuser] = par ;
    nopeuser++ ;
}

void System_of_eqs::add_ope (const char* nom, Term_eq (*pope) (const Term_eq&, const Term_eq&, Param*), Param* par) {

    names_opeuser_bin[nopeuser_bin] = new char [LMAX] ;
    trim_spaces (names_opeuser_bin[nopeuser_bin], nom) ;
    opeuser_bin[nopeuser_bin] = pope ;
    paruser_bin[nopeuser_bin] = par ;
    nopeuser_bin++ ;
}

// Add a  definition
void System_of_eqs::add_def (int dom, const char* nom) {

	// Récupère le = :
	char p1[LMAX] ;
	char p2[LMAX] ;
	bool indic = is_ope_bin(nom, p1, p2, '=') ;
	if (!indic) {
		cerr << "= needed for definitions" << endl ;
		abort() ;
	}

	//lhs is name of definition :
	names_def[ndef] = new char[LMAX] ;
	get_util (names_def[ndef], p1) ;
	
	int valence ;
	char* indices = nullptr ;
	Array<int>* ttype = nullptr ;
	bool auxi = is_tensor (p1, names_def[ndef], valence, indices, ttype) ;
	assert (auxi) ;
	
	def[ndef] = new Ope_def (this, give_ope (dom, p2), valence, indices, ttype) ;
	if (ttype!=nullptr)
		delete ttype ;
	if (indices!=nullptr)
		delete [] indices ;
	ndef++ ;
}

void System_of_eqs::add_def (const char* nom) {
  for (int dd=dom_min ; dd<=dom_max ; dd++)
    add_def(dd, nom) ;
}

// Add a definition invloving all the domains (i.e. integral)
void System_of_eqs::add_def_global (int dom, const char* nom) {

	// Récupère le = :
	char p1[LMAX] ;
	char p2[LMAX] ;
	bool indic = is_ope_bin(nom, p1, p2, '=') ;
	if (!indic) {
		cerr << "= needed for definitions" << endl ;
		abort() ;
	}
      
	//lhs is name of definition :
	names_def_glob[ndef_glob] = new char[LMAX] ;
	get_util (names_def_glob[ndef_glob], p1) ;
	def_glob[ndef_glob] = new Ope_def_global (this, dom, p2) ;
	ndef_glob++ ;
}

void System_of_eqs::add_def_global (const char* nom) {
  for (int dd=dom_min ; dd<=dom_max ; dd++)
    add_def_global(dd, nom) ;
}


void System_of_eqs::xx_to_ders (const Array<double>& xx) {

	assert (xx.get_ndim()==1) ;
	assert (xx.get_size(0) == nbr_unknowns) ;

	int conte = 0 ;
	int pos_term = 0 ;
	
	espace.xx_to_ders_variable_domains(xx, conte) ;
	
	for (int i=0 ; i<nvar_double ; i++) {
		for (int dd=dom_min ; dd<dom_max ; dd++) {
			term_double[pos_term]->set_der_d(xx(conte)) ;
			pos_term ++ ;
		}
		conte ++ ;
	}

	for (int i=0 ; i<nterm ; i++) {

		// Check for symetric tensor : 
		const Metric_tensor* pmet = dynamic_cast<const Metric_tensor*>(term[i]->val_t) ;
		if (pmet==nullptr) {
			Tensor auxi (term[i]->get_val_t(), false) ;
	
			int dom = term[i]->get_dom() ;
			espace.get_domain(dom)->affecte_tau(auxi, dom, xx, conte) ;
			
			term[i]->set_der_t(auxi) ;
			}
		else  {
			Metric_tensor auxi (*pmet, false) ;
	
			int dom = term[i]->get_dom() ;
			espace.get_domain(dom)->affecte_tau(auxi, dom, xx, conte) ;
			term[i]->set_der_t(auxi) ;
		}
	}
}

void System_of_eqs::xx_to_vars (const Array<double>& xx, int& conte) {

	assert (xx.get_ndim()==1) ;
	assert (xx.get_size(0) == nbr_unknowns) ;
	
	for (int i=0 ; i<nvar_double ; i++) {
		*var_double[i] = xx(conte) ; 
		conte ++ ;
	}
	
	for (int i=0 ; i<nvar ; i++) 
	  for (int dom=dom_min ; dom<=dom_max ; dom++) {
		espace.get_domain(dom)->affecte_tau(*var[i], dom, xx, conte) ;
	}
	
}

Tensor System_of_eqs::give_val_def (const char* so) const {
  
  char name[LMAX] ;
  trim_spaces(name, so) ;
  
  bool found = false ;
  int valence = -1;
  Array<int>* index ;
  Base_tensor* basis ;
  
  for (int i=0 ; i<ndef ; i++) 
    if (!found) {
      if (strcmp(names_def[i], name)==0) {
	found = true ;
	Tensor auxi (give_def(i)->get_res()->get_val_t()) ;
	valence = auxi.get_valence() ;
	index = new Array<int> (auxi.get_index_type()) ;
	basis = new Base_tensor (auxi.get_basis()) ;
      }
    }
    
  if (!found) {
    cerr << "Definition not found...." << endl ;
    abort() ;
  }
  
  Tensor res (espace, valence, *index, *basis) ;
  res = 0 ;
  
  delete index ;
  delete basis ;
  
  for (int i=0 ; i<ndef ; i++) 
    if (strcmp(names_def[i], name)==0) {
      int zedom = give_def(i)->get_res()->get_dom() ;
      Tensor auxi (give_def(i)->get_res()->get_val_t()) ;
      for (int n=0 ; n<res.get_n_comp() ; n++) {
	Array<int> ind (res.indices(n)) ;
	res.set(ind).set_domain(zedom) = auxi(ind)(zedom) ;
      }
	res.set_basis(zedom) = auxi.get_basis().get_basis(zedom) ;
    }

  return res ; 
}

void System_of_eqs::compute_matrix_cyclic(Array<double> &matrix, int n, int first_col, int n_col, int num_proc,
                                          bool transpose)
{
    assert(matrix.get_ndim()==2);
    n_col = (n_col == ALL_COLUMNS ? n : n_col);
    bool done = false;
    int current = 0;
    bool unable_to_compute{false};
    while (!done)
    {
        for (int i{0},k{first_col} ; i<n_col  && k < n; i++, k++) {
                Array<double> column{n};
                try {
                    column = do_col_J(k);
                }
                catch(std::exception const & e) {
                    std::cerr << "---> unable to compute column " << k << "/" << n << " (rank "
                             << mpi_proc_rank << "/" << mpi_world_size << ")" << std::endl;
                    unable_to_compute = true;
                    column = 0.;
                }
                if(transpose){
                    for (int j = 0; j < n; j++) matrix.set(current, j) = column(j);
                } else {
                    for (int j = 0; j < n; j++) matrix.set(j, current) = column(j);
                }
                current++;
        }
        first_col += num_proc * n_col;
        if (first_col>=n)
            done = true;
    }
    if(unable_to_compute) throw std::runtime_error{"Unable to compute the jacobian"};
}

void System_of_eqs::compute_matrix_adjacent(Array<double> &matrix, int n, int first_col, int n_col, int num_proc,
                                            bool transpose, std::vector<vector<std::size_t>> *dm)
{
    assert(matrix.get_ndim()==2);
    n_col = (n_col == ALL_COLUMNS ? n : n_col);
    for(int j{0};j<n_col;j++)
    {
        int const J{first_col+j};
        Array<double> column(do_col_J(J));
        if(transpose)
        {
            for(int i{0};i<n;i++) {matrix.set(J,i) = column(i); if(dm) (*dm)[i][J]++;}
        } else {
            for(int i{0};i<n;i++) {matrix.set(i,J) = column(i); if(dm) (*dm)[J][i]++;}
        }
    }

}

void System_of_eqs::newton_update_vars(Array<double> const &xx)
{
    int conte = 0;
    espace.xx_to_vars_variable_domains (this, xx, conte);

    double* old_var_double = (nvar_double==0) ? nullptr :new double[nvar_double];
    for (int i=0 ; i<nvar_double ; i++) old_var_double[i] = *var_double[i];
    Tensor** old_fields = new Tensor* [nvar];
    for (int i=0 ; i<nvar ; i++) old_fields[i] = new Tensor(*var[i]);
    xx_to_vars(xx, conte);
    for (int i=0 ; i<nvar ; i++)
    {
        *var[i] = *old_fields[i] - *var[i];
    }
    for (int i=0 ; i<nvar_double ; i++) *var_double[i] = old_var_double[i] - *var_double[i];
    if (old_var_double!=nullptr) delete [] old_var_double;
    for (int i=0 ; i<nvar ; i++) delete old_fields[i];
    delete [] old_fields;
}

}
