#include "kadath_spheric.hpp"

using namespace Kadath ;

int main(int argc, char** argv) {
		
	if (argc <2) {
       cout <<"File missing..." << endl ;
        abort() ;
    }

    char* name_fich = argv[1] ;
    FILE* fich = fopen(name_fich, "r") ;

    Space_spheric espace (fich) ;
    double n0 ;
    fread_be (&n0, sizeof(double), 1, fich) ;   
    double ome ;
    fread_be (&ome, sizeof(double), 1, fich) ;       
    double aa ;
    fread_be (&aa, sizeof(double), 1, fich) ;   
   
    Scalar conf (espace, fich) ;
    Scalar lapse (espace, fich) ;
    Tensor shift_lu (espace, fich) ;
    Vector shift(shift_lu) ;
    Metric_tensor gmet (espace, fich) ;
    fclose(fich) ;
    
    Base_tensor basis (shift.get_basis()) ;
    
    
    // Check system directly :
    int ndom = espace.get_nbr_domains() ;
    Metric_dirac met(gmet) ;
    Metric_tensor gfixed (espace, CON, basis) ;
    for (int i=1 ; i<=3 ; i++)
      for (int j=i ; j<=3 ; j++)
	if (i==j)
	    gfixed.set(i,j) = 1. ;
	else
	    gfixed.set(i,j).annule_hard() ;
    gfixed.std_base() ;

    Vector scov (espace, COV, basis) ;
    scov.set(1) = 1. ;
    scov.set(2) = 0. ;
    scov.set(3) = 0. ;
    scov.std_base() ;

    Vector mm (espace, CON, basis) ;
    for (int i=1 ; i<=3 ; i++)
	mm.set(i) = 0. ;
    Val_domain xx (espace.get_domain(1)->get_cart(1)) ;
    Val_domain yy (espace.get_domain(1)->get_cart(2)) ;
    mm.set(3).set_domain(1) = sqrt(xx*xx + yy*yy) ;
    mm.std_base() ;
   
    // Solve the equation in space outside the nucleus
	System_of_eqs syst (espace, 1, ndom-1) ;
	// Unknowns
	syst.add_var ("P", conf) ;
	syst.add_var ("N", lapse) ;
	syst.add_var ("bet", shift) ;
	met.set_system (syst, "g") ;

	// User defined constants
	syst.add_cst ("a", aa) ;
	syst.add_cst ("m", mm) ;
	syst.add_cst ("s", scov) ;	
	syst.add_cst ("n0", n0) ;
	syst.add_cst ("Ome", ome) ;
	syst.add_cst ("gf", gfixed) ;
	// definitions
	// For speed one stores derivatives of the CF fields :
	syst.add_def ("DN_i = D_i N") ;
	syst.add_def ("DP_i = D_i P") ;
	syst.add_def ("Dbet^ij = D^i bet^j") ;

	syst.add_def ("st^i = s^i / sqrt(s_i * s^i)") ;
	syst.add_def ("A^ij = (Dbet^ij + Dbet^ji - 2. / 3.* Dbet_k^k * g^ij)/2. / N ") ;
	syst.add_def ("LieK_ij = P^4 * (4 * A_ij * bet^k * DP_k / P + bet^k * D_k A_ij + A_ik * Dbet_j^k + A_jk * Dbet_i^k)") ;
	syst.add_def ("DDN_ij = D_i DN_j - 2 * DN_i * DP_j / P - 2 * DN_j * DP_i / P + 2 * g_ij * DN_k * DP^k / P") ;
        syst.add_def ("PartR_ij = R_ij + 6 * DP_i * DP_j / P^2 - 2 * D_i DP_j / P - 2 * g_ij * D_k DP^k / P - 2 * g_ij * DP_k * DP^k / P^2 - P^4 * 2 * A_ik * A_j^k") ;

	espace.add_inner_bc (syst, "N=n0") ; 
	espace.add_inner_bc (syst, "bet^i = n0 / P^2 * st^i + Ome * m^i * a") ;
	espace.add_inner_bc (syst, "4 * st^i * D_i P / P + D_i st^i + P^2 * A_ij * st^i * st^j = 0") ;
	
	// CFC Equations :
	espace.add_eq_full (syst, " D_i DN^i + 2 * DP_i * DN^i / P - N * P^4 * A_ij * A^ij = 0") ; 
	espace.add_eq_full (syst, "R - 8 * D_i DP^i / P - P^4 * A_ij * A^ij = 0") ;
	espace.add_eq_full (syst, "D^j A_ij + 6 * A_ij * DP^j / P =0") ;

	// Evolution
	espace.add_eq_full (syst, "DDN_ij - N * PartR_ij -LieK_ij =0") ;
	
	// Forgotten trace eq :
	espace.add_eq_full (syst, "D_k bet^k + 6 * bet^k * D_k P / P = 0") ;

	// Dirac gauge :
	espace.add_eq_full (syst, "dirac^i = 0") ;

	// Determinant
	espace.add_eq_full (syst, "determinant(g^ij) = 1") ;

	// Outer BC
	espace.add_outer_bc (syst, "N=1") ;
	espace.add_outer_bc (syst, "P=1") ;
	espace.add_outer_bc (syst, "bet^i=0") ;
	espace.add_outer_bc (syst, "g^ij=gf^ij") ;

	// Get Aij :
	Tensor aij (syst.give_def(5*(ndom-1)-1)->get_res()->get_val_t()) ;
	Val_domain a13 (aij(1,3)(ndom-1)) ;
	Val_domain integ_j (espace.get_domain(ndom-1)->mult_r(a13.mult_sin_theta())) ;
	double jadm = - espace.get_domain(ndom-1)->integ(integ_j, OUTER_BC)/8/M_PI ;

	// Computation adm mass :
	Val_domain integ_adm (conf(ndom-1).der_r()) ;
	double adm = -espace.get_domain(ndom-1)->integ(integ_adm, OUTER_BC)/2/M_PI ;

	// Computation of Komar :
	Val_domain integ_komar (lapse(ndom-1).der_r()) ;
	double komar = espace.get_domain(ndom-1)->integ(integ_komar, OUTER_BC)/4/M_PI ;


	Array<double> errors = syst.check_equations() ;
	double error_max = 0 ;
	for (int i=0 ; i<errors.get_size(0) ; i++)
	    if (fabs(errors(i))>error_max)
		error_max = fabs(errors(i)) ;
	cout << "Error max on the eqs " << endl ;
	cout << jadm/adm/adm << " " << error_max << endl ;
	cout << "Virial error " << endl ;
	cout << jadm/adm/adm << " " << fabs(adm-komar)/adm << endl ;
    
    Point M (3) ;
    double ta = 4 ;
    
    des_coupe(lapse, M, 1, -ta, ta, 3, -ta, ta, "N", "x", "z") ;
    des_coupe(conf, M, 1, -ta, ta, 3, -ta, ta, "P", "x", "z") ;
    for (int i=1 ; i<=3 ; i++)
        des_coupe(shift(1), M, 1, -ta, ta, 3, -ta, ta, "N\\ur", "x", "z") ;
    
    for (int i=1 ; i<=3 ; i++) {
	gmet.set(i,i) = gmet(i,i) - 1 ;
	gmet.set(i,i).set_domain(0) = 0 ;
    }
   gmet.std_base() ;
   des_coupe (gmet(1,1), M, 1, -ta, ta, 3, -ta, ta, "\\gg\\urr", "x", "z") ;


    return EXIT_SUCCESS ;
}

