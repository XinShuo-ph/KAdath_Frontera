#include "kadath_polar.hpp"


using namespace Kadath ;

int main(int argc, char** argv) {

  
	// Read the input data file
	int kk ;
	double omega ;
  
	FILE* fin = fopen (argv[1], "r") ;
	Space_polar space(fin) ;
	fread_be (&kk, sizeof(int), 1, fin) ;
	fread_be (&omega, sizeof(double), 1, fin) ;
	Scalar nu (space,fin) ;
	Scalar incA (space, fin) ;
	Scalar incB (space, fin) ;
	Scalar incNphi (space, fin) ;
	Scalar phi (space, fin) ;
	fclose(fin) ;
	phi.affect_parameters() ;
        phi.set_parameters()->set_m_quant() = kk ;
	
	int ndom = space.get_nbr_domains() ;
	int nr = space.get_domain(0)->get_nbr_coefs()(0) ;
	
	// Compute the "true" fields
	Scalar lapse (exp(nu)) ;
        lapse.std_base() ;
	
        Scalar bigA (exp(incA-nu)) ;
        bigA.std_base() ;

        Scalar bigB ((incB.div_rsint()+1)/lapse) ;
        bigB.std_base() ;
	
	// Check diff A and B
	int nbr = 200 ;
	double rmax = 20 ;
	double error_axe = 0 ;
	Point MM(2) ;
	for (int i=0 ; i<nbr ; i++) {
		double current = fabs(bigA.val_point(MM) - bigB.val_point(MM)) ;
		if (current>error_axe)
			error_axe = current ;
		MM.set(2) += rmax / nbr ;
	}

      	Scalar Np(incNphi.div_rsint()) ;

	// Madm
	Val_domain integadm (bigA(ndom-1).der_r()) ;
	double Madm = -space.get_domain(ndom-1)->integ(integadm, OUTER_BC) / 4 / M_PI ;
	
	// Mk 
	Val_domain integkomar (lapse(ndom-1).der_r()) ;
	double Mkomar = space.get_domain(ndom-1)->integ(integkomar, OUTER_BC) / 4 / M_PI ;
	
	
	// Get the momentum as integral at infinity :
	Val_domain auxiJs (Np(ndom-1).der_r_rtwo()) ;
	Val_domain integJs (auxiJs.mult_sin_theta().mult_sin_theta()) ;
	double Js = -space.get_domain(ndom-1)->integ(integJs, OUTER_BC) / 16 / M_PI ;
	
	// Moment as volume integral
	Scalar auxiJv (kk*(omega-Np*kk)*phi*phi/lapse*bigA*bigA*bigB) ;
	Scalar integJv (space) ;
	for (int d=0 ; d<ndom ; d++)
	  integJv.set_domain(d) = auxiJv(d).mult_r().mult_sin_theta() ;
	double Jv = 0 ;
	for (int d=0 ; d<ndom ; d++)
	  Jv += space.get_domain(d)->integ_volume(integJv(d))  ;
	Jv *= 2*M_PI ;
	

	cout << "Boson star with k = " << kk << " and omega = " << omega << endl ;
	cout << "Madm   = " << Madm << endl ;
	cout << "Mkomar = " << Mkomar << endl ; 
	cout << "Js     = " << Js << endl ;
	cout  << "Jv    = " << Jv << endl ;
	cout << "diff Komar ADM = " << fabs(Madm - Mkomar)/fabs(Madm+Mkomar) << endl ;
	cout << "diff Js and Jv = " << fabs(Js - Jv) / fabs(Js+Jv) << endl ;
	
	return EXIT_SUCCESS ;
}




