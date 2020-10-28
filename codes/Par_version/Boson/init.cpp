#include "kadath_polar.hpp"
#include "mpi.h"
#include "magma_interface.hpp"

using namespace Kadath ;

int main(int argc, char** argv) {

	int rc = MPI_Init(&argc, &argv) ;
	int rank = 0 ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;

#ifdef ENABLE_GPU_USE
	if(rank==0)
	{
		TESTING_CHECK(magma_init());
		magma_print_environment();
	}	
#endif

    // The new config :
    int dim = 2 ;
    // Number of points
    int type_coloc = CHEB_TYPE ;
    int resol = 13 ;
    Dim_array res (dim) ;
    res.set(0) = resol ; res.set(1) = resol ;

    Point center (2) ;
    for (int i=1 ; i<=dim ; i++)
    center.set(i) = 0 ;

    // Domains
    int ndom = 8 ;
    Array<double> bounds(ndom-1) ;
    bounds.set(0) = 1 ; bounds.set(1) = 2  ; bounds.set(2) = 4 ; bounds.set(3) = 8 ;
    bounds.set(4) = 16 ; bounds.set(5) = 32 ; bounds.set(6) = 48  ;
    Space_polar space(type_coloc, center, res, bounds) ;

    // Constants :
    int kk = 1 ;
    double omega = 0.99 ;
    double qpi = 4*M_PI ;

    // Fields :
    Scalar nu (space) ;
    nu.annule_hard() ;
    nu.std_base() ;

    Scalar incA (space) ;
    incA.annule_hard() ;
    incA.std_base() ;

    // Constant field
    Scalar one (space) ;
    one = 1 ;
    one.std_base() ;
    Scalar rsint (space) ;
    for (int d=0 ; d<ndom-1 ; d++)
    rsint.set_domain(d) = space.get_domain(d)->mult_r(space.get_domain(d)->mult_sin_theta(one(d))) ;
    rsint.set_domain(ndom-1) = space.get_domain(ndom-1)->mult_sin_theta(one(ndom-1)) ;

    double posmax = 20 ;
    double fmax = 0.002 ;
    double sigmax = 2*posmax*posmax/double(kk) ;
    double sigmaz = sigmax/4 ;
    double vmax = fmax / exp(0) / pow(posmax, kk) / exp(-posmax*posmax/sigmax) ;

    Scalar phi(space) ;
    for (int d=0 ; d<ndom-1 ; d++)
    phi.set_domain(d) = vmax * pow(rsint(d), kk) * exp(-pow(space.get_domain(d)->get_cart(1), 2)/sigmax) *  exp(-pow(space.get_domain(d)->get_cart(2),2)/ sigmaz);
    phi.set_domain(ndom-1) = 0 ;
    phi.set_parameters().set_m_quant() = kk ;

    phi.set_domain(ndom-1) = 0 ;
    phi.std_base() ;

    Scalar incB(rsint) ;
    incB.annule_hard() ;

    Scalar incNphi (rsint) ;
    incNphi.annule_hard() ;
     
    {
        Point Mc (2) ;
        Mc.set(1) = posmax ;
        double val = fmax ;

        // System of eqs
        System_of_eqs syst (space, 0, ndom-1) ;

        // variables
        syst.add_var ("nu", nu) ;
        syst.add_var ("phi", phi) ;
        syst.add_var ("ome", omega) ;

        // cst
        syst.add_cst ("rsint", rsint) ;
        syst.add_cst ("k", kk) ;
        syst.add_cst ("qpi", qpi) ;
        syst.add_cst ("val", val) ;

        // Definitions
        syst.add_def ("phisurrsint = divrsint(phi)") ;
        syst.add_def ("N = exp(nu)") ;
        syst.add_def ("E = 0.5* ome*ome/N^2*phi^2 + 0.5*scal(grad(phi),grad(phi)) + 0.5*phi^2 + 0.5*k*k*phisurrsint*phisurrsint") ;
        syst.add_def ("Pp = k/N* ome*phi^2") ;
        syst.add_def ("S = -0.5*scal(grad(phi), grad(phi)) - 0.5*k*k*phisurrsint*phisurrsint + 1.5*ome*ome/N^2 * phi^2- 1.5*phi^2") ;
        syst.add_def ("Spp = 0.5*ome*ome/N^2*phi^2 -0.5* scal(grad(phi),grad(phi))-0.5* phi^2 + 0.5* k*k*phisurrsint*phisurrsint") ;

        // Integral equation :
        space.add_eq_point (syst, Mc, "phi - val") ;

        // Equations
        // Equations
        for (int d=0 ; d<ndom-1 ; d++)
        syst.add_def (d, "eqnu = lap(nu) + scal(grad(nu), grad(nu)) - qpi*(E+S)") ;
        syst.add_def (ndom-1, "eqnu = lap(nu) + scal(grad(nu), grad(nu)) - qpi *(E+S)") ;

        for (int d=0 ; d<ndom-1 ; d++)
        syst.add_def (d, "eqphi = lap(phi) - (1-ome*ome/N^2)*phi + scal(grad(phi),grad(nu))") ;
        syst.add_def (ndom-1, "eqphi = lap(phi) - (1-ome*ome/N^2)*phi + scal(grad(phi),grad(nu))") ;


        space.add_eq (syst, "eqnu=0", "nu", "dn(nu)") ;
        space.add_eq (syst, "eqphi=0", "phi", "dn(phi)") ;


        syst.add_eq_bc (ndom-1, OUTER_BC, "nu=0") ;
        syst.add_eq_bc (ndom-1, OUTER_BC, "phi=0") ;

        double conv ;
        bool endloop = false ;
        int ite = 1 ;
        while (!endloop) {
            endloop = syst.do_newton(1e-6, conv, System_of_eqs::output_enabled) ;
            ite++ ;
        }

        syst.finalize_profiling();
        if (rank == 0)
            profiling_report(syst, std::cout);
    }

    {
        Point Mc(2);
        Mc.set(1) = posmax;
        double val = fmax;

        // System of eqs
        System_of_eqs syst(space, 0, ndom - 1);

        // variables
        syst.add_var("nu", nu);
        syst.add_var("incA", incA);
        syst.add_var("incB", incB);
        syst.add_var("phi", phi);
        syst.add_var("incNp", incNphi);
        syst.add_var("ome", omega);

        // cst
        syst.add_cst("rsint", rsint);
        syst.add_cst("k", kk);
        syst.add_cst("qpi", qpi);
        syst.add_cst("val", val);

        // Definitions
        syst.add_def("phisurrsint = divrsint(phi)");
        syst.add_def("N = exp(nu)");
        syst.add_def("Np = divrsint(incNp)");
        syst.add_def("B = (divrsint(incB) + 1)/N");
        syst.add_def("A = exp(incA - nu)");

        syst.add_def("E = 0.5*(ome-Np*k)^2/N^2*phi^2 + 0.5*scal(grad(phi),grad(phi))/A^2 + 0.5*phi^2 + 0.5*k*k*phisurrsint*phisurrsint/B^2");
        syst.add_def("Pp = k/N* (ome-Np*k)*phi^2");
        syst.add_def("S = -0.5*scal(grad(phi), grad(phi))/A^2 - 0.5*k*k*phisurrsint*phisurrsint/B^2 + 1.5*(ome-Np*k)^2/N^2 * phi^2- 1.5*phi^2");
        syst.add_def("Spp = 0.5*(ome-Np*k)^2/N^2*phi^2 -0.5* scal(grad(phi),grad(phi))/A^2 -0.5* phi^2 +0.5* k*k*phisurrsint*phisurrsint/B^2");

        // Integral equation :
        space.add_eq_point(syst, Mc, "phi - val");

        // Equations
        for (int d = 0; d < ndom - 1; d++)
            syst.add_def(d,"eqnu = lap(nu) - B^2*rsint^2/2/N^2 * scal(grad(Np) , grad(Np)) + scal(grad(nu), grad(nu+log(B))) - qpi*A^2*(E+S)");
        syst.add_def(ndom - 1,"eqnu = lap(nu) - B^2*rsint^2/2/N^2 * scal(multr(grad(Np)) , multr(grad(Np))) + scal(grad(nu), grad(nu+log(B))) - qpi*A^2*(E+S)");

        for (int d = 0; d < ndom - 1; d++)
            syst.add_def(d,"eqshift = lap(incNp) - divrsint(Np) - rsint*scal(grad(Np),grad(nu-3*log(B)))+4*qpi*N*A^2/B^2*divrsint(Pp)");
        syst.add_def(ndom - 1,"eqshift = lap(incNp) - divrsint(Np) - rsint*scal(multr(grad(Np)),grad(nu-3*log(B)))+4*qpi*N*A^2/B^2*divrsint(Pp)");


        for (int d = 0; d < ndom - 1; d++)
            syst.add_def(d, "eqB = lap2(incB) -2*qpi*N*A^2*B*rsint*(S-Spp)");
        syst.add_def(ndom - 1, "eqB = lap2(incB) -2*qpi*N*A^2*B*rsint*multr(S-Spp)");

        for (int d = 0; d < ndom - 1; d++)
            syst.add_def(d,"eqA = lap2(incA) - 2*qpi*A^2*Spp- 3 * B^2 * rsint^2 /4 / N^2 * scal(grad(Np) , grad(Np)) + scal(grad(nu),grad(nu))");
        syst.add_def(ndom - 1,"eqA =lap2(incA) - 2*qpi*A^2*Spp- 3 * B^2 * rsint^2 /4 / N^2 * scal(multr(grad(Np)) , multr(grad(Np))) + scal(grad(nu),grad(nu)) ");


        for (int d = 0; d < ndom - 1; d++)
            syst.add_def(d,"eqphi = lap(phi) - A^2*(1-ome*ome/N^2+2*Np/N^2*ome*k-Np^2/N^2*k*k)*phi + scal(grad(phi),grad(nu+log(B))) - divrsint(A^2/B^2 -1)*divrsint(phi)*k*k");
        syst.add_def(ndom - 1,
                     "eqphi = lap(phi) - A^2*(1-ome*ome/N^2+2*Np/N^2*ome*k-Np^2/N^2*k*k)*phi + scal(grad(phi),grad(nu+log(B))) - k*k*divrsint(A^2/B^2-1)*divrsint(phi)");


        space.add_eq(syst, "eqnu=0", "nu", "dn(nu)");
        space.add_eq(syst, "eqshift=0", "incNp", "dn(incNp)");
        space.add_eq(syst, "eqphi=0", "phi", "dn(phi)");
        space.add_eq(syst, "eqA=0", "incA", "dn(incA)");
        space.add_eq(syst, "eqB=0", "incB", "dn(incB)");

        syst.add_eq_bc(ndom - 1, OUTER_BC, "nu=0");
        syst.add_eq_bc(ndom - 1, OUTER_BC, "incNp=0");
        syst.add_eq_bc(ndom - 1, OUTER_BC, "phi=0");
        syst.add_eq_bc(ndom - 1, OUTER_BC, "incA=0");
        syst.add_eq_bc(ndom - 1, OUTER_BC, "incB=0");

        double conv;
        bool endloop = false;
        int ite = 1;
        while (!endloop) {
            endloop = syst.do_newton(1e-8, conv,System_of_eqs::output_enabled);
            ite++;
        }
        syst.finalize_profiling();
        if (rank == 0)
            profiling_report(syst, std::cout);
    }


    if (rank == 0) {

        char name[100];
        sprintf(name, "bosinit.dat");
        FILE *fiche = fopen(name, "w");
        space.save(fiche);
        fwrite_be(&kk, sizeof(int), 1, fiche);
        fwrite_be(&omega, sizeof(double), 1, fiche);
        nu.save(fiche);
        incA.save(fiche);
        incB.save(fiche);
        incNphi.save(fiche);
        phi.save(fiche);
        fclose(fiche);
    }
#ifdef ENABLE_GPU_USE
	if(rank==0)
	{
		TESTING_CHECK(magma_finalize());
	}
#endif

    MPI_Finalize();
	
    return EXIT_SUCCESS ;
}

