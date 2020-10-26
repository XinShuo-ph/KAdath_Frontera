/*
    Copyright 2020 sauliac

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


/**
 * \file Functionnal test using the Kerr sample code. Since no analytical solution
 * is known, the test simply consists in leading the computation to its end without
 * a scratch. Due to the relative complexity of the problem, it is only enalbled
 * for MPI parallel builds.
 */


#include "base_fftw.hpp"
#include "kadath_spheric.hpp"
#include "mpi.h"

using namespace Kadath ;

int main(int argc, char** argv) {

    int rc = MPI_Init (&argc, &argv) ;
    int rank = 0 ;
    MPI_Comm_rank (MPI_COMM_WORLD, &rank) ;

    // 3D :
    int dim = 3 ;

    // Number of points
    int nbr  = 11 ;
    int type_coloc = CHEB_TYPE ;
    Dim_array res (dim) ;
    res.set(0) = nbr ; res.set(1) = nbr ; res.set(2) = 1 ;

    // Center of the coordinates
    Point center (dim) ;

    // Number of domains and boundaries :
    int ndom = 3 ;
    Array<double> bounds (ndom-1) ;
    // Radius of the BH !
    double aa = 1. ;
    bounds.set(0) = aa ; bounds.set(1) = 2*aa ;
    // Spherical space :
    Space_spheric espace(type_coloc, center, res, bounds) ;

    // Spherical tensorial basis everywhere...
    Base_tensor basis (espace, SPHERICAL_BASIS) ;

    double n0 = 0.5 ; // Lapse on the horizon

    // Initial guess for the conformal factor :
    Scalar conf (espace) ;
    conf = 1  ;
    for (int i=1 ; i<ndom ; i++)
        conf.set_domain(i) = 1+aa/espace.get_domain(i)->get_radius() ;
    conf.std_base() ;

    // Lapse
    Scalar lapse (espace) ;
    lapse = 1 ;
    lapse.std_base() ;

    // Shift
    Vector shift (espace, CON, basis) ;
    for (int i=1 ; i<=3 ; i++)
        shift.set(i).annule_hard()  ;
    shift.std_base() ;

    Vector stilde (espace, CON, basis) ;
    for (int i=1 ; i<=3 ; i++)
        stilde.set(i) = 0. ;
    stilde.set(1).set_domain(1) = 1. ;
    stilde.set(2).set_domain(1) = 0. ;
    stilde.set(3).set_domain(1) = 0. ;
    stilde.std_base() ;

    // Flat metric :
    Metric_flat fmet (espace, basis) ;

    // Start to compute a Schwarzschild :
    double ome = 0. ;
    System_of_eqs syst_init (espace, 1, ndom-1) ;
    // Unknowns
    syst_init.add_var ("P", conf) ;
    syst_init.add_var ("N", lapse) ;
    syst_init.add_var ("bet", shift) ;

    // User defined constants
    syst_init.add_cst ("a", aa) ;
    syst_init.add_cst ("s", stilde) ;
    syst_init.add_cst ("n0", n0) ;

    // Metric :
    fmet.set_system(syst_init, "f") ;

    // definition of the extrinsic curvature :
    syst_init.add_def ("A^ij = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) /2. / N") ;

    // Inner BC :
    espace.add_inner_bc (syst_init, "N=n0") ;
    espace.add_inner_bc (syst_init, "bet^i = n0 / P^2 * s^i") ; // 1
    espace.add_inner_bc (syst_init, "dn(P) + 0.5 * P / a + P^3 * A_ij * s^i * s^j /4.= 0") ;

    // Equations :
    espace.add_eq (syst_init, "D_i D^i N + 2 * D_i P * D^i N / P - N * P^4 * A_ij *A^ij= 0", "N", "dn(N)") ;
    espace.add_eq (syst_init, "D_i D^i P + P^5 *A_ij * A^ij / 8= 0", "P", "dn(P)") ;
    espace.add_eq (syst_init, "D_j A^ij + 6 * A^ij * D_j P / P =0", "bet^i", "dn(bet^i)") ;

    // Outer BC
    espace.add_outer_bc (syst_init, "N=1") ;
    espace.add_outer_bc (syst_init, "P=1") ;
    espace.add_outer_bc (syst_init, "bet^i=0") ;

    {
        double conv ;
        bool endloop = false ;
        int ite = 1 ;
        if (rank==0) {
            cout << "Computation with omega =  0" << endl;
        }
        while (!endloop) {
            endloop = syst_init.do_newton(1e-6, conv) ;
            ite++ ;
        }
    }
    syst_init.finalize_profiling();

    if(rank==0) profiling_report(syst_init,std::cout);

    Metric_tensor gfixed (espace, CON, basis) ;
    for (int i=1 ; i<=3 ; i++)
        for (int j=i ; j<=3 ; j++)
            if (i==j)
                gfixed.set(i,j) = 1. ;
            else
                gfixed.set(i,j).annule_hard() ;
    gfixed.std_base() ;

    // Loop for omega :
    double step = 0.005 ;
    int nbr_ome = argc>1 ? std::atoi(argv[1]) : 1 ;

    // Associated metric
    Metric_tensor gmet(gfixed) ;
    Metric_dirac met (gmet) ;

    Vector scov (espace, COV, basis) ;
    scov.set(1) = 1. ;
    scov.set(2) = 0. ;
    scov.set(3) = 0. ;
    scov.std_base() ;

    Vector er (espace, CON, basis) ;
    er.set(1) = 1. ;
    er.set(2) = 0. ;
    er.set(3) = 0. ;

// Vector parallel to the sphere (needed only for inner BC)
    Vector mm (espace, CON, basis) ;
    for (int i=1 ; i<=3 ; i++)
        mm.set(i) = 0. ;
    Val_domain xx (espace.get_domain(1)->get_cart(1)) ;
    Val_domain yy (espace.get_domain(1)->get_cart(2)) ;
    mm.set(3).set_domain(1) = sqrt(xx*xx + yy*yy) ;
    mm.std_base() ;

    int n_evol_inner = 4 ;
    Array<int>** p_evol_inner = new Array<int>* [n_evol_inner] ;
    for (int i=0 ; i<n_evol_inner ; i++)
        p_evol_inner[i] = new Array<int> (2) ;
    p_evol_inner[0]->set(0) = 1 ; p_evol_inner[0]->set(1) = 1;
    p_evol_inner[1]->set(0) = 1 ; p_evol_inner[1]->set(1) = 3 ;
    p_evol_inner[2]->set(0) = 2 ; p_evol_inner[2]->set(1) = 2 ;
    p_evol_inner[3]->set(0) = 2 ; p_evol_inner[3]->set(1) = 3 ;

    int n_evol = 5 ;
    Array<int>** p_evol = new Array<int>* [n_evol] ;
    for (int i=0 ; i<n_evol ; i++)
        p_evol[i] = new Array<int> (2) ;
    p_evol[0]->set(0) = 1 ; p_evol[0]->set(1) = 1 ;
    p_evol[1]->set(0) = 1 ; p_evol[1]->set(1) = 2 ;
    p_evol[2]->set(0) = 1 ; p_evol[2]->set(1) = 3 ;
    p_evol[3]->set(0) = 2 ; p_evol[3]->set(1) = 2 ;
    p_evol[4]->set(0) = 2 ; p_evol[4]->set(1) = 3 ;

    int n_dirac = 1 ;
    Array<int>** p_dirac = new Array<int>* [n_dirac] ;
    for (int i=0 ; i<n_dirac ; i++)
        p_dirac[i] = new Array<int>(1) ;
    p_dirac[0]->set(0) = 2 ;

    for (int conte=0 ; conte<nbr_ome ; conte ++)
    {

        ome += step ;
        if (rank==0)
            cout << "Computation with omega = " << ome << endl ;

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
        syst.add_def ("LieK_ij = 4 * A_ij * bet^k * DP_k / P + bet^k * D_k A_ij + A_ik * Dbet_j^k + A_jk * Dbet_i^k") ;
        syst.add_def ("DDN_ij = D_i DN_j - 2 * DN_i * DP_j / P - 2 * DN_j * DP_i / P + 2 * g_ij * DN_k * DP^k / P") ;
        syst.add_def ("PartR_ij = R_ij + 6 * DP_i * DP_j / P^2 - 2 * D_i DP_j / P - 2 * g_ij * D_k DP^k / P - 2 * g_ij * DP_k * DP^k / P^2 - P^4 * 2 * A_ik * A_j^k") ;
        syst.add_def ("evol_ij = DDN_ij - N * PartR_ij - P^4 * LieK_ij") ;
        syst.add_def ("Pfourhor = N^2 / ( s_i * bet^i ) / ( s_i * bet^i )") ;

        espace.add_inner_bc (syst, "N=n0") ;
        espace.add_inner_bc (syst, "bet^i = n0 / P^2 * st^i + Ome * m^i * a") ;
        espace.add_inner_bc (syst, "4 * st^i * D_i P / P + D_i st^i + P^2 * A_ij * st^i * st^j = 0") ;
        espace.add_inner_bc (syst, "DDN_ij - N * PartR_ij - Pfourhor * LieK_ij=0", n_evol_inner, p_evol_inner) ;
        espace.add_inner_bc (syst, "dirac^i =0", n_dirac, p_dirac) ;

        // CFC Equations :
        espace.add_eq (syst, " D_i DN^i + 2 * DP_i * DN^i / P - N * P^4 * A_ij * A^ij = 0", "N", "dn(N)") ;
        espace.add_eq (syst, "R - 8 * D_i DP^i / P - P^4 * A_ij * A^ij = 0", "P", "dn(P)") ;
        espace.add_eq (syst, "D^j A_ij + 6 * A_ij * DP^j / P =0", "bet^i", "dn(bet^i)") ;

        // Evolution
        syst.add_eq_inside (1, "evol_ij =0", n_evol, p_evol) ;
        for (int d=2 ; d<ndom ; d++) {
            syst.add_eq_matching (d-1, OUTER_BC, "g^ij", n_evol, p_evol) ;
            syst.add_eq_matching (d-1, OUTER_BC, "dn(g^ij)", n_evol, p_evol) ;
            syst.add_eq_inside (d, "evol_ij=0", n_evol, p_evol) ;
        }


        espace.add_eq_full (syst, "determinant(g^ij) = 1") ;

        // Outer BC
        espace.add_outer_bc (syst, "N=1") ;
        espace.add_outer_bc (syst, "P=1") ;
        espace.add_outer_bc (syst, "bet^i=0") ;
        espace.add_outer_bc (syst, "g^ij=gf^ij", n_evol, p_evol) ;

        // Newton-Raphson
        double conv ;
        bool endloop = false ;
        int ite = 1 ;
        char name[100] ;
        sprintf(name, "kerr_%d_%f.dat", nbr, ome) ;

        while (!endloop) {
            endloop = syst.do_newton(1e-8, conv) ;
            ite++ ;
            // Save
            if (rank==0) {
                FILE* ff = fopen(name, "w") ;
                espace.save(ff) ;
                fwrite_be (&n0, sizeof(double), 1, ff) ;
                fwrite_be (&ome, sizeof(double), 1, ff) ;
                fwrite_be (&aa, sizeof(double), 1, ff) ;
                conf.save(ff) ;
                lapse.save(ff) ;
                shift.save(ff) ;
                gmet.save(ff) ;
                fclose(ff) ;
            }
        }
        syst.finalize_profiling();

        if(rank==0) profiling_report(syst_init,std::cout);
    }
    for (int i=0 ; i<n_evol_inner ; i++) delete p_evol_inner[i] ;
    delete [] p_evol_inner ;
    for (int i=0 ; i<n_evol ; i++) delete p_evol[i] ;
    delete [] p_evol ;
    for (int i=0 ; i<n_dirac ; i++) delete p_dirac[i] ;
    delete [] p_dirac ;

    MPI_Finalize() ;
    return EXIT_SUCCESS ;
}

