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
 * \file Functionnal test that performs the "schwarz" sample code and compare the
 * obtained approximation to the known analytical solution. The test fails if
 * the error is not small enough.
 */

#include <random>
#include "system_of_eqs.hpp"
#include "kadath_spheric.hpp"


constexpr double TESTS_TOLERANCE {1.e-8};

using namespace Kadath;

enum : bool
{
    DO_NOT_COMPUTE = false,
    COMPUTE = true
};

template<int R1=13,int R2=5,int R3=4> class Schwarz_test
{
public:
    static constexpr int dim {3};
    static constexpr int resolution[] = {R1,R2,R3};
    static constexpr double b_radius[] = {1.,1.7557,2.9861};
    static constexpr int ndom {4};
    static constexpr double aa {1.323};
    static constexpr int type_coloc {CHEB_TYPE};
    static constexpr double newton_tol {1.e-8};

    Dim_array res;
    Array<double> bounds;
    std::unique_ptr<Space_spheric> p_space;
    std::unique_ptr<Scalar> p_conf;
    std::unique_ptr<System_of_eqs> p_syst;
    bool converged;
    double newton_error;
    double error_l_infinity,error_l_2;
    bool initialize();
    bool init_syst();
    void do_newton();
    void check_solution();

public:
    double get_error_l_infinity() const {return error_l_infinity;}
    double get_error_l_2() const {return error_l_2;}

    Schwarz_test(bool compute = COMPUTE)
            : res{dim}, bounds{ndom-1}, p_space{nullptr}, p_conf{nullptr},
              p_syst{nullptr}, converged{false}, newton_error{-1.}
    {
        initialize();
        if(compute) {
            do_newton();
            check_solution();
        }
    }

    bool has_converged() const {return converged;}

};

template<int A,int B,int C> bool Schwarz_test<A,B,C>::initialize()
{
    res.set(0) = resolution[0] ; res.set(1) = resolution[1] ; res.set(2) = resolution[2] ;
    bounds.set(0) = b_radius[0]*aa ; bounds.set(1) = b_radius[1]*aa ; bounds.set(2) = b_radius[2]*aa ;
    p_space.reset(new Space_spheric{type_coloc, Point{dim}, res, bounds});
    p_conf.reset(new Scalar{*p_space});
    *p_conf = 1.;
    p_conf->std_base();
    bool const system_initialized {init_syst()};
    assert(system_initialized);
    p_syst->add_var ("P", *p_conf) ;
    p_syst->add_cst ("a", aa) ;
    p_space->add_inner_bc (*p_syst, "dn(P) + 0.5 / a * P = 0") ;
    p_space->add_eq (*p_syst, "Lap(P) = 0", "P", "dn(P)") ;
    p_space->add_outer_bc (*p_syst, "P=1") ;
    return system_initialized && p_space != nullptr && p_conf != nullptr ;
}

template<int A,int B,int C>
void Schwarz_test<A,B,C>::do_newton()
{
    unsigned short niter{0};
    while (!converged && niter<2) {
        converged = p_syst->do_newton(newton_tol, newton_error);
        niter++;
    }
}

template<int A,int B,int C> void Schwarz_test<A,B,C>::  check_solution()
{
    assert(p_conf);
    int resol = 100 ;
    double xxmin = bounds(0)*1.01 ;
    double xxmax = bounds(2)*5 ;
    double step = (xxmax-xxmin)/resol ;
    double xx = xxmin+step ;
    error_l_infinity = 0. ;
    error_l_2 = 0.;

    double tet = M_PI/2. ;
    double phi = -2.349 ;

    double xunit = sin(tet)*cos(phi) ;
    double yunit = sin(tet)*sin(phi) ;
    double zunit = cos(tet) ;

    Point M (3) ;
    for (int i=0 ; i<resol-1 ; i++) {

        M.set(1)=xx*xunit ;
        M.set(2)=xx*yunit ;
        M.set(3)=xx*zunit ;

        double ana = 1. + aa/xx ;
        double error_l1 = fabs (ana - p_conf->val_point(M)) ;
        error_l_2 += error_l1*error_l1;
        if (error_l1 > error_l_infinity)
            error_l_infinity = error_l1 ;
        xx+=step ;
    }
    error_l_2 = sqrt(error_l_2);
}

template<> bool Schwarz_test<>::init_syst()
{
    if(p_space) {
        p_syst.reset(new System_of_eqs{*p_space,1,ndom-1});
    }
    return p_syst != nullptr;
}

int main(int argc,char * argv[]) {
#ifdef PAR_VERSION
    MPI_Init(&argc,&argv);
#endif
    std::cout << "========================== t-schwarz_seq functional test ===========================\n\n";

    Schwarz_test<> test{};
    assert(test.converged);
    assert(test.error_l_infinity <= TESTS_TOLERANCE);
    assert(test.error_l_2 <= TESTS_TOLERANCE);
    std::cout << "\n\n====================================================================================";
#ifdef PAR_VERSION
    MPI_Finalize();
#endif
    return 0;
}