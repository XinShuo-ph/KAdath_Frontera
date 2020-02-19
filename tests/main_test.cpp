//
// Created by sauliac on 20/01/2020.
//

#include "gtest/gtest.h"
#include <utility>
#include "multi_thread.hpp"
#include "kadath_spheric.hpp"


constexpr double TESTS_TOLERANCE {1.e-8};

using namespace Kadath;


template<class System> class Schwarz_test
{
public:
    static constexpr int dim {3};
    static constexpr int resolution[] = {13,5,4};
    static constexpr double b_radius[] = {1.,1.7557,2.9861};
    static constexpr int ndom {4};
    static constexpr double aa {1.323};
    static constexpr int type_coloc {CHEB_TYPE};
    static constexpr double newton_tol {1.e-8};
protected:
    Dim_array res;
    Array<double> bounds;
    std::unique_ptr<Space_spheric> p_space;
    std::unique_ptr<Scalar> p_conf;
    std::unique_ptr<System> p_syst;
    bool converged;
    double newton_error;
    double error_l_infinity,error_l_2;
    bool initialize(std::size_t num_threads = 1);
    bool init_syst(std::size_t num_threads = 1);
    void do_newton();
    void check_solution();

public:
    double get_error_l_infinity() const {return error_l_infinity;}
    double get_error_l_2() const {return error_l_2;}

    Schwarz_test(std::size_t num_threads = 1) : res{dim}, bounds{ndom-1}, p_space{nullptr}, p_conf{nullptr},
                                                p_syst{nullptr}, converged{false}, newton_error{-1.}
    {
        initialize(num_threads);
        do_newton();
        check_solution();
    }

};

template<class System> bool Schwarz_test<System>::initialize(std::size_t num_threads)
{
    res.set(0) = resolution[0] ; res.set(1) = resolution[1] ; res.set(2) = resolution[2] ;
    bounds.set(0) = b_radius[0]*aa ; bounds.set(1) = b_radius[1]*aa ; bounds.set(2) = b_radius[2]*aa ;
    p_space.reset(new Space_spheric{type_coloc, Point{dim}, res, bounds});
    p_conf.reset(new Scalar{*p_space});
    *p_conf = 1.;
    p_conf->std_base();
    bool const system_initialized {init_syst(num_threads)};
    assert(system_initialized);
    p_syst->add_var ("P", *p_conf) ;
    p_syst->add_cst ("a", aa) ;
    p_space->add_inner_bc (*p_syst, "dn(P) + 0.5 / a * P = 0") ;
    p_space->add_eq (*p_syst, "Lap(P) = 0", "P", "dn(P)") ;
    p_space->add_outer_bc (*p_syst, "P=1") ;
    p_syst->disable_data_display();
    return system_initialized && p_space != nullptr && p_conf != nullptr ;
}

template<class System>
void Schwarz_test<System>::do_newton()
{
    while (!converged) converged = p_syst->do_newton(newton_tol, newton_error) ;

}

template<class System> void Schwarz_test<System>::check_solution()
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

template<> bool Schwarz_test<System_of_eqs>::init_syst(std::size_t num_threads)
{
    if(p_space) {
        p_syst.reset(new System_of_eqs{*p_space,1,ndom-1});
    }
    return p_syst != nullptr;
}

template<> bool Schwarz_test<System_of_eqs_threaded>::init_syst(std::size_t num_threads)
{
    if(p_space) {
        p_syst.reset(new System_of_eqs_threaded{num_threads,*p_space,1,ndom-1});
    }
    return p_syst != nullptr;
}

pair<double,double> const test_schwarz()
{
    // 3D :
    int dim = 3 ;

    // Number of points
    int nbr  = 13 ;
    Dim_array res (dim) ;
    res.set(0) = nbr ; res.set(1) = 5 ; res.set(2) = 4 ;

    // Center of the coordinates
    Point center (dim) ;
    for (int i=1 ; i<=dim ; i++)
        center.set(i) = 0 ;

    // Number of domains and boundaries :
    int ndom = 4 ;
    Array<double> bounds (ndom-1) ;
    // Radius of the BH !
    double aa = 1.323 ;
    bounds.set(0) = aa ; bounds.set(1) = 1.7557*aa ; bounds.set(2) = 2.9861*aa ;

    // Chebyshev or Legendre :
    int type_coloc = CHEB_TYPE ;

    // Sherical space :
    Space_spheric space(type_coloc, center, res, bounds) ;

    // Initial guess for the conformal factor :
    Scalar conf (space) ;
    conf = 1. ;
    conf.std_base() ;
/////

    // Solve the equation in space outside the nucleus
    System_of_eqs syst (space, 1, ndom-1) ;
    // Only one unknown
    syst.add_var ("P", conf) ;
    // One user defined constant
    syst.add_cst ("a", aa) ;

    // Inner BC
    space.add_inner_bc (syst, "dn(P) + 0.5 / a * P = 0") ;
    // Equation
    space.add_eq (syst, "Lap(P) = 0", "P", "dn(P)") ;
    // Outer BC
    space.add_outer_bc (syst, "P=1") ;
    syst.disable_data_display();

/////
    // Newton-Raphson
    double conv ;
    bool endloop = false ;
    int ite = 1 ;
    while (!endloop) {
        endloop = syst.do_newton(1e-8, conv) ;
        ite++ ;
    }


/////
    // Check of the solution
    int resol = 100 ;
    double xxmin = bounds(0)*1.01 ;
    double xxmax = bounds(2)*5 ;
    double step = (xxmax-xxmin)/resol ;
    double xx = xxmin+step ;
    double error_l_infinity {0.} ;
    double error_l_2 {0.};

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
        double error_l1 = fabs (ana - conf.val_point(M)) ;
        error_l_2 += error_l1*error_l1;
        if (error_l1 > error_l_infinity)
            error_l_infinity = error_l1 ;
        xx+=step ;
    }
    error_l_2 = sqrt(error_l_2);

    return {error_l_infinity,error_l_2};
}

class KadathSchwarzTest : public ::testing::Test, public Schwarz_test<System_of_eqs>
{
public:
    KadathSchwarzTest() : Test{} , Schwarz_test<System_of_eqs>{} {}
};

class KadathSchwarzMultiThreadTest : public ::testing::Test
{
public:
    static constexpr std::size_t num_threads {3};
protected:
    Schwarz_test<System_of_eqs_threaded> single_thread;
    Schwarz_test<System_of_eqs_threaded> multi_thread;
public:
    KadathSchwarzMultiThreadTest() : Test{}, single_thread{1}, multi_thread{num_threads} {}
};

TEST_F(KadathSchwarzTest,CheckErrorInLInfAndL2Norms)
{
    EXPECT_LE(error_l_infinity,TESTS_TOLERANCE);
    EXPECT_LE(error_l_2,TESTS_TOLERANCE);
}

TEST_F(KadathSchwarzMultiThreadTest,CheckErrorInLInfAndL2Norms)
{
    EXPECT_LE(single_thread.get_error_l_infinity(),TESTS_TOLERANCE);
    EXPECT_LE(single_thread.get_error_l_2(),TESTS_TOLERANCE);
    EXPECT_LE(multi_thread.get_error_l_infinity(),TESTS_TOLERANCE);
    EXPECT_LE(multi_thread.get_error_l_2(),TESTS_TOLERANCE);
}
