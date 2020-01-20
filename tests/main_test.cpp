//
// Created by sauliac on 20/01/2020.
//

#include "gtest/gtest.h"
#include <utility>
#include "kadath_spheric.hpp"

constexpr double TESTS_TOLERANCE {1.e-8};

using namespace Kadath;

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

    // Newton-Raphson
    double conv ;
    bool endloop = false ;
    int ite = 1 ;
    while (!endloop) {
        endloop = syst.do_newton(1e-8, conv) ;
        ite++ ;
    }

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

class KadathSchwarzTest : public ::testing::Test
{
protected:
    double error_l_infinity;
    double error_l_2;
public:
    KadathSchwarzTest();
};

TEST_F(KadathSchwarzTest,CheckErrorInLInfAndL2Norms)
{
    EXPECT_LE(error_l_infinity,TESTS_TOLERANCE);
    EXPECT_LE(error_l_2,TESTS_TOLERANCE);
}

KadathSchwarzTest::KadathSchwarzTest() : Test{}, error_l_infinity{}, error_l_2{}
{
    std::tie(error_l_infinity,error_l_2) = test_schwarz();
}