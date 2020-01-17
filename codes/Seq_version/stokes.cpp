//
// Created by sauliac on 16/01/2020.
//

#include "kadath_polar.hpp"

using namespace Kadath;

int main(int argc,char * argv[])
{
    constexpr int dim = 2;
    constexpr int nbr = 10;
    Dim_array res(dim);
    res.set(1) = res.set(0) = nbr;

    Point center(dim);
    center.set(1) = center.set(2) = 0;

    constexpr int ndom = 2;
    Array<double> bounds{ndom - 1};
    constexpr double R{1.};
    bounds.set(0) = R;

    Space_polar space{CHEB_TYPE, center, res, bounds};

    Scalar conf_ur{space};
    Scalar conf_utheta{space};
    Scalar conf_p{space};
    conf_ur = 1.;
    conf_utheta = 1.;
    conf_p = 0.;
    conf_ur.std_base();
    conf_utheta.std_base();
    conf_p.std_base();

    System_of_eqs system{space,1,ndom-1};
    system.add_var("ux",conf_ur);
    system.add_var("uy", conf_utheta);
    system.add_var("p",conf_p);

    return EXIT_SUCCESS;
}