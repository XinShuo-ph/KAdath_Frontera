//
// Created by sauliac on 29/05/2020.
//

#include "kadath_spheric.hpp"
#include "codes_utilities.hpp"

#ifndef __KADATH_CODES_SCHWARZ_HPP_
#define __KADATH_CODES_SCHWARZ_HPP_

using namespace Kadath;

//! Class demonstrating how to use Kadath to implement a solver for the Schwarz problem.
class Schwarz {
public:
    //! Overall dimension.
    static constexpr int dimension {3};

protected:
    //! Resolution for each coordinate.
    Dim_array number_of_points;
    //! Verbosity level.
    int verbosity;
public:
    //! Center of the coordinates
    Point center;
protected:
    //! Number of domains  :
    int number_of_domains; //special treatment, because setting it must resize the bounds array

public:
    //! boundaries
    Array<double> bounds ;
    //! Radius of the BH !
    double bh_radius;
    //! Coloc point.
    int type_coloc;

    //! Solving space.
    ptr_data_member(Space_spheric, space, unique);
    //! Solution in the conformal space.
    ptr_data_member(Scalar, conformal, unique);
    //! Pointer toward the system of equations object.
    ptr_data_member(System_of_eqs, system, unique);

    //! Current residue in the Newton-Rapthson algorithm
    internal_variable(double,newton_residue);
    //! Current number of iterations done in the NR algorithm.
    internal_variable(int,newton_nbr_iterations);
    //! Maximum allowed number of iterations (unlimited if zero or less).
    int newton_max_iterations;

    //! Difference between the approximated and analytical solutions in L-infinity norm.
    internal_variable(double,error_l_infinity);
    //! Difference between the approximated and analytical solutions in L-2 norm.
    internal_variable(double,error_l_2);

    //! Tolerance for error checking.
    double tolerance;

public:
    Dim_array const & get_number_of_points() const {number_of_points;}
    Dim_array & get_number_of_points() {number_of_points;}
    int get_verbosity() const {return verbosity;}
    Schwarz & set_verbosity(int new_value) {verbosity = new_value; return *this;}
    void set_number_of_points(int new_val) {number_of_points.set(0) = new_val;}
    int get_number_of_domains() const {return number_of_domains;}
    void set_number_of_domains(int new_val);

    // Simple constructor, just set the default values for all parameters.
    Schwarz(int nbr = 13,int ndom=4,double _bh_radius = 1.323,int _type_coloc=CHEB_TYPE) :
            number_of_points{dimension}, center{dimension}, number_of_domains{ndom},
            bounds{number_of_domains-1}, bh_radius{_bh_radius}, type_coloc{_type_coloc},
            space{nullptr}, conformal{nullptr}, system{nullptr},
            newton_residue{HUGE_VAL}, newton_nbr_iterations{0}, newton_max_iterations{-1},
            error_l_infinity{0.}, error_l_2{0.}, tolerance{1.e-8}
    {
        number_of_points.set(0) = nbr; number_of_points.set(1) = 5;
        number_of_points.set(2) = 4;
        for(int i{1};i<=dimension;i++) center.set(i) = 0;
        bounds.set(0) =          bh_radius;
        bounds.set(1) = 1.7557 * bh_radius;
        bounds.set(2) = 2.9861 * bh_radius;
    }

    /**
     * (Re)-set the space, solution scalar field and system of equations objects, based on
     * the parameter values (this should be re-called before solving if parameters are
     * changed through there mutators).
     */
    Schwarz & initialize() {
        // Sherical space :
        space.reset(new Space_spheric{type_coloc,center,number_of_points,bounds});

        // Initial guess for the conformal factor :
        conformal.reset(new Scalar{*space});
        *conformal = 1.;
        conformal->std_base();

        system.reset(new System_of_eqs{*space, 1, number_of_domains-1});
        // Only one unknown
        system->add_var ("P", *conformal) ;
        // One user defined constant
        system->add_cst ("a", bh_radius) ;

        // Inner BC
        space->add_inner_bc (*system, "dn(P) + 0.5 / a * P = 0") ;
        // Equation
        space->add_eq (*system, "Lap(P) = 0", "P", "dn(P)") ;
        // Outer BC
        space->add_outer_bc (*system, "P=1") ;

        newton_nbr_iterations = 0;
        newton_residue = HUGE_VAL;
        error_l_infinity = 0.;
        error_l_2 = 0.;

        this->set_verbosity(verbosity);//to call system settings without duplicating code. The verbosity assignment overhead is small.

        return *this;
    }

    /**
     * Performs the Newton-Raphson method.
     * @return \c true if \c newton_residue went lower than the \c tolerance value, \c false if the maximum
     * number of iterations is reached.
     */
    bool do_newton() {
        bool newton_success {false};
        bool const do_not_check_iter {newton_max_iterations < 0};
        while(!newton_success &&
                (do_not_check_iter || newton_nbr_iterations <= newton_max_iterations)) {
            newton_success = system->do_newton(tolerance,newton_residue, verbosity>0);
            newton_nbr_iterations++;
        }
    }

    /**
     * Compute the approximation error in the L-infinity and L2 norms.
     * @return \c true if the L-infinity error is lesser than the \c tolerance value.
     */
    bool check_solution() {
        int const resol = 100 ;
        double const xxmin = bounds(0)*1.01 ;
        double const xxmax = bounds(2)*5 ;
        double const step = (xxmax-xxmin)/resol ;
        double xx = xxmin+step ;

        double const tet = M_PI/2. ;
        double const phi = -2.349 ;

        double const xunit = sin(tet)*cos(phi) ;
        double const yunit = sin(tet)*sin(phi) ;
        double const zunit = cos(tet) ;

        Point M{3} ;
        for (int i=0 ; i<resol-1 ; i++) {

            M.set(1)=xx*xunit ;
            M.set(2)=xx*yunit ;
            M.set(3)=xx*zunit ;

            double ana = 1. + bh_radius/xx ;
            double error = fabs (ana - conformal->val_point(M)) ;
            if (error > error_l_infinity)
                error_l_infinity = error ;
            error_l_2 += error*error;
            xx+=step ;
        }
        return error_l_infinity <= tolerance;
    }

    //! Computes profiling datas (if enabled).
    Schwarz & finalize() {system->finalize_profiling(); return *this;}
    //! Sends the profiling datas in the passed output stream.
    void profiling_log(std::ostream & os) {
        if(verbosity >= 2) profiling_report(*system,os);
    }
};

inline void Schwarz::set_number_of_domains(int new_val) {
    assert(new_val > 0);
    number_of_domains = new_val;
    bounds.resize(number_of_domains-1);
}

#endif //__KADATH_CODES_SCHWARZ_HPP_
