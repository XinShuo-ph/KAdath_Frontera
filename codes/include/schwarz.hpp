//
// Created by sauliac on 29/05/2020.
//

#include "kadath_spheric.hpp"

#ifndef __KADATH_CODES_SCHWARZ_HPP_
#define __KADATH_CODES_SCHWARZ_HPP_

using namespace Kadath;

//! Class demonstrating how to use Kadath to implement a solver for the Schwarz problem.
class Schwarz {
public:
    //! Overall dimension.
    static constexpr int dimension {3};

private:
    //! Resolution for each coordinate.
    Dim_array number_of_points;
    //! Center of the coordinates
    Point center;
    //! Number of domains  :
    int number_of_domains;
    //! boundaries
    Array<double> bounds ;
    //! Radius of the BH !
    double bh_radius ;
    //! Coloc point.
    int type_coloc ;

    //! Solving space.
    std::unique_ptr<Space_spheric> space;
    //! Solution in the conformal space.
    std::unique_ptr<Scalar> conformal;
    //! Pointer toward the system of equations object.
    std::unique_ptr<System_of_eqs> system;

    //! Current residue in the Newton-Rapthson algorithm
    double newton_residue;
    //! Current number of iterations done in the NR algorithm.
    int newton_nbr_iterations;
    //! Maximum allowed number of iterations (unlimited if zero or less).
    int newton_max_iterations;

    //! Difference between the approximated and analytical solutions in L-infinity norm.
    double error_l_infinity;
    //! Difference between the approximated and analytical solutions in L-2 norm.
    double error_l_2;

    //! Tolerance for error checking.
    double tolerance;

public:
    // Accessors and mutators for parameters :
    int get_number_of_points(int i) const {return number_of_points(i);}
    Schwarz & set_number_of_points(int i,int new_val) {number_of_points.set(i) = new_val; return *this;}
    Point const & get_center() const {return center;}
    Schwarz & set_center(Point const & new_center) {center = new_center; return *this;}
    int get_number_of_domains() const {return number_of_domains;}
    Schwarz & set_number_of_domains(int new_value) ;
    Array<double> const & get_bounds() const {return bounds;}
    Array<double> & get_bounds() {return bounds;}
    double get_bh_radius() const {return bh_radius;}
    Schwarz & set_bh_radius(double new_val) {bh_radius = new_val;return *this;}
    int get_type_coloc() const {return type_coloc;}
    Schwarz & set_type_coloc(int new_val) {type_coloc = new_val; return *this;}
    int get_newton_max_iterations() const {return newton_max_iterations;}
    Schwarz & set_newton_max_iterations(int new_val) {newton_max_iterations = new_val; return *this;}
    double get_tolerance() const {return tolerance;}
    Schwarz & set_tolerance(double new_val) {tolerance = new_val; return *this;}

    // Accessor for internal read only values :
    double get_newton_residue() const {return newton_residue;}
    int get_newton_nbr_iterations() const {return newton_nbr_iterations;}
    double get_error_l_infinity() const {return error_l_infinity;}
    double get_error_l_2() const {return error_l_2;}

    // Accessors for Kadath space, system and solution objects :
    std::unique_ptr<Space_spheric> & get_space() {return space;}
    std::unique_ptr<Space_spheric> const & get_space() const {return space;}
    std::unique_ptr<Scalar> & get_conformal() {return conformal;}
    std::unique_ptr<Scalar> const & get_conformal() const {return conformal;}
    std::unique_ptr<System_of_eqs> & get_system() {return system;}
    std::unique_ptr<System_of_eqs> const & get_system() const {return system;}


    // Simple constructor, just set the default values for all parameters.
    Schwarz(int nbr = 13,int ndom=4,double _bh_radius = 1.323,int _type_coloc=CHEB_TYPE) :
            number_of_points{dimension}, center{dimension}, number_of_domains{ndom},
            bounds{number_of_domains-1}, bh_radius{_bh_radius}, type_coloc{_type_coloc},
            space{nullptr}, conformal{nullptr}, system{nullptr},
            newton_residue{HUGE_VAL}, newton_nbr_iterations{0}, newton_max_iterations{0},
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
    Schwarz & build_space_and_system() {
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

        return *this;
    }

    /**
     * Performs the Newton-Rapthson method.
     * @return \c true if \c newton_residue went lower than the \c tolerance value, \c false if the maximum
     * number of iterations is reached.
     */
    bool do_newton() {
        bool newton_success {false};
        bool const do_not_check_iter {newton_max_iterations <= 0};
        while(!newton_success &&
                (do_not_check_iter || newton_nbr_iterations <= newton_max_iterations)) {
            newton_success = system->do_newton(tolerance,newton_residue);
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
    void profiling_log(std::ostream & os) {profiling_report(*system,os);}
};

inline Schwarz & Schwarz::set_number_of_domains(int new_value) {
    assert(new_value > 0);
    number_of_domains = new_value;
    bounds.resize(number_of_domains-1);
    return *this;
}

#endif //__KADATH_CODES_SCHWARZ_HPP_
