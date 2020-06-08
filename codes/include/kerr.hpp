//
// Created by sauliac on 29/05/2020.
//
#include "codes_utilities.hpp"
#include "kadath_spheric.hpp"

#ifndef __KADATH_CODES_KERR_HPP_
#define __KADATH_CODES_KERR_HPP_

using namespace Kadath;

//! Class demonstrating how to use Kadath to implement a solver for the Kerr problem.
class Kerr_base {
public:
    //! Overall dimension.
    static constexpr int dimension {3};

protected:
    //! Resolution for each coordinate.
    Dim_array number_of_points;

public:
    //! Center of the coordinates
    Point center;
    //! Number of domains  :
    int number_of_domains;
    //! boundaries
    Array<double> bounds;
    //! Radius of the BH !
    double bh_radius ;
    //! Coloc point.
    int type_coloc ;
    //! Lapse value on horizon.
    double n0;
    double omega;

    //! Solving space.
    ptr_data_member(Space_spheric,space);
    //! Tensorial basis
    ptr_data_member(Base_tensor,basis);
    //! Solution in the conformal space.
    ptr_data_member(Scalar,conformal);
    //! Lapse.
    ptr_data_member(Scalar,lapse);
    //! Shift.
    ptr_data_member(Vector,shift);
    //! Pointer toward the system of equations object.
    ptr_data_member(System_of_eqs,system);

    //! Current residue in the Newton-Rapthson algorithm
    internal_variable(double,newton_residue);
    //! Current number of iterations done in the NR algorithm.
    internal_variable(int,newton_nbr_iterations);
    //! Maximum allowed number of iterations (unlimited if zero or less).
    int newton_max_iterations;
    //! Tolerance for error checking.
    double tolerance;
    //! MPI rank (0 if sequential).
    int mpi_rank;

public:
    Dim_array const & get_number_of_points() const {return number_of_points;}
    Dim_array & get_number_of_points() {return number_of_points;}
    void set_number_of_points(int new_value) {
        number_of_points.set(0) = new_value;
        number_of_points.set(1) = new_value;
    }

    // Simple constructor, just set the default values for all parameters.
    Kerr_base(int nbr = 17,int ndom=3,double _bh_radius = 1.,int _type_coloc=CHEB_TYPE) :
        number_of_points{dimension}, center{dimension}, number_of_domains{ndom},
        bounds{number_of_domains-1}, bh_radius{_bh_radius}, type_coloc{_type_coloc}, n0{0.5},
        omega{0.}, space{nullptr}, basis{nullptr}, conformal{nullptr}, system{nullptr},
        newton_residue{HUGE_VAL}, newton_nbr_iterations{0}, newton_max_iterations{-1},
        tolerance{1.e-6}, mpi_rank{0}
    {
        number_of_points.set(0) = nbr; number_of_points.set(1) = nbr;
        number_of_points.set(2) = 1;
        for(int i{1};i<=dimension;i++) center.set(i) = 0;
        bounds.set(0) = bh_radius;
        bounds.set(1) = 2 * bh_radius;
    }

    virtual ~Kerr_base() {};

    /** Resets the \c space data member (must be called when either \c type_coloc, the \c center
     * or \c number_of_points and \c bounds have been modified).
     * @return a reference to the current \c Kerr_base object.
     */
    virtual Kerr_base & reset_space() {
        space.reset(new Space_spheric{type_coloc,center,number_of_points,bounds});
        basis.reset(new Base_tensor{*space,SPHERICAL_BASIS});
        return *this;
    }

    /**
     * Reset the initial values for the unknowns of the problem. The \c space data member
     * must be non null when this method is called (either set it manually or through the
     * \c reset_space method).
     * @return a reference to the current object.
     */
    virtual Kerr_base & reset_initial_guess() {
        if(space && basis) {
            // Initial guess for the conformal factor :
            conformal.reset(new Scalar{*space});
            *conformal = 1;
            for (int i = 1; i < number_of_domains; i++)
                conformal->set_domain(i) = 1 + bh_radius / space->get_domain(i)->get_radius();
            conformal->std_base();
            // Initial guess for lapse :
            lapse.reset(new Scalar{*space});
            *lapse = 1;
            lapse->std_base();
            // initial guess for shift.
            shift.reset(new Vector{*space, CON, *basis});
            for (int i = 1; i <= 3; i++) shift->set(i).annule_hard();
            shift->std_base();
        }
        else initialization_order_error(std::cerr,__FILE__,__LINE__);
        return *this;
    }

    /**
     * Reset the system of equations (the \c space, tensorial \c basis, and all other
     * data members related to initial guess's must have been set before calling this method).
     * @return a reference toward the current object.
     */
    virtual Kerr_base & reset_system() {
        if(space && conformal && lapse && shift) {
            system.reset(new System_of_eqs{*space, 1, number_of_domains - 1});
            system->add_var("P", *conformal);
            system->add_var("N", *lapse);
            system->add_var("bet", *shift);
            system->add_cst ("a", bh_radius) ;
            system->add_cst ("n0", n0) ;
        } else initialization_order_error(std::cerr,__FILE__,__LINE__);
        return *this;
    }

    /**
    * (Re)-set the space, solution scalar field and system of equations objects, based on
    * the parameter values (this should be re-called before solving if parameters are
    * changed through there mutators).
    */
    virtual Kerr_base & build_space_and_system() {
        this->reset_space();
        this->reset_initial_guess();
        this->reset_system();

        newton_nbr_iterations = 0;
        newton_residue = HUGE_VAL;

        return *this;
    }
    /**
     * Performs the Newton-Rapthson method.
     * @return \c true if \c newton_residue went lower than the \c tolerance value, \c false if the maximum
     * number of iterations is reached.
     */
    virtual bool do_newton() {
        bool newton_success {false};
        bool const do_not_check_iter {newton_max_iterations < 0};
        if(mpi_rank==0) std::cout << "Computation with omega = " << omega << std::endl;
        while(!newton_success &&
              (do_not_check_iter || newton_nbr_iterations <= newton_max_iterations)) {
            newton_success = system->do_newton(tolerance,newton_residue);
            newton_nbr_iterations++;
        }
    }

    //! Computes profiling datas (if enabled).
    Kerr_base & finalize() {system->finalize_profiling(); return *this;}
    //! Sends the profiling datas in the passed output stream.
    void profiling_log(std::ostream & os) {if(mpi_rank==0) profiling_report(*system,os);}

    void initialization_order_error(std::ostream & os,std::string const & file,int line) const {
        if(mpi_rank==0) {
            os << "Error : in file " << file << " at line " << line << ", bad initialization order."
               << std::endl;
        }
        throw std::runtime_error{"bad initialization order"};
    }
};

class Kerr;

class Kerr_init : public Kerr_base {
    friend class Kerr;
protected:
    ptr_data_member(Vector,stilde);
    ptr_data_member(Metric_flat,fmet);

public:
    Kerr_init(int nbr = 17,int ndom=3,double _bh_radius = 1.,int _type_coloc=CHEB_TYPE) :
        Kerr_base{nbr,ndom,_bh_radius,_type_coloc}, stilde{nullptr}, fmet{nullptr} {}

    Kerr_init & reset_initial_guess() override {
        this->Kerr_base::reset_initial_guess();
        stilde.reset(new Vector{*space,CON,*basis});
        for (int i=1 ; i<=3 ; i++) stilde->set(i) = 0. ;
        stilde->set(1).set_domain(1) = 1. ;
        stilde->set(2).set_domain(1) = 0. ;
        stilde->set(3).set_domain(1) = 0. ;
        stilde->std_base() ;
        fmet.reset(new Metric_flat{*space,*basis});
        return *this;
    }

    Kerr_init & reset_system() override {
        this->Kerr_base::reset_system();
        if(stilde && fmet) {
            assert(system);
            system->add_cst("s", *stilde);
            fmet->set_system(*system, "f");

            // definition of the extrinsic curvature :
            system->add_def ("A^ij = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) /2. / N") ;

            // Inner BC :
            space->add_inner_bc (*system, "N=n0") ;
            space->add_inner_bc (*system, "bet^i = n0 / P^2 * s^i") ; // 1
            space->add_inner_bc (*system, "dn(P) + 0.5 * P / a + P^3 * A_ij * s^i * s^j /4.= 0") ;

            // Equations :
            space->add_eq (*system, "D_i D^i N + 2 * D_i P * D^i N / P - N * P^4 * A_ij *A^ij= 0", "N", "dn(N)") ;
            space->add_eq (*system, "D_i D^i P + P^5 *A_ij * A^ij / 8= 0", "P", "dn(P)") ;
            space->add_eq (*system, "D_j A^ij + 6 * A^ij * D_j P / P =0", "bet^i", "dn(bet^i)") ;

            // Outer BC
            space->add_outer_bc (*system, "N=1") ;
            space->add_outer_bc (*system, "P=1") ;
            space->add_outer_bc (*system, "bet^i=0") ;

        } else initialization_order_error(std::cerr,__FILE__,__LINE__);
    }
};


class Kerr : public Kerr_base {
public:
    // Some of these integers seems to be linked to the number of domains, so they should not
    // be static constant, but rather internal data members. That said, I can't spend time
    // to investigate on this, so if someone wants to adapt, he is welcome.
    static constexpr int n_evol_inner {4};
    static constexpr int n_evol {5};
    static constexpr int n_dirac {1};

public:
    //! Steps for the increment of \c omega.
    double omega_step;
    //! Maximum number of increment for \c omega.
    int nbr_max_omega_val;
    internal_variable(int,count_omega_val);
    ptr_data_member(Metric_tensor,gfixed);
    bool save_to_file;

    ptr_data_member(Metric_tensor,gmet);
    ptr_data_member(Metric_dirac,met);
    ptr_data_member(Vector,scov);
    ptr_data_member(Vector,er);
    ptr_data_member(Vector,mm);

    Array<int> ** p_evol_inner;
    Array<int> ** p_evol;
    Array<int> ** p_dirac;

public:

    Kerr(Kerr_init & kerr_init) : Kerr_base{kerr_init.number_of_points(0),
                                            kerr_init.number_of_domains,
                                            kerr_init.bh_radius,kerr_init.type_coloc},
                                  omega_step{0.005}, nbr_max_omega_val{40},
                                  count_omega_val{0}, gfixed{nullptr}, save_to_file{false},
                                  gmet{nullptr}, met{nullptr}, scov{nullptr}, er{nullptr},
                                  mm{nullptr}, p_evol_inner{new Array<int>* [n_evol_inner]},
                                  p_evol{new Array<int>* [n_evol] },
                                  p_dirac{new Array<int>* [n_dirac] }
    {
        number_of_points = kerr_init.number_of_points;
        center = kerr_init.center;
        number_of_domains = kerr_init.number_of_domains;
        bounds = kerr_init.bounds;
        bh_radius = kerr_init.bh_radius;
        type_coloc = kerr_init.type_coloc;
        n0 = kerr_init.n0;
        omega = kerr_init.omega;
        // pointers are swapped, since kerr_init is now useless.
        space.swap(kerr_init.space);
        basis.swap(kerr_init.basis);
        conformal.swap(kerr_init.conformal);
        lapse.swap(kerr_init.lapse);
        shift.swap(kerr_init.shift);
        system.swap(kerr_init.system);

        newton_nbr_iterations = kerr_init.newton_max_iterations;
        tolerance = kerr_init.tolerance;

        // here again, the values may be number of domains dependants...
        for (int i=0 ; i<n_evol_inner ; i++) p_evol_inner[i] = new Array<int>{2} ;
        p_evol_inner[0]->set(0) = 1 ; p_evol_inner[0]->set(1) = 1;
        p_evol_inner[1]->set(0) = 1 ; p_evol_inner[1]->set(1) = 3 ;
        p_evol_inner[2]->set(0) = 2 ; p_evol_inner[2]->set(1) = 2 ;
        p_evol_inner[3]->set(0) = 2 ; p_evol_inner[3]->set(1) = 3 ;

        for (int i=0 ; i<n_evol ; i++) p_evol[i] = new Array<int>{2};
        p_evol[0]->set(0) = 1 ; p_evol[0]->set(1) = 1 ;
        p_evol[1]->set(0) = 1 ; p_evol[1]->set(1) = 2 ;
        p_evol[2]->set(0) = 1 ; p_evol[2]->set(1) = 3 ;
        p_evol[3]->set(0) = 2 ; p_evol[3]->set(1) = 2 ;
        p_evol[4]->set(0) = 2 ; p_evol[4]->set(1) = 3 ;

        for (int i=0 ; i<n_dirac ; i++) p_dirac[i] = new Array<int>{1} ;
        p_dirac[0]->set(0) = 2 ;

        this->reset_initial_guess();
    }

    ~Kerr() {
        if(p_evol_inner) {
            for (int i = 0; i < n_evol_inner; i++) safe_delete(p_evol_inner[i]);
            delete[] p_evol_inner;
        }
        if(p_evol) {
            for (int i = 0; i < n_evol; i++) safe_delete(p_evol[i]);
            delete[] p_evol;
        }
        if(p_dirac) {
            for (int i = 0; i < n_dirac; i++) safe_delete(p_dirac[i]);
            delete[] p_dirac;
        }
    }

    Kerr & reset_initial_guess() override {
        gfixed.reset(new Metric_tensor{*space, CON, *basis}) ;
        for (int i=1 ; i<=3 ; i++)
            for (int j=i ; j<=3 ; j++)
                if (i==j) gfixed->set(i,j) = 1. ;
                else gfixed->set(i,j).annule_hard() ;
        gfixed->std_base() ;
        gmet.reset(new Metric_tensor{*gfixed});
        met.reset(new Metric_dirac{*gmet});
        scov.reset(new Vector{*space, COV, *basis}) ;
        scov->set(1) = 1. ;
        scov->set(2) = 0. ;
        scov->set(3) = 0. ;
        scov->std_base() ;

        er.reset(new Vector{*space, CON, *basis}) ;
        er->set(1) = 1. ;
        er->set(2) = 0. ;
        er->set(3) = 0. ;

        mm.reset(new Vector{*space, CON, *basis}) ;
        for (int i=1 ; i<=3 ; i++)
            mm->set(i) = 0. ;
        Val_domain xx (space->get_domain(1)->get_cart(1)) ;
        Val_domain yy (space->get_domain(1)->get_cart(2)) ;
        mm->set(3).set_domain(1) = sqrt(xx*xx + yy*yy) ;
        mm->std_base() ;
        return *this;
    }

    Kerr & reset_system() override {
        // Solve the equation in space outside the nucleus
        system.reset(new System_of_eqs{*space, 1, number_of_domains-1}) ;
        // Unknowns
        system->add_var ("P", *conformal) ;
        system->add_var ("N", *lapse) ;
        system->add_var ("bet", *shift) ;
        met->set_system (*system, "g") ;

        // User defined constants
        system->add_cst ("a", bh_radius) ;
        system->add_cst ("m", *mm) ;
        system->add_cst ("s", *scov) ;
        system->add_cst ("n0", n0) ;
        system->add_cst ("Ome", omega) ;
        system->add_cst ("gf", *gfixed) ;

        // definitions
        // For speed one stores derivatives of the CF fields :
        system->add_def ("DN_i = D_i N") ;
        system->add_def ("DP_i = D_i P") ;
        system->add_def ("Dbet^ij = D^i bet^j") ;

        system->add_def ("st^i = s^i / sqrt(s_i * s^i)") ;
        system->add_def ("A^ij = (Dbet^ij + Dbet^ji - 2. / 3.* Dbet_k^k * g^ij)/2. / N ") ;
        system->add_def ("LieK_ij = 4 * A_ij * bet^k * DP_k / P + bet^k * D_k A_ij + A_ik * Dbet_j^k + A_jk * Dbet_i^k") ;
        system->add_def ("DDN_ij = D_i DN_j - 2 * DN_i * DP_j / P - 2 * DN_j * DP_i / P + 2 * g_ij * DN_k * DP^k / P") ;
        system->add_def ("PartR_ij = R_ij + 6 * DP_i * DP_j / P^2 - 2 * D_i DP_j / P - 2 * g_ij * D_k DP^k / P - 2 * g_ij * DP_k * DP^k / P^2 - P^4 * 2 * A_ik * A_j^k") ;
        system->add_def ("evol_ij = DDN_ij - N * PartR_ij - P^4 * LieK_ij") ;
        system->add_def ("Pfourhor = N^2 / ( s_i * bet^i ) / ( s_i * bet^i )") ;

        space->add_inner_bc (*system, "N=n0") ;
        space->add_inner_bc (*system, "bet^i = n0 / P^2 * st^i + Ome * m^i * a") ;
        space->add_inner_bc (*system, "4 * st^i * D_i P / P + D_i st^i + P^2 * A_ij * st^i * st^j = 0") ;
        space->add_inner_bc (*system, "DDN_ij - N * PartR_ij - Pfourhor * LieK_ij=0", n_evol_inner, p_evol_inner) ;
        space->add_inner_bc (*system, "dirac^i =0", n_dirac, p_dirac) ;

        // CFC Equations :
        space->add_eq (*system, " D_i DN^i + 2 * DP_i * DN^i / P - N * P^4 * A_ij * A^ij = 0", "N", "dn(N)") ;
        space->add_eq (*system, "R - 8 * D_i DP^i / P - P^4 * A_ij * A^ij = 0", "P", "dn(P)") ;
        space->add_eq (*system, "D^j A_ij + 6 * A_ij * DP^j / P =0", "bet^i", "dn(bet^i)") ;

        // Evolution
        system->add_eq_inside (1, "evol_ij =0", n_evol, p_evol) ;
        for (int d=2 ; d<number_of_domains ; d++) {
            system->add_eq_matching (d-1, OUTER_BC, "g^ij", n_evol, p_evol) ;
            system->add_eq_matching (d-1, OUTER_BC, "dn(g^ij)", n_evol, p_evol) ;
            system->add_eq_inside (d, "evol_ij=0", n_evol, p_evol) ;
        }


        space->add_eq_full (*system, "determinant(g^ij) = 1") ;

        // Outer BC
        space->add_outer_bc (*system, "N=1") ;
        space->add_outer_bc (*system, "P=1") ;
        space->add_outer_bc (*system, "bet^i=0") ;
        space->add_outer_bc (*system, "g^ij=gf^ij", n_evol, p_evol) ;

        newton_nbr_iterations = 0;
        newton_residue = HUGE_VAL;
        return *this;
    }

    Kerr & build_space_and_system() override {
        this->reset_initial_guess();
        this->reset_system();
        return *this;
    }

    //made virtual to add rank condition for the MPI version.
    virtual void save(char * file_name) {
        FILE* ff = fopen(file_name, "w") ;
        space->save(ff) ;
        fwrite_be (&n0, sizeof(double), 1, ff) ;
        fwrite_be (&omega, sizeof(double), 1, ff) ;
        fwrite_be (&bh_radius, sizeof(double), 1, ff) ;
        conformal->save(ff) ;
        lapse->save(ff) ;
        shift->save(ff) ;
        gmet->save(ff) ;
        fclose(ff) ;
    }

    //! Performs Newton's method for the current value of \c omega.
    bool do_newton() override {
        bool newton_success {false};
        char name[100] ;
        sprintf(name, "kerr_%d_%f.dat", number_of_points(0), omega) ;
        bool const do_not_check_iter {newton_max_iterations < 0};
        if(mpi_rank==0) std::cout << "Computation with omega = " << omega << std::endl;
        while(!newton_success &&
              (do_not_check_iter || newton_nbr_iterations <= newton_max_iterations)) {
            newton_success = system->do_newton(tolerance,newton_residue);
            newton_nbr_iterations++;
            if(save_to_file && mpi_rank==0) this->save(name);
        }
    }

    /**
     * Increment the value of \c omega if possible.
     * @return true if the increment is done, false if the max number of iterations with
     * respect to \c omega has been reached.
     */
    bool increment_omega() {
        if (count_omega_val < nbr_max_omega_val) {
            omega += omega_step;
            count_omega_val++;
            return true;
        } else return false;
    }


};



#endif //__KADATH_CODES_KERR_HPP_
