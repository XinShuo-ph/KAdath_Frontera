//
// Created by sauliac on 20/01/2020.
//

#include "gtest/gtest.h"
#include <utility>
#include <random>
#include "multi_thread.hpp"
#include "kadath_spheric.hpp"


constexpr double TESTS_TOLERANCE {1.e-8};

using namespace Kadath;

enum : bool
{
    DO_NOT_COMPUTE = false,
    COMPUTE = true
};

template<class System,int R1=13,int R2=5,int R3=4> class Schwarz_test
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

    Schwarz_test(std::size_t num_threads = 1,bool compute = COMPUTE)
        : res{dim}, bounds{ndom-1}, p_space{nullptr}, p_conf{nullptr},
                                                p_syst{nullptr}, converged{false}, newton_error{-1.}
    {
        initialize(num_threads);
        if(compute) {
            do_newton();
            check_solution();
        }
    }

    bool has_converged() const {return converged;}

};

template<class System,int A,int B,int C> bool Schwarz_test<System,A,B,C>::initialize(std::size_t num_threads)
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

template<class System,int A,int B,int C>
void Schwarz_test<System,A,B,C>::do_newton()
{
    unsigned short niter{0};
    while (!converged && niter<2) {
        converged = p_syst->do_newton(newton_tol, newton_error);
        niter++;
    }
}

template<class System,int A,int B,int C> void Schwarz_test<System,A,B,C>::check_solution()
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
    EXPECT_TRUE(converged);
    EXPECT_LE(error_l_infinity,TESTS_TOLERANCE);
    EXPECT_LE(error_l_2,TESTS_TOLERANCE);
}
//
//TEST_F(KadathSchwarzMultiThreadTest,CheckErrorInLInfAndL2Norms)
//{
//    double matrix_diff_l2 {0.};
//    double matrix_t_diff_l2 {0.};
//    Array<double> const * const matrix_single_thread {single_thread.get_matrix()};
//    Array<double> const * const matrix_multi_thread {multi_thread.get_matrix()};
//    if(matrix_single_thread && matrix_multi_thread)
//    {
//        EXPECT_EQ(matrix_single_thread->get_dimensions()(0),matrix_multi_thread->get_dimensions()(0));
//        EXPECT_EQ(matrix_single_thread->get_dimensions()(1),matrix_multi_thread->get_dimensions()(1));
//        auto N = matrix_single_thread->get_dimensions()(0);
//        auto M = matrix_multi_thread->get_dimensions()(1);
//        EXPECT_EQ(N,M);
//        for(int i=0;i<N;i++) {
//            for (int j = 0; j < M; j++) {
//                double const diff{(*matrix_single_thread)(i, j) - (*matrix_multi_thread)(i, j)};
//                double const diff_t{(*matrix_single_thread)(i, j) - (*matrix_multi_thread)(j, i)};
//                matrix_diff_l2 += (diff * diff);
//                matrix_t_diff_l2 += (diff_t * diff_t);
//            }
//        }
//        EXPECT_NEAR(matrix_diff_l2,0.,TESTS_TOLERANCE);
//        if(matrix_t_diff_l2 <= TESTS_TOLERANCE)
//        {
//            std::cout << "A matrix has been transposed." << std::endl;
//        }
//    }
//    EXPECT_TRUE(single_thread.has_converged());
//    EXPECT_LE(single_thread.get_error_l_infinity(),TESTS_TOLERANCE);
//    EXPECT_LE(single_thread.get_error_l_2(),TESTS_TOLERANCE);
//    EXPECT_TRUE(multi_thread.has_converged());
//    EXPECT_LE(multi_thread.get_error_l_infinity(),TESTS_TOLERANCE);
//    EXPECT_LE(multi_thread.get_error_l_2(),TESTS_TOLERANCE);
//}


class KadathMultiThreadTests : public ::testing::Test
{
public:
    static constexpr std::size_t num_threads {5};
    static constexpr std::size_t num_col_checks {25};
protected:
    Schwarz_test<System_of_eqs> sequential;
    Schwarz_test<System_of_eqs_threaded> multi_thread;
public:
    KadathMultiThreadTests() : Test{}, sequential{1,DO_NOT_COMPUTE},
                               multi_thread{num_threads,DO_NOT_COMPUTE} {}

    bool check_second_members_seq();
    bool check_second_members_par();

    bool check_do_col_j_seq();
    bool check_do_col_j_par();
};


bool KadathMultiThreadTests::check_second_members_seq() {
    //First checks what happens when second members is called twice in a row...
    Array<double> sm_from_sequential{sequential.p_syst->sec_member()},
            sm_from_sequential_b{sequential.p_syst->sec_member()};
    Array<double> diff_sequential{fabs(sm_from_sequential - sm_from_sequential_b)};
    double const error_sequential{max(diff_sequential)};
    EXPECT_EQ(error_sequential, 0.);

    Array<double> sm_from_multi_thread{multi_thread.p_syst->sec_member()};
    Array<double> sm_from_multi_thread_b{multi_thread.p_syst->sec_member()};
    auto diff_seq_mth = fabs(sm_from_sequential - sm_from_multi_thread);
    auto diff_seq_mth_b = fabs(sm_from_sequential - sm_from_multi_thread_b);
    double const err_seq_mth{max(diff_seq_mth)};
    double const err_seq_mth_b{max(diff_seq_mth_b)};
    EXPECT_EQ(err_seq_mth, 0.);
    EXPECT_EQ(err_seq_mth_b, 0.);

    System_of_eqs_threaded::container &systems = multi_thread.p_syst->get_thread_specific_systems();
    std::vector<std::unique_ptr<Array<double>>> second_members(num_threads - 1);
    for (std::size_t i{0}; i < num_threads - 1; i++) {
        second_members.at(i).reset(new Array<double>{systems.at(i)->sec_member()});
    }
    for (auto const &sm : second_members) {
        Array<double> diff_seq_sm{fabs(*sm - sm_from_sequential)};
        double const err{max(diff_seq_sm)};
        EXPECT_EQ(err, 0.);
    }
}

bool KadathMultiThreadTests::check_second_members_par()
{
    Array<double> sm_from_sequential{sequential.p_syst->sec_member()};
    Array<double> sm_from_multi_thread{multi_thread.p_syst->sec_member()};
    EXPECT_EQ(max(fabs(sm_from_sequential-sm_from_multi_thread)),0.);
    System_of_eqs_threaded::container &systems = multi_thread.p_syst->get_thread_specific_systems();
    std::vector<std::unique_ptr<Array<double>>> second_members(num_threads - 1);
    std::deque<std::thread> pool;
    for(std::size_t i{0};i<num_threads-1;i++)
    {
        pool.emplace_back([&second_members,&systems,i]()
                            {second_members.at(i).reset(new Array<double>{systems.at(i)->sec_member()});});
    }
    for(auto & p : pool) p.join();
    for(std::size_t i{0};i<num_threads-1;i++)
    {
        Array<double> diff_seq_sm {fabs(*second_members.at(i) - sm_from_sequential)};
        double const err {max(diff_seq_sm)};
        EXPECT_EQ(err,0.);
    }
}

bool KadathMultiThreadTests::check_do_col_j_seq()
{
    using VUPAD = std::vector<std::unique_ptr<Array<double>>>;
    Array<double> second_member{multi_thread.p_syst->sec_member()};
    Array<double> smc{sequential.p_syst->sec_member()};
    EXPECT_EQ(max(fabs(second_member-smc)),0.);

    std::random_device rd;
    std::mt19937 gen{rd()};
    std::uniform_int_distribution<> dis{0,second_member.get_size(0)};

    std::vector<int> idx_to_check(num_col_checks);
    for(auto & x : idx_to_check) x = dis(gen);

    VUPAD correct_values(num_col_checks);
    std::vector<VUPAD> sequentially_obtained_values(num_threads);
    std::vector<VUPAD> simultaneously_obtained_values(num_threads);
    for(std::size_t j{0};j<num_col_checks;j++)
    {
        correct_values.at(j).reset(new Array<double>{sequential.p_syst->do_col_J(idx_to_check.at(j))});
    }

    //Sequential call :
    for(std::size_t i{0};i<num_threads;i++)
    {
        System_of_eqs & syst {i == 0 ? *multi_thread.p_syst
                : *multi_thread.p_syst->get_thread_specific_systems().at(i-1) };
        sequentially_obtained_values.at(i).resize(num_col_checks);
        for(std::size_t j{0};j<num_col_checks;j++)
        {
            sequentially_obtained_values.at(i).at(j).reset(new Array<double>{syst.do_col_J(idx_to_check.at(j))});
        }
    }
    for(std::size_t i{0};i<num_threads;i++)
    {
        for(std::size_t j{0};j<num_col_checks;j++)
        {
            EXPECT_EQ(max(fabs(*sequentially_obtained_values.at(i).at(j) - *correct_values.at(j) )),0.);
        }
    }
}

bool KadathMultiThreadTests::check_do_col_j_par()
{
    using VUPAD = std::vector<std::unique_ptr<Array<double>>>;
    Array<double> second_member{multi_thread.p_syst->sec_member()};
    Array<double> smc{sequential.p_syst->sec_member()};
    EXPECT_EQ(max(fabs(second_member-smc)),0.);

    std::random_device rd;
    std::mt19937 gen{rd()};
    std::uniform_int_distribution<> dis{0,second_member.get_size(0)};

    std::vector<int> idx_to_check(num_col_checks);
    for(auto & x : idx_to_check) x = dis(gen);

    VUPAD correct_values(num_col_checks);
    std::vector<VUPAD> sequentially_obtained_values(num_threads);
    std::vector<VUPAD> simultaneously_obtained_values(num_threads);
    for(std::size_t j{0};j<num_col_checks;j++)
    {
        correct_values.at(j).reset(new Array<double>{sequential.p_syst->do_col_J(idx_to_check.at(j))});
    }

    //Parallel call.
    std::deque<std::thread> pool;
    for(std::size_t i{0};i<num_threads;i++)
    {
        System_of_eqs & syst {i == 0 ? *multi_thread.p_syst
                                     : *multi_thread.p_syst->get_thread_specific_systems().at(i-1) };
        pool.emplace_back([i,&syst,&simultaneously_obtained_values,&idx_to_check](){
            simultaneously_obtained_values.at(i).resize(num_col_checks);
            for(std::size_t j{0};j<num_col_checks;j++)
            {
                simultaneously_obtained_values.at(i).at(j).reset(
                        new Array<double>{syst.do_col_J(idx_to_check.at(j))});
            }});
    }
    for(auto & p : pool) p.join();
    for(std::size_t i{0};i<num_threads;i++)
    {
        for(std::size_t j{0};j<num_col_checks;j++)
        {
            EXPECT_EQ(max(fabs(*simultaneously_obtained_values.at(i).at(j) - *correct_values.at(j) )),0.);
        }
    }

}


TEST_F(KadathMultiThreadTests,CheckSecondMemberSequential)
{
    this->check_second_members_seq();
}
//TEST_F(KadathMultiThreadTests,CheckSecondMemberParallel)
//{
//    this->check_second_members_par();
//}

TEST_F(KadathMultiThreadTests,CheckDoColJSequential)
{
    this->check_do_col_j_seq();
}

//TEST_F(KadathMultiThreadTests,CheckDoColJParallel)
//{
//    this->check_do_col_j_par();
//}