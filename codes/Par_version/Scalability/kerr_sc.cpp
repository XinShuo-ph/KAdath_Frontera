//
// Created by sauliac on 29/05/2020.
//

#include "kerr.hpp"
#include "mpi.h"
#include "magma_interface.hpp"

//if the profiler does not keep a record of the computational times, nothing can be done...
#ifdef ENABLE_INTERNAL_PROFILER
constexpr char matrix_computation_str[] = "matrix_computation";

using Profiler_stat = System_of_eqs::Statistics;

struct Profiler_data {
    int problem_size;
    Profiler_stat statistics;
    Profiler_data() = default;
    Profiler_data(int ps,Profiler_stat const & stat) :
        problem_size{ps}, statistics{stat} {}
    Profiler_data & operator+=(Profiler_stat const & another_data) {
        unsigned long const n1{statistics.n_samples}, n2{another_data.n_samples};
        unsigned long const ntot{n1+n2};
        double const T1 {statistics.total_duration}, T2{another_data.total_duration};
        double const mu1{statistics.average_duration}, mu2{another_data.average_duration};
        double const sigma1_2{statistics.std_deviation}, sigma2_2{another_data.std_deviation};
        statistics.total_duration = T1 + T2;
        statistics.n_samples = ntot;
        statistics.average_duration = (n1*mu1 + n2*mu2) / ntot;
        statistics.std_deviation = (n1*sigma1_2 + n2*sigma2_2 + n1*n2*(mu1-mu2)*(mu1-mu2)/ntot)/ntot;
        return *this;
    }
};

class Stat_extractor {
public:
    static constexpr int all_sizes {-1};
    static constexpr bool do_not_print_data_header {false}, print_data_header{true};
    using MPI_comm_size_t = int;
    using Kadath_system_size_t = int;
    using Nproc_data_map = std::map<MPI_comm_size_t,Profiler_data>;
    using Complete_data_map = std::map<Kadath_system_size_t ,Nproc_data_map>;
private:
    MPI_comm_size_t nb_mpi_process;
    std::string extraction_key;
    Complete_data_map extracted_stats;
    std::string data_file_name_prefix;
    std::map<int,std::ofstream> data_files;

public:
    int get_nb_mpi_process() const {return nb_mpi_process;}
    Stat_extractor & set_nb_mpi_process(MPI_comm_size_t new_value) {nb_mpi_process = new_value; return *this;}
    std::string const & get_extraction_key() const {return extraction_key;}
    Stat_extractor & set_extraction_key(std::string const & new_key) {extraction_key = new_key; return *this;}
    std::string const & get_data_file_name_prefix() const {return data_file_name_prefix;}
    Stat_extractor & set_data_file_name_prefix(std::string const & new_prefix) {data_file_name_prefix = new_prefix; return *this;}


    Stat_extractor(MPI_comm_size_t mpi_comm_size, std::string const &key,
                    std::string const & file_name_prefix = "kerr_scal") :
            nb_mpi_process{mpi_comm_size}, extraction_key{key} ,
            extracted_stats{}, data_file_name_prefix{file_name_prefix}, data_files{} {}

    void print(std::ostream & os,bool print_header=print_data_header,int size = all_sizes);
    static void spwan_header(std::ostream & os);
    void save();
    void extract();
};



int main(int argc,char ** argv) {
    MPI_Init(&argc, &argv) ;
    int mpi_world_size{0};
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_world_size);
    int mpi_proc_rank{0};
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_proc_rank) ;
    Arguments_parser arg_parser{argc, argv};

//    auto min_nb_mpi_proc = arg_parser.get_option_value<int>("-npmin","Sets the starting number of mpi process",1);
    auto data_file_prefix = arg_parser.get_option_value<std::string>("-f", "Sets the output file name prefix.", "kerr_scal");
    auto max_iterations = arg_parser.get_option_value<int>("-niter",
                                                           "Sets the maximum number of iteration for Newton-Raphson method. A null or negative sets this limit to infinity.",
                                                           0);
    auto max_nb_omega = arg_parser.get_option_value<int>("-nomega",
                                                         "Sets the number of increments toward Omega to perform.",
                                                         40);
    auto nb_points = arg_parser.get_option_value<int>("-npts",
                                                      "Sets the number of collocation points (note that this value is constraint by the spectral method).",
                                                      17);
    auto save_to_file = arg_parser.get_option_value<bool>("-s","Enables / disable data save to file.",true);
    auto verbosity_level = arg_parser.get_option_value<int>("-v", "Sets the verbosity level", 2);
    bool const show_help{arg_parser.find_option("-h", "Display this help message.")};
    if (show_help) {
        if (mpi_proc_rank == 0) arg_parser.display(std::cout);
        return 0;
    }
    Stat_extractor matrix_computation_time_extractor{mpi_world_size, matrix_computation_str,
                                                     data_file_prefix.first};


    Kerr_init kerr_init{nb_points.first};
    kerr_init.mpi_rank = mpi_proc_rank;
    kerr_init.set_verbosity(verbosity_level.first);
    if (max_iterations.second) kerr_init.newton_max_iterations = max_iterations.first;

    // build all internal data.
    kerr_init.build_space_and_system();

    kerr_init.do_newton();

    kerr_init.finalize();

    Kerr kerr{kerr_init};

    kerr.mpi_rank = mpi_proc_rank;
    kerr.nbr_max_omega_val = max_nb_omega.second ? max_nb_omega.first : 40;
    kerr.reset_initial_guess();

    while (kerr.increment_omega()) {
        //re-build the system with the new value of omega
        kerr.reset_system();
        // perform Newton-Raphson method for this value
        kerr.do_newton();
    }
    kerr.finalize();


    if(mpi_proc_rank == 0) System_of_eqs::display(std::cout);
    matrix_computation_time_extractor.extract();
    if(mpi_proc_rank == 0) {
        std::cout << std::endl << std::endl;
        matrix_computation_time_extractor.print(std::cout);
        if(save_to_file.first) matrix_computation_time_extractor.save();
    }

    MPI_Finalize() ;

    return 0;
}

void Stat_extractor::save() {
    for(auto const & data_block : extracted_stats) {
        int const problem_size {data_block.first};
        std::string file_name{data_file_name_prefix + "_"};
        file_name += std::to_string(problem_size) + ".txt";
        std::ifstream file{file_name};
        if(file.is_open()) {
            bool end_of_file{false};
            while (!end_of_file) {
                Profiler_stat stat;
                stat.user_key = "matrix_computation";
                int size, nproc;
                file >> size;
                file >> nproc;
                file >> stat.average_duration;
                file >> stat.std_deviation;
                file >> stat.n_samples;
                end_of_file = file.eof();
                if(!end_of_file) {
                    assert(size == problem_size);
                    stat.std_deviation = stat.std_deviation*stat.std_deviation;
                    stat.total_duration = stat.average_duration * stat.n_samples;
                    auto &data_map = extracted_stats[size];
                    auto entry = data_map.emplace(nproc, Profiler_data{size, stat});
                    if (!entry.second) {
                        assert(entry.first->first == nproc);
                        entry.first->second += stat;
                    }
                }
            }
            file.close();
        }
        auto pos = data_files.try_emplace(problem_size,std::ofstream{file_name});
        this->print(pos.first->second,do_not_print_data_header,problem_size);
    }
}

void Stat_extractor::spwan_header(std::ostream &os) {
    //                                                               <--------20-------->1<--------20-------->1<------16------>
    //     <------16------>   +  1   +  <------16------>   +  1  +   <------------------|-|---58-------------|-|-------------->
    os << "     System_size" << ' ' << "  nb_mpi_process" << ' ' << "                    matrix_computation                    " << std::endl;
    os << "                                  "                   << "            mean (s)          std_dev (s)       nb_samples" << std::endl;

}

void Stat_extractor::print(std::ostream &os,bool print_header,int size)  {
    if(print_header) spwan_header(os);
    for(auto const & data_map : extracted_stats) {
        if(print_header) os << std::endl;
        os << std::right;
        auto const psize {data_map.first};
        if(size == all_sizes || psize == size) {
            for (auto const &e : data_map.second) {
                assert(e.second.problem_size == psize);
                os << std::setw(16) << e.second.problem_size << ' ';
                os << std::setw(16) << e.first << ' ';
                os << std::setw(20) << e.second.statistics.average_duration << ' ';
                os << std::setw(20) << std::sqrt(e.second.statistics.std_deviation) << ' ';
                os << std::setw(16) << e.second.statistics.n_samples << " \n";
            }
        }
        if(print_header) os << std::endl;
   }
}

void Stat_extractor::extract() {
    auto const & stat_map = System_of_eqs::get_statistic_map();
    for(auto const & e : stat_map) {
        auto const & entry_stats = e.second;
        std::string const &key = entry_stats.user_key;
        auto pos = key.find(extraction_key);
        if(pos != std::string::npos) {
            std::stringstream ss{key};
            Kadath_system_size_t problem_size{};
            std::string operation,buffer;
            ss >> buffer >> buffer >> problem_size >> operation;
            assert(operation == extraction_key);
            auto & data_map = extracted_stats[problem_size];
            auto entry = data_map.emplace(nb_mpi_process,Profiler_data{problem_size,entry_stats});
            if(!entry.second) {
                assert(entry.first->first == nb_mpi_process);
                entry.first->second += entry_stats;
            }
        }
    }
}

#else
int main(int argc,char ** argv) {
    MPI_Init(&argc, &argv) ;
    int this_proc_world_rank{0};
    MPI_Comm_rank(MPI_COMM_WORLD, &this_proc_world_rank) ;
    if(this_proc_world_rank == 0) {
        std::cout << "Error : this executable must be build with a version of Kadath compiled with the"
                     " ENABLE_INTERNAL_PROFILER option turned On." << std::endl;
    }
    MPI_Finalize() ;
    return 0;
}
#endif