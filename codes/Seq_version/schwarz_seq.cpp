#include "codes_utilities.hpp"
#include "schwarz.hpp"

using namespace Kadath ;

int main(int argc,char** argv) {

    Arguments_parser arg_parser{argc,argv};
    auto max_iterations = arg_parser.get_option_value<int>("-niter","Sets the maximum number of iteration for Newton-Raphson method. A null or negative sets this limit to infinity.",0);
    auto nb_points = arg_parser.get_option_value<int>("-npts","Sets the number of collocation points (note that this value is constraint by the spectral method).",13);
    bool const show_help {arg_parser.find_option("-h","Display this help message.")};
    if(show_help) {
        arg_parser.display(std::cout);
        return 0;
    }
    Schwarz schwarz_solver{nb_points.first};
    if(max_iterations.second) schwarz_solver.newton_max_iterations = max_iterations.first;

    schwarz_solver.initialize();

    schwarz_solver.do_newton();

    schwarz_solver.check_solution();

    cout << "Error max : " << schwarz_solver.get_error_l_infinity() << endl ;
    cout << "Error L2  : " << schwarz_solver.get_error_l_2() << endl;

    schwarz_solver.finalize();

    return EXIT_SUCCESS ;
}

