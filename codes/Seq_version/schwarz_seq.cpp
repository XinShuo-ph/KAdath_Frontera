#include "utility.hpp"
#include "schwarz.hpp"

using namespace Kadath ;

int main(int argc,char** argv) {

    Arguments_parser arg_parser{argc,argv};
    auto max_iterations = arg_parser.get_option_value<int>("-iter");
    auto nb_points = arg_parser.get_option_value<int>("-npts");
    if(!nb_points.second) nb_points.first = 13;

    Schwarz schwarz_solver{nb_points.first};
    if(max_iterations.second) schwarz_solver.set_newton_max_iterations(max_iterations.first);

    schwarz_solver.build_space_and_system();

    schwarz_solver.do_newton();

    schwarz_solver.check_solution();

    cout << "Error max : " << schwarz_solver.get_error_l_infinity() << endl ;
    cout << "Error L2  : " << schwarz_solver.get_error_l_2() << endl;

    schwarz_solver.finalize();

    return EXIT_SUCCESS ;
}

