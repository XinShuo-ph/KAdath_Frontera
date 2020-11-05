/*
    Copyright 2020 sauliac

    This file is part of Kadath.

    Kadath is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Kadath is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Kadath.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "solvers.hpp"

using namespace Kadath;

int main(int argc, char * argv[]) {

    std::cout << "============================= t-solvers unit-tests set =============================\n\n";
    Solver solver;
    solver.display(std::cout);
    assert(solver.get_target_error() == Tolerance.default_value);
    assert(solver.get_target_nb_iteration() == MaxNbIter.default_value);
    assert(solver.get_target_duration() == MaxElapsedTime.default_value);
    assert(solver.get_target_error_decrease() == MinImprovement.default_value);
    assert(solver.get_output());
    assert(!solver.get_enable_gpu());

    solver.set(MaxNbIter = 10, DisableOutput);
    solver.display(std::cout);
    assert(solver.get_target_error() == Tolerance.default_value);
    assert(solver.get_target_nb_iteration() == 10);
    assert(solver.get_target_duration() == MaxElapsedTime.default_value);
    assert(solver.get_target_error_decrease() == MinImprovement.default_value);
    assert(!solver.get_output());
    assert(!solver.get_enable_gpu());

    Solver another_solver{EnableGPU,MaxElapsedTime = 666.,Tolerance = 0.01,MinImprovement = 0.0001,EnableOutput};
    another_solver.display(std::cout);
    assert(another_solver.get_target_error() == 0.01);
    assert(another_solver.get_target_nb_iteration() == MaxNbIter.default_value);
    assert(another_solver.get_target_duration() == 666.);
    assert(another_solver.get_target_error_decrease() == 0.0001);
    assert(another_solver.get_output());
    assert(another_solver.get_enable_gpu());
    std::cout << "\n\n====================================================================================";

}