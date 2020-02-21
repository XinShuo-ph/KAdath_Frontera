
/*
    Copyright 2017 Philippe Grandclement

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
#include <config.h>

#ifdef PAR_VERSION
#include "mpi.h"
#endif

#include "system_of_eqs.hpp"
#include "matrice.hpp"
#include "scalar.hpp"
#include "array_math.cpp"




namespace Kadath {
    template<>
    bool System_of_eqs::do_newton<Computational_model::sequential>(double precision, double& error,std::ostream &os)
    {
#ifdef PAR_VERSION
        int rank;
        MPI_Comm_rank (MPI_COMM_WORLD, &rank);
        if(rank==0) {
#endif
            niter++;
            if(niter==1 && display_newton_data)
            {
                display_do_newton_report_header(os,precision);
            }
            clock_t begin, end;
            vars_to_terms();
            Array<double> second(sec_member());
            error = max(fabs(second));
            if (error < precision)
            {
                if(display_newton_data)
                {
                    display_do_newton_ending_line(os,precision,error);
                    os  << endl;
                }
                return true;
            } else {
                int nn(second.get_size(0));
                if (nbr_unknowns != nn) {
                    cerr << "N unknowns  = " << nbr_unknowns << endl;
                    cerr << "N equations = " << nn << endl;
                    abort();
                }
                cout.flush();

                Hash_key chrono_key = this->start_chrono("do_newton | problem size = ", nn, " | matrix computation ");
                Matrice ope(nn, nn);
                compute_matrix(ope.get_array(), nn);
                Duration const
                        t_load_matrix{this->stop_chrono(chrono_key)};


                chrono_key = this->start_chrono("do_newton | problem size = ", nn, " | matrix inversion ");
                ope.set_lu();
                Array<double> xx(ope.solve(second));
                Duration const
                        t_inv_matrix{this->stop_chrono(chrono_key)};

                chrono_key = this->start_chrono("do_newton | problem size = ", nn, " | Newton update ");
                newton_update_vars(xx);
                Duration const t_newton_update
                        {this->stop_chrono(chrono_key)};
                return false;
            }
#ifdef PAR_VERSION
        }
        else {
            return false;
        }
#endif
    }

}

