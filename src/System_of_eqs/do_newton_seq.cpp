
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
            clock_t begin, end;
            vars_to_terms();
            Array<double> second(sec_member());
            error = max(fabs(second));
            os << "Error init = " << error << endl;
            if (error < precision) return true;
            int nn(second.get_size(0));
            if (nbr_unknowns != nn) {
                cerr << "N unknowns  = " << nbr_unknowns << endl;
                cerr << "N equations = " << nn << endl;
                abort();
            }
            os << "Size = " << nn << endl;
            os << "Progression computation : ";
            cout.flush();

            Hash_key chrono_key = this->start_chrono("do_newton | problem size = ", nn, " | matrix computation ");
            Array<double> jx(nn);
            Matrice ope(nn, nn);
            for (int col(nn - 1); col >= 0; col--) {
                if ((col != 0) && (col % int(nn / 10) == 0)) {
                    os << "*";
                    os.flush();
                }
                jx = do_col_J(col);
                for (int line(0); line < nn; line++) ope.set(line, col) = jx(line);
            }
            os << endl;
            Duration const
                    t_load_matrix{this->stop_chrono(chrono_key)};
            os << "Loading the matrix : " << to_seconds(t_load_matrix) << " seconds" << endl;

            chrono_key = this->start_chrono("do_newton | problem size = ", nn, " | matrix inversion ");
            ope.set_lu();
            Array<double> xx(ope.solve(second));
            Duration const
                    t_inv_matrix{this->stop_chrono(chrono_key)};
            os << "Inverting the matrix : " << to_seconds(t_inv_matrix) << " seconds" << endl;
            int conte(0);

            chrono_key = this->start_chrono("do_newton | problem size = ", nn, " | Newton update ");
            espace.xx_to_vars_variable_domains(this, xx, conte);
            double *old_var_double(new double[nvar_double]);
            for (int i(0); i < nvar_double; ++i) old_var_double[i] = *var_double[i];
            Tensor **old_fields(new Tensor *[nvar]);
            for (int i(0); i < nvar; i++) old_fields[i] = new Tensor(*var[i]);
            xx_to_vars(xx, conte);
            for (int i(0); i < nvar_double; i++) *var_double[i] = old_var_double[i] - (*var_double[i]);
            for (int i(0); i < nvar; i++) *var[i] = *old_fields[i] - (*var[i]);
            delete[] old_var_double;
            for (int i(0); i < nvar; i++) delete old_fields[i];
            delete[] old_fields;
            Duration const t_newton_update
                    {this->stop_chrono(chrono_key)};
            os << "Newton update : " << to_seconds(t_newton_update) << " seconds" << endl;
#ifdef PAR_VERSION
        }
#endif
        return false;
    }

}

