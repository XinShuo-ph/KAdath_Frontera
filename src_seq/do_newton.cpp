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

#include "system_of_eqs.hpp"
#include "matrice.hpp"
#include "scalar.hpp"
#include "metric.hpp"
#include "array_math.cpp"

#include <ctime>
namespace Kadath {
bool System_of_eqs::do_newton(double precision, double& error) 
{
   clock_t begin, end;
   vars_to_terms();
   Array<double> second(sec_member());
   error = max(fabs(second));
   cout << "Error init = " << error << endl;
   if (error<precision) return true;
   int nn(second.get_size(0));
   if (nbr_unknowns!=nn) 
   {
      cerr << "N unknowns  = " << nbr_unknowns << endl;
      cerr << "N equations = " << nn << endl;
      abort();
   }
   cout << "Size = " << nn << endl;
   cout << "Progression computation : ";
   cout.flush();
   begin = clock();
   Array<double> jx(nn);
   Matrice ope(nn,nn);
   for (int col(nn-1) ; col >= 0 ; col--) 
   {
      if ((col != 0) && (col%int(nn/10) == 0)) 
      {
         cout << "*";
         cout.flush();
      }
      jx = do_col_J(col);
      for (int line(0) ; line < nn ; line++) ope.set(line,col) = jx(line);
   }
   cout << endl;
   end = clock();
   cout << "Loading the matrix : " << static_cast<double>(end - begin)/CLOCKS_PER_SEC << " seconds" << endl;
   begin = clock();
   ope.set_lu();
   Array<double> xx(ope.solve(second));

   end = clock();
   cout << "Inverting the matrix : " << static_cast<double>(end - begin)/CLOCKS_PER_SEC << " seconds" << endl;
   int conte(0);
   espace.xx_to_vars_variable_domains(this, xx, conte);
   double* old_var_double(new double[nvar_double]);
   for (int i(0) ; i < nvar_double ; ++i) old_var_double[i] = *var_double[i];
   Tensor** old_fields(new Tensor* [nvar]);
   for (int i(0) ; i < nvar ; i++) old_fields[i] = new Tensor(*var[i]);
   xx_to_vars(xx, conte);
   for (int i(0) ; i < nvar_double ; i++) *var_double[i] = old_var_double[i] - (*var_double[i]);
   for (int i(0) ; i < nvar ; i++) *var[i] = *old_fields[i] - (*var[i]);
   delete [] old_var_double;
   for (int i(0) ; i<nvar ; i++) delete old_fields[i];
   delete [] old_fields;
   return false;
}


bool System_of_eqs::do_newton_with_linesearch(double, double& , int, double ) 
{
  cerr << "Newton-Raphson with lineseach Not implemented yet in the sequential version" << endl ;
  abort() ;
}

void System_of_eqs::check_positive(double)
{  cerr << "Newton-Raphson with lineseach Not implemented yet in the sequential version" << endl ;
  abort() ;
}

void System_of_eqs::check_negative(double )
{
   cerr << "Newton-Raphson with lineseach Not implemented yet in the sequential version" << endl ;
  abort() ;
}

double System_of_eqs::compute_f(Array<double> const&)
{  cerr << "Newton-Raphson with lineseach Not implemented yet in the sequential version" << endl ;
  abort() ;
}

void System_of_eqs::check_size_VS_unknowns(int)
{
    cerr << "Newton-Raphson with lineseach Not implemented yet in the sequential version" << endl ;
  abort() ;
}

void System_of_eqs::check_bsize(int )
{
    cerr << "Newton-Raphson with lineseach Not implemented yet in the sequential version" << endl ;
  abort() ;
}

void System_of_eqs::compute_matloc(Array<double>& , int , int )
{
     cerr << "Newton-Raphson with lineseach Not implemented yet in the sequential version" << endl ;
  abort() ;
}

void System_of_eqs::translate_second_member(Array<double>& , Array<double> const& , int , int , int , int , int )
{
     cerr << "Newton-Raphson with lineseach Not implemented yet in the sequential version" << endl ;
  abort() ;
}

void System_of_eqs::get_global_solution(Array<double>&, Array<double> const& , int, int , int , int , int )
{  cerr << "Newton-Raphson with lineseach Not implemented yet in the sequential version" << endl ;
  abort() ;
}

void System_of_eqs::update_fields(double, vector<double> const& , vector<Tensor> const& , vector<double> const& , vector<Tensor> const& )
{  cerr << "Newton-Raphson with lineseach Not implemented yet in the sequential version" << endl ;
  abort() ;
}

void System_of_eqs::compute_old_and_var(Array<double> const& , vector<double>& , vector<Tensor>& , vector<double>& , vector<Tensor>& )
{  cerr << "Newton-Raphson with lineseach Not implemented yet in the sequential version" << endl ;
  abort() ;
  
}

void System_of_eqs::compute_p(Array<double>&, Array<double> const& , int )
{
    cerr << "Newton-Raphson with lineseach Not implemented yet in the sequential version" << endl ;
  abort() ;
}
}


