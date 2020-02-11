//
// Created by sauliac on 09/01/2020.
//
#include <config.h>

#ifdef PAR_VERSION
/*
    Copyright 2017 Philippe Grandclement & Gregoire Martinon

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

#include "mpi.h"
#include "system_of_eqs.hpp"
#include "matrice.hpp"
#include "scalar.hpp"
#include "metric.hpp"
#include "array_math.cpp"

#include <ctime>
namespace Kadath {
extern "C" {
    void sl_init_ (int*, int*, int*);
    void blacs_gridexit_ (int*);
    void blacs_gridinfo_ (int*, int*, int*, int*, int*);
    int numroc_ (int*, int*, int*, int*, int*);
    void descinit_ (int*, int*, int*, int*, int* , int*, int*, int*, int*, int*);
    void pdgesv_ (int*, int*, double*, int*, int*, int*, int*, double*, int*, int*, int*, int*);
    int blacs_pnum_ (int*, int*, int*);
    void Cpdgemr2d (int, int, double*, int, int, int*, double*, int, int, int*, int);
}

void split_two_d (int target, int& low, int& up) {
      int start = int(sqrt(target));
      low = start;
      up = start;

      bool endloop = (low*up==target) ? true : false;

      bool first = true;
      while (!endloop) {
	    if (!first)
	      low--;
	    else
	      first= false;

	    up = int(target/low);

	if (low*up==target)
	    endloop = true;
      }

}


bool System_of_eqs::do_newton(double precision, double& error) {

   clock_t begin, end;
	bool res;
	int bsize = 64;

	// rank and nproc from MPI :
	int rank, nproc;
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Comm_size (MPI_COMM_WORLD, &nproc);

	vars_to_terms();
	Array<double> second (sec_member());
	error = max(fabs(second));
	if (rank==0)
		cout << "Error init = " << error << endl;
	if (error<precision) {
		res = true;
	}
	else {

	int nn = second.get_size(0);
	if (nn!=nbr_unknowns) {
		cerr << "Number of unknowns is different from the number of equations" << endl;
		cerr << "nbr equations = " << nn << endl;
		cerr << "nbr unknowns  = " << nbr_unknowns << endl;
		abort();
	}

	if (rank==0)
		cout << "Size of the system " << nn << endl;

	// Computation in a 1D distributed matrice
	int zero = 0;
	int ictxt_in;
        int nprow_in = 1;
        int npcol_in = nproc;
        sl_init_ (&ictxt_in, &nprow_in, &npcol_in);

	// Get my row and mycol
        int myrow_in, mycol_in;
        blacs_gridinfo_ (&ictxt_in, &nprow_in, &npcol_in, &myrow_in, &mycol_in);

	while (bsize*nproc>nn) {
		bsize = div(bsize,2).quot;
	}
	if (bsize<1) {
		cerr << "Too many processrs in do_newton" << endl;
		abort();
	}

	int nrowloc_in = numroc_ (&nn, &bsize, &myrow_in, &zero, &nprow_in);
	int ncolloc_in = numroc_ (&nn, &bsize, &mycol_in, &zero, &npcol_in);

	Array<double> matloc_in (ncolloc_in, nrowloc_in);
	int start = bsize*rank;
	bool done = false;
	int current = 0;

	hash_key chrono_key = this->start_chrono("MPI parallel do_newton | problem size = ",
	        nn," | matrix computation ");

	while (!done) {
		for (int i=0 ; i<bsize ; i++)
			if (start+i<nn) {
			Array<double> column (do_col_J(start+i));
			for (int j=0 ; j<nn ; j++)
				matloc_in.set(current,j) = column(j);
			current++;
		}
		start += nproc*bsize;
		if (start>=nn)
			done = true;
	}

	 // Descriptor of the matrix :
        Array<int> descamat_in(9);
        int info;
        descinit_ (descamat_in.set_data(), &nn, &nn, &bsize, &bsize, &zero, &zero, &ictxt_in, &nrowloc_in, &info);

	// Wait for everybody
	MPI_Barrier(MPI_COMM_WORLD);

    duration const t_load_matrix {this->stop_chrono(chrono_key)};

    if (rank == 0) cout << "Loading the matrix : " << to_seconds(t_load_matrix) << " seconds" << endl;

    chrono_key = this->start_chrono("MPI parallel do_newton | problem size = ",
                                        nn," | matrix translation ");
	// Now translate to a 2D cyclic decomposition
	int npcol, nprow;
	split_two_d (nproc, npcol, nprow);
	int ictxt;
	sl_init_ (&ictxt, &nprow, &npcol);

        int myrow, mycol;
        blacs_gridinfo_ (&ictxt, &nprow, &npcol, &myrow, &mycol);

        int nrowloc = numroc_ (&nn, &bsize, &myrow, &zero, &nprow);
        int ncolloc = numroc_ (&nn, &bsize, &mycol, &zero, &npcol);

	Array<double> matloc (ncolloc, nrowloc);
        Array<int> descamat(9);
        descinit_ (descamat.set_data(), &nn, &nn, &bsize, &bsize, &zero, &zero, &ictxt, &nrowloc, &info);

        Cpdgemr2d (nn, nn, matloc_in.set_data(), 1, 1, descamat_in.set_data(), matloc.set_data(), 1, 1, descamat.set_data(), ictxt);

	matloc_in.delete_data();

	// Translate the second member :
	Array<double> secloc (nrowloc);
        for (int row=0 ; row<nn ; row++) {
                int pi = div(row/bsize, nprow).rem;
                int li = int(row/(nprow*bsize));
                int xi = div(row, bsize).rem;
                if ((pi==myrow) && (mycol==0))
                        secloc.set(li*bsize+xi) = second(row);
        }

        // Descriptor of the second member
        Array<int> descsec(9);
        int one = 1;
        descinit_ (descsec.set_data(), &nn, &one, &bsize, &bsize, &zero, &zero, &ictxt, &nrowloc, &info);
        duration const t_trans_matrix {this->stop_chrono(chrono_key)};

        if (rank == 0) cout << "Translating the matrix : " << to_seconds(t_trans_matrix) << " seconds" << endl;

	// Inversion
        Array<int> ipiv (nn);
        chrono_key = this->start_chrono("MPI parallel do_newton | problem size = ",
                                          nn," | matrix inversion ");
        pdgesv_ (&nn, &one, matloc.set_data(), &one, &one, descamat.set_data(), ipiv.set_data(), secloc.set_data(), &one, &one, descsec.set_data(), &info);
        duration const t_inv_matrix {this->stop_chrono(chrono_key)};
        if (rank == 0) cout << "Inverting the matrix : " << to_seconds(t_inv_matrix) << " seconds" << endl;

        chrono_key = this->start_chrono("MPI parallel do newton | problem size = ", nn, " | update ");
	// Get the global solution
        Array<double> auxi(nn);
        auxi = 0.;
        for (int row=0 ; row<nn ; row++) {
                int pi = div(row/bsize, nprow).rem;
                int li = int(row/(nprow*bsize));
                int xi = div(row, bsize).rem;
                if ((pi==myrow) && (mycol==0))
                        auxi.set(row) = secloc(li*bsize+xi);
        }

        Array<double> xx (nn);
        MPI_Allreduce (auxi.set_data(), xx.set_data(), nn, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	blacs_gridexit_ (&ictxt);
	blacs_gridexit_ (&ictxt_in);

	int conte = 0;
	espace.xx_to_vars_variable_domains (this, xx, conte);

	double* old_var_double = (nvar_double==0) ? 0x0 :
					new double[nvar_double];
	for (int i=0 ; i<nvar_double ; i++)
		old_var_double[i] = *var_double[i];

	Tensor** old_fields = new Tensor* [nvar];
	for (int i=0 ; i<nvar ; i++)
		old_fields[i] = new Tensor(*var[i]);

	xx_to_vars(xx, conte);

	for (int i=0 ; i<nvar ; i++)
		*var[i] = *old_fields[i] - *var[i];

	for (int i=0 ; i<nvar_double ; i++)
		*var_double[i] = old_var_double[i] - *var_double[i];

	if (old_var_double!=0x0)
	    delete [] old_var_double;
	for (int i=0 ; i<nvar ; i++)
		delete old_fields[i];

	delete [] old_fields;
	duration const t_newton_update {this->stop_chrono(chrono_key)};
    if (rank == 0) cout << "Newton update : " << to_seconds(t_newton_update) << " seconds" << endl << endl;
	res = false;
      }
	return res;
}


bool System_of_eqs::do_newton_with_linesearch(double precision, double& error, int ntrymax, double stepmax)
{
   // Numerical recipes 2007, section 9.7.1
   bool res(false);
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   Array<double> second(sec_member());  // rhs of J.x = F
   error = max(fabs(second));
   if (rank == 0) cout << "Error init = " << error << endl;
   if (error < precision) res = true;
   else
   {
      int nn(second.get_size(0));
      check_size_VS_unknowns(nn);
      if (rank == 0) cout << "Size of the system " << nn << endl;

      Array<double> p(nn);             // solution of J.x = F
      double fold, f, slope, f2(0.0);
      double lambda(1.0), lambdatmp(0.0), lambda2(0.0);
      const double lambdamin(1.0e-10);
      const double alpha(1.0e-4);

      fold = compute_f(second);    // old value of f = 0.5*F.F
      slope = -2.0*fold;
      check_negative(slope);
      compute_p(p, second, nn);    // all the parallel inversion of the jacobian is crammed into this function
      double pmax(max(p));
      if (rank == 0) cout << "max(newton step) = " << pmax << endl;

      if (pmax > stepmax) p *= stepmax/pmax;      // rescale in case p is too large
      vector<double> old_var_double, p_var_double;
      vector<Tensor> old_var_fields, p_var_fields;
      compute_old_and_var(p, old_var_double, old_var_fields, p_var_double, p_var_fields);
      for (int itry(1) ; itry <= ntrymax ; ++itry)
      {
         update_fields(lambda, old_var_double, old_var_fields, p_var_double, p_var_fields);    // var now updated
         second = sec_member();
         f = compute_f(second);
         if (f <= 5.0e-13) break; // in double precision, linesearch in spoiled by numerical error as soon as the error is close to 1.0e-7 (f \propto error^2)
         if (f <= fold + alpha*lambda*slope or lambda <= lambdamin) break;
         if (rank == 0) cout << "line search " << itry << " ====> ";
         if (lambda == 1.0) lambdatmp = -0.5*slope/(f - fold - slope);
         else
         {
            double rhs1(f - fold - lambda*slope);
            double rhs2(f2 - fold - lambda2*slope);
            double a((rhs1/pow(lambda,2) - rhs2/pow(lambda2,2))/(lambda-lambda2));
            double b((-lambda2*rhs1/pow(lambda,2) + lambda*rhs2/pow(lambda2,2))/(lambda-lambda2));
            if (a == 0.0) lambdatmp = -0.5*slope/b;
            double delta(pow(b,2) - 3.0*a*slope);
            check_positive(delta);
            lambdatmp = (-b + sqrt(delta))/(3.0*a);
            if (lambdatmp > 0.5*lambda) lambdatmp = 0.5*lambda;
         }
         lambda2 = lambda;
         f2 = f;
         lambda = max(lambdatmp, 0.1*lambda);
         if (rank == 0) cout << "lambda = " << lambda << endl;
         if (itry == ntrymax) update_fields(lambda, old_var_double, old_var_fields, p_var_double, p_var_fields);
      }
   }
   return res;
}

void System_of_eqs::check_positive(double delta)
{
   if (delta < 0.0)
   {
      cerr << "roundoff error in do_newton_with_linesearch" << endl;
      abort();
   }
}

void System_of_eqs::check_negative(double delta)
{
   if (delta >= 0.0)
   {
      cerr << "roundoff error in do_newton_with_linesearch" << endl;
      abort();
   }
}

double System_of_eqs::compute_f(Array<double> const& second)
{
   double f(0.0);
   for (int i(0) ; i < second.get_nbr() ; ++i) f += pow(second.get_data()[i],2) ;
   return f;
}

void System_of_eqs::check_size_VS_unknowns(int n)
{
   if (n != nbr_unknowns)
   {
      cerr << "Number of unknowns is different from the number of equations" << endl;
      cerr << "nbr equations = " << n << endl;
      cerr << "nbr unknowns  = " << nbr_unknowns << endl;
      abort();
   }
}

void System_of_eqs::check_bsize(int bsize)
{
   if (bsize < 1)
   {
      cerr << "Too many processrs in do_newton" << endl;
      abort();
   }
}

void System_of_eqs::compute_matloc(Array<double>& matloc_in, int nn, int bsize)
{
   int rank, nproc;
   MPI_Comm_rank (MPI_COMM_WORLD, &rank);
   MPI_Comm_size (MPI_COMM_WORLD, &nproc);

   int start(bsize*rank);
   bool done(false);
   int current(0);
   while (!done)
   {
      for (int i(0) ; i < bsize ; ++i)
         if (start+i < nn)
         {
            Array<double> column(do_col_J(start+i));
            for (int j(0) ; j < nn ; ++j) matloc_in.set(current,j) = column(j);
            current++;
         }
      start += nproc*bsize;
      if (start >= nn) done = true;
   }
}

void System_of_eqs::translate_second_member(Array<double>& secloc, Array<double> const& second, int nn, int bsize, int nprow, int myrow, int mycol)
{
   for (int row(0) ; row < nn ; ++row)
   {
      int pi(div(row/bsize, nprow).rem);
      int li(int(row/(nprow*bsize)));
      int xi(div(row, bsize).rem);
      if ((pi == myrow) and (mycol == 0)) secloc.set(li*bsize+xi) = second(row);
   }
}

void System_of_eqs::get_global_solution(Array<double>& auxi, Array<double> const& secloc, int nn, int bsize, int nprow, int myrow, int mycol)
{
   auxi = 0.0;
   for (int row(0) ; row < nn ; ++row)
   {
      int pi(div(row/bsize, nprow).rem);
      int li(int(row/(nprow*bsize)));
      int xi(div(row, bsize).rem);
      if ((pi == myrow) and (mycol == 0)) auxi.set(row) = secloc(li*bsize+xi);
   }
}

void System_of_eqs::update_fields(double lambda, vector<double> const& old_var_double, vector<Tensor> const& old_var_fields, vector<double> const& p_var_double, vector<Tensor> const& p_var_fields)
{
   int rank;
   MPI_Comm_rank (MPI_COMM_WORLD, &rank);
   if (nvar_double > 0) for (int i(0) ; i < nvar_double ; ++i) *var_double[i] = old_var_double[i] - lambda*p_var_double[i];
   for (int i(0) ; i < nvar ; ++i) *var[i] = old_var_fields[i] - lambda*p_var_fields[i];
}

void System_of_eqs::compute_old_and_var(Array<double> const& xx, vector<double>& old_var_double, vector<Tensor>& old_var_fields, vector<double>& p_var_double, vector<Tensor>& p_var_fields)
{
   int rank;
   MPI_Comm_rank (MPI_COMM_WORLD, &rank);
   int conte(0);
   espace.xx_to_vars_variable_domains(this, xx, conte);
   if (nvar_double > 0) for (int i(0) ; i < nvar_double ; ++i) old_var_double.push_back(*var_double[i]);
   for (int i(0) ; i < nvar ; ++i) old_var_fields.push_back(Tensor(*var[i]));
   xx_to_vars(xx, conte);
   if (nvar_double > 0) for (int i(0) ; i < nvar_double ; ++i) p_var_double.push_back(*var_double[i]);
   for (int i(0) ; i < nvar ; ++i) p_var_fields.push_back(Tensor(*var[i]));
}

void System_of_eqs::compute_p(Array<double>& xx, Array<double> const& second, int nn)
{
   int bsize(64);
   // rank and nproc from MPI :
   int rank, nproc;
   MPI_Comm_rank (MPI_COMM_WORLD, &rank);
   MPI_Comm_size (MPI_COMM_WORLD, &nproc);

   // Computation in a 1D distributed matrice
   int zero(0);
   int ictxt_in;
   int nprow_in(1);
   int npcol_in(nproc);
   sl_init_(&ictxt_in, &nprow_in, &npcol_in);

   // Get my row and mycol
   int myrow_in, mycol_in;
   blacs_gridinfo_ (&ictxt_in, &nprow_in, &npcol_in, &myrow_in, &mycol_in);
   while (bsize*nproc > nn) bsize = div(bsize,2).quot;
   check_bsize(bsize);
   int nrowloc_in(numroc_(&nn, &bsize, &myrow_in, &zero, &nprow_in));
   int ncolloc_in(numroc_(&nn, &bsize, &mycol_in, &zero, &npcol_in));
   Array<double> matloc_in(ncolloc_in, nrowloc_in);
   compute_matloc(matloc_in, nn, bsize);

    // Descriptor of the matrix :
   Array<int> descamat_in(9);
   int info;
   descinit_(descamat_in.set_data(), &nn, &nn, &bsize, &bsize, &zero, &zero, &ictxt_in, &nrowloc_in, &info);

   // Wait for everybody
   MPI_Barrier(MPI_COMM_WORLD);

   // Now translate to a 2D cyclic decomposition
   int npcol, nprow;
   int ictxt;
   int myrow, mycol;
   split_two_d(nproc, npcol, nprow);
   sl_init_(&ictxt, &nprow, &npcol);
   blacs_gridinfo_(&ictxt, &nprow, &npcol, &myrow, &mycol);
   int nrowloc(numroc_(&nn, &bsize, &myrow, &zero, &nprow));
   int ncolloc(numroc_(&nn, &bsize, &mycol, &zero, &npcol));
   Array<double> matloc(ncolloc, nrowloc);
   Array<int> descamat(9);
   descinit_(descamat.set_data(), &nn, &nn, &bsize, &bsize, &zero, &zero, &ictxt, &nrowloc, &info);
   Cpdgemr2d(nn, nn, matloc_in.set_data(), 1, 1, descamat_in.set_data(), matloc.set_data(), 1, 1, descamat.set_data(), ictxt);
   matloc_in.delete_data();

   // Translate the second member :
   Array<double> secloc(nrowloc);
   translate_second_member(secloc, second, nn, bsize, nprow, myrow, mycol);

   // Descriptor of the second member
   Array<int> descsec(9);
   int one(1);
   descinit_(descsec.set_data(), &nn, &one, &bsize, &bsize, &zero, &zero, &ictxt, &nrowloc, &info);

   // Inversion
   Array<int> ipiv(nn);
   pdgesv_(&nn, &one, matloc.set_data(), &one, &one, descamat.set_data(), ipiv.set_data(), secloc.set_data(), &one, &one, descsec.set_data(), &info);

   // Get the global solution
   Array<double> auxi(nn);
   get_global_solution(auxi, secloc, nn, bsize, nprow, myrow, mycol);
   MPI_Allreduce(auxi.set_data(), xx.set_data(), nn, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   blacs_gridexit_(&ictxt);
   blacs_gridexit_(&ictxt_in);
}
}


#else //ifdef PAR_VERSION

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

        hash_key chrono_key = this->start_chrono("do_newton | problem size = ",nn," | matrix computation ");
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
        duration const
            t_load_matrix {this->stop_chrono(chrono_key)};
        cout << "Loading the matrix : " << to_seconds(t_load_matrix) << " seconds" << endl;

        chrono_key = this->start_chrono("do_newton | problem size = ",nn," | matrix inversion ");
        ope.set_lu();
        Array<double> xx(ope.solve(second));
        duration const
                t_inv_matrix {this->stop_chrono(chrono_key)};
        cout << "Inverting the matrix : " << to_seconds(t_inv_matrix) << " seconds" << endl;
        int conte(0);

        chrono_key = this->start_chrono("do_newton | problem size = ",nn," | Newton update ");
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
        duration const t_newton_update
            {this->stop_chrono(chrono_key)};
        cout << "Newton update : " << to_seconds(t_newton_update) << " seconds" << endl;
        return false;
    }


    bool System_of_eqs::do_newton_with_linesearch(double, double& , int, double )
    {
        cerr << "Newton-Raphson with line-search Not implemented yet in the sequential version" << endl ;
        abort() ;
    }

    void System_of_eqs::check_positive(double)
    {  cerr << "Newton-Raphson with line-search Not implemented yet in the sequential version" << endl ;
        abort() ;
    }

    void System_of_eqs::check_negative(double )
    {
        cerr << "Newton-Raphson with line-search Not implemented yet in the sequential version" << endl ;
        abort() ;
    }

    double System_of_eqs::compute_f(Array<double> const&)
    {  cerr << "Newton-Raphson with line-search Not implemented yet in the sequential version" << endl ;
        abort() ;
    }

    void System_of_eqs::check_size_VS_unknowns(int)
    {
        cerr << "Newton-Raphson with line-search Not implemented yet in the sequential version" << endl ;
        abort() ;
    }

    void System_of_eqs::check_bsize(int )
    {
        cerr << "Newton-Raphson with line-search Not implemented yet in the sequential version" << endl ;
        abort() ;
    }

    void System_of_eqs::compute_matloc(Array<double>& , int , int )
    {
        cerr << "Newton-Raphson with line-search Not implemented yet in the sequential version" << endl ;
        abort() ;
    }

    void System_of_eqs::translate_second_member(Array<double>& , Array<double> const& , int , int , int , int , int )
    {
        cerr << "Newton-Raphson with line-search Not implemented yet in the sequential version" << endl ;
        abort() ;
    }

    void System_of_eqs::get_global_solution(Array<double>&, Array<double> const& , int, int , int , int , int )
    {  cerr << "Newton-Raphson with line-search Not implemented yet in the sequential version" << endl ;
        abort() ;
    }

    void System_of_eqs::update_fields(double, vector<double> const& , vector<Tensor> const& , vector<double> const& , vector<Tensor> const& )
    {  cerr << "Newton-Raphson with line-search Not implemented yet in the sequential version" << endl ;
        abort() ;
    }

    void System_of_eqs::compute_old_and_var(Array<double> const& , vector<double>& , vector<Tensor>& , vector<double>& , vector<Tensor>& )
    {  cerr << "Newton-Raphson with line-search Not implemented yet in the sequential version" << endl ;
        abort() ;

    }

    void System_of_eqs::compute_p(Array<double>&, Array<double> const& , int )
    {
        cerr << "Newton-Raphson with line-search Not implemented yet in the sequential version" << endl ;
        abort() ;
    }
}



#endif //ifdef PAR_VERSION
