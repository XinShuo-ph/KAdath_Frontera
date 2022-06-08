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
#include "array.hpp"

#include <chrono>
#include <fstream>

namespace Kadath {
extern "C" {
    void sl_init_ (int*, int*, int*);
    void blacs_gridexit_ (int*);
    void blacs_gridinfo_ (int*, int*, int*, int*, int*);
    void blacs_gridinit_ (int*, char*, int*, int*);
    int numroc_ (int*, int*, int*, int*, int*);
    int indxl2g_ (int*, int*, int*, int*, int*);
    void descinit_ (int*, int*, int*, int*, int* , int*, int*, int*, int*, int*);

    void pdgesv_ (int*, int*, double*, int*, int*, int*, int*, double*, int*, int*, int*, int*);
    void pdtran_ (int*, int*, double*, double*, int*, int*, int*, double*, double*, int*, int*, int*);
    void pdgemv_ (char*, int*, int*, double*, double*, int*, int*, int*, double*, int*, int*, int*, int*, double*, double*, int*, int*, int*, int*);
    void pdgemm_ (char* transa, char* transb, int* m, int* n, int* k, double* alpha, double* a, int* ia, int* ja, int* desc_a, double* b, int* ib, int* jb, int* desc_b, double* beta, double* c, int* ic, int* jc, int* desc_c);

    int blacs_pnum_ (int*, int*, int*);
    void pdgemr2d_ (int *m , int *n , double *a , int *ia , int *ja , int *desca , double *b , int *ib , int *jb , int *descb , int *ictxt );
    void Cpdgemr2d (int, int, double*, int, int, int*, double*, int, int, int*, int);
}

double elapsed_time(std::chrono::time_point<std::chrono::system_clock> const & begin) {
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> diff = end - begin;
  return diff.count();
}

void split_2d (int target, int& low, int& up) {
  up = std::ceil(std::sqrt(target));

  while(target % up)
    ++up;

  low = target / up;
}

bool System_of_eqs::do_newton(double precision, double& error, System_of_eqs::SOLVER solver) {
	// rank and nproc from MPI :
	int rank, nproc;
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Comm_size (MPI_COMM_WORLD, &nproc);

	vars_to_terms();
	Array<double> second (sec_member());
	Array<double> abs_errs(check_equations());
  error = max(abs_errs);

	if (rank==0) {
    if(abs_errs.nbr != eq_list.size() + eq_int_list.size()) {
      std::cout << "There is a mismatch between the number of eqs and the registered eq strings."
                   "Please check the add_eq routines of your space." << std::endl;
  		//abort();
    }

    double max = abs_errs.data[0];
    int i_max = 0;

  	for (int i=0; i<abs_errs.nbr; i++)
  		if (abs_errs.data[i] > max) {
  			max = abs_errs.data[i];
  			i_max = i;
  		}

		std::cout << "Error init = " << error << ",";

		if(i_max < neq_int && neq_int > 0) {
      auto [ eq, dom, boundary ] = eq_int_list[i_max];
      std::cout << " Eq: " << eq << " Dom: " << dom << " Bounday: " << boundary << std::endl;
		}
		else {
      i_max -= neq_int;
      auto [ eq, dom, boundary ] = eq_list[i_max];
      std::cout << " Eq: " << eq << " Dom: " << dom << " Bounday: " << boundary << std::endl;
		}
	}

	if (error < precision) {
		return true;
	}
	else {

	int m = second.get_size(0);
  int n = nbr_unknowns;

	if (solver == NEWTON_RAPHSON && n != m) {
		cerr << "Number of unknowns is different from the number of equations" << endl;
		cerr << "nbr equations = " << m << endl;
		cerr << "nbr unknowns  = " << n << endl;
		abort();
	}

	if (rank==0) {
    double size = 1.0 * n * m;
		std::cout << "DOF / rank = " << n << " x " << m << " / " << nproc <<  " = " << size / nproc << endl;
	}

	// Computation in a 1D distributed matrice
	int zero_i = 0;
	int one_i = 1;

	double zero_d = 0.;
	double one_d = 1.;

  // create an input process grid to generate J
  // processor grid context
	int ictxt_in;

  // processor grid extents
  // we want to fill full cols locally, so generate a 1d process grid
  int nprow_in = 1;
  int npcol_in = nproc;

  // init context
  sl_init_ (&ictxt_in, &nprow_in, &npcol_in);

	// Get processor row and col
  int prow_in, pcol_in;
  blacs_gridinfo_ (&ictxt_in, &nprow_in, &npcol_in, &prow_in, &pcol_in);

  int bsize = 256;

  while(nproc * bsize > m)
    bsize = bsize / 2;

  if(rank == 0) {
    if(bsize < 64) {
      std::cout << "Block size is " << bsize << ", consider lowering the number of cores." << std::endl;
    }
  }

  // local row and col extent
	int nrow_in = numroc_ (&m, &bsize, &prow_in, &zero_i, &nprow_in);
	int ncol_in = numroc_ (&n, &bsize, &pcol_in, &zero_i, &npcol_in);

  // collect global columns assigned to this proc, covers full rows because of the linear process grid above
  std::vector<int> globcol_in;
  for(int j = 1; j <= ncol_in; ++j)
    globcol_in.push_back(indxl2g_(&j,&bsize,&pcol_in,&zero_i,&npcol_in) - 1);

  auto begin = std::chrono::system_clock::now();

  // actually compute the local cols of J
  std::vector<Array<double>> col_J;
  for(int j = 0; j < ncol_in; ++j) {
    col_J.push_back(do_col_J(globcol_in[j]));
  }

  // distributed storage for local columns of J
  // we need col major order!
	Array<double> J_in (ncol_in, nrow_in);

  // every proc spans all rows in this case because of the 1d processor grid
  for(int j = 0; j < ncol_in; ++j) {
    for(int i = 0; i < nrow_in; ++i) {
      J_in.set(j,i) = col_J[j](i);
    }
  }

  // catching errors
  int info;

  // descriptor of J
  Array<int> J_desc_in(9);
  descinit_ (J_desc_in.set_data(), &m, &n, &bsize, &bsize, &zero_i, &zero_i, &ictxt_in, &nrow_in, &info);

  // create the compute processor grid
	int ictxt;

  // processor grid extents
  int nprow;
  int npcol;
  split_2d(nproc, nprow, npcol);

  // init context
  sl_init_ (&ictxt, &nprow, &npcol);

	// Get processor row and col
  int prow, pcol;
  blacs_gridinfo_ (&ictxt, &nprow, &npcol, &prow, &pcol);

  MPI_Barrier(MPI_COMM_WORLD);

  // local row and col extent
	int nrow = numroc_ (&m, &bsize, &prow, &zero_i, &nprow);
	int ncol = numroc_ (&n, &bsize, &pcol, &zero_i, &npcol);

  // local copy of J for the computation
  // we need col major order!
	Array<double> J (ncol, nrow);

  // descriptor of J
  Array<int> J_desc(9);
  descinit_ (J_desc.set_data(), &m, &n, &bsize, &bsize, &zero_i, &zero_i, &ictxt, &nrow, &info);

  // now distribute J_in over the compute grid to J
  pdgemr2d_ (&m, &n,
             J_in.set_data(), &one_i , &one_i, J_desc_in.set_data(),
             J.set_data(), &one_i , &one_i, J_desc.set_data(),
             &ictxt);

  J_in.delete_data();

  // local rows of B
	int nrow_B = numroc_ (&m, &bsize, &prow, &zero_i, &nprow);
	Array<double> B (nrow_B);

  for(int i = 1; i <= nrow_B; ++i) {
    int glob_row = indxl2g_(&i,&bsize,&prow,&zero_i,&nprow) - 1;
    B.set(i-1) = second(glob_row);
  }

  // descriptor of B
  Array<int> B_desc(9);
  descinit_ (B_desc.set_data(), &m, &one_i, &bsize, &one_i, &zero_i, &zero_i, &ictxt, &nrow_B, &info);

  double time = elapsed_time(begin);

  double maxtime;
  double mintime;
  double avgtime;

  MPI_Reduce(&time, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&time, &mintime, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&time, &avgtime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    avgtime /= nproc;
    std::cout << "Load and redistrib Jacobian: " << avgtime << " / " << mintime << " / " << maxtime << " seconds (avg/min/max)" << endl;
  }

  begin = std::chrono::system_clock::now();

  // local copy of the solution
    Array<double>& X_sol = B;

    if(solver == GAUSS_NEWTON) {
    // local rows of J^T B
  	int nrow_X = numroc_ (&n, &bsize, &prow, &zero_i, &nprow);
  	Array<double> JTB (nrow_X);

    // descriptor of J^T B
    Array<int> JTB_desc(9);
    descinit_ (JTB_desc.set_data(), &n, &one_i, &bsize, &one_i, &zero_i, &zero_i, &ictxt, &nrow_X, &info);

    // local copy of J^T J
  	int nrow_sym = numroc_ (&n, &bsize, &prow, &zero_i, &nprow);
  	int ncol_sym = numroc_ (&n, &bsize, &pcol, &zero_i, &npcol);
  	Array<double> JTJ (ncol_sym, nrow_sym);

    // descriptor of J^T J
    Array<int> JTJ_desc(9);
    descinit_ (JTJ_desc.set_data(), &n, &n, &bsize, &bsize, &zero_i, &zero_i, &ictxt, &nrow_sym, &info);

    // FIXME precomputing J^T and storing it could be a bit faster, since we need it after this again
    // compute J^T B
    char trans = 'T';
    pdgemv_(&trans, &m, &n,
            &one_d,
            J.set_data(), &one_i, &one_i, J_desc.set_data(),
            B.set_data(), &one_i, &one_i, B_desc.set_data(), &one_i,
            &zero_d,
            JTB.set_data(), &one_i, &one_i, JTB_desc.set_data(), &one_i);

    // compute J^T J
    char no_trans = 'N';
    pdgemm_ (&trans, &no_trans,
             &n, &n, &m, &one_d,
             J.set_data(), &one_i, &one_i, J_desc.set_data(),
             J.set_data(), &one_i, &one_i, J_desc.set_data(),
             &zero_d,
             JTJ.set_data(), &one_i, &one_i, JTJ_desc.set_data());

  	// solve the linear system
    Array<int> ipiv (n);

    pdgesv_ (&n, &one_i,
             JTJ.set_data(), &one_i, &one_i, JTJ_desc.set_data(),
             ipiv.set_data(),
             JTB.set_data(), &one_i, &one_i, JTB_desc.set_data(), &info);

    // local copy of the solution
    X_sol = JTB;
  }

  if(solver == NEWTON_RAPHSON) {
  	// solve the linear system
    Array<int> ipiv (n);

    pdgesv_ (&n, &one_i,
             J.set_data(), &one_i, &one_i, J_desc.set_data(),
             ipiv.set_data(),
             B.set_data(), &one_i, &one_i, B_desc.set_data(), &info);
  }

  if(info != 0) {
    if (rank == 0) {
      std::cout << "Jacobian is singular! Stopping." << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    abort();
  }

  if (rank == 0) {
    std::cout << "Inverting the matrix: " << elapsed_time(begin) << " seconds" << endl;
  }

  begin = std::chrono::system_clock::now();

  // map number of elements to rank
  std::vector<int> n_elements(nproc);

  for(int i = 0; i < nprow; ++i) {
    for(int j = 0; j < npcol; ++j) {
      // rank as function of row and col proc number (default ordering)
      int r = i * npcol + j;
      // vectors are only distributed over first proc column
    	int n_elem = j == 0 ? numroc_(&n, &bsize, &i, &zero_i, &nprow) : 0;

      n_elements[r] = n_elem;
    }
  }

  // no displacement
  std::vector<int> disp(nproc, 0);

  int offset = 0;
  for(int r = 0; r < nproc; ++r) {
    disp[r] = offset;
    offset += n_elements[r];
  }

  // recieve buffer
  std::vector<double> recv(n);

  // distribute full solution to all processors
  MPI_Allgatherv(X_sol.set_data(), n_elements[rank], MPI_DOUBLE,
                 recv.data(), n_elements.data(), disp.data(), MPI_DOUBLE,
                 MPI_COMM_WORLD);

  // rearrange data to correct ordering
  Array<double>  X(n);

  for(int i = 0; i < nprow; ++i) {
    for(int j = 0; j < npcol; ++j) {
      // rank as function of row and col proc number (default ordering)
      int r = i * npcol + j;

      if(n_elements[r] > 0) {
        // data is ordered by ranks, so get correct offset
        int pos = 0;
        for(int k = 0; k < r; ++k)
          pos += n_elements[k];

        for(int k = 0; k < n_elements[r]; ++k) {
          int elem = k + 1;
          int glob_row = indxl2g_(&elem,&bsize,&i,&zero_i,&nprow) - 1;

          X.set(glob_row) = recv[pos + k];
        }
      }
    }
  }

	blacs_gridexit_ (&ictxt_in);
	blacs_gridexit_ (&ictxt);

	int conte = 0;
	espace.xx_to_vars_variable_domains (this, X, conte);

	double* old_var_double = (nvar_double==0) ? 0x0 : new double[nvar_double];

	for (int i=0 ; i<nvar_double ; i++)
		old_var_double[i] = *var_double[i];

	Tensor** old_fields = new Tensor* [nvar];
	for (int i=0 ; i<nvar ; i++)
		old_fields[i] = new Tensor(*var[i]);

	xx_to_vars(X, conte);

	for (int i=0 ; i<nvar ; i++)
		*var[i] = *old_fields[i] - *var[i];

	for (int i=0 ; i<nvar_double ; i++)
		*var_double[i] = old_var_double[i] - *var_double[i];

	if (old_var_double!=0x0)
	    delete [] old_var_double;

	for (int i=0 ; i<nvar ; i++)
		delete old_fields[i];

	delete [] old_fields;

  if (rank == 0) {
    std::cout << "Distributing solution: " << elapsed_time(begin) << " seconds" << endl;
  }

	return false;
  }
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

bool System_of_eqs::do_newton_seq(double precision, double &error)
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

