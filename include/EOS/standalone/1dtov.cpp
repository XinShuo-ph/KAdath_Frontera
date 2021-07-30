/*
 * =====================================================================================
 *
 *       Filename:  1dtov.cpp
 *
 *    Description:  Simple TOV test code to test tov.hh
 *
 *        Version:  1.0
 *        Created:  1/07/2021
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Samuel David Tootle, tootle@itp.uni-frankfurt.de
 *   Organization:  Goethe University Frankfurt
 *   Notes: this is a minimal test for the tov.hh code
 *
 * =====================================================================================
 */


#define PWPOLY_SETUP
#define DEBUG
#include "cold_table.hh"
#include "cold_table_implementation.hh"
#include "setup_cold_table.cc"
#include "tov.hh"
#include <memory>

int main() {
  using namespace Kadath::Margherita;
  setup_Cold_Table("togashi.lorene", 2000);
  auto tov = std::make_unique<MargheritaTOV<Cold_Table>>();
//  tov->adaptive = false;
//  tov->rk45 = false;
  //tov->solve(1.37e-3);
  tov->solve_for_MADM(2.225);
  auto& state = tov->state;
  std::cout << *tov << std::endl;
  return 0;
}
