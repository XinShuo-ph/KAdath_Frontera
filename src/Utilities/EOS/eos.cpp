//
// This file is part of Margherita, the light-weight EOS framework
//
//  Copyright (C) 2017, Elias Roland Most
//                      <emost@th.physik.uni-frankfurt.de>
//  Copyright (C) 2020, Ludwig Jens Papenfort
//                      <@th.physik.uni-frankfurt.de>
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <EOS/EOS.hh>

using namespace Kadath::Margherita;

// storage for static members
std::array<double, Cold_PWPoly::max_num_pieces> Cold_PWPoly::k_tab;
std::array<double, Cold_PWPoly::max_num_pieces> Cold_PWPoly::gamma_tab;
std::array<double, Cold_PWPoly::max_num_pieces> Cold_PWPoly::rho_tab;
std::array<double, Cold_PWPoly::max_num_pieces> Cold_PWPoly::eps_tab;
std::array<double, Cold_PWPoly::max_num_pieces> Cold_PWPoly::P_tab;
std::array<double, Cold_PWPoly::max_num_pieces> Cold_PWPoly::h_tab;
int Cold_PWPoly::num_pieces = 1;

double Cold_PWPoly::rhomin;
double Cold_PWPoly::rhomax;
