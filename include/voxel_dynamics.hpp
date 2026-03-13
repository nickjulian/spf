/* ----------------------------------------------------------------------
    SPF - Stochastic Phase Field
    Copyright (C) 2025
    Peng Geng <penggeng@g.ucla.edu>
    Nicholas Huebner Julian <njulian@ucla.edu>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the README file in the top-level SPF directory.
---------------------------------------------------------------------- */
// File: voxel_dynamics.hpp

#ifndef VOXEL_DYNAMICS_HPP
#define VOXEL_DYNAMICS_HPP

#include <cstdlib>   // EXIT_SUCCESS, EXIT_FAILURE
#include <vector>
#include <iostream>  // cout, cin, cerr, endl
#include <algorithm> // sort
#include <math.h>    // sqrt, round // maybe try <cmath>

#include "rand.hpp"
#include "flags.hpp"
#include "macheps.hpp"  // determines the local machine's relative
                        //  rounding error

namespace SPF_NS
{
   int enforce_bounds_int_outward(
         std::vector<double>& phi_local_flux,   // integers
         const std::vector<double>& phi_local,  // integers
         const std::vector<double>& phi_local_rates,  // doubles
         SPF_NS::random& rr,
         //const std::vector<size_t>& neigh_order,
         const size_t& Nvoxel_neighbors,
         const double& phi_lower_limit,
         const double& phi_upper_limit,
         const int& Nx_local,
         const int& Ny,
         const int& Nz,
         const epsilon& eps,
         int_flags& flags
         );

   int enforce_bounds_int_inward(
         std::vector<double>& phi_local_flux,   // integers
         const std::vector<double>& phi_local,  // integers
         const std::vector<double>& phi_local_rates,  // doubles
         SPF_NS::random& rr,
         //const std::vector<size_t>& neigh_order,
         const size_t& Nvoxel_neighbors,
         const double& phi_lower_limit,
         const double& phi_upper_limit,
         const int& Nx_local,
         const int& Ny,
         const int& Nz,
         const epsilon& eps,
         int_flags& flags
         );

   double laplacian1d(
         const double& hh_x,
         const std::vector<double>& local_field,
         const size_t& idx,
         const size_t& neigh_idx_x_a,
         const size_t& neigh_idx_x_b
         );

   int identify_local_neighbors(
      size_t& neigh_idx_x_a,
      size_t& neigh_idx_x_b,
      size_t& neigh_idx_y_a,
      size_t& neigh_idx_y_b,
      size_t& neigh_idx_z_a,
      size_t& neigh_idx_z_b,
      //size_t* const neigh_x_idx,
      //size_t* const neigh_y_idx,
      //size_t* const neigh_z_idx,
      const size_t& ii,
      const size_t& jj,
      const size_t& kk,
      //const int& Nx_total,
      const int& Ny,
      const int& Nz
      );

    double laplacian1dT(
         const double& hh,
         const std::vector<double>& T_field,
         const size_t& idx,
         const size_t& neigh_idx_x_a,
         const size_t& neigh_idx_x_b
         );

} // SPF_NS
#endif
