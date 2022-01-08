/* ----------------------------------------------------------------------
    SPF - Stochastic Phase Field
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>

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
#include <math.h>    // sqrt

#include "rand.hpp"
#include "flags.hpp"
#include "macheps.hpp"  // determines the local machine's relative 
                        //  rounding error 

namespace SPF_NS
{
   int randomize_neighbor_order(
         std::vector<size_t>& neigh_order,  // 6 elements
         SPF_NS::random& rr, // random generator
         std::uniform_real_distribution<double>& rand_decimal,
         std::vector<double>& rand_decimals1, // reused vector, 6 elements
         std::vector<double>& rand_decimals2, // reused vector, 6 elements
         const int_flags& flags
         );

   int enforce_bounds_generic(
         std::vector<double>& phi_local_flux,
         const std::vector<double>& phi_local,
         const std::vector<double>& phi_local_rates,
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

   int enforce_bounds_pairwise_int_inward(
         std::vector<double>& phi_local_flux,   // integers
         const std::vector<double>& phi_local,  // integers
         const std::vector<double>& phi_local_rates, // doubles
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

   int enforce_bounds_pairwise_int_outward(
         std::vector<double>& phi_local_flux,   // integers
         const std::vector<double>& phi_local,  // integers
         const std::vector<double>& phi_local_rates, // doubles
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

   int enforce_bounds_pairwise_dbl_inward(
         std::vector<double>& phi_local_flux,   // integers
         const std::vector<double>& phi_local,  // integers
         const std::vector<double>& phi_local_rates, // doubles
         SPF_NS::random& rr,
         const size_t& Nvoxel_neighbors,
         const double& phi_lower_limit,
         const double& phi_upper_limit,
         const int& Nx_local,
         const int& Ny,
         const int& Nz,
         const epsilon& eps,
         int_flags& flags
         );

   int enforce_bounds_pairwise_dbl_outward(
         std::vector<double>& phi_local_flux,   // integers
         const std::vector<double>& phi_local,  // integers
         const std::vector<double>& phi_local_rates, // doubles
         SPF_NS::random& rr,
         const size_t& Nvoxel_neighbors,
         const double& phi_lower_limit,
         const double& phi_upper_limit,
         const int& Nx_local,
         const int& Ny,
         const int& Nz,
         const epsilon& eps,
         int_flags& flags
         );

   int enforce_bounds_dbl_outward(
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

   int enforce_bounds_dbl_inward(
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

   //int conserved_gaussian_flux( 
   //int conserved_gaussian_flux_single_distribution_stratonovich(
   //int conserved_gaussian_flux_separate_distributions(
   //int conserved_gaussian_flux_separate_distributions_stratonovich(
   //int conserved_gaussian_flux_separate_distributions_milstein(

   //int conserved_gaussian_flux_single_distribution(
   //   std::vector<double>& local_change,// must be same size as local_field
   //   const std::vector<double>& local_field,
   //   SPF_NS::random& rr,
   //   //const double& rate_scale_factor,
   //   const double& dt,
   //   const size_t& idx,
   //   const size_t& neigh_idx_x_a,
   //   const size_t& neigh_idx_x_b,
   //   const size_t& neigh_idx_y_a,
   //   const size_t& neigh_idx_y_b,
   //   const size_t& neigh_idx_z_a,
   //   const size_t& neigh_idx_z_b,
   //   //const int& Nx_total,
   //   const int& Ny,
   //   const int& Nz
   //   );

   //int conserved_gaussian_flux_single_distribution_milstein(
   //   std::vector<double>& local_change,// must be same size as local_field
   //   const std::vector<double>& local_field,
   //   SPF_NS::random& rr,
   //   //const double& rate_scale_factor,
   //   const double& dt,
   //   const size_t& idx,
   //   const size_t& neigh_idx_x_a,
   //   const size_t& neigh_idx_x_b,
   //   const size_t& neigh_idx_y_a,
   //   const size_t& neigh_idx_y_b,
   //   const size_t& neigh_idx_z_a,
   //   const size_t& neigh_idx_z_b,
   //   //const int& Nx_total,
   //   const int& Ny,
   //   const int& Nz
   //   );

   //int conserved_gaussian_flux_separate_distributions(
   //   std::vector<double>& local_change, // must be same size as local_field
   //   const std::vector<double>& local_field,
   //   SPF_NS::random& rr,
   //   //const double& rate_scale_factor,
   //   const double& dt,
   //   const size_t& idx,
   //   const size_t& neigh_idx_x_a,
   //   const size_t& neigh_idx_x_b,
   //   const size_t& neigh_idx_y_a,
   //   const size_t& neigh_idx_y_b,
   //   const size_t& neigh_idx_z_a,
   //   const size_t& neigh_idx_z_b,
   //   //const int& Nx_total,
   //   const int& Ny,
   //   const int& Nz
   //   );

   //int conserved_gaussian_flux_separate_distributions_gradient_milstein(
   //   std::vector<double>& local_change, // must be same size as local_field
   //   const std::vector<double>& local_field,
   //   SPF_NS::random& rr,
   //   //const double& rate_scale_factor,
   //   const double& dt,
   //   const size_t& idx,
   //   const size_t& neigh_idx_x_a,
   //   const size_t& neigh_idx_x_b,
   //   const size_t& neigh_idx_y_a,
   //   const size_t& neigh_idx_y_b,
   //   const size_t& neigh_idx_z_a,
   //   const size_t& neigh_idx_z_b,
   //   //const int& Nx_total,
   //   const int& Ny,
   //   const int& Nz
   //   );

   //int conserved_jump_flux_single_distribution( 
   //   std::vector<double>& local_change, // must be same size as local_field
   //   const std::vector<double>& local_field,
   //   SPF_NS::random& rr,
   //   //const double& rate_scale_factor,
   //   const double& dt,
   //   const size_t& idx,
   //   const size_t& neigh_idx_x_a,
   //   const size_t& neigh_idx_x_b,
   //   const size_t& neigh_idx_y_a,
   //   const size_t& neigh_idx_y_b,
   //   const size_t& neigh_idx_z_a,
   //   const size_t& neigh_idx_z_b,
   //   //const int& Nx_total,
   //   const int& Ny,
   //   const int& Nz);

   int conserved_jump_flux_separate_distributions( 
      std::vector<double>& pairwise_flux, //( 6 * local_field size )
      SPF_NS::random& rr,
      const std::vector<double>& phi_local,  // integers
      const std::vector<double>& jump_rates,  // 6 elements
      const double& dt,
      const size_t& Nvoxel_neighbors,
      const std::vector<size_t>& neigh_idxs,
      const double& phi_upper_limit,
      const double& phi_lower_limit,
      const size_t& idx);

   int conserved_jump_flux_pairwise_distributions( 
      std::vector<double>& pairwise_flux, // 6* local_field size
      SPF_NS::random& rr,
      const std::vector<double>& phi_local,  // integers
      const std::vector<double>& jump_rates, // 6* local_field size
      const double& dt,
      const size_t& Nvoxel_neighbors,
      const std::vector<size_t>& neigh_idxs,
      const double& phi_upper_limit,
      const double& phi_lower_limit,
      const int& Nx_local, const int& Ny, const int& Nz,
      const size_t& ii, const size_t& jj, const size_t& kk
      //const size_t& idx
      );

   int conserved_jump_flux_pairwise_drift_distributions( 
      std::vector<double>& pairwise_flux, // 6* local_field size
      SPF_NS::random& rr,
      const std::vector<double>& phi_local,  // integers
      const std::vector<double>& jump_rates, // 6* local_field size
      const double& dt,
      const size_t& Nvoxel_neighbors,
      const std::vector<size_t>& neigh_idxs,
      const double& phi_upper_limit,
      const double& phi_lower_limit,
      const int& Nx_local, const int& Ny, const int& Nz,
      const size_t& ii, const size_t& jj, const size_t& kk
      //const size_t& idx
      );

   int conserved_gaussian_flux_separate_distributions( 
      std::vector<double>& pairwise_flux, // 6* local_field size
      SPF_NS::random& rr,
      const std::vector<double>& phi_local,  // integers
      const std::vector<double>& jump_rates_sqrt,  // 6 elements
      const std::vector<double>& jump_rate_sqrt_derivatives,  // 6 elements
      const double& dt,
      const size_t& Nvoxel_neighbors,
      const std::vector<size_t>& neigh_idxs,
      const double& phi_upper_limit,
      const double& phi_lower_limit,
      const size_t& idx);

   int conserved_gaussian_flux_pairwise_distributions( 
      std::vector<double>& pairwise_flux, // 6* local_field size
      SPF_NS::random& rr,
      const std::vector<double>& phi_local,  // integers
      const std::vector<double>& jump_rates_sqrt,  // 6 elements
      const std::vector<double>& jump_rate_sqrt_derivatives,  // 6 elements
      const double& dt,
      const size_t& Nvoxel_neighbors,
      const std::vector<size_t>& neigh_idxs,
      const double& phi_upper_limit,
      const double& phi_lower_limit,
      const size_t& idx);

   //int conserved_gaussian_flux_separate_distributions( 
   //   std::vector<double>& local_change,// must be same size as local_field
   //   const std::vector<double>& local_field,
   //   SPF_NS::random& rr,
   //   const std::vector<double>& jump_rates_sqrt,  // 6 elements
   //   const std::vector<double>& jump_rate_sqrt_derivatives,  // 6 elements
   //   const double& dt,
   //   const size_t& idx,
   //   const std::vector<size_t>& neigh_idxs,  // 6 elements
   //   const std::vector<size_t>& neigh_order,  // 6 elements
   //   const int& Ny,
   //   const int& Nz);

   int conserved_gaussian_flux_separate_distributions_ito( 
      std::vector<double>& local_change,// must be same size as local_field
      const std::vector<double>& local_field,
      SPF_NS::random& rr,
      const std::vector<double>& jump_rates,  // 6 elements
      const std::vector<double>& jump_rate_derivatives,  // 6 elements
      const double& dt,
      const size_t& idx,
      const std::vector<size_t>& neigh_idxs,  // 6 elements
      const std::vector<size_t>& neigh_order,  // 6 elements
      const int& Ny,
      const int& Nz);

   double laplacian1d(
         const double& hh_x,
         const std::vector<double>& local_field,
         const size_t& idx,
         const size_t& neigh_idx_x_a,
         const size_t& neigh_idx_x_b
         );

   double laplacian3d(
         const double& hh_x,
         const double& hh_y,
         const double& hh_z,
         const size_t& Nz,
         const std::vector<double>& local_field,
         const size_t& idx,
         const size_t& neigh_idx_x_a,
         const size_t& neigh_idx_x_b,
         const size_t& neigh_idx_y_a,
         const size_t& neigh_idx_y_b,
         const size_t& neigh_idx_z_a,
         const size_t& neigh_idx_z_b
         //const size_t* const neigh_x_idx,
         //const size_t* const neigh_y_idx,
         //const size_t* const neigh_z_idx
         );

   int laplacian_flux( 
      std::vector<double>& local_change, // must be same size as local_field
      const std::vector<double>& local_field,
      const double& hh_x, const double& hh_y, const double& hh_z,
      const double& dt,
      const double& diffusivity,
      const size_t& idx,
      const size_t& neigh_idx_x_a,
      const size_t& neigh_idx_x_b,
      const size_t& neigh_idx_y_a,
      const size_t& neigh_idx_y_b,
      const size_t& neigh_idx_z_a,
      const size_t& neigh_idx_z_b,
      //const size_t& ii, const size_t& jj, const size_t& kk,
      //const int& Nx_total,
      const int& Ny, const int& Nz
      );

   int flux_test( 
      std::vector<double>& local_change,
      const std::vector<double>& local_field,
      const size_t& ii,
      const size_t& jj,
      const size_t& kk,
      const int& Nx_total,
      const int& Ny,
      const int& Nz
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
} // SPF_NS
#endif
