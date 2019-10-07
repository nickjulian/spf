/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: voxel_dynamics.hpp
// Purpose:

#ifndef VOXEL_DYNAMICS_HPP
#define VOXEL_DYNAMICS_HPP

#include <cstdlib>   // EXIT_SUCCESS, EXIT_FAILURE
#include <vector>
#include <iostream> // cout, cin, cerr, endl

#include "rand.hpp"

namespace SPF_NS
{
   int conserved_gaussian_flux( 
      std::vector<double>& local_change, // must be same size as local_field
      const std::vector<double>& local_field,
      SPF_NS::random& rr,
      const double& rate_scale_factor,
      const double& dt,
      const size_t& idx,
      const size_t& neigh_idx_x_a,
      const size_t& neigh_idx_x_b,
      const size_t& neigh_idx_y_a,
      const size_t& neigh_idx_y_b,
      const size_t& neigh_idx_z_a,
      const size_t& neigh_idx_z_b,
      //const int& Nx_total,
      const int& Ny,
      const int& Nz
      );

   int conserved_jump_flux( 
      std::vector<double>& local_change, // must be same size as local_field
      const std::vector<double>& local_field,
      SPF_NS::random& rr,
      const double& rate_scale_factor,
      const double& dt,
      const size_t& idx,
      const size_t& neigh_idx_x_a,
      const size_t& neigh_idx_x_b,
      const size_t& neigh_idx_y_a,
      const size_t& neigh_idx_y_b,
      const size_t& neigh_idx_z_a,
      const size_t& neigh_idx_z_b,
      //const int& Nx_total,
      const int& Ny,
      const int& Nz);

   double laplacian(
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
