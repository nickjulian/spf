/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: voxel_fluxes.hpp
// Purpose:

#ifndef VOXEL_FLUXES_HPP
#define VOXEL_FLUXES_HPP

#include <cstdlib>   // EXIT_SUCCESS, EXIT_FAILURE
#include <vector>

namespace SPF_NS
{
   int flux_test( 
      std::vector<double>& local_change,
      const std::vector<double>& local_field,
      const size_t& ii,
      const size_t& jj,
      const size_t& kk,
      const int& Nx,
      const int& Ny,
      const int& Nz
      );

} // SPF_NS
#endif
