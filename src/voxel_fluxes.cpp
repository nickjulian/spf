/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: voxel_fluxes.cpp
// Purpose:

#ifndef VOXEL_FLUXES_CPP
#define VOXEL_FLUXES_CPP

#include "voxel_fluxes.hpp"


int SPF_NS::laplacian_flux( 
   std::vector<double>& local_change, // must be same size as local_field
   const std::vector<double>& local_field,
   const double& diffusivity,
   const size_t& ii,
   const size_t& jj,
   const size_t& kk,
   const int& Nx_total,
   const int& Ny,
   const int& Nz
   )
{
   // assuming periodic boundary conditions
   size_t idx; 
   size_t neigh_x_idx[2];
   size_t neigh_y_idx[2];
   size_t neigh_z_idx[2];

   // index of local_field, non-ghost elements
   idx = kk + Nz*(jj + Ny*ii);

   // indices wrt local_change
   identify_local_neighbors(
         neigh_x_idx, neigh_y_idx, neigh_z_idx,
         ii, jj, kk,
         //Nx_total, 
         Ny, Nz
         );

   /////////////////////////////////////////////////////////
   // evaluate local fluxes and changes
   // using local_field[idx]
   // and its neighbors indexed by neigh_x/y/z_idx[] 
   // and accumulate field changes to 
   //  local_change[idx]
   //  and its neighbors local_change[neigh_x/y/z_idx[]]]

   // change in local_field[idx] may be estimated by a 
   //  function of its neighboring values
   
   //   The second spatial derivative along x may be estimated by 
   //    (1/h^2)*( u_{i+1,j} -2*u_{i,j}) + u_{i-1,j} )
   //    where h = 1/(n+1), and n+2 is the number of points along x
   //   some function of 
   //    local_field[neigh_x_[0]], local_field[neigh_x_[1]]
   //    local_field[neigh_y_[0]], local_field[neigh_y_[1]]
   //    and if Nz >1, 
   //    local_field[neigh_z_[0]] and local_field[neigh_z_[1]]
   
   local_change[idx] += //local_field[ idx ];
       diffusivity*(Nx_total -1)*(Nx_total -1)*    // 1/h^2
          (   
             local_field[neigh_x_idx[0]]
             -2.0* local_field[idx]
             + local_field[neigh_x_idx[1]]
             );

   local_change[idx] +=
       diffusivity*(Ny -1)*(Ny -1)*(   // 1/h^2
             local_field[neigh_y_idx[0]]
             -2.0* local_field[idx]
             + local_field[neigh_y_idx[1]]
             );
   if ( Nz >1 )
   {
      local_change[idx] +=
          diffusivity*(Nz -1)*(Nz -1)*(   // 1/h^2
                local_field[neigh_z_idx[0]]
                -2.0* local_field[idx]
                + local_field[neigh_z_idx[1]]
                );
   }
   /////////////////////////////////////////////////////////

   return EXIT_SUCCESS;
}

int SPF_NS::flux_test( 
   std::vector<double>& local_change,
   const std::vector<double>& local_field,
   const size_t& ii,
   const size_t& jj,
   const size_t& kk,
   const int& Nx_total,
   const int& Ny,
   const int& Nz
   )
{
   // assuming periodic boundary conditions
   size_t idx; 
   size_t neigh_x_idx[2];
   size_t neigh_y_idx[2];
   size_t neigh_z_idx[2];

   // index of local_field, non-ghost elements
   idx = kk + Nz*(jj + Ny*ii);

   
   // indices wrt local_change
   identify_local_neighbors(
         neigh_x_idx, neigh_y_idx, neigh_z_idx,
         ii, jj, kk,
         //Nx_total, 
         Ny, Nz
         );

   /////////////////////////////////////////////////////////
   // evaluate local fluxes and changes
   // using local_field[idx]
   // and its neighbors indexed by neigh_x/y/z_idx[] 
   // and accumulate field changes to 
   //  local_change[idx]
   //  and its neighbors local_change[neigh_x/y/z_idx[]]]

   // testing with constant increments
   local_change[neigh_x_idx[0]] += 0.01 * local_field[idx];
   local_change[neigh_x_idx[1]] += 0.01 * local_field[idx];
   local_change[neigh_y_idx[0]] += 0.0001 * local_field[idx];
   local_change[neigh_y_idx[1]] += 0.0001 * local_field[idx];
   if ( Nz > 1)
   {
      local_change[neigh_z_idx[0]] 
         += 0.000001 * local_field[idx];
      local_change[neigh_z_idx[1]] 
         += 0.000001 * local_field[idx];
   }
   /////////////////////////////////////////////////////////

   return EXIT_SUCCESS;
}

int SPF_NS::identify_local_neighbors( 
   size_t* const neigh_x_idx,
   size_t* const neigh_y_idx,
   size_t* const neigh_z_idx,
   const size_t& ii,
   const size_t& jj,
   const size_t& kk,
   //const int& Nx_total,
   const int& Ny,
   const int& Nz
   )
{
   neigh_x_idx[0] = kk + Nz*(jj + Ny*(ii -1));
   neigh_x_idx[1] = kk + Nz*(jj + Ny*(ii +1));

   if ( jj == 0 ) // below along y, periodic boundary
      neigh_y_idx[0] = kk + Nz*((Ny -1) + Ny*ii);
   else
      neigh_y_idx[0] = kk + Nz*((jj -1) + Ny*ii);

   if ( jj == Ny -1 ) // above along y, periodic boundary
      neigh_y_idx[1] = kk + Nz*(0 + Ny*ii);
   else
      neigh_y_idx[1] = kk + Nz*((jj +1) + Ny*ii);

   if ( kk == 0 ) // below along z, periodic boundary
      neigh_z_idx[0] = (Nz -1) + Nz*(jj + Ny*ii);
   else
      neigh_z_idx[0] = (kk -1) + Nz*(jj + Ny*ii);

   if ( kk == Nz -1 ) // above along z, periodic boundary
      neigh_z_idx[1] = 0 + Nz*(jj + Ny*ii);
   else
      neigh_z_idx[1] = (kk +1) + Nz*(jj + Ny*ii);

   return EXIT_SUCCESS;
}

/////////////////////////////////////////////////////////
// TODO
// outward flux - Brownian motion
// outward flux - jump process
// non-conserved changes
/////////////////////////////////////////////////////////

#endif
