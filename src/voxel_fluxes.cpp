/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: voxel_fluxes.cpp
// Purpose:

#ifndef VOXEL_FLUXES_CPP
#define VOXEL_FLUXES_CPP

#include "voxel_fluxes.hpp"

int SPF_NS::flux_test( 
   std::vector<double>& local_change,
   const std::vector<double>& local_field,
   const size_t& ii,
   const size_t& jj,
   const size_t& kk,
   const int& Nx,
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

   //if ( (ii==0) && (jj==0) && (kk==0) ) // debug
   //{
   //   cout 
   //      << "node " << mynode 
   //      << " Nx_local+2,Ny,Nz : " 
   //      << Nx_local + 2 << ", " << Ny << ", " << Nz << endl
   //      << "node " << mynode << "i,j,k " 
   //      << ii << ", " << jj << ", " << kk << endl
   //      << "node " << mynode 
   //      << " neigh_x_idx[0] : " << neigh_x_idx[0] << endl
   //      << "node " << mynode 
   //      << " neigh_x_idx[1] : " << neigh_x_idx[1] << endl
   //      << "node " << mynode 
   //      << " neigh_y_idx[0] : " << neigh_y_idx[0] << endl
   //      << "node " << mynode 
   //      << " neigh_y_idx[1] : " << neigh_y_idx[1] << endl
   //      << "node " << mynode 
   //      << " neigh_z_idx[0] : " << neigh_z_idx[0] << endl
   //      << "node " << mynode 
   //      << " neigh_z_idx[1] : " << neigh_z_idx[1] << endl;
   //}

   /////////////////////////////////////////////////////////
   // evaluate local fluxes and changes
   // using local_field[idx]
   // and its neighbors indexed by neigh_x/y/z_idx[] 
   // and accumulate field changes to 
   //  local_change[idx]
   //  and its neighbors local_change[neigh_x/y/z_idx[]]]

   // debug
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
   // end debug
   /////////////////////////////////////////////////////////

   return EXIT_SUCCESS;
}

/////////////////////////////////////////////////////////
// TODO
// outward flux - Brownian motion
// outward flux - jump process
// non-conserved changes
/////////////////////////////////////////////////////////

#endif
