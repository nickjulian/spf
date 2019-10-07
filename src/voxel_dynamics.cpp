/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: voxel_dynamics.cpp
// Purpose:
//  The functions implemented in this file evaluate changes to the fields 
//   which will be summed to realize their evolution in time.

#ifndef VOXEL_DYNAMICS_CPP
#define VOXEL_DYNAMICS_CPP

#include "voxel_dynamics.hpp"

int SPF_NS::conserved_gaussian_flux( 
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
   )
{
   // assuming periodic boundary conditions

   /////////////////////////////////////////////////////////
   // evaluate local fluxes and changes
   // using local_field[idx]
   // and its neighbors indexed by neigh_x/y/z_idx[] 
   // and accumulate field changes to 
   //  local_change[idx]
   //  and its neighbors local_change[neigh_x/y/z_idx[]]]

   double jump_magnitudes[6]; // assumes ndims == 3
   for (size_t ii=0; ii < 6; ++ii) jump_magnitudes[ii] = 0;

   int exiting_current_voxel; exiting_current_voxel = 0;

   size_t dest_idx;

   // evaluate jump rates in each direction
   // TODO: establish reasoning for these choices ... !
   double jump_rate; 

   // evaluate jump magnitude in each direction

   jump_rate = rate_scale_factor 
                  * (local_field[idx]);// + change_to_current_voxel);
   //std::cout << "jump rate : " << jump_rate << std::endl; //debug
   if ( jump_rate > 0.0 ) 
   {
      std::normal_distribution<double> gd( jump_rate * dt, jump_rate *dt);
      // multiply by 0.5 to compensate for the absolute value 
      // since Gaussian sigma_{2*x}^2 == 2*sigma_x^2 for Gaussian x.
      //exiting_current_voxel = round(abs(0.5* gd( rr.generator )));
      double jump; jump = gd( rr.generator );
      if ( jump < 0 ) jump =0;
      exiting_current_voxel = round( jump );

      //std::list<int> dest_idx_randomized;
      if ( exiting_current_voxel >= local_field[idx] )
      {
         exiting_current_voxel = local_field[idx];
      }

      if ( exiting_current_voxel > 0) 
      {
         for (int ii=0; ii < exiting_current_voxel; ++ii)
         {
            // choose which neighbor to send this walker to
            std::uniform_int_distribution<int> ud(0,5);
            dest_idx = ud( rr.generator );
            jump_magnitudes[dest_idx] += 1.0;
         }
      }
      else 
         if ( exiting_current_voxel < 0 )
         {
            std::cout  << "somehow 'exiting_current_voxel' < 0" 
               << std::endl;
         }

      local_change[neigh_idx_x_a] += jump_magnitudes[0];
      local_change[neigh_idx_x_b] += jump_magnitudes[1];
      local_change[neigh_idx_y_a] += jump_magnitudes[2];
      local_change[neigh_idx_y_b] += jump_magnitudes[3];
      local_change[neigh_idx_z_a] += jump_magnitudes[4];
      local_change[neigh_idx_z_b] += jump_magnitudes[5];
   }
   local_change[idx] -= exiting_current_voxel;
   
   /////////////////////////////////////////////////////////

   return EXIT_SUCCESS;
}


int SPF_NS::conserved_jump_flux( 
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
   )
{
   // assuming periodic boundary conditions

   /////////////////////////////////////////////////////////
   // evaluate local fluxes and changes
   // using local_field[idx]
   // and its neighbors indexed by neigh_x/y/z_idx[] 
   // and accumulate field changes to 
   //  local_change[idx]
   //  and its neighbors local_change[neigh_x/y/z_idx[]]]

   //double jump_magnitude_x_a; // jump_magnitude_x = 0.0;
   //double jump_magnitude_x_b; // jump_magnitude_x = 0.0;
   //double jump_magnitude_y_a; // jump_magnitude_y = 0.0;
   //double jump_magnitude_y_b; // jump_magnitude_y = 0.0;
   //double jump_magnitude_z_a; // jump_magnitude_z = 0.0;
   //double jump_magnitude_z_b; // jump_magnitude_z = 0.0;

   double jump_magnitudes[6]; // assumes ndims == 3
   for (size_t ii=0; ii < 6; ++ii) jump_magnitudes[ii] = 0;

   //double change_to_current_voxel; change_to_current_voxel = 0;
   int exiting_current_voxel; exiting_current_voxel = 0;

   size_t dest_idx;

   // evaluate jump rates in each direction
   // TODO: establish reasoning for these choices ... !
   double jump_rate; 
   //jump_rate = rate_scale_factor * local_field[idx];
   //rr.update_poisson_rate( jump_rate, dt );

   //if ( rr.get_subinterval() != jump_rate * dt * exp( -1.0*jump_rate*dt) )
   //   std::cout << "warning: subinterval not updated : "
   //      << rr.get_subinterval() << " != " <<  jump_rate * dt * exp( -1.0*jump_rate*dt)  << std::endl; // debug
   // evaluate jump magnitude in each direction

   //SPF_NS::random rr2( jump_rate, dt);

   //jump_magnitude_x_a = rr.poisson_event_count( rr.generator );
   //jump_rate = rate_scale_factor * (local_field[idx] - local_change[idx]);
   
   jump_rate = rate_scale_factor 
                  * (local_field[idx]);// + change_to_current_voxel);
   //std::cout << "jump rate : " << jump_rate << std::endl; //debug
   if ( jump_rate > 0.0 ) 
   {
      std::poisson_distribution<int> pd( jump_rate * dt);
      exiting_current_voxel = pd( rr.generator );
      //std::cout << " exiting_current_voxel : " 
      //   << exiting_current_voxel 
      //   << std::endl;
      //std::cout << " random integer: " << pd(rr.generator )
      //   << std::endl;
      //local_change[idx] -= jump_magnitude_x_a;

      //std::list<int> dest_idx_randomized;
      if ( exiting_current_voxel >= local_field[idx] )
      {
         exiting_current_voxel = local_field[idx];
      }

      if ( exiting_current_voxel > 0) 
      {
         for (int ii=0; ii < exiting_current_voxel; ++ii)
         {
            // choose which neighbor to send this walker to
            std::uniform_int_distribution<int> ud(0,5);
            dest_idx = ud( rr.generator );
            jump_magnitudes[dest_idx] += 1.0;
         }
      }
      else 
         if ( exiting_current_voxel < 0 )
         {
            std::cout  << "somehow 'exiting_current_voxel' < 0" 
               << std::endl;
         }

      local_change[neigh_idx_x_a] += jump_magnitudes[0];
      local_change[neigh_idx_x_b] += jump_magnitudes[1];
      local_change[neigh_idx_y_a] += jump_magnitudes[2];
      local_change[neigh_idx_y_b] += jump_magnitudes[3];
      local_change[neigh_idx_z_a] += jump_magnitudes[4];
      local_change[neigh_idx_z_b] += jump_magnitudes[5];
   }
   local_change[idx] -= exiting_current_voxel;
   
   //jump_rate = rate_scale_factor 
   //               * (local_field[idx] + change_to_current_voxel);
   //if ( jump_rate > 0.0 ) 
   //{
   //   //std::cout << "jump rate > 0" << std::endl; //debug
   //   std::poisson_distribution<int> pd1( jump_rate * dt);
   //   jump_magnitude_x_b = pd1( rr.generator );
   //   //local_change[idx] -= jump_magnitude_x_a;
   //   change_to_current_voxel -= jump_magnitude_x_a;
   //   local_change[neigh_idx_x_a] += jump_magnitude_x_a;
   //}

   //jump_rate = rate_scale_factor 
   //               * (local_field[idx] + change_to_current_voxel);
   //if ( jump_rate > 0.0 ) 
   //{
   //   //std::cout << "jump rate > 0" << std::endl; //debug
   //   std::poisson_distribution<int> pd2( jump_rate * dt);
   //   jump_magnitude_x_b = pd2( rr.generator );
   //   //local_change[idx] -= jump_magnitude_x_b;
   //   change_to_current_voxel -= jump_magnitude_x_b;
   //   local_change[neigh_idx_x_b] += jump_magnitude_x_b;
   //}

   //jump_rate = rate_scale_factor 
   //               * (local_field[idx] + change_to_current_voxel);
   //if ( jump_rate > 0.0 ) 
   //{
   //   //std::cout << "jump rate > 0" << std::endl; //debug
   //   std::poisson_distribution<int> pd3( jump_rate * dt);
   //   jump_magnitude_y_a = pd3( rr.generator );
   //   //local_change[idx] -= jump_magnitude_y_a;
   //   change_to_current_voxel -= jump_magnitude_y_a;
   //   local_change[neigh_idx_y_a] += jump_magnitude_y_a;
   //}

   //jump_rate = rate_scale_factor 
   //               * (local_field[idx] + change_to_current_voxel);
   //if ( jump_rate > 0.0 ) 
   //{
   //   //std::cout << "jump rate > 0" << std::endl; //debug
   //   std::poisson_distribution<int> pd4( jump_rate * dt);
   //   jump_magnitude_y_b = pd4( rr.generator );
   //   //local_change[idx] -= jump_magnitude_y_b;
   //   change_to_current_voxel -= jump_magnitude_y_b;
   //   local_change[neigh_idx_y_b] += jump_magnitude_y_b;
   //}

   //jump_rate = rate_scale_factor 
   //               * (local_field[idx] + change_to_current_voxel);
   //if ( jump_rate > 0.0 ) 
   //{
   //   //std::cout << "jump rate > 0" << std::endl; //debug
   //   std::poisson_distribution<int> pd5( jump_rate * dt);
   //   jump_magnitude_z_a = pd5( rr.generator );
   //   //local_change[idx] -= jump_magnitude_z_a;
   //   change_to_current_voxel -= jump_magnitude_z_a;
   //   local_change[neigh_idx_z_a] += jump_magnitude_z_a;
   //}

   //jump_rate = rate_scale_factor 
   //               * (local_field[idx] + change_to_current_voxel);
   //if ( jump_rate > 0.0 ) 
   //{
   //   //std::cout << "jump rate > 0" << std::endl; //debug
   //   std::poisson_distribution<int> pd6( jump_rate * dt);
   //   jump_magnitude_z_b = pd6( rr.generator );
   //   //local_change[idx] -= jump_magnitude_z_b;
   //   change_to_current_voxel -= jump_magnitude_z_b;
   //   local_change[neigh_idx_z_b] += jump_magnitude_z_b;
   //}

   //jump_magnitude_x_b = rr.poisson_event_count( rr.generator );
   //jump_magnitude_y_a = rr.poisson_event_count( rr.generator );
   //jump_magnitude_y_b = rr.poisson_event_count( rr.generator );
   //jump_magnitude_z_a = rr.poisson_event_count( rr.generator );
   //jump_magnitude_z_b = rr.poisson_event_count( rr.generator );

   // TODO: check that the following won't cause a jump to below 0
   //local_change[idx] -= jump_magnitude_x_a;
                        // + jump_magnitude_x_b 
                        //+ jump_magnitude_y_a + jump_magnitude_y_b
                        //+ jump_magnitude_z_a + jump_magnitude_z_b;

   //if ( local_field[idx] + change_to_current_voxel < 0 )
   //{
   //   std::cout << "warning: jump quantity exiting voxel " << idx 
   //      << " is larger than the quantity present." << std::endl;

   //   local_change[idx] = -1.0 * local_field[idx];
   //}
   //else local_change[idx] += change_to_current_voxel;

   //local_change[neigh_idx_x_a] += jump_magnitude_x_a;
   //local_change[neigh_idx_x_b] += jump_magnitude_x_b;
   //local_change[neigh_idx_y_a] += jump_magnitude_y_a;
   //local_change[neigh_idx_y_b] += jump_magnitude_y_b;
   //local_change[neigh_idx_z_a] += jump_magnitude_z_a;
   //local_change[neigh_idx_z_b] += jump_magnitude_z_b;

   /////////////////////////////////////////////////////////

   return EXIT_SUCCESS;
}


double SPF_NS::laplacian(
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
         )
{
   double ll; ll = 0.0;
   ll +=
     (1.0/(hh_x*hh_x))*
     ( 
           //local_field[neigh_x_idx[0]]
           local_field[neigh_idx_x_a]
           -2.0* local_field[idx]
           + local_field[neigh_idx_x_b]
           //+ local_field[neigh_x_idx[1]]
     );

   ll +=
     (1.0/(hh_y*hh_y))*
     (
           //local_field[neigh_y_idx[0]]
           local_field[neigh_idx_y_a]
           -2.0* local_field[idx]
           + local_field[neigh_idx_y_b]
      );

   if (Nz > 1)
   {
      ll +=
        (1.0/(hh_z*hh_z))*
        (
              local_field[neigh_idx_z_a]
              -2.0* local_field[idx]
              + local_field[neigh_idx_z_b]
         );
   }

   return ll;
}

int SPF_NS::laplacian_flux( 
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
   )
{
   // assuming periodic boundary conditions
   //size_t idx; 
   //size_t neigh_x_idx[2];
   //size_t neigh_y_idx[2];
   //size_t neigh_z_idx[2];

   // index of local_field, non-ghost elements
   //idx = kk + Nz*(jj + Ny*ii);

   // indices wrt local_change
   //identify_local_neighbors(
   //      neigh_x_idx, neigh_y_idx, neigh_z_idx,
   //      ii, jj, kk,
   //      //Nx_total, 
   //      Ny, Nz
   //      );

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

   // For the dynamics to be stable, require
   //  \delta t < (h^2)/(6*diffusivity)
   
   local_change[ idx ] 
      += diffusivity * dt *
         laplacian(
                  hh_x, hh_y, hh_z,
                  Nz,
                  local_field,
                  idx,
                  neigh_idx_x_a, 
                  neigh_idx_x_b, 
                  neigh_idx_y_a, 
                  neigh_idx_y_b, 
                  neigh_idx_z_a,
                  neigh_idx_z_b
                  //neigh_x_idx, neigh_y_idx, neigh_z_idx
                  );

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
         //neigh_x_idx, neigh_y_idx, neigh_z_idx,
         neigh_x_idx[0], 
         neigh_x_idx[1], 
         neigh_y_idx[0], 
         neigh_y_idx[1], 
         neigh_z_idx[0],
         neigh_z_idx[1],
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
   //size_t* const neigh_x_idx,
   //size_t* const neigh_y_idx,
   //size_t* const neigh_z_idx,
   size_t& neigh_idx_x_a,
   size_t& neigh_idx_x_b,
   size_t& neigh_idx_y_a,
   size_t& neigh_idx_y_b,
   size_t& neigh_idx_z_a,
   size_t& neigh_idx_z_b,
   const size_t& ii,
   const size_t& jj,
   const size_t& kk,
   //const int& Nx_total,
   const int& Ny,
   const int& Nz
   )
{
   neigh_idx_x_a = kk + Nz*(jj + Ny*(ii -1));
   neigh_idx_x_b = kk + Nz*(jj + Ny*(ii +1));

   if ( jj == 0 ) // below along y, periodic boundary
      neigh_idx_y_a = kk + Nz*((Ny -1) + Ny*ii);
   else
      neigh_idx_y_a = kk + Nz*((jj -1) + Ny*ii);

   if ( jj == Ny -1 ) // above along y, periodic boundary
      neigh_idx_y_b = kk + Nz*(0 + Ny*ii);
   else
      neigh_idx_y_b = kk + Nz*((jj +1) + Ny*ii);

   if ( kk == 0 ) // below along z, periodic boundary
      neigh_idx_z_a = (Nz -1) + Nz*(jj + Ny*ii);
   else
      neigh_idx_z_a = (kk -1) + Nz*(jj + Ny*ii);

   if ( kk == Nz -1 ) // above along z, periodic boundary
      neigh_idx_z_b = 0 + Nz*(jj + Ny*ii);
   else
      neigh_idx_z_b = (kk +1) + Nz*(jj + Ny*ii);
   //neigh_x_idx[0] = kk + Nz*(jj + Ny*(ii -1));
   //neigh_x_idx[1] = kk + Nz*(jj + Ny*(ii +1));

   //if ( jj == 0 ) // below along y, periodic boundary
   //   neigh_y_idx[0] = kk + Nz*((Ny -1) + Ny*ii);
   //else
   //   neigh_y_idx[0] = kk + Nz*((jj -1) + Ny*ii);

   //if ( jj == Ny -1 ) // above along y, periodic boundary
   //   neigh_y_idx[1] = kk + Nz*(0 + Ny*ii);
   //else
   //   neigh_y_idx[1] = kk + Nz*((jj +1) + Ny*ii);

   //if ( kk == 0 ) // below along z, periodic boundary
   //   neigh_z_idx[0] = (Nz -1) + Nz*(jj + Ny*ii);
   //else
   //   neigh_z_idx[0] = (kk -1) + Nz*(jj + Ny*ii);

   //if ( kk == Nz -1 ) // above along z, periodic boundary
   //   neigh_z_idx[1] = 0 + Nz*(jj + Ny*ii);
   //else
   //   neigh_z_idx[1] = (kk +1) + Nz*(jj + Ny*ii);

   return EXIT_SUCCESS;
}

/////////////////////////////////////////////////////////
// TODO
// outward flux - Brownian motion
// outward flux - jump process
// non-conserved changes
/////////////////////////////////////////////////////////

#endif
