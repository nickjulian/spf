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

#ifndef ONESIXTH 
#define ONESIXTH 0.16666666666666666666666666666666666666666666666666666667
#endif

int SPF_NS::randomize_neighbor_order(
      std::vector<size_t>& neigh_order,  // 6 elements
      SPF_NS::random& rr,
      std::uniform_real_distribution<double>& rand_decimal,
      std::vector<double>& rand_decimals1, // reused, 6 elements
      std::vector<double>& rand_decimals2, // reused, 6 elements
      const int_flags& flags
      )
{
   //std::vector<double> rand_decimals1(6,0);
   //std::vector<double> rand_decimals2(6,0);

   bool randomized_flag; randomized_flag = true;
   //std::uniform_real_distribution<double> rand_decimal(0,1);
   do {
      for( size_t ii=0; ii < 6; ++ii) 
      {
         rand_decimals1[ii] = rand_decimal( rr.generator );
         rand_decimals2[ii] = rand_decimals1[ii];
      }
      std::sort( rand_decimals1.begin(), rand_decimals1.end());

      // check that no number was drawn twice by rand_decimal()
      std::vector<double>::iterator decimals_itr;
      decimals_itr = std::unique(rand_decimals1.begin(), 
                                 rand_decimals1.end() );
      rand_decimals1.resize( 
               std::distance( rand_decimals1.begin(), decimals_itr) );
      if ( rand_decimals1.size() != 6 )
      {
         randomized_flag = false;
         if ( flags.debug )
         {
            std::cout << "note: rand_decimals1 not unique;" 
                      << " running randomization loop again" << std::endl;
         }
      }
   }
   while ( ! randomized_flag );

   for( size_t ii=0; ii < 6; ++ii)
   {
      for( size_t jj=0; jj < 6; ++jj)
      {
            if ( rand_decimals1[jj] == rand_decimals2[ii]) 
            {
               neigh_order[jj] = ii;
               continue;
            }
      }
   }

   return EXIT_SUCCESS;
}

int SPF_NS::enforce_bounds_generic(
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
         )
{
   // updates phi_local_flux with acceptable flux values

   // Considered as an outward flux, neighbor orders are 
   //  equally balanced (assuming equal barrier heights).
   // But when considered as an inward flux, the neighbor
   //  orders may be balanced using first passage 
   //  distributions.

   double outward_flux, inward_flux;
   //std::vector<size_t> neigh_idxs(Nvoxel_neighbors, 0);
   std::vector<size_t> neigh_idxs(6, 0);  // TODO: make 2-D compatable
   //std::vector<size_t> neigh_order(Nvoxel_neighbors, 0);
   std::uniform_real_distribution<double> rand_decimal(0,1);// for order
   std::vector<double> rand_decimals1(Nvoxel_neighbors, 0);
   std::vector<double> rand_decimals2(Nvoxel_neighbors, 0);
   //std::vector<double> flux_reduction_factors(6,1);
   double flux_reduction_factor; flux_reduction_factor =1;
   double total_flux_rate;
   std::vector<size_t> neigh_pairs(Nvoxel_neighbors, 0);
   neigh_pairs[0] = 1;  // x upward
   neigh_pairs[1] = 0;  // x downward
   size_t idx; idx = 0;

   if ( Nvoxel_neighbors >= 2 )
   {
      neigh_pairs[2] = 3;  // y upward
      neigh_pairs[3] = 2;  // y downward
   }
   if ( Nvoxel_neighbors >= 6 )
   {
      neigh_pairs[4] = 5;  // z upward
      neigh_pairs[5] = 4;  // z downward
   }

   if ( Nvoxel_neighbors != 6)
   {
      std::cout << "Error: 2-D enforce_bounds_...() not yet"
        << " compatible with 2-D" << std::endl;
      return EXIT_FAILURE;
   }

   // iterate over local voxels
   for (size_t ii=1; ii < Nx_local +1; ++ii) // loop over non-ghosts
      for ( size_t jj=0; jj < Ny; ++jj)
         for ( size_t kk=0; kk < Nz; ++kk)
         {
            idx = kk + Nz*(jj + Ny*ii);
            identify_local_neighbors(
                  neigh_idxs[0], 
                  neigh_idxs[1], 
                  neigh_idxs[2], 
                  neigh_idxs[3], 
                  neigh_idxs[4],
                  neigh_idxs[5],
                  ii, jj, kk,
                  Ny, Nz
                  );

            //randomize_neighbor_order(
            //      neigh_order,
            //      rr,   // random generator
            //      rand_decimal,  // uniform_distribution<int>
            //      rand_decimals1,// reused vector, not useful outside
            //      rand_decimals2, // reused vector, not useful outside
            //      flags
            //      );

            outward_flux = 0;
            for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
            {
               outward_flux += phi_local_flux[nn + Nvoxel_neighbors*idx];
            }
            if ( phi_local[idx] - outward_flux < phi_lower_limit )
            {
               // debug
               //std::cout << "reducing outward flux " 
               //   << "; phi_local[" << idx << "] - outward_flux : "
               //   << phi_local[idx] << " - " << outward_flux 
               //   << std::endl;
               // end debug
               total_flux_rate = 0;
               //for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
               //{
               //   total_flux_rate 
               //      += phi_local_rates[nn  + Nvoxel_neighbors*idx];
               //}
               for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
               {
                  // balance outward fluxes to not exceed phi_local[idx]
                  //flux_reduction_factor
                  phi_local_flux[nn + Nvoxel_neighbors*idx]
                     *= //phi_local_rates[nn + Nvoxel_neighbors*idx]
                        (phi_local[idx] - phi_lower_limit)
                        //* phi_local_flux[nn + Nvoxel_neighbors*idx]
                        / outward_flux;
                     //   / total_flux_rate;
                     // = phi_local_flux[nn + Nvoxel_neighbors*idx]
                        // / outward_flux;
                  
                  //phi_local_flux[nn + Nvoxel_neighbors*idx]
                  //   *= flux_reduction_factor;
                  //phi_local_flux[nn + Nvoxel_neighbors*idx]
                  //   = 
                  //      (phi_local[idx] - phi_lower_limit)
                  //         * phi_local_rates[nn + Nvoxel_neighbors*idx]
                  //         / total_flux_rate;
               }
               // debug
               //std::cout << "Check to see if outward_flux was appropriately normalized; total_flux_rate " << total_flux_rate
               //   << ", phi_local[" << idx << "] "
               //   << phi_local[idx]
               //   << ", initial outward_flux "
               //   << outward_flux
               //   //<< ", flux_reduction_factor "
               //   //<< flux_reduction_factor
               //   << ", phi_local_rates[]:";
               //for (size_t nn=0; nn < Nvoxel_neighbors; ++nn)
               //      std::cout << phi_local_rates[nn +Nvoxel_neighbors*idx]
               //         << ", ";
               //std::cout << std::endl;
               //outward_flux = 0;
               //for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
               //{
               //   outward_flux += phi_local_flux[nn +Nvoxel_neighbors*idx];
               //}
               //std::cout << " renormalized total outward_flux "
               //   << outward_flux 
               //   << std::endl;
               // end debug
            }
            // debug
            //outward_flux = 0;
            //for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
            //{
            //   outward_flux += phi_local_flux[nn + Nvoxel_neighbors*idx];
            //}
            //if ( phi_local[idx] - outward_flux < phi_lower_limit )
            //{
            //   std::cout << "Error: outward_flux was not appropriately normalized; total_flux_rate " << total_flux_rate
            //      << ", phi_local[" << idx << "] "
            //      << phi_local[idx]
            //      << ", outward_flux "
            //      << outward_flux
            //      << ", phi_local_rates[]:";
            //   for (size_t nn=0; nn < Nvoxel_neighbors; ++nn)
            //      std::cout << phi_local_rates[nn + Nvoxel_neighbors*idx]
            //            << ", ";
            //   std::cout << std::endl;
            //}
            // end debug

            // Iterate over neighbors of local voxels.
            // Ensure inward flux + local population doesn't exceed upper 
            //  limit.

            inward_flux = 0;
            //size_t mm; mm=0;
            for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
            {
               inward_flux += phi_local_flux[ 
                                    neigh_pairs[nn]   // mm 
                                    + Nvoxel_neighbors * neigh_idxs[nn]
                                 ];
            }
            if ( inward_flux + phi_local[idx] > phi_upper_limit )
            {  // normalize the inward fluxes by their relative rates
               total_flux_rate = 0;
               for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
               {
                  total_flux_rate 
                     += phi_local_rates[
                                 neigh_pairs[nn]   // mm 
                                    + Nvoxel_neighbors * neigh_idxs[nn]];
               }
               for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
               {
                  phi_local_flux[ 
                                 neigh_pairs[nn]   // mm 
                                    + Nvoxel_neighbors * neigh_idxs[nn]]
                     *= (phi_upper_limit - phi_local[idx])
                        * phi_local_rates[
                                 neigh_pairs[nn]   // mm 
                                    + Nvoxel_neighbors * neigh_idxs[nn]
                           ] / total_flux_rate;
               }
            }

            if (flags.debug != 0)
            { // debug
               outward_flux = 0;
               inward_flux = 0;
               for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
               {
                  outward_flux += phi_local_flux[
                                    nn + Nvoxel_neighbors*idx];
                  inward_flux += phi_local_flux[ 
                                    neigh_pairs[nn]   // mm 
                                    + Nvoxel_neighbors * neigh_idxs[nn]
                                 ];
               }

               if ( phi_local[idx] - outward_flux < phi_lower_limit )
               {
                  // TODO: the following is just a guess of error limit
                  if (abs(phi_local[idx] - phi_lower_limit - outward_flux)
                       > (eps.dblsqrt* (phi_local[idx] - phi_lower_limit))
                     )
                     std::cout << "Warning: reduced outward flux still greater than voxel contents. (outward_flux, phi_local[" << idx << "]): (" << outward_flux << ", " << phi_local[idx] << ")" << std::endl;
                  // else phi_local[idx] =0 if <0 in parent function
               }

               if ( inward_flux + phi_local[idx] > phi_upper_limit )
               {
                  // TODO: the following is just a guess of error limit
                  if ( (inward_flux + phi_local[idx] - phi_upper_limit)
                       > (eps.dblsqrt* (phi_local[idx] - phi_upper_limit))
                     )
                     std::cout << "Warning: reduced inward flux + previous population still greater than voxel upper limit. (inward_flux, phi_local[" << idx << "]): (" << inward_flux << ", " << phi_local[idx] << ")"
                        << std::endl;
               }
            } // end debug
         }
   return EXIT_SUCCESS;
}

//int SPF_NS::conserved_gaussian_flux_single_distribution_milstein(
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
//   )
//{
//   // assuming periodic boundary conditions
//
//   /////////////////////////////////////////////////////////
//   // evaluate local fluxes and changes
//   // using local_field[idx]
//   // and its neighbors indexed by neigh_x/y/z_idx[] 
//   // and accumulate field changes to 
//   //  local_change[idx]
//   //  and its neighbors local_change[neigh_x/y/z_idx[]]]
//
//   double jump_magnitudes[6]; // assumes ndims == 3
//   for (size_t ii=0; ii < 6; ++ii) jump_magnitudes[ii] = 0;
//
//   int exiting_current_voxel; exiting_current_voxel = 0;
//
//   //size_t dest_idx;
//
//   // evaluate jump rates in each direction
//   // TODO: establish reasoning for these choices ... !
//   double jump_rate; 
//
//   jump_rate = //rate_scale_factor * 
//                  (local_field[idx]);// + change_to_current_voxel);
//   
//   // evaluate jump magnitude in each direction
//   double total_exiting_current_voxel, dw, drift;
//   total_exiting_current_voxel  = 0;
//   drift = jump_rate;
//
//   //std::cout << "jump rate : " << jump_rate << std::endl; //debug
//   if ( jump_rate > 0.0 ) 
//   {
//      //std::normal_distribution<double> gd(jump_rate *dt, jump_rate *dt);
//      std::normal_distribution<double> gd( 0.0, 1.0 );
//      dw = sqrt(dt) * gd( rr.generator );
//
//      total_exiting_current_voxel
//         = (drift * dt)               // f(x) dt
//            + (jump_rate * dw)        // g(x) dw
//            + (0.5 * jump_rate //* rate_scale_factor // + 0.5 g(x) g'(x)
//                * (dw * dw - dt));    // *( (dw)^2 - dt )
//
//      if ( total_exiting_current_voxel > local_field[idx])
//      {
//         total_exiting_current_voxel = local_field[idx];
//      }
//         
//      if ( total_exiting_current_voxel < 0)
//      {
//         total_exiting_current_voxel  = 0;
//      }
//      else
//      {
//         // choose how to distribute the flux among neighbors
//         std::vector<double> random_points;
//         std::uniform_real_distribution<double> ud(0,1);
//         for ( size_t ii=0; ii < 5; ++ii)
//         {
//            random_points.push_back( ud( rr.generator ) );
//         }
//
//         // sort the random points into ascending order
//         std::sort( random_points.begin(), random_points.end() );
//         jump_magnitudes[0] = total_exiting_current_voxel
//                                 * random_points[0];
//         jump_magnitudes[1] = total_exiting_current_voxel
//                                 * (random_points[1] - random_points[0]);
//         jump_magnitudes[2] = total_exiting_current_voxel
//                                 * (random_points[2] - random_points[1]);
//         jump_magnitudes[3] = total_exiting_current_voxel
//                                 * (random_points[3] - random_points[2]);
//         jump_magnitudes[4] = total_exiting_current_voxel
//                                 * (random_points[4] - random_points[3]);
//         jump_magnitudes[5] = total_exiting_current_voxel
//                                 * (1.0 - random_points[4]);
//      }
//
//      local_change[neigh_idx_x_a] += jump_magnitudes[0];
//      local_change[neigh_idx_x_b] += jump_magnitudes[1];
//      local_change[neigh_idx_y_a] += jump_magnitudes[2];
//      local_change[neigh_idx_y_b] += jump_magnitudes[3];
//      local_change[neigh_idx_z_a] += jump_magnitudes[4];
//      local_change[neigh_idx_z_b] += jump_magnitudes[5];
//
//   }
//   local_change[idx] -= total_exiting_current_voxel;
//   
//   /////////////////////////////////////////////////////////
//
//   return EXIT_SUCCESS;
//}

//int SPF_NS::conserved_gaussian_flux_single_distribution(
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
//   )
//{
//   // assuming periodic boundary conditions
//
//   /////////////////////////////////////////////////////////
//   // evaluate local fluxes and changes
//   // using local_field[idx]
//   // and its neighbors indexed by neigh_x/y/z_idx[] 
//   // and accumulate field changes to 
//   //  local_change[idx]
//   //  and its neighbors local_change[neigh_x/y/z_idx[]]]
//
//   double jump_magnitudes[6]; // assumes ndims == 3
//   for (size_t ii=0; ii < 6; ++ii) jump_magnitudes[ii] = 0;
//
//   int exiting_current_voxel; exiting_current_voxel = 0;
//
//   //size_t dest_idx;
//
//   // evaluate jump rates in each direction
//   // TODO: establish reasoning for these choices ... !
//   double jump_rate; 
//
//   jump_rate = //rate_scale_factor * 
//                  (local_field[idx]);// + change_to_current_voxel);
//   
//   // evaluate jump magnitude in each direction
//   double total_exiting_current_voxel, dw, drift;
//   drift = jump_rate;
//
//   //std::cout << "jump rate : " << jump_rate << std::endl; //debug
//   if ( jump_rate > 0.0 ) 
//   {
//      //std::normal_distribution<double> gd(jump_rate *dt, jump_rate *dt);
//      //std::normal_distribution<double> gd( drift, jump_rate );
//      std::normal_distribution<double> gd( 0.0, 1.0);
//      dw = sqrt(dt) * gd( rr.generator );
//
//      //total_exiting_current_voxel = dw;
//      total_exiting_current_voxel = (drift * dt) + jump_rate * dw;
//
//      if ( total_exiting_current_voxel > local_field[idx])
//      {
//         total_exiting_current_voxel = local_field[idx];
//      }
//         
//      if ( total_exiting_current_voxel < 0)
//      {
//         total_exiting_current_voxel  = 0;
//      }
//      else
//      {
//         // choose how to distribute the flux among neighbors
//         std::vector<double> random_points;
//         std::uniform_real_distribution<double> ud(0,1);
//         for ( size_t ii=0; ii < 5; ++ii)
//         {
//            random_points.push_back( ud( rr.generator ) );
//         }
//
//         // sort the random points into ascending order
//         std::sort( random_points.begin(), random_points.end() );
//         jump_magnitudes[0] = total_exiting_current_voxel
//                                 * random_points[0];
//         jump_magnitudes[1] = total_exiting_current_voxel
//                                 * (random_points[1] - random_points[0]);
//         jump_magnitudes[2] = total_exiting_current_voxel
//                                 * (random_points[2] - random_points[1]);
//         jump_magnitudes[3] = total_exiting_current_voxel
//                                 * (random_points[3] - random_points[2]);
//         jump_magnitudes[4] = total_exiting_current_voxel
//                                 * (random_points[4] - random_points[3]);
//         jump_magnitudes[5] = total_exiting_current_voxel
//                                 * (1.0 - random_points[4]);
//      }
//
//      local_change[neigh_idx_x_a] += jump_magnitudes[0];
//      local_change[neigh_idx_x_b] += jump_magnitudes[1];
//      local_change[neigh_idx_y_a] += jump_magnitudes[2];
//      local_change[neigh_idx_y_b] += jump_magnitudes[3];
//      local_change[neigh_idx_z_a] += jump_magnitudes[4];
//      local_change[neigh_idx_z_b] += jump_magnitudes[5];
//
//   }
//   local_change[idx] -= total_exiting_current_voxel;
//   
//   /////////////////////////////////////////////////////////
//
//   return EXIT_SUCCESS;
//}

//int SPF_NS::conserved_gaussian_flux_separate_distributions_gradient_milstein(
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
//   )
//{
//   // assuming periodic boundary conditions
//   // WARNING: this is assymetrical and could lead to excess flow
//   //          along the negative direction along each axis
//   // TODO: maybe find a better solution
//
//   /////////////////////////////////////////////////////////
//   // evaluate local fluxes and changes
//   // using local_field[idx]
//   // and its neighbors indexed by neigh_x/y/z_idx[] 
//   // and accumulate field changes to 
//   //  local_change[idx]
//   //  and its neighbors local_change[neigh_x/y/z_idx[]]]
//
//   //double jump_magnitudes[6]; // assumes ndims == 3
//   double drift[6]; // assumes ndims == 3
//   double jump_variance[6]; // assumes ndims == 3
//   double jump_rate[6]; // assumes ndims == 3
//   double dw[6]; // assumes ndims == 3
//   double exiting_current_voxel_x_a; exiting_current_voxel_x_a = 0;
//   //double exiting_current_voxel_x_b; exiting_current_voxel_x_b = 0;
//   double exiting_current_voxel_y_a; exiting_current_voxel_y_a = 0;
//   //double exiting_current_voxel_y_b; exiting_current_voxel_y_b = 0;
//   double exiting_current_voxel_z_a; exiting_current_voxel_z_a = 0;
//   //double exiting_current_voxel_z_b; exiting_current_voxel_z_b = 0;
//
//   //for (size_t ii=0; ii < 6; ++ii) jump_magnitudes[ii] = 0;
//
//   // evaluate jump rates in each direction
//   // TODO: establish reasoning for these choices ... !
//
//   drift[0] = //rate_scale_factor * 
//                  ONESIXTH 
//                  * (local_field[idx] - local_field[neigh_idx_x_a]);
//   //drift[1] = rate_scale_factor 
//   //               * (local_field[idx] - local_field[neigh_idx_x_b]);
//   drift[2] = //rate_scale_factor * 
//                  ONESIXTH 
//                  * (local_field[idx] - local_field[neigh_idx_y_a]);
//   //drift[3] = rate_scale_factor 
//   //               * (local_field[idx] - local_field[neigh_idx_y_b]);
//   drift[4] = //rate_scale_factor * 
//                  ONESIXTH 
//                  * (local_field[idx] - local_field[neigh_idx_z_a]);
//   //drift[5] = rate_scale_factor 
//   //               * (local_field[idx] - local_field[neigh_idx_z_b]);
//
//   jump_rate[0] = ONESIXTH //* rate_scale_factor 
//                     * local_field[idx];
//   //jump_rate[1] = ONESIXTH * rate_scale_factor * local_field[idx];
//   jump_rate[2] = ONESIXTH //* rate_scale_factor 
//                     * local_field[idx];
//   //jump_rate[3] = ONESIXTH * rate_scale_factor * local_field[idx];
//   jump_rate[4] = ONESIXTH //* rate_scale_factor 
//                     * local_field[idx];
//   //jump_rate[5] = ONESIXTH * rate_scale_factor * local_field[idx];
//
//   jump_variance[0] = jump_rate[0];
//   //jump_variance[1] = jump_rate[1];
//   jump_variance[2] = jump_rate[2];
//   //jump_variance[3] = jump_rate[3];
//   jump_variance[4] = jump_rate[4];
//   //jump_variance[5] = jump_rate[5];
//
//   //total_exiting_current_voxel
//   //      = drift * dt               // f(x) dt
//   //         + jump_rate * dw        // g(x) dw
//   //         + 0.5 * jump_rate * rate_scale_factor // + 0.5 g(x) g'(x)
//   //            * (dw * dw - dt);    // *( (dw)^2 - dt )
//
//   std::normal_distribution<double> gd( 0.0, 1.0 );
//
//   dw[0] = sqrt(dt) * gd( rr.generator );
//   //dw[1] = sqrt(dt) * gd( rr.generator );
//   dw[2] = sqrt(dt) * gd( rr.generator );
//   //dw[3] = sqrt(dt) * gd( rr.generator );
//   dw[4] = sqrt(dt) * gd( rr.generator );
//   //dw[5] = sqrt(dt) * gd( rr.generator );
//
//   exiting_current_voxel_x_a 
//      = drift[0] * dt                              // f(x) dt
//         + jump_rate[0] * dw[0]                    // + g(x) dw
//         + 0.5 * jump_rate[0]                      // + 0.5 g(x) 
//           * ONESIXTH //* rate_scale_factor           // * g'(x)
//           * ( dw[0] * dw[0] - dt);                // *( (dw)^2 - dt )
//   //exiting_current_voxel_x_b
//   exiting_current_voxel_y_a 
//      = drift[2] * dt                              // f(x) dt
//         + jump_rate[2] * dw[2]                    // + g(x) dw
//         + 0.5 * jump_rate[2]                      // + 0.5 g(x)
//           * ONESIXTH //* rate_scale_factor           // *  g'(x)
//           * ( dw[2] * dw[2] - dt);                // *( (dw)^2 - dt )
//   //exiting_current_voxel_y_b
//   exiting_current_voxel_z_a 
//      = drift[4] * dt                              // f(x) dt
//         + jump_rate[4] * dw[4]                    // + g(x) dw
//         + 0.5 * jump_rate[4]                      // + 0.5 g(x) 
//           * ONESIXTH //* rate_scale_factor           // * g'(x)
//           * ( dw[4] * dw[4] - dt);                // *( (dw)^2 - dt )
//   //exiting_current_voxel_z_b
//   
//   double total_exiting_current_voxel;
//   total_exiting_current_voxel
//         = exiting_current_voxel_x_a //+ exiting_current_voxel_x_b 
//            + exiting_current_voxel_y_a //+ exiting_current_voxel_y_b 
//            + exiting_current_voxel_z_a; //+ exiting_current_voxel_z_b;
//
//   // if total amount leaving is greater than is present, rescale them all
//   if (total_exiting_current_voxel >local_field[idx])//+local_change[idx])
//   {
//      double rescaling_factor; 
//      rescaling_factor = local_field[idx]/total_exiting_current_voxel;
//
//      exiting_current_voxel_x_a 
//         = exiting_current_voxel_x_a * rescaling_factor;
//      //exiting_current_voxel_x_b
//      //   = exiting_current_voxel_x_b * rescaling_factor;
//      exiting_current_voxel_y_a 
//         = exiting_current_voxel_y_a * rescaling_factor;
//      //exiting_current_voxel_y_b
//      //   = exiting_current_voxel_y_b * rescaling_factor;
//      exiting_current_voxel_z_a 
//         = exiting_current_voxel_z_a * rescaling_factor;
//      //exiting_current_voxel_z_b
//      //   = exiting_current_voxel_z_b * rescaling_factor;
//      
//      total_exiting_current_voxel = local_field[idx];
//   }
//
//   // If the neighbors have enough to supply the change, apply it,
//   //  otherwise restrict the change to the amount they currently have.
//   // WARNING: this is assymetrical and could lead to excess flow
//   //          along the negative direction along each axis
//   // TODO: maybe find a better solution
//   //  The preceding block rescales all outward flux of a voxel,
//   //   but in this inward flux case we don't know what other fluxes
//   //   the neighbor has so we can't rescale them.
//   if ( local_field[neigh_idx_x_a] + exiting_current_voxel_x_a 
//         + local_change[neigh_idx_x_a] < 0)
//   {
//      exiting_current_voxel_x_a 
//        = -1.0*(local_change[neigh_idx_x_a] + local_field[neigh_idx_x_a]);
//   }
//   if ( local_field[neigh_idx_y_a] + exiting_current_voxel_y_a 
//         + local_change[neigh_idx_y_a] < 0)
//   {
//      exiting_current_voxel_y_a 
//        = -1.0*(local_change[neigh_idx_y_a] + local_field[neigh_idx_y_a]);
//   }
//   if ( local_field[neigh_idx_z_a] + exiting_current_voxel_z_a 
//         + local_change[neigh_idx_z_a] < 0)
//   {
//      exiting_current_voxel_z_a 
//        = -1.0*(local_change[neigh_idx_z_a] + local_field[neigh_idx_z_a]);
//   }
//
//   local_change[neigh_idx_x_a] += exiting_current_voxel_x_a;
//   local_change[neigh_idx_y_a] += exiting_current_voxel_y_a;
//   local_change[neigh_idx_z_a] += exiting_current_voxel_z_a;
//
//
//   local_change[idx] -= total_exiting_current_voxel;
//   
//   /////////////////////////////////////////////////////////
//
//   return EXIT_SUCCESS;
//}


//int SPF_NS::conserved_gaussian_flux_separate_distributions(
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
//   )
//{
//   // assuming periodic boundary conditions
//   // WARNING: this is assymetrical and could lead to excess flow
//   //          along the negative direction along each axis
//   // TODO: maybe find a better solution
//
//   /////////////////////////////////////////////////////////
//   // evaluate local fluxes and changes
//   // using local_field[idx]
//   // and its neighbors indexed by neigh_x/y/z_idx[] 
//   // and accumulate field changes to 
//   //  local_change[idx]
//   //  and its neighbors local_change[neigh_x/y/z_idx[]]]
//
//   //double jump_magnitudes[6]; // assumes ndims == 3
//   double jump_mean[6]; // assumes ndims == 3
//   double jump_variance[6]; // assumes ndims == 3
//   double exiting_current_voxel_x_a; exiting_current_voxel_x_a = 0;
//   //double exiting_current_voxel_x_b; exiting_current_voxel_x_b = 0;
//   double exiting_current_voxel_y_a; exiting_current_voxel_y_a = 0;
//   //double exiting_current_voxel_y_b; exiting_current_voxel_y_b = 0;
//   double exiting_current_voxel_z_a; exiting_current_voxel_z_a = 0;
//   //double exiting_current_voxel_z_b; exiting_current_voxel_z_b = 0;
//
//   //for (size_t ii=0; ii < 6; ++ii) jump_magnitudes[ii] = 0;
//
//   // evaluate jump rates in each direction
//   // TODO: establish reasoning for these choices ... !
//
//   jump_mean[0] = //rate_scale_factor * 
//                  (local_field[idx] - local_field[neigh_idx_x_a]);
//   //jump_mean[1] = rate_scale_factor 
//   //               * (local_field[idx] - local_field[neigh_idx_x_b]);
//   jump_mean[2] = //rate_scale_factor * 
//                  (local_field[idx] - local_field[neigh_idx_y_a]);
//   //jump_mean[3] = rate_scale_factor 
//   //               * (local_field[idx] - local_field[neigh_idx_y_b]);
//   jump_mean[4] = //rate_scale_factor * 
//                  (local_field[idx] - local_field[neigh_idx_z_a]);
//   //jump_mean[5] = rate_scale_factor 
//   //               * (local_field[idx] - local_field[neigh_idx_z_b]);
//
//   jump_variance[0] = ONESIXTH //* rate_scale_factor 
//                        * local_field[idx];
//   //jump_variance[1] = ONESIXTH * rate_scale_factor * local_field[idx];
//   jump_variance[2] = ONESIXTH //* rate_scale_factor 
//                        * local_field[idx];
//   //jump_variance[3] = ONESIXTH * rate_scale_factor * local_field[idx];
//   jump_variance[4] = ONESIXTH //* rate_scale_factor 
//                        * local_field[idx];
//   //jump_variance[5] = ONESIXTH * rate_scale_factor * local_field[idx];
//
//   std::normal_distribution<double>
//         gd_x_a( jump_mean[0], jump_variance[0] );
//   exiting_current_voxel_x_a = sqrt(dt) * gd_x_a( rr.generator );
//   //std::normal_distribution<double>
//   //      gd_x_b( jump_mean[1], jump_variance[1] );
//   //exiting_current_voxel_x_b = sqrt(dt) * gd_x_b( rr.generator );
//   std::normal_distribution<double>
//         gd_y_a( jump_mean[2], jump_variance[2] );
//   exiting_current_voxel_y_a = sqrt(dt) * gd_y_a( rr.generator );
//   //std::normal_distribution<double>
//   //      gd_y_b( jump_mean[3], jump_variance[3] );
//   //exiting_current_voxel_y_b = sqrt(dt) * gd_y_b( rr.generator );
//
//   std::normal_distribution<double>
//         gd_z_a( jump_mean[4], jump_variance[4] );
//   exiting_current_voxel_z_a = sqrt(dt) * gd_z_a( rr.generator );
//   //std::normal_distribution<double>
//   //      gd_z_b( jump_mean[5], jump_variance[5] );
//   //exiting_current_voxel_z_b = sqrt(dt) * gd_z_b( rr.generator );
//   
//   double total_exiting_current_voxel;
//   total_exiting_current_voxel
//         = exiting_current_voxel_x_a //+ exiting_current_voxel_x_b 
//            + exiting_current_voxel_y_a //+ exiting_current_voxel_y_b 
//            + exiting_current_voxel_z_a; //+ exiting_current_voxel_z_b;
//
//   // if total amount leaving is greater than is present, rescale them all
//   if (total_exiting_current_voxel >local_field[idx])//+local_change[idx])
//   {
//      double rescaling_factor; 
//      rescaling_factor = local_field[idx]/total_exiting_current_voxel;
//
//      exiting_current_voxel_x_a 
//         = exiting_current_voxel_x_a * rescaling_factor;
//      //exiting_current_voxel_x_b
//      //   = exiting_current_voxel_x_b * rescaling_factor;
//      exiting_current_voxel_y_a 
//         = exiting_current_voxel_y_a * rescaling_factor;
//      //exiting_current_voxel_y_b
//      //   = exiting_current_voxel_y_b * rescaling_factor;
//      exiting_current_voxel_z_a 
//         = exiting_current_voxel_z_a * rescaling_factor;
//      //exiting_current_voxel_z_b
//      //   = exiting_current_voxel_z_b * rescaling_factor;
//      
//      total_exiting_current_voxel = local_field[idx];
//   }
//
//   // If the neighbors have enough to supply the change, apply it,
//   //  otherwise restrict the change to the amount they currently have.
//   // WARNING: this is assymetrical and could lead to erroneous flow
//   //          along the negative direction along each axis
//   // TODO: maybe find a better solution
//   //  The preceding block rescales all outward flux of a voxel,
//   //   but in this inward flux case we don't know what other fluxes
//   //   the neighbor has so we can't rescale them.
//   if ( local_field[neigh_idx_x_a] + exiting_current_voxel_x_a 
//         + local_change[neigh_idx_x_a] < 0)
//   {
//      exiting_current_voxel_x_a 
//        = -1.0*(local_change[neigh_idx_x_a] + local_field[neigh_idx_x_a]);
//   }
//   if ( local_field[neigh_idx_y_a] + exiting_current_voxel_y_a 
//         + local_change[neigh_idx_y_a] < 0)
//   {
//      exiting_current_voxel_y_a 
//        = -1.0*(local_change[neigh_idx_y_a] + local_field[neigh_idx_y_a]);
//   }
//   if ( local_field[neigh_idx_z_a] + exiting_current_voxel_z_a 
//         + local_change[neigh_idx_z_a] < 0)
//   {
//      exiting_current_voxel_z_a 
//        = -1.0*(local_change[neigh_idx_z_a] + local_field[neigh_idx_z_a]);
//   }
//
//   local_change[neigh_idx_x_a] += exiting_current_voxel_x_a;
//   local_change[neigh_idx_y_a] += exiting_current_voxel_y_a;
//   local_change[neigh_idx_z_a] += exiting_current_voxel_z_a;
//
//
//   local_change[idx] -= total_exiting_current_voxel;
//   
//   /////////////////////////////////////////////////////////
//
//   return EXIT_SUCCESS;
//}


//int SPF_NS::conserved_jump_flux_single_distribution( 
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
//   )
//{
//   // assuming periodic boundary conditions
//
//   /////////////////////////////////////////////////////////
//   // evaluate local fluxes and changes
//   // using local_field[idx]
//   // and its neighbors indexed by neigh_x/y/z_idx[] 
//   // and accumulate field changes to 
//   //  local_change[idx]
//   //  and its neighbors local_change[neigh_x/y/z_idx[]]]
//
//   //double jump_magnitude_x_a; // jump_magnitude_x = 0.0;
//   //double jump_magnitude_x_b; // jump_magnitude_x = 0.0;
//   //double jump_magnitude_y_a; // jump_magnitude_y = 0.0;
//   //double jump_magnitude_y_b; // jump_magnitude_y = 0.0;
//   //double jump_magnitude_z_a; // jump_magnitude_z = 0.0;
//   //double jump_magnitude_z_b; // jump_magnitude_z = 0.0;
//
//   double jump_magnitudes[6]; // assumes ndims == 3
//   for (size_t ii=0; ii < 6; ++ii) jump_magnitudes[ii] = 0;
//
//   //double change_to_current_voxel; change_to_current_voxel = 0;
//   int exiting_current_voxel; exiting_current_voxel = 0;
//
//   size_t dest_idx;
//
//   // evaluate jump rates in each direction
//   // TODO: establish reasoning for these choices ... !
//   double jump_rate; 
//   //jump_rate = rate_scale_factor * local_field[idx];
//   //rr.update_poisson_rate( jump_rate, dt );
//
//   //if ( rr.get_subinterval() != jump_rate * dt * exp( -1.0*jump_rate*dt) )
//   //   std::cout << "warning: subinterval not updated : "
//   //      << rr.get_subinterval() << " != " <<  jump_rate * dt * exp( -1.0*jump_rate*dt)  << std::endl; // debug
//   // evaluate jump magnitude in each direction
//
//   //jump_magnitude_x_a = rr.poisson_event_count( rr.generator );
//   //jump_rate = rate_scale_factor * (local_field[idx] - local_change[idx]);
//   
//   jump_rate = //rate_scale_factor * 
//                  (local_field[idx]);// + change_to_current_voxel);
//   //std::cout << "jump rate : " << jump_rate << std::endl; //debug
//   if ( jump_rate > 0.0 ) 
//   {
//      std::poisson_distribution<int> pd( jump_rate * dt);
//      exiting_current_voxel = pd( rr.generator );
//      //std::cout << " exiting_current_voxel : " 
//      //   << exiting_current_voxel 
//      //   << std::endl;
//      //std::cout << " random integer: " << pd(rr.generator )
//      //   << std::endl;
//      //local_change[idx] -= jump_magnitude_x_a;
//
//      //std::list<int> dest_idx_randomized;
//      if ( exiting_current_voxel >= local_field[idx] )
//      {
//         exiting_current_voxel = local_field[idx];
//      }
//
//      if ( exiting_current_voxel > 0) 
//      {
//         for (int ii=0; ii < exiting_current_voxel; ++ii)
//         {
//            // choose which neighbor to send this walker to
//            std::uniform_int_distribution<int> ud(0,5);
//            dest_idx = ud( rr.generator );
//            jump_magnitudes[dest_idx] += 1.0;
//         }
//      }
//      else 
//         if ( exiting_current_voxel < 0 )
//         {
//            std::cout  << "somehow 'exiting_current_voxel' < 0" 
//               << std::endl;
//         }
//
//      local_change[neigh_idx_x_a] += jump_magnitudes[0];
//      local_change[neigh_idx_x_b] += jump_magnitudes[1];
//      local_change[neigh_idx_y_a] += jump_magnitudes[2];
//      local_change[neigh_idx_y_b] += jump_magnitudes[3];
//      local_change[neigh_idx_z_a] += jump_magnitudes[4];
//      local_change[neigh_idx_z_b] += jump_magnitudes[5];
//   }
//   local_change[idx] -= exiting_current_voxel;
//   
//   //jump_rate = rate_scale_factor 
//   //               * (local_field[idx] + change_to_current_voxel);
//   //if ( jump_rate > 0.0 ) 
//   //{
//   //   //std::cout << "jump rate > 0" << std::endl; //debug
//   //   std::poisson_distribution<int> pd1( jump_rate * dt);
//   //   jump_magnitude_x_b = pd1( rr.generator );
//   //   //local_change[idx] -= jump_magnitude_x_a;
//   //   change_to_current_voxel -= jump_magnitude_x_a;
//   //   local_change[neigh_idx_x_a] += jump_magnitude_x_a;
//   //}
//
//   //jump_rate = rate_scale_factor 
//   //               * (local_field[idx] + change_to_current_voxel);
//   //if ( jump_rate > 0.0 ) 
//   //{
//   //   //std::cout << "jump rate > 0" << std::endl; //debug
//   //   std::poisson_distribution<int> pd2( jump_rate * dt);
//   //   jump_magnitude_x_b = pd2( rr.generator );
//   //   //local_change[idx] -= jump_magnitude_x_b;
//   //   change_to_current_voxel -= jump_magnitude_x_b;
//   //   local_change[neigh_idx_x_b] += jump_magnitude_x_b;
//   //}
//
//   //jump_rate = rate_scale_factor 
//   //               * (local_field[idx] + change_to_current_voxel);
//   //if ( jump_rate > 0.0 ) 
//   //{
//   //   //std::cout << "jump rate > 0" << std::endl; //debug
//   //   std::poisson_distribution<int> pd3( jump_rate * dt);
//   //   jump_magnitude_y_a = pd3( rr.generator );
//   //   //local_change[idx] -= jump_magnitude_y_a;
//   //   change_to_current_voxel -= jump_magnitude_y_a;
//   //   local_change[neigh_idx_y_a] += jump_magnitude_y_a;
//   //}
//
//   //jump_rate = rate_scale_factor 
//   //               * (local_field[idx] + change_to_current_voxel);
//   //if ( jump_rate > 0.0 ) 
//   //{
//   //   //std::cout << "jump rate > 0" << std::endl; //debug
//   //   std::poisson_distribution<int> pd4( jump_rate * dt);
//   //   jump_magnitude_y_b = pd4( rr.generator );
//   //   //local_change[idx] -= jump_magnitude_y_b;
//   //   change_to_current_voxel -= jump_magnitude_y_b;
//   //   local_change[neigh_idx_y_b] += jump_magnitude_y_b;
//   //}
//
//   //jump_rate = rate_scale_factor 
//   //               * (local_field[idx] + change_to_current_voxel);
//   //if ( jump_rate > 0.0 ) 
//   //{
//   //   //std::cout << "jump rate > 0" << std::endl; //debug
//   //   std::poisson_distribution<int> pd5( jump_rate * dt);
//   //   jump_magnitude_z_a = pd5( rr.generator );
//   //   //local_change[idx] -= jump_magnitude_z_a;
//   //   change_to_current_voxel -= jump_magnitude_z_a;
//   //   local_change[neigh_idx_z_a] += jump_magnitude_z_a;
//   //}
//
//   //jump_rate = rate_scale_factor 
//   //               * (local_field[idx] + change_to_current_voxel);
//   //if ( jump_rate > 0.0 ) 
//   //{
//   //   //std::cout << "jump rate > 0" << std::endl; //debug
//   //   std::poisson_distribution<int> pd6( jump_rate * dt);
//   //   jump_magnitude_z_b = pd6( rr.generator );
//   //   //local_change[idx] -= jump_magnitude_z_b;
//   //   change_to_current_voxel -= jump_magnitude_z_b;
//   //   local_change[neigh_idx_z_b] += jump_magnitude_z_b;
//   //}
//
//   //jump_magnitude_x_b = rr.poisson_event_count( rr.generator );
//   //jump_magnitude_y_a = rr.poisson_event_count( rr.generator );
//   //jump_magnitude_y_b = rr.poisson_event_count( rr.generator );
//   //jump_magnitude_z_a = rr.poisson_event_count( rr.generator );
//   //jump_magnitude_z_b = rr.poisson_event_count( rr.generator );
//
//   // TODO: check that the following won't cause a jump to below 0
//   //local_change[idx] -= jump_magnitude_x_a;
//                        // + jump_magnitude_x_b 
//                        //+ jump_magnitude_y_a + jump_magnitude_y_b
//                        //+ jump_magnitude_z_a + jump_magnitude_z_b;
//
//   //if ( local_field[idx] + change_to_current_voxel < 0 )
//   //{
//   //   std::cout << "warning: jump quantity exiting voxel " << idx 
//   //      << " is larger than the quantity present." << std::endl;
//
//   //   local_change[idx] = -1.0 * local_field[idx];
//   //}
//   //else local_change[idx] += change_to_current_voxel;
//
//   //local_change[neigh_idx_x_a] += jump_magnitude_x_a;
//   //local_change[neigh_idx_x_b] += jump_magnitude_x_b;
//   //local_change[neigh_idx_y_a] += jump_magnitude_y_a;
//   //local_change[neigh_idx_y_b] += jump_magnitude_y_b;
//   //local_change[neigh_idx_z_a] += jump_magnitude_z_a;
//   //local_change[neigh_idx_z_b] += jump_magnitude_z_b;
//
//   /////////////////////////////////////////////////////////
//
//   return EXIT_SUCCESS;
//}

int SPF_NS::conserved_jump_flux_separate_distributions( 
            std::vector<double>& pairwise_flux, // 6* local_field size
            SPF_NS::random& rr,
            const std::vector<double>& jump_rates, // 6* local_field size
            const double& dt,
            const size_t& Nvoxel_neighbors,
            const size_t& idx
            )
{
   // assuming periodic boundary conditions
   // jump_rates must be evaluated before calling 

   /////////////////////////////////////////////////////////
   // Evaluate local fluxes and changes due to jump processes
   //  using jump_rates. Enforce boundary conditions elsewhere.

   if ( pairwise_flux.size() // prevent assignments out of bounds
         >= ((Nvoxel_neighbors -1) + Nvoxel_neighbors*idx))
   {
      for( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
      {
         if( jump_rates[nn + Nvoxel_neighbors * idx] > 0.0)
         {
            std::poisson_distribution<int> 
               pd( dt * jump_rates[nn + Nvoxel_neighbors*idx]);

            pairwise_flux[nn + Nvoxel_neighbors*idx] = pd( rr.generator);
         }
         else
            pairwise_flux[nn + Nvoxel_neighbors*idx ] = 0.0;
      }
   }
   else
   {
      std::cout << "Error: pairwise_flux.size() not large enough"
         << std::endl;
   }

   return EXIT_SUCCESS;
}

int SPF_NS::conserved_gaussian_flux_separate_distributions( 
   std::vector<double>& pairwise_flux, // 6* local_field size
   SPF_NS::random& rr,
   const std::vector<double>& jump_rates_sqrt,  // 6 elements
   const std::vector<double>& jump_rate_sqrt_derivatives,  // 6 elements
   const double& dt,
   const size_t& Nvoxel_neighbors,
   const size_t& idx
   )
{
   // assuming periodic boundary conditions
   // jump_rates must be evaluated before calling 
   double sqrtdt; sqrtdt = sqrt(dt);

   /////////////////////////////////////////////////////////
   // evaluate local fluxes and changes due to gaussian processes
   //  using jump_rates. Enforce boundary conditions elsewhere.

   std::normal_distribution<double> gd( 0.0, 1.0);

   if ( pairwise_flux.size() // prevent assignments out of bounds
         >= ((Nvoxel_neighbors -1) + Nvoxel_neighbors*idx))
   {
      for( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
      {
         if( jump_rates_sqrt[ nn + Nvoxel_neighbors * idx] > 0.0)
         {
            // Ito
            //jump_magnitudes[*neigh_num_itr]
            //  = jump_rates[*neigh_num_itr] * dt
            //   + jump_rates[*neigh_num_itr] * sqrtdt * gd(rr.generator);

            // Stratonovich   (a(x) = mean, b(x) = variance of B(t))
            // current_state  + a(x)*dt + 0.5* b'(x)*b(x)*dt + b(x)dB(s,x)
            // b(x) = (1/6)*rate_scale_factor*x
            pairwise_flux[
                  nn + Nvoxel_neighbors * idx
               ] // stratonovich
               = ( //drift is sigma^2
                     jump_rates_sqrt[nn + Nvoxel_neighbors * idx] 
                     * jump_rates_sqrt[nn + Nvoxel_neighbors * idx]) *dt
                 // W-Z correction
                 + (0.5 
                     * jump_rate_sqrt_derivatives[
                           nn + Nvoxel_neighbors * idx] 
                     * jump_rates_sqrt[nn + Nvoxel_neighbors * idx] * dt)
                 // noise contribution
                 + (jump_rates_sqrt[nn + Nvoxel_neighbors * idx] 
                        * sqrtdt * gd(rr.generator)); // Ito of dB

            //if ( jump_magnitudes[*neigh_num_itr] < 0.0)
            //{
            //      jump_magnitudes[*neigh_num_itr] = 0.0;
            //}
            //if( (exiting_current_voxel + jump_magnitudes[*neigh_num_itr])
            //      > local_field[idx] )
            //{
            //   jump_magnitudes[*neigh_num_itr] =
            //      local_field[idx] - exiting_current_voxel ;
            //   if ( jump_magnitudes[*neigh_num_itr] < 0.0 ) 
            //   {
            //      jump_magnitudes[*neigh_num_itr] = 0.0;
            //      exiting_current_voxel = local_field[idx];
            //   }
            //}
            //if ( jump_magnitudes[*neigh_num_itr] + 
            //      local_change[neigh_idxs[*neigh_num_itr]] > 1.0)
            //{
            //   jump_magnitudes[*neigh_num_itr] 
            //      = 1.0 - local_change[neigh_idxs[*neigh_num_itr]];
            //}
            //exiting_current_voxel += jump_magnitudes[*neigh_num_itr];
            //local_change[neigh_idxs[*neigh_num_itr]] 
            //               += jump_magnitudes[*neigh_num_itr];
         }
         else
         {
            pairwise_flux[nn + Nvoxel_neighbors * idx] = 0;
         }
      } // end loop over neighbors
   }
   else
   {
      std::cout << "Error: pairwise_flux.size() not large enough"
         << std::endl;
   }

   return EXIT_SUCCESS;
}

int SPF_NS::conserved_gaussian_flux_separate_distributions_ito( 
   std::vector<double>& local_change, // must be same size as local_field
   const std::vector<double>& local_field,
   SPF_NS::random& rr,
   const std::vector<double>& jump_rates,  // 6 elements
   const std::vector<double>& jump_rate_derivatives,  // 6 elements
   const double& dt,
   const size_t& idx,
   const std::vector<size_t>& neigh_idxs,  // 6 elements
   const std::vector<size_t>& neigh_order,  // 6 elements
   const int& Ny,
   const int& Nz
   )
{
   // assuming periodic boundary conditions
   // jump_rates must be evaluated before calling 
   double sqrtdt; sqrtdt = sqrt(dt);

   /////////////////////////////////////////////////////////
   // evaluate local fluxes and changes due to jump processes
   // using local_field[idx]
   // and its neighbors indexed by neigh_x/y/z_idx[] 
   // and accumulate field changes to 
   //  local_change[idx]
   //  and its neighbors local_change[neigh_x/y/z_idx[]]]


   double jump_magnitudes[6]; // assumes ndims == 3
   for (size_t ii=0; ii < 6; ++ii) jump_magnitudes[ii] = 0;

   //double change_to_current_voxel; change_to_current_voxel = 0;
   double exiting_current_voxel; exiting_current_voxel = 0;

   size_t dest_idx;

   // evaluate jump rates in each direction
   //for( size_t ii=0; ii < 6; ++ii)
   //{
   //   jump_rate[ii] = rate_scale_factor   // TODO: change to be specific to physical phenomenon; pass rates into this function
   //                  * (local_field[idx]);// + change_to_current_voxel);
   //}
   
   // randomize the order in which walkers are distributed to neighbors
   //std::vector<size_t> neigh_num_scrambled(6,0);
   //std::vector<double> rand_decimals1(6,0);
   //std::vector<double> rand_decimals2(6,0);

   //bool randomized_flag; randomized_flag = true;
   //std::uniform_real_distribution<double> rand_decimal(0,1);
   //do {
   //   for( size_t ii=0; ii < 6; ++ii) 
   //   {
   //      rand_decimals1[ii] = rand_decimal( rr.generator );
   //      rand_decimals2[ii] = rand_decimals1[ii];
   //   }
   //   std::sort( rand_decimals1.begin(), rand_decimals1.end());

   //   // check that no number was drawn twice by rand_decimal()
   //   std::vector<double>::iterator decimals_itr;
   //   decimals_itr = std::unique(rand_decimals1.begin(), 
   //                              rand_decimals1.end() );
   //   rand_decimals1.resize( 
   //            std::distance( rand_decimals1.begin(), decimals_itr) );
   //   if ( rand_decimals1.size() != 6 )
   //   {
   //      randomized_flag = false;
   //   }
   //}
   //while ( ! randomized_flag );

   //for( size_t ii=0; ii < 6; ++ii)
   //{
   //   for( size_t jj=0; jj < 6; ++jj)
   //   {
   //         if ( rand_decimals1[jj] == rand_decimals2[ii]) 
   //         {
   //            neigh_num_scrambled[jj] = ii;
   //            continue;
   //         }
   //   }
   //}

   // debug
   //std::cout << "neigh_idx_scramble[] : " ;
   //for( size_t ii=0; ii < 6; ++ii) 
   //{
   //   std::cout << neigh_num_scrambled[ii] << ", ";
   //}
   //std::cout << std::endl;
   // end debug

   std::normal_distribution<double> gd( 0.0, 1.0);
   for( std::vector<size_t>::const_iterator 
         neigh_num_itr = neigh_order.begin(); 
         neigh_num_itr != neigh_order.end(); 
         ++neigh_num_itr)
         //size_t ii=0; ii < neigh_num_scrambled.size(); ++ii)
   {
      if( jump_rates[ *neigh_num_itr ] > 0.0)
      {

         // Ito
         jump_magnitudes[*neigh_num_itr]
            = jump_rates[*neigh_num_itr] * dt
               + jump_rates[*neigh_num_itr] * sqrtdt * gd(rr.generator);

         // Stratonovich   (a(x) = mean, b(x) = variance of B(t))
         //    current_state  + a(x)*dt + 0.5* b'(x)*b(x)*dt + b(x)dB(s,x)
         // b(x) = (1/6)*rate_scale_factor*x
         //jump_magnitudes[*neigh_num_itr] // stratonovich
         //   = jump_rates[*neigh_num_itr] * dt // drift
         //      + (0.5 *jump_rate_derivatives[*neigh_num_itr] 
         //         *jump_rates[*neigh_num_itr] * dt)   // W-Z correction
         //      + (jump_rates[*neigh_num_itr] 
         //            * sqrtdt * gd(rr.generator)); // Ito of dB

         if ( jump_magnitudes[*neigh_num_itr] < 0.0)
         {
               jump_magnitudes[*neigh_num_itr] = 0.0;
         }
         if( (exiting_current_voxel + jump_magnitudes[*neigh_num_itr])
               > local_field[idx] )
         {
            jump_magnitudes[*neigh_num_itr] =
               local_field[idx] - exiting_current_voxel ;
            if ( jump_magnitudes[*neigh_num_itr] < 0.0 ) 
            {
               jump_magnitudes[*neigh_num_itr] = 0.0;
               exiting_current_voxel = local_field[idx];
            }
         }
         if ( jump_magnitudes[*neigh_num_itr] + 
               local_change[neigh_idxs[*neigh_num_itr]] > 1.0)
         {
            jump_magnitudes[*neigh_num_itr] 
               = 1.0 - local_change[neigh_idxs[*neigh_num_itr]];
         }
         exiting_current_voxel += jump_magnitudes[*neigh_num_itr];
         local_change[neigh_idxs[*neigh_num_itr]] 
                        += jump_magnitudes[*neigh_num_itr];
      }
      //if( jump_rates[ *neigh_num_itr ] <= 0.0)
   }
   local_change[idx] -= exiting_current_voxel;

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

int SPF_NS::identify_local_neighbors(  // TODO: make 2-D compatable
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
   //const int& Nx_local,
   const int& Ny,
   const int& Nz
   )
{
   // periodic boundary along x is enforced during ghost assignments
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
