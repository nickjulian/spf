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

int SPF_NS::enforce_bounds_int_outward(
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
   //std::uniform_real_distribution<double> rand_decimal(0,1);// for order
   std::uniform_int_distribution<int> ud(0, Nvoxel_neighbors -1);
   //std::vector<double> rand_decimals1(Nvoxel_neighbors, 0);
   //std::vector<double> rand_decimals2(Nvoxel_neighbors, 0);
   //std::vector<double> flux_reduction_factors(6,1);
   double flux_reduction_factor; flux_reduction_factor =1;
   double total_flux_rate;
   std::vector<size_t> neigh_pairs(Nvoxel_neighbors, 0);
   neigh_pairs[0] = 1;  // x upward
   neigh_pairs[1] = 0;  // x downward
   size_t idx; idx = 0;
   int rounding_error, dest_idx; rounding_error = 0; dest_idx = 0;

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

   bool dest_flag; dest_flag = false;
   double current_flux; current_flux = 0;

   // iterate over local voxels
   
   // Ensure outward flux isn't too high
   for (size_t ii=1; ii < Nx_local +1; ++ii){
      size_t ioffset = Ny * ii;                  // loop over non-ghosts
      for ( size_t jj=0; jj < Ny; ++jj){
         size_t joffset = Nz * (jj + ioffset);
         for ( size_t kk=0; kk < Nz; ++kk)
         {
            idx = kk + joffset;
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
               //total_flux_rate = 0;
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
                  
                  phi_local_flux[nn + Nvoxel_neighbors*idx]
                     = round(phi_local_flux[nn + Nvoxel_neighbors*idx]);
               }
               // check to see if any walkers were lost due to rounding
               rounding_error = 0;
               for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
               {
                  rounding_error 
                     += phi_local_flux[nn + Nvoxel_neighbors*idx];
               }
               rounding_error = int(phi_local[idx] 
                                          - phi_lower_limit 
                                          - rounding_error );
               while ( rounding_error < 0 )
               {  // outward_flux is too large
                  // find a neighbor to reduce flux of
                  dest_flag = false;
                  dest_idx = ud(rr.generator);// random initial neigh
                  if (phi_local_flux[dest_idx + Nvoxel_neighbors*idx] >0)
                  {
                     dest_flag = true;
                  }
                  //while ( 
                  //   phi_local_flux[dest_idx + Nvoxel_neighbors*idx] <=0)
                  //{
                  //   dest_idx = ud(rr.generator);
                  //   dest_flag = true;
                  //}

                  // randomize the order of index choice
                  std::vector<int>
                     neighIdxsToSample( Nvoxel_neighbors, 0);
                  for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
                  {
                     neighIdxsToSample[nn] = (int) nn;
                  }

                  std::vector<int>
                     randomNeighIdxs( Nvoxel_neighbors, 0);
                  for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
                  {
                     std::uniform_int_distribution<int>
                        uid( 0, neighIdxsToSample.size() -1);
                     int chosenIdx = uid( rr.generator);

                     randomNeighIdxs[nn]
                        = neighIdxsToSample[ chosenIdx];

                     // erase the chosen idx from neighIdxsToSample
                     neighIdxsToSample.erase(
                           neighIdxsToSample.begin() + chosenIdx
                           );
                  }

                  for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
                  {  // choose neigh with the lowest outward flux rate
                     if((phi_local_flux[ 
                              randomNeighIdxs[ nn] 
                              + Nvoxel_neighbors*idx] > 0)
                          &&
                        ((phi_local_rates[dest_idx + Nvoxel_neighbors*idx]
                           > 
                          phi_local_rates[
                           randomNeighIdxs[ nn] + Nvoxel_neighbors*idx]
                         )
                        || (dest_flag == false)
                        )) // this will prioritize lower index neighbors...
                     {
                        dest_idx = randomNeighIdxs[nn];
                        dest_flag = true;
                     }
                  }

                  if ( dest_flag )
                  {
                     phi_local_flux[dest_idx + Nvoxel_neighbors*idx] -= 1;
                     rounding_error += 1;
                  }
                  else
                  {
                     break;
                  }
               } // while ( rounding_error < 0 )
               while ( rounding_error > 0 )
               {  // outward_flux too small
                  // find the neighbor flux having highest rate 
                  //   and increment its flux 
                  dest_flag = false;
                  dest_idx = ud(rr.generator);// initially random neigh
                  if ( phi_local[ neigh_idxs[dest_idx]]
                           < phi_upper_limit)
                  {
                     dest_flag = true;
                  }

                  // randomize the choice of neighbor indices
                  std::vector<int> 
                     neighIdxsToSample( Nvoxel_neighbors, 0);
                  for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
                  {
                     neighIdxsToSample[nn] = (int) nn;
                  }

                  std::vector<int> 
                     randomNeighIdxs( Nvoxel_neighbors, 0);
                  for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
                  {
                     std::uniform_int_distribution<int> 
                        uid( 0, neighIdxsToSample.size() -1);
                     int chosenIdx = uid( rr.generator);

                     randomNeighIdxs[nn] 
                        = neighIdxsToSample[ chosenIdx];

                     // erase the chosen idx from neighIdxsToSample
                     neighIdxsToSample.erase(
                           neighIdxsToSample.begin() + chosenIdx
                           );
                  }

                  for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
                  {  // choose neigh with the greatest outward flux rate
                     if ((
                         (phi_local_rates[dest_idx + Nvoxel_neighbors*idx]
                           < 
                          phi_local_rates[
                           randomNeighIdxs[ nn] + Nvoxel_neighbors*idx]
                         )
                              || (dest_flag == false)
                         ) && (
                           // ensure neighbor isn't full
                           phi_local[ neigh_idxs[
                              randomNeighIdxs[nn]]]
                              < phi_upper_limit
                         )
                        )
                     {
                           dest_idx = randomNeighIdxs[nn];
                           dest_flag = true;
                     }
                  }
                  if ( dest_flag )
                  { // only fix the rounding error if a flux can be
                    //  modified without overfilling a neighbor
                     phi_local_flux[ 
                        dest_idx + Nvoxel_neighbors*idx
                                    ] += 1;
                     rounding_error -= 1;
                  }
                  else
                  {
                     break;
                  }
               } // while ( rounding_error > 0)
               //if( rounding_error == 1)
               //{
               //   //std::uniform_int_distribution<int> 
               //   //   ud(0, Nvoxel_neighbors -1);
               //   dest_idx = ud(rr.generator); 
               //   for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
               //   {  // find the greatest outward flux rate
               //      if (phi_local_rates[dest_idx + Nvoxel_neighbors*idx]
               //            < phi_local_rates[nn + Nvoxel_neighbors*idx])
               //         dest_idx = nn;
               //   }
               //   // add the lost walker to the chosen flux
               //   phi_local_flux[dest_idx + Nvoxel_neighbors*idx] += 1;
               //   rounding_error -= 1;
               //}
               // debug
               //if ( rounding_error != 0)
               //{
               //   std::cout 
               //      << "Somehow rounding_error was neither 0 or 1: "
               //      << rounding_error
               //      << " line 495, ";
               //   for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
               //   {  // find the greatest outward flux rate
               //      std::cout 
               //         << phi_local_flux[ nn + Nvoxel_neighbors*idx]
               //                  << ", ";
               //   }
               //   std::cout << std::endl;
               //}
               // end debug
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

         }// outward flux boundary check
      }
   }
   return EXIT_SUCCESS;
}

int SPF_NS::enforce_bounds_int_inward(
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
   //std::uniform_real_distribution<double> rand_decimal(0,1);// for order
   std::uniform_int_distribution<int> ud(0, Nvoxel_neighbors -1);
   //std::vector<double> rand_decimals1(Nvoxel_neighbors, 0);
   //std::vector<double> rand_decimals2(Nvoxel_neighbors, 0);
   //std::vector<double> flux_reduction_factors(6,1);
   double flux_reduction_factor; flux_reduction_factor =1;
   double total_flux_rate;
   std::vector<size_t> neigh_pairs(Nvoxel_neighbors, 0);
   neigh_pairs[0] = 1;  // x upward
   neigh_pairs[1] = 0;  // x downward
   size_t idx; idx = 0;
   int rounding_error, dest_idx; rounding_error = 0; dest_idx = 0;

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

   bool dest_flag; dest_flag = false;
   double current_flux; current_flux = 0;

   // iterate over local voxels
   
   // Ensure inward flux isn't too high
   for (size_t ii=1; ii < Nx_local +1; ++ii){
      size_t ioffset = Ny * ii;                 // loop over non-ghosts
      for ( size_t jj=0; jj < Ny; ++jj){
         size_t joffset = Nz * (jj + ioffset);
         for ( size_t kk=0; kk < Nz; ++kk)
         {
            idx = kk + joffset;
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
            ///////////////////////////////////////////////////////////

            // Iterate over neighbors of local voxels.
            // Ensure inward flux + local population doesn't exceed upper 
            //  limit.

            inward_flux = 0;
            for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
            {
               inward_flux += phi_local_flux[ 
                                    neigh_pairs[nn]   // mm 
                                    + Nvoxel_neighbors * neigh_idxs[nn]
                                 ];
            }
            if ( inward_flux + phi_local[idx] > phi_upper_limit )
            {  // renormalize the inward fluxes
               //total_flux_rate = 0;
               //total_flux = 0;
               //for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
               //{
               //   total_flux_rate 
               //      += phi_local_rates[
               //                  neigh_pairs[nn]   // mm 
               //                    + Nvoxel_neighbors * neigh_idxs[nn]];
               //}
               //if ( flags.debug != 0)
               //{
               //   std::cout << "renormalizing phi_local_flux " ;
               //}
               for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
               {
                  phi_local_flux[
                                 neigh_pairs[nn]   // mm 
                                    + Nvoxel_neighbors * neigh_idxs[nn]]
                     *= (phi_upper_limit - phi_local[idx])
                        //* phi_local_flux[
                        //         neigh_pairs[nn]   // mm 
                        //            + Nvoxel_neighbors * neigh_idxs[nn]]
                                    / inward_flux;
                                    // / total_flux
                     //(phi_upper_limit - phi_local[idx])
                     //   * phi_local_rates[
                     //            neigh_pairs[nn]   // mm 
                     //               + Nvoxel_neighbors * neigh_idxs[nn]
                     //      ] / total_flux_rate;
                  // round the flux to integers
                  //renormalized_flux[nn]

                  phi_local_flux[ neigh_pairs[nn]
                                    + Nvoxel_neighbors * neigh_idxs[nn]]
                     = round( phi_local_flux[ 
                                 neigh_pairs[nn]
                                    + Nvoxel_neighbors * neigh_idxs[nn]]);
                  //         renormalized_flux[nn]);
                  //if ( flags.debug != 0)
                  //{
                  //   std::cout << ", phi_local_flux["
                  //               << neigh_pairs[nn]
                  //                  + Nvoxel_neighbors * neigh_idxs[nn]
                  //                  << "] "
                  //      << phi_local_flux[
                  //               neigh_pairs[nn]
                  //                  + Nvoxel_neighbors * neigh_idxs[nn]];
                  //}
               }
               //if ( flags.debug != 0) std::cout << std::endl;
               // check rounding error
               inward_flux = 0;
               for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
               {
                  inward_flux += phi_local_flux[ 
                                       neigh_pairs[nn]   // mm 
                                       + Nvoxel_neighbors * neigh_idxs[nn]
                                    ];
               }
               rounding_error = int( phi_upper_limit - phi_local[idx] 
                                       - inward_flux);
               while ( rounding_error > 0 )
               {  // inward_flux is too small
                  // find a neighbor to pull walkers from
                  dest_flag = false;
                  dest_idx = ud(rr.generator); // initially random neigh
                  // ensure initial choice has walkers to spare
                  current_flux = 0;
                  for ( size_t mm=0; mm < Nvoxel_neighbors; ++mm)
                  {
                     current_flux += 
                        phi_local_flux[
                           mm + Nvoxel_neighbors* neigh_idxs[dest_idx]];
                  }
                  if (( current_flux <
                                 phi_local[ neigh_idxs[dest_idx]]
                              )
                        &&
                              // disallow taking from ghosts
                              //  because nodes aren't updating phi_local
                              //  ghosts from which they take from
                              ( neigh_idxs[dest_idx] >= Ny*Nz )
                              &&
                              (neigh_idxs[dest_idx] < (Nx_local+1)*Ny*Nz)
                              // TODO: resolve this assymmetry
                        )
                  {
                     dest_flag = true;
                  }

                  // ensure order of redistribution is random
                  std::vector<int> 
                     neighIdxsToSample( Nvoxel_neighbors, 0);
                  for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
                  {
                     neighIdxsToSample[nn] = (int) nn;
                  }

                  std::vector<int> 
                     randomNeighIdxs( Nvoxel_neighbors, 0);
                  for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
                  {
                     std::uniform_int_distribution<int> 
                        uid( 0, neighIdxsToSample.size() -1);
                     int chosenIdx = uid( rr.generator);

                     randomNeighIdxs[nn] 
                        = neighIdxsToSample[ chosenIdx];

                     // erase the chosen idx from neighIdxsToSample
                     neighIdxsToSample.erase(
                           neighIdxsToSample.begin() + chosenIdx
                           );
                  }

                  for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
                  {  // choose neigh with the greatest outward flux rate
                     if ((phi_local_rates[
                           neigh_pairs[dest_idx] 
                              + Nvoxel_neighbors*neigh_idxs[dest_idx]]
                           < 
                         phi_local_rates[
                              neigh_pairs[ randomNeighIdxs[nn]] 
                                 + Nvoxel_neighbors
                                    *neigh_idxs[ randomNeighIdxs[nn]]]
                        )
                           || (dest_flag == false)
                        )
                     {
                        // ensure voxel nn has walkers to spare
                        current_flux =0;
                        for ( size_t mm=0; mm < Nvoxel_neighbors; ++mm)
                        {
                           current_flux += 
                              phi_local_flux[
                                 mm + Nvoxel_neighbors 
                                    * neigh_idxs[ randomNeighIdxs[nn]]];
                        }
                        if (( current_flux
                                   <
                                 phi_local[ 
                                    neigh_idxs[ randomNeighIdxs[nn]]]
                              )&&
                              // disallow taking from ghosts
                              // TODO: resolve this assymmetry
                              ( neigh_idxs[
                                 randomNeighIdxs[nn]] >= Ny*Nz )
                              &&
                              ( neigh_idxs[
                                 randomNeighIdxs[ nn]]
                                < (Nx_local+1)*Ny*Nz)
                           )
                        {
                           dest_flag = true;
                           dest_idx = randomNeighIdxs[ nn];
                        }
                     }
                  }
                  // add the lost walker to the chosen flux
                  if ( dest_flag )
                  {
                     phi_local_flux[
                              neigh_pairs[dest_idx] 
                                 + Nvoxel_neighbors*neigh_idxs[dest_idx]
                           ] += 1;
                     rounding_error -= 1;
                  }
                  else
                  {
                     break;
                  }
               }
               while ( rounding_error < 0 )
               {  // inward_flux is too large
                  // find the neighbor flux having smallest rate 
                  //  and non-zero flux and decrement its flux 
                  //
                  dest_flag = false;
                  dest_idx = ud(rr.generator); // initially random neigh
                  if ( phi_local_flux[ // ensure the neighs flux isn't 0
                              neigh_pairs[dest_idx] 
                                 + Nvoxel_neighbors* neigh_idxs[dest_idx]]
                                    > 0)
                  {
                     dest_flag = true;
                  }

                  // ensure order of redistribution is random
                  std::vector<int> 
                     neighIdxsToSample( Nvoxel_neighbors, 0);
                  for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
                  {
                     neighIdxsToSample[nn] = (int) nn;
                  }

                  std::vector<int> 
                     randomNeighIdxs( Nvoxel_neighbors, 0);
                  for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
                  {
                     std::uniform_int_distribution<int> 
                        uid( 0, neighIdxsToSample.size() -1);
                     int chosenIdx = uid( rr.generator);

                     randomNeighIdxs[nn] 
                        = neighIdxsToSample[ chosenIdx];

                     // erase the chosen idx from neighIdxsToSample
                     neighIdxsToSample.erase(
                           neighIdxsToSample.begin() + chosenIdx
                           );
                  }
                  for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
                  {  // choose neigh with the lowest outward flux rate
                     if (((phi_local_rates[
                           neigh_pairs[dest_idx] 
                              + Nvoxel_neighbors*neigh_idxs[dest_idx]]
                           >
                         phi_local_rates[
                              neigh_pairs[ randomNeighIdxs[nn]] 
                                 + Nvoxel_neighbors
                                    *neigh_idxs[ randomNeighIdxs[nn]]]
                         )
                              || (dest_flag == false)
                         )
                           && 
                        ( phi_local_flux[
                                 neigh_pairs[ randomNeighIdxs[nn]] 
                                    + Nvoxel_neighbors
                                      * neigh_idxs[ randomNeighIdxs[nn]]]
                                    > 0)
                        )
                     {
                        dest_flag = true;
                        dest_idx = randomNeighIdxs[nn];
                     }
                  }
                  // remove the extra walker to the chosen flux
                  if ( dest_flag )
                  {
                     phi_local_flux[
                              neigh_pairs[dest_idx] 
                                 + Nvoxel_neighbors*neigh_idxs[dest_idx]
                           ] -= 1;
                     rounding_error += 1;
                  }
                  else
                  {
                     break;
                  }
               }
               //if ( rounding_error != 0)
               //{
               //   std::cout 
               //      << "Somehow rounding_error was neither 0 or 1: "
               //      << rounding_error
               //      << ", inward_flux " << inward_flux
               //      << ", phi_local[" << idx << "] " << phi_local[idx]
               //      << ", phi_upper_limit " << phi_upper_limit
               //      //<< ", total_flux_rate " << total_flux_rate;
               //      << ", phi_loca[" << dest_idx << "] " 
               //      << phi_local[dest_idx];
               //   for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
               //   {
               //      std::cout << ", phi_local_flux[" 
               //      << neigh_pairs[nn]+ Nvoxel_neighbors * neigh_idxs[nn]
               //      << "] "
               //      << phi_local_flux[ 
               //                     neigh_pairs[nn] 
               //                       + Nvoxel_neighbors * neigh_idxs[nn]]
               //         << std::endl;
               //   }
               //}
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
                       > 10*eps.dbl
                     )
                     std::cout << "Warning: reduced outward flux still greater than voxel contents. (outward_flux, phi_local[" << idx << "], phi_local[] - phi_lower_limit - outward_flux): (" << outward_flux << ", " << phi_local[idx] << ", " << phi_local[idx] - phi_lower_limit - outward_flux << ")" << std::endl;
                  // else phi_local[idx] =0 if <0 in parent function
               }

               if ( inward_flux + phi_local[idx] > phi_upper_limit )
               {
                  // TODO: the following is just a guess of error limit
                  if ( (inward_flux + phi_local[idx] - phi_upper_limit)
                       > 10*eps.dbl
                     )
                     std::cout << "Warning: reduced inward flux + previous population still greater than voxel upper limit. (inward_flux, phi_local[" << idx << "], " << "inward_flux + phi_local[] - phi_upper_limit" << "): (" << inward_flux << ", " << phi_local[idx] << ", " << inward_flux + phi_local[idx] - phi_upper_limit << ")"
                        << std::endl;
               }
            } // end debug
         }
      }
   }
   return EXIT_SUCCESS;
}

double SPF_NS::laplacian1d(
         const double& hh,
         const std::vector<double>& local_field,
         const size_t& idx,
         const size_t& neigh_idx_x_a,
         const size_t& neigh_idx_x_b
      )
{
   return (1.0/(hh*hh))*
            (
               local_field[neigh_idx_x_a]
               -2.0* local_field[idx]
               + local_field[neigh_idx_x_b]
            );
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

double SPF_NS::laplacian1dT(
         const double& hh,
         const std::vector<double>& T_field,
         const size_t& idx,
         const size_t& neigh_idx_x_a,
         const size_t& neigh_idx_x_b
      )
{
   return (1.0/(hh*hh))*
            (
               T_field[neigh_idx_x_a]
               -2.0* T_field[idx]
               + T_field[neigh_idx_x_b]
            );
}
#endif
