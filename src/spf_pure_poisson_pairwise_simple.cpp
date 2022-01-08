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
// File: spf_pure_poisson_pairwise_simple.cpp

#include <iostream>  // cout, cin, cerr, endl
#include <iomanip>   // setw, setprecision
#include <fstream>   // ifstream, ofstream
#include <cstdlib>   // EXIT_SUCCESS, EXIT_FAILURE
#include <vector>
#include <string>
#include <math.h> // floor, isnan
#include <time.h> // time_t, time, ctime

#include <mpi.h>
#include "../include/hdf5.h"

#include "check_for_failure.hpp"
#include "read_cmdline.hpp"
#include "readHDF5c.hpp"
#include "writeHDF5c.hpp"
#include "spf_communication.hpp"

#include "stochastic_rates.hpp"
#include "voxel_dynamics.hpp"

#include "flags.hpp" // int_flags flags;

#include "rand.hpp"
//#include "jumps.hpp"
//#include "rungekutta.hpp"
#include "macheps.hpp"  // determines the local machine's relative 
                        //  rounding error 

using std::cout;
using std::endl;
using std::cerr;
using std::setw;
using namespace std;
using namespace SPF_NS;

int main( int argc, char* argv[])
{
   // initialize MPI
   int rootnode; rootnode = 0; // arbitrary selection identifying rootnode
   int mynode, totalnodes;
   MPI_Comm world_comm, neighbors_comm; //, workers_comm; // communicators
   MPI_Status mpi_status;

   MPI_Init( &argc, &argv);
   world_comm = MPI_COMM_WORLD;
   MPI_Comm_size(world_comm, &totalnodes);//totalnodes= # of nodes in comm
   MPI_Comm_rank(world_comm, &mynode ); // mynode = rank of current node

   hid_t fa_plist_id, dx_plist_id;
   // file access property list
   fa_plist_id = H5Pcreate( H5P_FILE_ACCESS );
   H5Pset_fapl_mpio(fa_plist_id, world_comm, MPI_INFO_NULL);//info); 
   // dataset transfer property list
   dx_plist_id = H5Pcreate( H5P_DATASET_XFER );
   H5Pset_dxpl_mpio( dx_plist_id, H5FD_MPIO_COLLECTIVE);

   ////////////////////////////////////////////////////////////////////
   // In case it becomes useful to dedicate a group to a type of process:
   //MPI_Group world_group, worker_group;
   //int ranks[1];  // array of ranks
   //MPI_Comm_group( world_comm, &world_group); extract world group
   //int server; 
   //server = totalnodes - 1; // identify the node to act as a server
   //// form a new group by excluding nodes specified by ranks
   //ranks[0] = server; // ranks of nodes to exclude from a group
   //MPI_Group_excl( world_group, 1, ranks, &worker_group);
   //// create a new communicator (workers) from the group worker_group
   //MPI_Comm_create( world_comm, worker_group, &workers_comm); 
   //// release the communicator assigned to the group worker_group
   //MPI_Comm_free( &workers_comm); 
   //MPI_Group_free(&worker_group);// may be freed before or after comm
   ////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////
   // dimensional parameters
   std::vector<hsize_t> dims;

   // TODO: why are these array dimensions integers rather than size_t?
   int Nx_total, Nx_local, Ny, Nz, Nt;
   Nx_total = 1; Nx_local = 1; Ny = 1; Nz = 1;
   int time_step; time_step = 0;
   int write_period; write_period = 1;
   double time, dt; time = 0; dt = 1;

   // TODO: erase this and read Nt from the cmdline
   //double diffusivityT; diffusivityT = 3.0E-4;
   int Nv; Nv = 1; // number of walkers possible in a voxel

   ////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////
   // read input parameters
   std::string inputFileName;
   std::vector<string> args( argv, argv + argc );
   std::string output_prefix;
   int_flags flags;
   flags.fail = 0;

   // NOTE: there might not be a way to safely Bcast strings without  
   //       assuming they're ASCII, so read cmdline options on all nodes.
   if ( read_cmdline_options(
            args,
            dt,
            Nt,
            Nv,
            write_period,
            flags,
            output_prefix,
            inputFileName,
            mynode,
            rootnode,
            MPI_COMM_WORLD
            ) != EXIT_SUCCESS )
   {
      if ( mynode == rootnode ) PRINT_USAGE;

      //MPI_Comm_free( &neighbors_comm); 
      MPI_Finalize();
      return EXIT_FAILURE;
   }
   if ( (flags.debug != 0 ) && (mynode == rootnode ))
   {
      std::cout << "running with parameters: "
                << "  -o " << output_prefix << std::endl
                << "  -i " << inputFileName << std::endl
                << "  -Nt " << Nt << std::endl
                << "  -Nv " << Nv << std::endl
                << "  -dt " << dt << std::endl
                << "  -wp " << write_period << std::endl;
   }

   // establish field limits
   double phi_upper_limit, phi_lower_limit;
   phi_upper_limit = Nv;
   phi_lower_limit = 0.0;
   ////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////
   // read input data
   hid_t inFile_id, outFile_id, phi_dataset_id, phi_dataspace_id;

   // open the file containing the initial state
   inFile_id = H5Fopen( inputFileName.c_str(), 
                        H5F_ACC_RDONLY, fa_plist_id);

   if (inFile_id < 0) flags.fail = -1;
   if ( check_for_failure( flags, world_comm) )
   {
      MPI_Finalize();
      if ( mynode == rootnode )
      {
         cout << "Error, failed to open input file: " 
            << inputFileName << endl;
      }
      return EXIT_FAILURE;
   }

   if ( read_dims_from_hdf5(
         inFile_id,
         dims,
         flags,
         mynode,
         rootnode,
         totalnodes,
         world_comm
         ) == EXIT_FAILURE )
   {
      flags.fail = -1;
   }

   if ( totalnodes > dims[0] )
   {
      flags.fail = -1;
      if ( mynode == rootnode )
         cout << "Error, number of compute nodes "
            "> number of points on phi domain" << endl;
   }

   if ( check_for_failure( flags, world_comm) )
   {
      MPI_Finalize();
      if ( mynode == rootnode )
      {
         cout << "Error, failed to read dimensionality from input file: " 
            << inputFileName << endl;
      }
      return EXIT_FAILURE;
   }

   int ndims; 
   ndims = dims.size();

   // machine relative error
   epsilon eps;
   if ( flags.debug != 0)
   {
      std::cout << "node " << mynode 
         << " machine epsilon: " << eps.dbl << std::endl;
   }

   // indices of local (non-ghost) data for use in pre-split data array
   std::vector<size_t> idx_start(ndims,0);
   std::vector<size_t> idx_end(ndims,0);
   for (size_t ii=1; ii < ndims; ++ii) idx_end[ii] = dims[ii] -1;
   
   if ( determine_local_idxs(
                        dims,
                        mynode,
                        rootnode,
                        totalnodes,
                        Nx_local,
                        idx_start,
                        idx_end
                        ) == EXIT_FAILURE)
   {
      flags.fail = -1;
   }

   if ( check_for_failure( flags, world_comm) )
   {
      H5Fclose( inFile_id );
      MPI_Finalize();
      if ( mynode == rootnode )
      {
         cout << 
            "Error, failed to determine local indices w.r.t. global data: "
            << endl;
      }
      return EXIT_FAILURE;
   }
   
   if ( ndims == 1 ) 
   {
      Nx_total = dims[0];
      Ny = 1;
      Nz = 1;
   }
   else if ( ndims == 2 ) 
   {
      Nx_total = dims[0];
      Ny = dims[1];
      Nz = 1;
   }
   else if ( ndims == 3 ) 
   {
      Nx_total = dims[0];
      Ny = dims[1];
      Nz = dims[2];
   }

   int periodic; periodic = 1; // assume periodic boundary conditions
   std::vector<int> periodicity(ndims, periodic); // all identical for now

   size_t Nvoxel_neighbors; Nvoxel_neighbors = 2; 
   std::vector<double> phi_local(Nx_local + 2, 0);
   //if ( ndims == 1 ) phi_local.resize(Nx_local + 2, 0);
   if ( ndims == 2 ) 
   {
      phi_local.resize((Nx_local + 2)*Ny, 0);
      Nvoxel_neighbors = 4;
      if ( mynode == rootnode )
         std::cout << "Error: program not yet capable of 2-D, only 3-D."
            << std::endl;
      return EXIT_FAILURE;
   }
   else if ( ndims == 3 ) 
   {
      phi_local.resize((Nx_local + 2)*Ny*Nz, 0);
      Nvoxel_neighbors = 6;
   }
   else
   {
      if ( mynode == rootnode )
         std::cout << "Error: program only currently capable of 3-D"
            << std::endl;
      return EXIT_FAILURE;
   }

   if ( read_phi_from_hdf5( 
                        inFile_id,
                        phi_local, 
                        idx_start, 
                        idx_end,
                        flags,
                        //periodicity,
                        mynode, rootnode, totalnodes, world_comm
                        ) == EXIT_FAILURE)
   {
      flags.fail = -1;
   }

   if ( check_for_failure( flags, world_comm) )
   {
      H5Fclose( inFile_id );
      MPI_Finalize();
      if ( mynode == rootnode )
      {
         cout << 
            "Error, failed to read phi from file: "
            << inputFileName
            << endl;
      }
      return EXIT_FAILURE;
   }

   // close the initial state HDF5 file
   H5Fclose( inFile_id );

   // TODO: make voxel separation distance part of the input file
   //double hh_x; double hh_y; double hh_z;
   //hh_x = 1.0/(Nx_total ); 
   //hh_y = 1.0/(Ny ); 
   //hh_z = 1.0/(Nz );

   // TODO: make stability checks dynamic when diffusivity becomes variable
   // For the dynamics to be stable, require
   //  \delta t <= (h^2)/(2*D_{T}) on each independent axis.
   //if ( dt > hh_x*hh_x/(6.0*diffusivityT) ) failflag = -1;
   //if ( check_for_failure( flags.fail, world_comm) )
   //{
   //   MPI_Finalize();
   //   if ( mynode == rootnode )
   //   {
   //      cout << "Error, discretization is unstable: dt > dx/(2 D_T)" 
   //         << inputFileName
   //         << endl;
   //   }
   //   return EXIT_FAILURE;
   //}
   //if ( dt > hh_y*hh_y/(6.0*diffusivityT) ) failflag = -1;
   //if ( check_for_failure( flags.fail, world_comm) )
   //{
   //   MPI_Finalize();
   //   if ( mynode == rootnode )
   //   {
   //      cout << "Error, discretization is unstable: dt > dy/(2 D_T)" << inputFileName
   //         << endl;
   //   }
   //   return EXIT_FAILURE;
   //}
   //if ( dt > hh_z*hh_z/(6.0*diffusivityT) ) failflag = -1;
   //if ( check_for_failure( flags.fail, world_comm) )
   //{
   //   MPI_Finalize();
   //   if ( mynode == rootnode )
   //   {
   //      cout << "Error, discretization is unstable: dt > dz/(2 D_T)" << inputFileName
   //         << endl;
   //   }
   //   return EXIT_FAILURE;
   //}
   //size_t Nphi_local; Nphi_local = phi_local.size();
   ////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////
   // write run parameters to a text file
   ofstream log_file;
   if ((flags.debug != 0) && ( mynode == rootnode ))
   {
      cout << "Saving run parameters to file: " 
         << output_prefix + ".log" << endl;
      log_file.open(output_prefix + ".log", 
                        ios::app | ios::ate);

      if ( ! log_file.good() )
      {
         cout << "warning: could not open log file: "
            << output_prefix + ".log"
            << endl;
         log_file.close();
      }
      time_t start_time;
      std::time(&start_time);
      log_file << "Start " << std::ctime(&start_time) << endl;
      log_file << "-o " << output_prefix << endl;
      log_file << "-i " << inputFileName << endl;
      log_file << "-Nt " << Nt << endl;
      log_file << "-Nv " << Nv << endl;
      log_file << "-dt " << dt << endl;
      log_file << "-wp " << write_period << endl;
      log_file << endl;
      if ( flags.calcstat != 0 ) log_file << "-stat " << endl;
      log_file.close();
   }
   ////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////
   // Create communicators between cartesian neighbors
   int reorder; reorder = 0;
   int nprocs[ndims]; // number of processes per dimension

   // splitting only one dimension among processors
   nprocs[0] = totalnodes; 
   for (size_t i=1; i<ndims; ++i) nprocs[i] = 1;

   //cout << "node: " << mynode << ", nprocs[]: "; // debug
   //for (size_t i=0; i<ndims; ++i) cout << nprocs[i] << " ";// debug
   //cout << endl;// debug

   MPI_Cart_create( 
         world_comm, ndims, nprocs, 
         &periodicity[0], 
         reorder, &neighbors_comm);

   // identify the ranks of neighboring nodes
   int neighbor_x_lower, neighbor_x_higher;
   MPI_Cart_shift( neighbors_comm, 
         0, // direction (index of dimension) of shift
         1, // displacement, > 0 or < 0
         &neighbor_x_lower,
         &neighbor_x_higher
         );
   ////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////
   // open a file to write the fields evolutions to
   string outputFileName; outputFileName = output_prefix + ".h5";
   outFile_id = H5Fcreate(outputFileName.c_str(), 
                        //H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
                        H5F_ACC_TRUNC, H5P_DEFAULT, fa_plist_id);
   // open file
   if (outFile_id < 0) flags.fail= -1;
   if ( check_for_failure( flags, world_comm) )
   {
      MPI_Finalize();
      if ( mynode == rootnode )
      {
         cout << "Error, failed to open file for writing: " 
            << outputFileName << endl;
      }
      H5Fclose( outFile_id );
      return EXIT_FAILURE;
   }
   ////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////
   // open a file to store statistical information
   ofstream stat_file;
   if ( (flags.calcstat != 0) && (mynode == rootnode) )
   {
      cout << "opening file: " << output_prefix + "_phi_variance.txt" 
        << endl; // debug 
      stat_file.open(output_prefix + "_phi_variance.txt", 
                        ios::app | ios::ate);

      if ( ! stat_file.good() )
      {
         cout << "Error: could not open output file: "
            << output_prefix + "_phi_variance.txt"
            << endl;
         stat_file.close();
      }
   }
   ////////////////////////////////////////////////////////////////////
   
   ////////////////////////////////////////////////////////////////////
   // instantiate containers of local field changes
   //std::vector<double> phi_local_change(Nx_local*Ny*Nz,0);
   //std::vector<double> phi_local_change(phi_local.size(),0);
   //std::vector<double> phi_local_change(Nphi_local,0);

   //std::vector<std::vector<double>> // flux[i,j,k][2]
   //   phi_local_flux(phi_local.size(),std::vector<double>(2,0)); 
   std::vector<double> // (0:x-, 1:x+, 2:y-, 3:y+, 4:z-, 5:z+)
      phi_local_flux( Nvoxel_neighbors * phi_local.size() );//[n,i,j,k]

   std::vector<double> // (0:x-, 1:x+, 2:y-, 3:y+, 4:z-, 5:z+)
      phi_local_rates( 
            Nvoxel_neighbors * phi_local.size(), 0.0);//[n,i,j,k]
   ////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////
   // instantiate reusable variables
   
   // fluxes, assuming each element has only a single source voxel
   std::vector<double> phi_flux_upward( Ny * Nz, 0);
   std::vector<double> phi_flux_downward( Ny * Nz, 0);
   std::vector<double> phi_flux_from_above( Ny * Nz, 0);
   std::vector<double> phi_flux_from_below( Ny * Nz, 0);

   // flux rates of first passage to determine voxel filling order
   std::vector<double> phi_flux_upward_rates( Ny * Nz, 0);
   std::vector<double> phi_flux_downward_rates( Ny * Nz, 0);
   std::vector<double> phi_flux_from_above_rates( Ny * Nz, 0);
   std::vector<double> phi_flux_from_below_rates( Ny * Nz, 0);

   // remainders returned to neighbors after fillling
   //std::vector<double> phi_flux_upward_remainder( Ny * Nz, 0);
   //std::vector<double> phi_flux_downward_remainder( Ny * Nz, 0);
   //std::vector<double> phi_flux_from_above_remainder( Ny * Nz, 0);
   //std::vector<double> phi_flux_from_below_remainder( Ny * Nz, 0);

      
   std::vector<size_t> neigh_idxs(Nvoxel_neighbors, 0);

   //std::vector<double> jump_rates(6,0); // overwritten at every voxel
   // (0:x-, 1:x+, 2:y-, 3:y+, 4:z-, 5:z+)

   //std::vector<size_t> neigh_order(Nvoxel_neighbors, 0);
   //std::uniform_real_distribution<double> rand_decimal(0,1);// for order
   //std::vector<double> rand_decimals1(Nvoxel_neighbors, 0);
   //std::vector<double> rand_decimals2(Nvoxel_neighbors, 0);
   
   MPI_Request halo_flux_recv_requests[4]; // four Irecv per halo
   MPI_Request halo_flux_send_requests[4]; // four Isend per halo

   MPI_Request halo_accepted_flux_recv_requests1[2]; // two Irecv per halo
   MPI_Request halo_accepted_flux_send_requests1[2]; // two Isend per halo
   MPI_Request halo_accepted_flux_recv_requests2[2]; // two Irecv per halo
   MPI_Request halo_accepted_flux_send_requests2[2]; // two Isend per halo

   // TODO: see if the random class may be instantiated only once
   SPF_NS::random rr;//( 0.01, dt);

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

   // variables for statistical analysis
   double phi_local_sum; phi_local_sum = 0;
   double phi_local_sqr_sum; phi_local_sqr_sum = 0;
   double phi_total_sum; phi_total_sum = 0;
   double phi_total_sqr_sum; phi_total_sqr_sum = 0;
   double phi_mean; phi_mean = 0;
   double phi_variance; phi_variance = 0;

   ////////////////////////////////////////////////////////////////////

   /*-----------------------------------------------------------------*/
   /* begin loop over time -------------------------------------------*/
   /*-----------------------------------------------------------------*/
   for (time_step = 0; time_step <= Nt; ++time_step)
   {

      //////////////////////////////////////////////////////////////////
      // append fields to a file
      if ( time_step % write_period == 0 )
      {
         // write the local subset of phi to the file
         if ( append_phi_to_hdf5_multinode(
                  outFile_id,
                  time_step,
                  time,
                  phi_local,
                  Nx_local,
                  dims, 
                  idx_start,
                  idx_end,
                  dx_plist_id,
                  mynode, rootnode, totalnodes, world_comm
                  ) == EXIT_FAILURE)
         {
            flags.fail = -1;
         }
         if ( check_for_failure( flags, world_comm) )
         {
            if ( mynode == rootnode )
            {
               cout << "Error, failure while writing file: " 
                  << outputFileName << endl;
            }
            H5Fclose( outFile_id );
            MPI_Comm_free( &neighbors_comm); 
            MPI_Finalize();
            return EXIT_FAILURE;
         }
      }
      //////////////////////////////////////////////////////////////////
   
      //////////////////////////////////////////////////////////////////
      // increment time and exit if it exceeds the requested number
      time += dt;
      if ( time_step > Nt )
      {
         H5Fclose( outFile_id );
         MPI_Comm_free( &neighbors_comm); 
         MPI_Finalize();
         return EXIT_SUCCESS;
      }
      //////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////
      // update ghosts between neighbors
      
      //cout << "node " << mynode << " local data before comms:" // debug
      //    << endl;//debug
      //for (size_t i=0; i < (Nx_local +2); ++i) // debug
      //{ // debug
      //   for (size_t j=0; j < Ny; ++j) // debug
      //   {
      //      cout << "node " << mynode << " ["; // debug
      //      for (size_t k=0; k < Nz; ++k) // debug
      //         cout << setw(5) << phi_local[ k + Nz*(j + i*Ny) ];//debug
      //      cout << "] " << endl; // debug
      //   }
      //   cout << endl; // debug
      //} // debug
      //cout << endl;// debug
      
      update_ghosts(
                     phi_local,
                     Nx_local,
                     Ny,
                     Nz,
                     neighbor_x_higher,
                     neighbor_x_lower, 
                     neighbors_comm
                     );
      
      // cout << "node " << mynode << " local data after comms:" // debug
      //   << endl;//debug
      //for (size_t i=0; i < (Nx_local +2); ++i) // debug
      //{ // debug
      //   for (size_t j=0; j < Ny; ++j) // debug
      //   {
      //      cout << "node " << mynode << " ["; // debug
      //      for (size_t k=0; k < Nz; ++k) // debug
      //         cout << setw(8) << setprecision(10) // debug
      //          << phi_local[ k + Nz*(j + i*Ny) ]; // debug
      //      cout << "] " << endl; // debug
      //   }
      //   cout << endl; // debug
      //} // debug
      //cout << endl;// debug
      //////////////////////////////////////////////////////////////////


      //////////////////////////////////////////////////////////////////
      // receive inward flux from neighboring nodes
      flux_exchange_irecv(
            phi_flux_from_above, // store data from neighbor_x_higher here
            phi_flux_from_above_rates,
            phi_flux_from_below, // store data from neighbor_x_lower here
            phi_flux_from_below_rates,
            Ny, Nz,
            neighbor_x_higher, 
            neighbor_x_lower, 
            halo_flux_recv_requests, neighbors_comm
            );

      // NOTE:
      // Filling is a process of varying inward flux, in which
      //  the first passage determines which neighbor 
      //  contributes.
      // Emptying is a process of outward flux and distributes
      //  among neighbors according to either the respective 
      //  barrier height or gradient if there is no barrier.
      // Temporal ordering of flux at both upper and lower 
      //  bounds of voxel population requires simultaneous 
      //  knowledge of all of the fluxes to and from a voxel, 
      //  so a pairwise flux variable should be used.
      // Combining outward fluxes without attention to the
      //  upper filling limit disregards their order of 
      //  arrival and prioritizes flux from those voxels that
      //  happen to have their outward fluxes evaluated first.
      // Knowledge of inward flux is required to enforce the
      //  upper bound of a voxel's population, while knowledge
      //  of outward flux is required to enforce the lower
      //  bound.

      /****************************************************************/
      /* loop over voxels in phi_local, skipping ghosts ***************/

      // Local field doesn't change until after evaluating flux for 
      //  every voxel.
      
      size_t idx; 
      size_t mm;
      for (size_t ii=1; ii < Nx_local +1; ++ii) // loop over non-ghosts
         for ( size_t jj=0; jj < Ny; ++jj)
            for ( size_t kk=0; kk < Nz; ++kk)
            {
               // indices wrt local field
               //size_t idx; 
               idx = kk + Nz*(jj + Ny*ii);
               //std::cout << "phi_local[" << idx << "]: " // debug
               //   << phi_local[idx] << std::endl;// debug

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

               // assign values to jump_rates[]
               for( size_t nn=0; nn < (Nvoxel_neighbors/2); ++nn)
               {  // performing this once for each neighbor pair achieved
                  // by only evaluating the neighbors in the positive 
                  // x, y, and z directions and assuming periodic 
                  // boundary conditions.

                  mm = (2*nn) +1; // 1:x+, 3:y+, 5:z+
                  simple_identity_rate_gradient(
                        phi_local_rates[
                                    mm + Nvoxel_neighbors * idx],
                        phi_local,
                        neigh_idxs[mm],
                        idx
                        );

                  if ( isnan( phi_local_rates[
                                 mm + Nvoxel_neighbors*idx ] ))
                  {
                     if ( flags.debug != 0)
                        std::cout << "Error, node "
                           << mynode
                           << ": phi_local_rates["
                           << mm + Nvoxel_neighbors*idx 
                           << "] is a NaN." 
                           << " phi_local[idx] : "
                           << phi_local[idx]
                           << std::endl;
                        flags.fail = 1;
                  }

               }

               // evaluate stochastic fluxes to neighbor cells
               conserved_jump_flux_pairwise_distributions(
                     phi_local_flux,
                     rr,
                     phi_local,
                     phi_local_rates,
                     dt,
                     Nvoxel_neighbors,
                     neigh_idxs,
                     phi_upper_limit,
                     phi_lower_limit,
                     Nx_local, Ny, Nz,
                     ii, jj, kk
                     );

               for(size_t nn=0; nn < (Nvoxel_neighbors/2); ++nn)
               {
                  mm = (2*nn) +1; // 1:x+, 3:y+, 5:z+
                  //if( phi_local_flux[mm + Nvoxel_neighbors*idx] < 0)
                  //{ 
                  //   // debug
                  //   //std::cout << "Warning: step "
                  //   //   << time_step << " phi_local_flux["
                  //   //   << nn + Nvoxel_neighbors*idx << "] < 0;"
                  //   //   << " setting it to 0." << std::endl;
                  //   //// end debug
                  //   phi_local_flux[mm + Nvoxel_neighbors*idx] = 0;
                  //}

                  if ( isnan( phi_local_flux[mm + Nvoxel_neighbors*idx]))
                  {
                     if ( flags.debug != 0)
                     {
                        std::cout << "Error, node " << mynode 
                           << ": phi_local_flux["
                           << mm + Nvoxel_neighbors*idx << "] is a NaN."
                           << " phi_local["
                           << mm + Nvoxel_neighbors*idx << "]: "
                           << phi_local[mm + Nvoxel_neighbors*idx]
                           << std::endl;
                     }
                     flags.fail = 1;
                  }
                  //if ( (phi_local[idx] == 0) &&
                  //      (phi_local_flux[mm + Nvoxel_neighbors*idx] > 0)
                  //   )
                  //{
                  //   phi_local_flux[mm + Nvoxel_neighbors*idx] = 0.0;
                  //}
                  //mm = 2*nn;
                  //if ( (phi_local[neigh_idxs[mm]] == 0) &&
                  //      (phi_local_flux[
                  //            mm + Nvoxel_neighbors*idx
                  //            ] < 0)
                  //   )
                  //{
                  //   phi_local_flux[ mm + Nvoxel_neighbors*idx] = 0.0;
                  //}
                  //if( phi_local_flux[nn + Nvoxel_neighbors*idx] 
                  //      > phi_local[idx])
                  //{
                  //   // debug
                  //   std::cout << "Warning: step "
                  //      << time_step 
                  //      << " phi_local_flux["
                  //      << nn + Nvoxel_neighbors*idx 
                  //      << "] "
                  //      << phi_local_flux[nn + Nvoxel_neighbors*idx]
                  //      << " > phi_local[" << idx 
                  //      << "]  " << phi_local[idx]
                  //      << "; phi_local_rates_sqrt[" 
                  //      << nn + Nvoxel_neighbors*idx 
                  //      << "] : " 
                  //      << phi_local_rates_sqrt[
                  //            nn + Nvoxel_neighbors*idx ]
                  //      << std::endl;
                  //   // end debug
                  //}
               }

               //conserved_jump_flux_singl_distribution( 
               //      phi_local_change,
               //      phi_local,
               //      rr,
               //      dt,
               //      idx,
               //      neigh_idx_x_a,
               //      neigh_idx_x_b,
               //      neigh_idx_y_a,
               //      neigh_idx_y_b,
               //      neigh_idx_z_a,
               //      neigh_idx_z_b,
               //      //Nx_total,
               //      Ny,
               //      Nz
               //      );

               //// heat equation
               ////  \frac{\partial T}{\partial t} = D_{T} \nabla^{2} T
               //laplacian_flux(
               //      phi_local_change, 
               //      phi_local, 
               //      hh_x, hh_y, hh_z,
               //      dt,
               //      diffusivityT,
               //      idx,
               //      neigh_idx_x_a, 
               //      neigh_idx_x_b, 
               //      neigh_idx_y_a, 
               //      neigh_idx_y_b, 
               //      neigh_idx_z_a,
               //      neigh_idx_z_b,
               //      //ii, jj, kk, 
               //      //Nx_total, 
               //      Ny, Nz);
               
               //if ( phi_local_change[idx] > 2.0 ) // debug
               //   cout  // debug
               //      << "hh_x " << hh_x << ", "
               //      << "phi_local[ " << idx << "] " // debug
               //      << phi_local[idx]  // debug
               //      << ", phi_local_change[ " << idx << "] " // debug
               //      << phi_local_change[idx]  // debug
               //      << ", neighs "  // debug
               //      << phi_local[neigh_idx_x_a] << ", " // debug
               //      << phi_local[neigh_idx_x_b] << ", " // debug
               //      << phi_local[neigh_idx_y_a] << ", " // debug
               //      << phi_local[neigh_idx_y_b] << ", " // debug
               //      << phi_local[neigh_idx_z_a] << ", " // debug
               //      << phi_local[neigh_idx_z_b] << ", "
               //      << "formula gives: "
               //      << 
               //      (
               //       phi_local[neigh_idx_x_a]
               //       +phi_local[neigh_idx_x_b]
               //       -2*phi_local[idx]
               //      )/(hh_x*hh_x)
               //      << endl; // debug

               // jump process
               //conserved_jump_flux(
               //      phi_local_change, 
               //      phi_local, 
               //      scale_factor,
               //      ii, jj, kk, 
               //      //Nx_total, 
               //      Ny, Nz);

               //flux_test(
               //      phi_local_change, 
               //      phi_local, 
               //      ii, jj, kk, Nx_total, Ny, Nz);


               //// Copy outward flux to be sent to neighboring nodes.
               //// Also copy flux rates, to order inward flux to 
               ////  neighboring node voxels.

               //if (ii == 1 )  // lower x-axis boundary of non-ghosts
               //{
               //   phi_flux_downward[  kk + Nz*jj] 
               //      = phi_local_flux[ 0+ Nvoxel_neighbors * idx ]; 
               //   // (0:x-, 1:x+, 2:y-, 3:y+, 4:z-, 5:z+)

               //   phi_flux_downward_rates[kk + Nz*jj] 
               //      = phi_local_rates[0 + Nvoxel_neighbors * idx];
               //}
               //if (ii == Nx_local) // upper x-axis boundary of non-ghosts
               //{
               //   phi_flux_upward[ kk + Nz*jj] 
               //      = phi_local_flux[ 1+ Nvoxel_neighbors * idx];

               //   phi_flux_upward_rates[kk + Nz*jj] 
               //      = phi_local_rates[1 + Nvoxel_neighbors * idx];
               //}
            }

      /* end loop over voxels *****************************************/
      /****************************************************************/
      
      for ( size_t jj=0; jj < Ny; ++jj)
         for ( size_t kk=0; kk < Nz; ++kk)
         {
            // lower x-axis boundary of non-ghosts
            idx = kk + Nz*(jj + Ny*1);

            // Copy outward flux to be sent to neighboring nodes.
            // Also copy flux rates, to order inward flux to 
            //  neighboring node voxels.

            phi_flux_downward[  kk + Nz*jj] 
               = phi_local_flux[ 0+ Nvoxel_neighbors * idx ]; 
            // (0:x-, 1:x+, 2:y-, 3:y+, 4:z-, 5:z+)

            phi_flux_downward_rates[kk + Nz*jj] 
               = phi_local_rates[0 + Nvoxel_neighbors * idx];

            // upper x-axis boundary of non-ghosts
            idx = kk + Nz*(jj + Ny*Nx_local);

            phi_flux_upward[ kk + Nz*jj] 
               = phi_local_flux[ 1+ Nvoxel_neighbors * idx];

            phi_flux_upward_rates[kk + Nz*jj] 
               = phi_local_rates[1 + Nvoxel_neighbors * idx];
         }
      
      // send local flux info to neighboring nodes
      flux_exchange_isend(
         phi_flux_upward, // Ny*Nz
         phi_flux_upward_rates, // Ny*Nz
         phi_flux_downward, // Ny*Nz
         phi_flux_downward_rates, // Ny*Nz
         Ny, Nz,
         neighbor_x_higher, neighbor_x_lower, 
         halo_flux_send_requests, neighbors_comm
         );

      MPI_Waitall(4, halo_flux_recv_requests, MPI_STATUSES_IGNORE);
      // copy received values into local phi_local_flux
      for (size_t jj=0; jj < Ny; ++jj)
         for( size_t kk=0; kk < Nz; ++kk)
         {
            // flux to lower neighbor of upper ghost along x-axis
            // ghost x_idx = Nx_local+1; ghost neigh_idx = 0 (x-)
            idx = 0 + Nvoxel_neighbors* (kk + Nz*(jj + Ny*(Nx_local+1)));
            phi_local_flux[ idx ] 
               = phi_flux_from_above[kk + Nz*jj];

            phi_local_rates[ idx ]
               = phi_flux_from_above_rates[kk + Nz*jj];

            // flux to upper neighbor of lower ghost along x-axis
            // ghost x_idx = 0; ghost neigh_idx = 1 (x+)
            idx = 1 + Nvoxel_neighbors* (kk + Nz*(jj + Ny*0));
            phi_local_flux[ idx ]
               = phi_flux_from_below[kk + Nz*jj];

            phi_local_rates[ idx ]
               = phi_flux_from_below_rates[kk + Nz*jj];
         }

      // all ghost fluxes should now contain the iniital fluxes
      
      flux_accepted_irecv(
            phi_flux_from_below, // store data from neighbor_x_lower here
            phi_flux_from_above, // store data from neighbor_x_higher here
            Ny, Nz,
            neighbor_x_higher,
            neighbor_x_lower,
            halo_accepted_flux_recv_requests1,
            neighbors_comm
            );

      // debug
      //std::cout 
      //   << "boundary fluxes before enforcing of outward limits: ";
      //std::cout << "from below" << std::endl;
      //for ( size_t jj=0; jj < Ny; ++jj)
      //{
      //   for ( size_t kk=0; kk < Nz; ++kk)
      //   {
      //      idx = kk + Nz*(jj + Ny*0);
      //      std::cout << phi_local_flux[1 + Nvoxel_neighbors * idx]
      //         << ", ";
      //   }
      //   std::cout << std::endl;
      //}
      //std::cout << "upward" << std::endl;
      //for ( size_t jj=0; jj < Ny; ++jj)
      //{
      //   for ( size_t kk=0; kk < Nz; ++kk)
      //   {
      //      idx = kk + Nz*(jj + Ny*Nx_local);
      //      std::cout << phi_local_flux[1 + Nvoxel_neighbors * idx]
      //         << ", ";
      //   }
      //   std::cout << std::endl;
      //}
      // end debug
      /****************************************************************/
      /* enforce voxel value bounds ***********************************/
      enforce_bounds_pairwise_int_outward(
            // updates phi_local_flux with acceptable flux values
            phi_local_flux,
            phi_local,
            phi_local_rates,  // not necessarily integers
            rr,
            // neigh_order,
            Nvoxel_neighbors,
            phi_lower_limit,
            phi_upper_limit,
            Nx_local, Ny, Nz,
            eps,
            flags
            );

      MPI_Waitall(4, halo_flux_send_requests, MPI_STATUSES_IGNORE);

      // debug
      //std::cout 
      //   << "boundary fluxes before communication of outward flux: ";
      //std::cout << "from below" << std::endl;
      //for ( size_t jj=0; jj < Ny; ++jj)
      //{
      //   for ( size_t kk=0; kk < Nz; ++kk)
      //   {
      //      idx = kk + Nz*(jj + Ny*0);
      //      std::cout << phi_local_flux[1 + Nvoxel_neighbors * idx]
      //         << ", ";
      //   }
      //   std::cout << std::endl;
      //}
      //std::cout << "upward" << std::endl;
      //for ( size_t jj=0; jj < Ny; ++jj)
      //{
      //   for ( size_t kk=0; kk < Nz; ++kk)
      //   {
      //      idx = kk + Nz*(jj + Ny*Nx_local);
      //      std::cout << phi_local_flux[1 + Nvoxel_neighbors * idx]
      //         << ", ";
      //   }
      //   std::cout << std::endl;
      //}
      // end debug
      // TODO: check that accepted outward fluxes are communicated
      // update neighbor boundary voxels with approved outward flux
      for ( size_t jj=0; jj < Ny; ++jj)
         for ( size_t kk=0; kk < Nz; ++kk)
         { // Copy flux value to be sent to neighboring nodes.
            
            // acceptable flux at lower ghosts
            idx = kk + Nz*(jj + Ny*0);

            if ( phi_local_flux[ 1+ Nvoxel_neighbors * idx ] <= 0)
            {
               phi_flux_downward[  kk + Nz*jj] 
                  = phi_local_flux[ 1+ Nvoxel_neighbors * idx ]; 
            }

            // (0:x-, 1:x+, 2:y-, 3:y+, 4:z-, 5:z+)

            // upper x-axis boundary of non-ghosts
            // acceptable flux into upper ghosts from local voxels
            idx = kk + Nz*(jj + Ny*Nx_local);

            if ( phi_local_flux[ 1+ Nvoxel_neighbors * idx] >= 0) 
            {
               phi_flux_upward[ kk + Nz*jj] 
                  = phi_local_flux[ 1+ Nvoxel_neighbors * idx];
            }
         }
      
      //////////////////////////////////////////////////////////////////
      // copy outward flux to be sent to neighboring nodes 
      //for (size_t nn=0; nn < Nvoxel_neighbors; ++nn)

      //for (size_t jj=0; jj < Ny; ++jj)
      //   for (size_t kk=0; kk < Nz; ++kk)
      //{// one flux up and one down for each boundary voxel
      //   phi_flux_downward[  kk + Nz*jj] 
      //      = phi_local_flux[ // (0:x-, 1:x+, 2:y-, 3:y+, 4:z-, 5:z+)
      //            0+ Nvoxel_neighbors*(kk + Nz*(jj+ Ny*0))
      //           ];

      //   phi_flux_upward[ kk + Nz*jj] 
      //      = phi_local_flux[ 
      //            1+ Nvoxel_neighbors*(kk + Nz*(jj+ Ny*(Nx_local+1)))
      //           ];
      //}

      //for (size_t jj=0; jj < Ny; ++jj) // debug
      //{ // debug
      //   cout << "node " << mynode << " sending upward[ "; // debug
      //   for (size_t kk=0; kk < Nz; ++kk) // debug
      //   { // debug
      //  cout << setw(5) << phi_flux_upward[kk + Nz*jj] << ", "; // debug
      //   } // debug
      //   cout << "]" << endl; // debug
      //} // debug
      //for (size_t jj=0; jj < Ny; ++jj) // debug
      //{ // debug
      //   cout << "node " << mynode << " sending downward [ "; // debug
      //   for (size_t kk=0; kk < Nz; ++kk) // debug
      //   { // debug
      //  //cout << setw(5) << phi_flux_upward[kk+Nz*jj] << ", "// debug
      //  cout << setw(5) << phi_flux_downward[kk+Nz*jj] << ", ";// debug
      //   } // debug
      //   cout << "]" << endl; // debug
      //} // debug
      //////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////
      // send locally acceptable outward flux to neighboring nodes
      flux_accepted_isend(
         phi_flux_upward, // data to send to neighbor_x_higher
         phi_flux_downward, // data to send to neighbor_x_lower
         Ny, Nz,
         neighbor_x_higher, 
         neighbor_x_lower, 
         halo_accepted_flux_send_requests1, neighbors_comm
         );

      //////////////////////////////////////////////////////////////////

      MPI_Waitall(2, halo_accepted_flux_recv_requests1, 
                     MPI_STATUSES_IGNORE);

      //for (size_t jj=0; jj < Ny; ++jj) // debug
      //{ // debug
      //   cout << "node " << mynode << " received from above[ "; // debug
      //   for (size_t kk=0; kk < Nz; ++kk) // debug
      //   { // debug
      //      cout << setw(5) // debug
      //       << phi_flux_from_above[kk + Nz*jj] << ", "; // debug
      //   } // debug
      //   cout << "]" << endl; // debug
      //} // debug
      //for (size_t jj=0; jj < Ny; ++jj) // debug
      //{ // debug
      //   cout << "node " << mynode << " received from below[ "; // debug
      //   for (size_t kk=0; kk < Nz; ++kk) // debug
      //   { // debug
      //  cout << setw(5) << phi_flux_from_below[kk +Nz*jj] << ", ";// debug
      //   } // debug
      //   cout << "]" << endl; // debug
      //} // debug

      //if ( check_for_failure( flags, world_comm) )
      //{
      //   if ( mynode == rootnode )
      //   {
      //      cout << "Error, failure while writing file: " 
      //         << outputFileName << endl;
      //   }
      //   H5Fclose( outFile_id );
      //   MPI_Comm_free( &neighbors_comm); 
      //   MPI_Finalize();
      //   return EXIT_FAILURE;
      //}


      // copy flux from other nodes into phi_local_flux
      for (size_t jj=0; jj < Ny; ++jj)
         for( size_t kk=0; kk < Nz; ++kk)
         {
            idx = 1 + Nvoxel_neighbors* (kk + Nz*(jj + Ny*Nx_local));
            if ( phi_local_flux[ idx ] < 0) // outward wrt neighbor
            {
               if ( phi_flux_from_above[kk + Nz*jj]
                     > phi_local_flux[ idx ] )
               {
                  phi_local_flux[ idx ] 
                     = phi_flux_from_above[kk + Nz*jj];
               }
            }

            idx = 1 + Nvoxel_neighbors* (kk + Nz*(jj + Ny*0));
            if ( phi_local_flux[ idx ] > 0) // outward wrt neighbor
            {
               if ( phi_flux_from_below[kk + Nz*jj]
                     < phi_local_flux[ idx ])
               {
                  phi_local_flux[ idx ]
                     = phi_flux_from_below[kk + Nz*jj];
               }
            }


            //debug
            //if ((phi_local[kk + Nz*(jj + Ny*0)] == 0) 
            //      && (phi_local_flux[1+ Nvoxel_neighbors * (
            //            kk + Nz*(jj+ Ny*0)
            //            )] > 0))
            //{
            //   std::cout << "Warning: phi_local[" << 
            //            kk + Nz*(jj+ Ny*0)
            //      << "] == 0, but phi_local_flux["
            //      << 1+ Nvoxel_neighbors * idx << "] > 0 : "
            //      << phi_local_flux[1+ Nvoxel_neighbors * 
            //            ( kk + Nz*(jj+ Ny*0)) ] 
            //      << ", phi_flux_from_below["
            //      <<      kk + Nz*(jj+ Ny*0)
            //      << "] " 
            //      << phi_flux_from_below[ kk + Nz*(jj + Ny*0) ]
            //      << std::endl;
            //}
            //end debug
         }

      flux_accepted_irecv(
            phi_flux_from_below, // store data from neighbor_x_lower here
            phi_flux_from_above, // store data from neighbor_x_higher here
            Ny, Nz,
            neighbor_x_higher,
            neighbor_x_lower,
            halo_accepted_flux_recv_requests2,
            neighbors_comm
            );

      // debug
      //std::cout 
      //   << "boundary fluxes after communication of outward flux: ";
      //std::cout << "from below" << std::endl;
      //for ( size_t jj=0; jj < Ny; ++jj)
      //{
      //   for ( size_t kk=0; kk < Nz; ++kk)
      //   {
      //      idx = kk + Nz*(jj + Ny*0);
      //      std::cout << phi_local_flux[1 + Nvoxel_neighbors * idx]
      //         << ", ";
      //   }
      //   std::cout << std::endl;
      //}
      //std::cout << "upward" << std::endl;
      //for ( size_t jj=0; jj < Ny; ++jj)
      //{
      //   for ( size_t kk=0; kk < Nz; ++kk)
      //   {
      //      idx = kk + Nz*(jj + Ny*Nx_local);
      //      std::cout << phi_local_flux[1 + Nvoxel_neighbors * idx]
      //         << ", ";
      //   }
      //   std::cout << std::endl;
      //}
      // end debug

      // enforce_field_bounds( field, flux )
      //        using phi_local_flux[] for flux magnitudes
      //        and phi_local_rates[] for balancing those magnitudes
      //        saving acceptable fluxes in phi_local_flux
      // Considered as an outward flux, neighbor orders are 
      //  equally balanced (assuming equal barrier heights).
      // But when considered as an inward flux, the neighbor
      //  orders may be balanced using first passage 
      //  distributions.

      //randomize_neighbor_order(
      //      neigh_order,
      //      rr,   // random generator
      //      rand_decimal,  // uniform_distribution<int>
      //      rand_decimals1,// reused vector, not useful outside
      //      rand_decimals2, // reused vector, not useful outside
      //      flags
      //      );

      // debug
      //std::cout 
      //   << "boundary fluxes before enforcing inward flux bounds: ";
      //std::cout << "from below" << std::endl;
      //for ( size_t jj=0; jj < Ny; ++jj)
      //{
      //   for ( size_t kk=0; kk < Nz; ++kk)
      //   {
      //      idx = kk + Nz*(jj + Ny*0);
      //      std::cout << phi_local_flux[1 + Nvoxel_neighbors * idx]
      //         << ", ";
      //   }
      //   std::cout << std::endl;
      //}
      //std::cout << "upward" << std::endl;
      //for ( size_t jj=0; jj < Ny; ++jj)
      //{
      //   for ( size_t kk=0; kk < Nz; ++kk)
      //   {
      //      idx = kk + Nz*(jj + Ny*Nx_local);
      //      std::cout << phi_local_flux[1 + Nvoxel_neighbors * idx]
      //         << ", ";
      //   }
      //   std::cout << std::endl;
      //}
      // end debug

      // enforce_bounds_generic renormalizes excessive loss to available
      //  walkers, and excessive gain flux to their rates.
      //  Neither method guarantees that rounding the results will 
      //  yield the same total number of walkers lost or gained.

      // check that inward fluxes don't exceed local bounds
      enforce_bounds_pairwise_int_inward(
            // updates phi_local_flux with acceptable flux values
            phi_local_flux,
            phi_local,
            phi_local_rates,   // maybe use \sigma^2 ?
            rr,
            // neigh_order,
            Nvoxel_neighbors,
            phi_lower_limit,
            phi_upper_limit,
            Nx_local, Ny, Nz,
            eps,
            flags
            );

      MPI_Waitall(2, halo_accepted_flux_send_requests1, 
                     MPI_STATUSES_IGNORE);

      // debug
      //std::cout 
      //   << "boundary fluxes before communication of inward flux: ";
      //std::cout << "from below" << std::endl;
      //for ( size_t jj=0; jj < Ny; ++jj)
      //{
      //   for ( size_t kk=0; kk < Nz; ++kk)
      //   {
      //      idx = kk + Nz*(jj + Ny*0);
      //      std::cout << phi_local_flux[1 + Nvoxel_neighbors * idx]
      //         << ", ";
      //   }
      //   std::cout << std::endl;
      //}
      //std::cout << "upward" << std::endl;
      //for ( size_t jj=0; jj < Ny; ++jj)
      //{
      //   for ( size_t kk=0; kk < Nz; ++kk)
      //   {
      //      idx = kk + Nz*(jj + Ny*Nx_local);
      //      std::cout << phi_local_flux[1 + Nvoxel_neighbors * idx]
      //         << ", ";
      //   }
      //   std::cout << std::endl;
      //}
      // end debug

      // Update phi_flux_from_below / above with accepted inward fluxes 
      for (size_t jj=0; jj < Ny; ++jj)
         for (size_t kk=0; kk < Nz; ++kk)
         {
            idx = 1 + Nvoxel_neighbors* (kk + Nz*(jj + Ny*Nx_local));
            if ( phi_local_flux[ idx ] <= 0)
            {
               phi_flux_upward[kk + Nz*jj]
                  = phi_local_flux[ idx ];
            }

            idx = 1 + Nvoxel_neighbors* (kk + Nz*(jj + Ny*0));
            if ( phi_local_flux[ idx ] >= 0)
            {
               phi_flux_downward[kk + Nz*jj]
                  = phi_local_flux[ idx ];
            }
         }

      flux_accepted_isend(
         phi_flux_upward, // data to send to neighbor_x_higher
         phi_flux_downward, // data to send to neighbor_x_lower
         Ny, Nz,
         neighbor_x_higher, 
         neighbor_x_lower, 
         halo_accepted_flux_send_requests2, neighbors_comm
         );


      // wait for accepted fluxes to be received
      MPI_Waitall(2, halo_accepted_flux_recv_requests2, 
                  MPI_STATUSES_IGNORE);

      // return the accepted inward fluxes to their orgin flux variable
      for (size_t jj=0; jj < Ny; ++jj)
         for (size_t kk=0; kk < Nz; ++kk)
         {
            // TODO: double check that this is correct
            idx = 1 + Nvoxel_neighbors* (kk + Nz*(jj + Ny*Nx_local));
            if ( phi_local_flux[idx] > 0) // inward wrt neighbor
            {
               if ( phi_flux_from_above[kk + Nz*jj]
                     <
                     phi_local_flux[ idx ])
               {
                  phi_local_flux[ idx ] 
                     = phi_flux_from_above[kk + Nz*jj];
               }
            }

            idx = 1 + Nvoxel_neighbors* (kk + Nz*(jj + Ny*0));
            if ( phi_local_flux[ idx ] < 0) // inward wrt neighbor
            {
               if ( phi_flux_from_below[kk + Nz*jj]
                     > phi_local_flux[ idx ])
               {
                  phi_local_flux[ idx ]
                        = phi_flux_from_below[kk + Nz*jj];
               }
            }
         }
      // debug
      //std::cout 
      //   << "boundary fluxes after communication of inward flux: ";
      //std::cout << "from below" << std::endl;
      //for ( size_t jj=0; jj < Ny; ++jj)
      //{
      //   for ( size_t kk=0; kk < Nz; ++kk)
      //   {
      //      idx = kk + Nz*(jj + Ny*0);
      //      std::cout << phi_local_flux[1 + Nvoxel_neighbors * idx]
      //         << ", ";
      //   }
      //   std::cout << std::endl;
      //}
      //std::cout << "upward" << std::endl;
      //for ( size_t jj=0; jj < Ny; ++jj)
      //{
      //   for ( size_t kk=0; kk < Nz; ++kk)
      //   {
      //      idx = kk + Nz*(jj + Ny*Nx_local);
      //      std::cout << phi_local_flux[1 + Nvoxel_neighbors * idx]
      //         << ", ";
      //   }
      //   std::cout << std::endl;
      //}
      // end debug


      // Update phi_local by applying the accepted fluxes
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

               // Apply flux to and from each voxel (phi_local)
               for ( size_t nn=0; nn < (Nvoxel_neighbors/2); ++nn)
               {
                  mm = (2*nn) +1;

                  //if ( (phi_local[idx] 
                  //         - phi_local_flux[mm + Nvoxel_neighbors*idx])
                  //      < phi_lower_limit )
                  //{
                  //   phi_local_flux[mm + Nvoxel_neighbors*idx]
                  //      = phi_local[idx] - phi_lower_limit;

                  //   phi_local[idx] = phi_lower_limit;
                  //}

                  // Assuming all fluxes have been approved by boundary
                  //  checks.
                  phi_local[idx] 
                     -= phi_local_flux[mm + Nvoxel_neighbors*idx];

                  //phi_local[neigh_idxs[mm]] 
                  //   += phi_local_flux[mm + Nvoxel_neighbors*idx];
                  //   ^ this won't work if neigh_idx[mm] is a ghost.

                  mm = 2*nn;
                  phi_local[idx]
                     += phi_local_flux[
                           neigh_pairs[mm] 
                              + Nvoxel_neighbors*neigh_idxs[mm]];
               }

               // check to ensure enforce_bounds_generic worked
               if ( phi_local[idx] < phi_lower_limit )
               {
                  if ((abs(phi_local[idx] - phi_lower_limit )
                        > 2*eps.dblsqrt )
                     && (flags.debug != 0))// && (mynode == rootnode))
                  {
                     std::cout << "Warning: step " 
                        << time_step 
                        << " flux out of voxel caused lower bound breach"
                        << " phi_local[" << idx << "]: "
                        << phi_local[idx] 
                        << ", eps.dbl " 
                        << eps.dbl 
                        << ", (phi_local[]-phi_lower_limit) "
                        << (phi_local[idx]-phi_lower_limit)
                        << ", setting phi_local to its lower limit"
                        << std::endl;
                  }

                  // debug
                  //std::cout << "outward fluxes: ";
                  //for(size_t nn=0; nn < Nvoxel_neighbors; ++nn)
                  //{
                  //   std::cout 
                  //      << phi_local_flux[nn + Nvoxel_neighbors*idx]
                  //      << ", ";
                  //}
                  //std::cout << std::endl;
                  //if ( ii == 1 )
                  //{
                  //   std::cout << "phi_flux_downward[" << kk + Nz*jj 
                  //      << "] : " << phi_flux_downward[kk + Nz*jj ]
                  //      << std::endl;
                  //   std::cout << "phi_flux_from_above[" << kk + Nz*jj 
                  //      << "] : " << phi_flux_from_above[kk + Nz*jj ]
                  //      << std::endl;
                  //}
                  //if ( ii == Nx_local )
                  //{
                  //   std::cout << "phi_flux_upward[" << kk + Nz*jj 
                  //      << "] : " << phi_flux_upward[kk + Nz*jj ]
                  //      << std::endl;
                  //   std::cout << "phi_flux_from_below[" << kk + Nz*jj 
                  //      << "] : " << phi_flux_from_below[kk + Nz*jj ]
                  //      << std::endl;
                  //}
                  // end debug

                  phi_local[idx] = phi_lower_limit;
               }

               // Add inward flux to phi_local
               // TODO: adjust for pairwise flux
               //for ( size_t nn=0; nn < 0.5*Nvoxel_neighbors; ++nn)
               //{
               //   mm = (2*nn) +1;
               //   phi_local[idx]
               //      += phi_local_flux[
               //            neigh_pairs[nn] 
               //               + Nvoxel_neighbors * neigh_idxs[nn]
               //         ];
               //}

               if ( phi_local[idx] > phi_upper_limit )
               {
                  if ((abs(phi_local[idx] - phi_upper_limit )
                           > 
                           2*eps.dbl
                      )
                           // ^ guess of error to be acceptably lost
                     && (flags.debug != 0))// && (mynode == rootnode))
                  {
                     std::cout << "Warning: step " 
                        << time_step 
                        << " flux into voxel caused upper bound breach"
                        << " phi_local[" << idx << "] - phi_upper_limit: "
                        << phi_local[idx] - phi_upper_limit
                        << std::endl;
                     for ( size_t mm=0; mm < (Nvoxel_neighbors/2); ++mm)
                     {
                        // << " setting phi_local to its upper limit"
                        size_t oo;
                        oo = 2*mm +1; 
                        std::cout << "phi_local_flux[" << 
                           oo + Nvoxel_neighbors*idx
                           << "] : "
                        << phi_local_flux[ oo + Nvoxel_neighbors*idx ]
                        << std::endl;

                        size_t ee;
                        ee = 2*mm;
                        std::cout << "phi_local_flux[" << 
                              neigh_pairs[ee] 
                              + Nvoxel_neighbors*neigh_idxs[ee] 
                              << "] : " << phi_local_flux[ 
                                 neigh_pairs[ee] 
                                    + Nvoxel_neighbors*neigh_idxs[ee] ] 
                              << std::endl;
                     }
                  }
                  phi_local[idx] = phi_upper_limit;
               }
            }
      // end of loop over non-ghosts

      //if ( flags.debug != 0 )
      //{
      //   cout << "node " << mynode // debug
      //    << " local data after adding neighbor flux:"
      //      << endl;//debug
      //   for (size_t i=0; i < (Nx_local +2); ++i)
      //   {
      //      for (size_t j=0; j < Ny; ++j)
      //      {
      //         cout << "node " << mynode << " [";
      //         for (size_t k=0; k < Nz; ++k)
      //            cout << setw(9) << phi_local[ k + Nz*(j + i*Ny) ];
      //         cout << "] " << endl;
      //      }
      //      cout << endl;
      //   }
      //   cout << endl;
      //}

      /* end field value boundary enforcement *************************/
      /****************************************************************/

      // wait for remaining flux sends to complete
      MPI_Waitall(2, halo_accepted_flux_send_requests2, 
                  MPI_STATUSES_IGNORE);


      /////////////////////////////////////////////////////////////////
      // calculate mean and variance of all voxel populations
      phi_local_sum  = 0.0;
      phi_local_sqr_sum  = 0.0;
      if ( flags.calcstat != 0 )
      {
         for (size_t ii=1; ii < Nx_local+1; ++ii)
         {
            for (size_t jj=0; jj < Ny; ++jj)
            {
               for (size_t kk=0; kk < Nz; ++kk)
               {
                  // sum local populations
                  phi_local_sum += phi_local[ kk + Nz*(jj + Ny*(ii)) ];
                  // sum square of local populations
                  phi_local_sqr_sum 
                     += (phi_local[ kk + Nz*(jj + Ny*(ii))])
                        *(phi_local[ kk + Nz*(jj + Ny*(ii))]);
               }
            }
         }
         // reduce local population sums and sums of squares to root node
         MPI_Allreduce(&phi_local_sum, &phi_total_sum, 1, 
                        MPI_DOUBLE, MPI_SUM, world_comm );
         MPI_Allreduce(&phi_local_sqr_sum, &phi_total_sqr_sum, 1, 
                        MPI_DOUBLE, MPI_SUM, world_comm );
         // calculate mean and variance 
         phi_mean = phi_total_sum / (Nx_total * Ny * Nz);
         phi_variance = phi_total_sqr_sum / (Nx_total * Ny * Nz)
                        - phi_mean * phi_mean;
         /////////////////////////////////////////////////////////////////

         /////////////////////////////////////////////////////////////////
         // write mean and variation to a file

         if (mynode == rootnode) // debug
         {
            //cout << "time_step: " << time_step  // debug
            //   << ", (phi_mean, phi_variance): (" // debug
            //  << phi_mean << ", " << phi_variance << ")" << endl;//debug
         
            if (stat_file.good())
            {
               stat_file << setw(10) << setprecision(8) 
                  << time_step << " "
                  << setw(10) << setprecision(8) 
                  << time << " "
                  << setw(10) << setprecision(8) 
                  << phi_mean << " "
                  << setw(10) << setprecision(8) 
                  << phi_variance << endl;
            }

            // debug
            //std::cout << "step " << time_step << " , population: "
            //   << setw(12) << setprecision(10) << phi_total_sum 
            //   << std::endl;
            // end debug

            phi_total_sum = 0.0;
            phi_total_sqr_sum = 0.0;
         }

         phi_local_sum = 0.0;
         phi_local_sqr_sum = 0.0;
         /////////////////////////////////////////////////////////////////
      }

   }
   /*-----------------------------------------------------------------*/
   /* end loop over time ---------------------------------------------*/
   /*-----------------------------------------------------------------*/

   if ( (mynode == rootnode) && stat_file.is_open()) stat_file.close();
   H5Fclose( outFile_id );

   ////////////////////////////////////////////////////////////////////
   // write time to a log file
   if ( mynode == rootnode )
   {
      log_file.open(output_prefix + ".log", 
                        ios::app | ios::ate);

      if ( ! log_file.good() )
      {
         cout << "warning: could not open log file: "
            << output_prefix + ".log"
            << endl;
         log_file.close();
      }
      time_t end_time;
      std::time(&end_time);
      log_file << "End " << std::ctime(&end_time) << endl;
      //log_file << endl;
      log_file.close();
   }
   ////////////////////////////////////////////////////////////////////
   
   ////////////////////////////////////////////////////////////////////
   MPI_Comm_free( &neighbors_comm); 
   MPI_Finalize();
   return EXIT_SUCCESS;
}

