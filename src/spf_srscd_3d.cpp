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
// File: spf_srscd_3d.cpp

#include <iostream>  // cout, cin, cerr, endl
#include <iomanip>   // setw, setprecision
#include <fstream>   // ifstream, ofstream
#include <cstdlib>   // EXIT_SUCCESS, EXIT_FAILURE
#include <vector>
#include <string>
#include <math.h> // floor
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

   //double diffusivityT; diffusivityT = 3.0E-4;
   int Nv; Nv = 1; // number of walkers possible in a voxel
   //double tilt_alpha; tilt_alpha = 0.999; // when 1, potential isn't tilted
   //double ww; ww = -0.1; // order energy
   //double TT; TT = 540; // temperature

   // TODO: erase this and read from from a file containing parameters
   double ww; ww = -1.6; // order energy
   double shape_constant; shape_constant = 5.0;
   double mobility; mobility = 5.0;
   double kappa; kappa = 3.0;
   double c_alpha; c_alpha = 0.2;
   double c_beta; c_beta = 0.8;
   double hh_x; hh_x = 0.0;
   double hh_y; hh_y = 0.0;
   double hh_z; hh_z = 0.0;

   ////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////
   // read input parameters
   std::string inputFileName;
   std::vector<string> args( argv, argv + argc );
   std::string output_prefix;
   std::string datasetPath;
   int_flags flags;
   flags.fail = 0;

   // NOTE: there might not be a way to safely Bcast strings without  
   //       assuming they're ASCII, so read cmdline options on all nodes.
   if ( read_cmdline_options(
            args,
            dt,
            Nt,
            Nv,
            hh_x,
            ww,
            shape_constant,
            mobility,
            kappa,
            c_alpha,
            c_beta,
            write_period,
            flags,
            output_prefix,
            inputFileName,
            datasetPath,
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
                << "  -datasetPath " << datasetPath << std::endl
                << "  -Nt " << Nt << std::endl
                << "  -Nv " << Nv << std::endl
                << "  -mesh-size " << hh_x << std::endl
                << "  -order-energy " << ww << std::endl
                << "  -shape-constant " << shape_constant << std::endl
                << "  -mobility " << mobility << std::endl
                << "  -kappa " << kappa << std::endl
                << "  -c-alpha " << c_alpha << std::endl
                << "  -c-beta " << c_beta << std::endl
                << "  -dt " << dt << std::endl
                << "  -wp " << write_period << std::endl;
   }

   // establish field limits
   double phi_upper_limit, phi_lower_limit;
   phi_upper_limit = Nv;
   phi_lower_limit = 0.0;
   // copy mesh size (assuming qubic)
   hh_y = hh_x;
   hh_z = hh_x;
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
         datasetPath,
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

   size_t Nvoxel_neighbors; 
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
   if ( ndims == 3 ) 
   {
      phi_local.resize((Nx_local + 2)*Ny*Nz, 0);
      Nvoxel_neighbors = 6;

   }

   if ( read_dataset_from_hdf5( 
                        inFile_id,
                        phi_local, 
                        datasetPath,
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

   // If mesh size wasn't specified by input, assume system side length =1
   if ( hh_x == 0.0 ) // note Nx_total, Ny, Nz !=0 since initialized to 1
   {
      hh_x = 1.0/(Nx_total ); 
   }
   if ( hh_y == 0.0 )
   {
      hh_y = 1.0/(Ny ); 
   }
   if ( hh_z == 0.0 )
   {
      hh_z = 1.0/(Nz ); 
   }

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
   if ( mynode == rootnode )
   {
      if (flags.debug != 0) 
      {
         cout << "Saving run parameters to file: " 
            << output_prefix + ".log" << endl;
         log_file.open(output_prefix + ".log", 
                           ios::app | ios::ate);
      }

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
      phi_local_rates( Nvoxel_neighbors * phi_local.size() );//[n,i,j,k]
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

   MPI_Request halo_accepted_flux_recv_requests[2]; // two Irecv per halo
   MPI_Request halo_accepted_flux_send_requests[2]; // two Isend per halo

   // TODO: see if the random class may be instantiated only once
   SPF_NS::random rr;//( 0.01, dt);

   std::vector<size_t> neigh_pairs(Nvoxel_neighbors, 0);
   neigh_pairs[0] = 1;  // x upward (neighbor above along x)
   neigh_pairs[1] = 0;  // x downward (neighbor below along x)
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

      double flux_total; flux_total = 0.0;// debug
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
            phi_flux_from_above,
            phi_flux_from_above_rates,
            phi_flux_from_below,
            phi_flux_from_below_rates,
            Ny, Nz,
            neighbor_x_higher, neighbor_x_lower, 
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

               // assign values to jump rates

               //double upward_shift;
               //upward_shift = -0.5*0.00008617*TT*ww;

               /* Implementation for using 3D Laplacian */
               //double lap3d;  // 3-D Laplacian
               //lap3d = laplacian3d( hh_x, hh_y, hh_z,
               //                     Nz,
               //                     phi_local,
               //                     idx,
               //                     neigh_idxs[neigh_pairs[0]], //x+
               //                     neigh_idxs[neigh_pairs[1]], //x-
               //                     neigh_idxs[neigh_pairs[2]], //y+
               //                     neigh_idxs[neigh_pairs[3]], //y-
               //                     neigh_idxs[neigh_pairs[4]], //z+
               //                     neigh_idxs[neigh_pairs[5]]  //z-
               //         );

               //for( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
               //{
               //   //std::cout << "phi_local[" << idx << "]/Nv : "
               //   //   << phi_local[idx]/Nv << std::endl;// debug
               //   if ( phi_local[idx] > 0)
               //   {
               //      phi_local_rates[nn + Nvoxel_neighbors*idx]
               //         = double_well_srscd(
               //            *(double_well_srscd(
               //               phi_local[idx],
               //               c_alpha,
               //               c_beta,
               //               lap3d,
               //               shape_constant,
               //               kappa
               //               ) + ww); // order energy ww here is < 0
               //
               //         //= double_well_tilted(
               //         //   phi_local[idx]/Nv, // convert to conc
               //         //   upward_shift,
               //         //   ww,
               //         //   TT,
               //         //   tilt_alpha
               //         //   );// - 0.5*0.00008617*TT*ww;
               //   //// debug
               //   //std::cout << "node " << mynode 
               //   //   << " phi_local_rates[" 
               //   //   <<  nn + Nvoxel_neighbors*idx 
               //   //   << "] "
               //   //   << phi_local_rates[ 
               //   //   nn + Nvoxel_neighbors*idx ]
               //   //   << ", phi_local[" << idx << "]"
               //   //   << phi_local[idx]
               //   //   << "phi_local[" << idx << "]/" << Nv << " "
               //   //   << phi_local[idx ]/Nv
               //   //   //<< ", upward_shift " << upward_shift 
               //   //   //<< ", ww " << ww
               //   //   //<< ", TT " << TT
               //   //   //<< ", tilt_alpha " << tilt_alpha
               //   //   << std::endl; 
               //   //// end debug

               //      // convert to # walkers
               //      phi_local_rates[ nn + Nvoxel_neighbors*idx ] *= Nv;
               //   }
               //   else
               //   {
               //      phi_local_rates[nn + Nvoxel_neighbors*idx] = 0;
               //   }
               //}

               /* Implementation for using 1D 2nd derivative*/
               double lap1d;  // 1-D Laplacian
               double hh[3]; hh[0] = hh_x; hh[1] = hh_y; hh[2] = hh_z;
               for( size_t nn=0; nn < 0.5*Nvoxel_neighbors; ++nn)
               {
                  //std::cout << "phi_local[" << idx << "]/Nv : "
                  //   << phi_local[idx]/Nv << std::endl;// debug
                  if ( phi_local[idx] > 0)
                  {
                     lap1d = laplacian1d(
                              hh[nn],
                              phi_local,
                              idx,
                              neigh_idxs[neigh_pairs[2*nn+1]],//neigh below
                              neigh_idxs[neigh_pairs[2*nn]] //neigh above
                              ) / Nv;  // convert from pop. to conc.

                     // downward rate
                     phi_local_rates[ 2*nn +Nvoxel_neighbors*idx]
                        = (mobility/hh[nn])
                           *(double_well_srscd(
                              phi_local[idx]/Nv,
                              c_alpha,
                              c_beta,
                              lap1d,
                              shape_constant,
                              kappa
                              ) - ww); // order energy ww here is < 0

                     // copy to upward rate
                     phi_local_rates[ (2*nn+1) + Nvoxel_neighbors*idx]
                             = phi_local_rates[
                                 2*nn + Nvoxel_neighbors*idx];
                  //// debug
                  //std::cout << "node " << mynode 
                  //   << " phi_local_rates[" 
                  //   <<  nn + Nvoxel_neighbors*idx 
                  //   << "] "
                  //   << phi_local_rates[ nn + Nvoxel_neighbors*idx ]
                  //   << ", phi_local[" << idx << "] "
                  //   << phi_local[idx]
                  //   << ", phi_local[" << idx << "]/" << Nv << " "
                  //   << phi_local[idx ]/Nv
                  //   //<< ", ww " << ww
                  //   << std::endl; 
                  //// end debug

                     // convert to # walkers
                     phi_local_rates[ nn + Nvoxel_neighbors*idx ] *= Nv;
                  }
                  else
                  {
                     phi_local_rates[nn + Nvoxel_neighbors*idx] = 0;
                  }
               }

               // evaluate stochastic fluxes to neighbor cells
               conserved_jump_flux_separate_distributions( 
                     phi_local_flux,
                     rr,
                     phi_local,
                     phi_local_rates,
                     dt,
                     Nvoxel_neighbors,
                     neigh_idxs,
                     phi_upper_limit,
                     phi_lower_limit,
                     idx
                     );
               //for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
               //{
               //   if ( flags.debug != 0)
               //   {
               //      //std::cout 
               //      //   << time_step << " phi_local_flux["
               //      //   << nn + Nvoxel_neighbors*idx << "]: "
               //      //   << phi_local_flux[ nn + Nvoxel_neighbors*idx]
               //      //   << ", phi_local_rates[" 
               //      //   << nn +Nvoxel_neighbors*idx
               //      //   << "] : " 
               //      //   << phi_local_rates[ nn +Nvoxel_neighbors*idx] 
               //      //   << std::endl;
               //      flux_total 
               //         += phi_local_flux[ nn + Nvoxel_neighbors*idx]/Nv;
               //   }

               //}

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

      //// debug
      //std::cout << time_step << "flux_total : " << flux_total 
      //   << std::endl;
      //// end debug
      /* end loop over voxels *****************************************/
      /****************************************************************/
      
      enforce_bounds_int_outward(
            // updates phi_local_flux with acceptable flux values
            phi_local_flux,   // integers
            phi_local,        // integers
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

      //for (size_t ii=1; ii < Nx_local +1; ++ii) // loop over non-ghosts
      size_t iii;
      for ( size_t jj=0; jj < Ny; ++jj)
         for ( size_t kk=0; kk < Nz; ++kk)
         {
            iii = 1; // lower x-axis boundary of non-ghosts
            idx = kk + Nz*(jj + Ny*iii);

            // Copy outward flux to be sent to neighboring nodes.
            // Also copy flux rates, to order inward flux to 
            //  neighboring node voxels.

            phi_flux_downward[  kk + Nz*jj] 
               = phi_local_flux[ 0+ Nvoxel_neighbors * idx ]; 
            // (0:x-, 1:x+, 2:y-, 3:y+, 4:z-, 5:z+)

            phi_flux_downward_rates[kk + Nz*jj] 
               = phi_local_rates[0 + Nvoxel_neighbors * idx];

            iii = Nx_local; // upper x-axis boundary of non-ghosts
            idx = kk + Nz*(jj + Ny*iii);

            phi_flux_upward[ kk + Nz*jj] 
               = phi_local_flux[ 1+ Nvoxel_neighbors * idx];

            phi_flux_upward_rates[kk + Nz*jj] 
               = phi_local_rates[1 + Nvoxel_neighbors * idx];
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
      flux_exchange_isend(
         phi_flux_upward, // Ny*Nz
         phi_flux_upward_rates, // Ny*Nz
         phi_flux_downward, // Ny*Nz
         phi_flux_downward_rates, // Ny*Nz
         Ny, Nz,
         neighbor_x_higher, neighbor_x_lower, 
         halo_flux_send_requests, neighbors_comm
         );

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
      //  cout << setw(5) << flux_from_below[kk + Nz*jj] << ", "; // debug
      //   } // debug
      //   cout << "]" << endl; // debug
      //} // debug
      //////////////////////////////////////////////////////////////////

      MPI_Waitall(4, halo_flux_recv_requests, MPI_STATUSES_IGNORE);

      //////////////////////////////////////////////////////////////////
      // combine received flux with local values
      //cout << "node " << mynode // debug
      // << " local data before adding neighbor flux:" // debug
      //   << endl;//debug
      //for (size_t i=0; i < (Nx_local +2); ++i) // debug
      //{ // debug
      //   for (size_t j=0; j < Ny; ++j) // debug
      //   {
      //      cout << "node " << mynode << " ["; // debug
      //      for (size_t k=0; k < Nz; ++k) // debug
      //         cout << setw(5) << phi_local[ k + Nz*(j + i*Ny) ]
      //            << " "; // debug
      //      cout << "] " << endl; // debug
      //   }
      //   cout << endl; // debug
      //} // debug
      //cout << endl;// debug

      // copy flux from other nodes into phi_local_flux
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

      // enforce_bounds_generic renormalizes excessive loss to available
      //  walkers, and excessive gain flux to their rates.
      //  Neither method guarantees that rounding the results will 
      //  yield the same total number of walkers lost or gained.

      // check that inward fluxes don't exceed local bounds
      enforce_bounds_int_inward(
            // updates phi_local_flux with acceptable flux values
            phi_local_flux,   // integers
            phi_local,        // integers
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

      // Update phi_flux_from_below / above with accepted inward fluxes 
      //        from ghosts residing in phi_local_flux
      for (size_t jj=0; jj < Ny; ++jj)
         for (size_t kk=0; kk < Nz; ++kk)
         {
            idx = 0 + Nvoxel_neighbors* (kk + Nz*(jj + Ny*(Nx_local+1)));
            phi_flux_from_above[kk + Nz*jj]
               = phi_local_flux[ idx ];

            idx = 1 + Nvoxel_neighbors* (kk + Nz*(jj + Ny*0));
            phi_flux_from_below[kk + Nz*jj]
               = phi_local_flux[ idx ];
         }

      MPI_Waitall(4, halo_flux_send_requests, MPI_STATUSES_IGNORE);

      flux_accepted_irecv(
            phi_flux_downward,
            phi_flux_upward,
            Ny,
            Nz,
            neighbor_x_higher,
            neighbor_x_lower,
            halo_accepted_flux_recv_requests,
            neighbors_comm
            );

      // send accepted flux values back to sources
      flux_accepted_isend(
            phi_flux_from_above, // Ny*Nz
            phi_flux_from_below, // Ny*Nz
            Ny,
            Nz,
            neighbor_x_higher,
            neighbor_x_lower, 
            halo_accepted_flux_send_requests, // two Isend per halo
            neighbors_comm
            );

      // wait for accepted fluxes to be received
      MPI_Waitall(2, halo_accepted_flux_recv_requests, 
                  MPI_STATUSES_IGNORE);

      // return the accepted outward fluxes to their orgin flux variable
      for (size_t jj=0; jj < Ny; ++jj)
         for (size_t kk=0; kk < Nz; ++kk)
         {
            // Reduce the upward or downward flux to ghosts only if
            //  the neighboring node requests smaller values than 
            //  determined by the local node's boundary enforcement.
            idx = 0 + Nvoxel_neighbors* (kk + Nz*(jj + Ny*(1)));
            //if ( phi_flux_downward[kk + Nz*jj] < phi_local_flux[ idx ] )
            //{
            phi_local_flux[ idx ] = phi_flux_downward[kk + Nz*jj];
            //}

            idx = 1 + Nvoxel_neighbors* (kk + Nz*(jj + Ny*(Nx_local)));
            //if ( phi_flux_upward[kk + Nz*jj] < phi_local_flux[ idx ] )
            //{
            phi_local_flux[ idx ]
                  = phi_flux_upward[kk + Nz*jj];
            //}
            // After this loop, both local voxels and neighbor's ghosts
            //  should contain the smaller of the fluxes determined by
            //  both nodes, to prevent both over filling and 
            //  over drawing.
         }

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

               // Subtract outward flux from each voxel (phi_local)
               for(size_t nn=0; nn < Nvoxel_neighbors; ++nn)
               {
                  phi_local[idx] 
                     -= phi_local_flux[nn + Nvoxel_neighbors*idx];
               }

               // check to ensure enforce_bounds_generic worked
               if ( phi_local[idx] < phi_lower_limit )
               {
                  //if ( abs(phi_local[idx] - phi_lower_limit - outward_flux)
                  //      > (eps.dblsqrt 
                  //            * (phi_local[idx] - phi_lower_limit)))
                  if ( flags.debug != 0 )
                  {
                     std::cout << "node " << mynode << " Warning: step " 
                         << time_step 
                         << " flux out of voxel caused lower bound breach"
                         << " phi_local[" << idx << "]: "
                         << phi_local[idx] 
                         << ", (phi_local[]-phi_lower_limit) "
                         << (phi_local[idx]-phi_lower_limit)
                         << ", setting phi_local to its lower limit"
                         << ", phi_local_flux[nn + Nvoxel_neighbors*" 
                         << idx << "]: ";
                     for (size_t mm=0; mm < Nvoxel_neighbors; ++mm)
                     {
                        std::cout 
                           << phi_local_flux[mm + Nvoxel_neighbors*idx]
                           << ", ";
                     }
                     std::cout << std::endl;
                  }
                  flags.fail = 1;

                  phi_local[idx] = phi_lower_limit;
               }

               // Add inward flux to phi_local
               for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
               {
                  phi_local[idx]
                     += phi_local_flux[
                           neigh_pairs[nn] 
                              + Nvoxel_neighbors * neigh_idxs[nn]
                        ];
               }

               if ( phi_local[idx] > phi_upper_limit )
               {
                  if ( flags.debug != 0 )
                  {
                     std::cout << "node " << mynode << " Warning: step " 
                        << time_step 
                        << " flux into voxel caused upper bound breach"
                        << " phi_local[" << idx << "]: "
                        << phi_local[idx] 
                        // << " setting phi_local to its upper limit"
                        << std::endl;
                  }
                  phi_local[idx] = phi_upper_limit;
                  flags.fail = 1;
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
      //////////////////////////////////////////////////////////////////

      /////////////////////////////////////////////////////////////////
      // calculate mean and variance of all voxel populations
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

            phi_total_sum = 0.0;
            phi_total_sqr_sum = 0.0;
         }

         phi_local_sum = 0.0;
         phi_local_sqr_sum = 0.0;
         /////////////////////////////////////////////////////////////////
      }

      // wait for remainders to be sent
      MPI_Waitall(2, halo_accepted_flux_send_requests, 
                     MPI_STATUSES_IGNORE);
      //// end debug

      if ( flags.debug != 0)
      {
         if ( check_for_failure( flags, world_comm) )
         {
            H5Fclose( outFile_id );
            MPI_Comm_free( &neighbors_comm); 
            MPI_Finalize();
            if ( mynode == rootnode )
            {
               std::cout << "Error : failed to enforce voxel value limits"
                  << std::endl;
            }
            return EXIT_FAILURE;
         }
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

