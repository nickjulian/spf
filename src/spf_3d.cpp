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
// File: spf_3d.cpp

#include <iostream>  // cout, cin, cerr, endl
#include <iomanip>   // setw, setprecision
#include <fstream>   // ifstream, ofstream
#include <cstdlib>   // EXIT_SUCCESS, EXIT_FAILURE
#include <vector>
#include <string>
#include <math.h> // floor
#include <chrono>

#include <mpi.h>
#include "../include/hdf5.h"

#include "check_for_failure.hpp"
#include "read_parameter_file.hpp"
#include "readHDF5c.hpp"
#include "writeHDF5c.hpp"
#include "spf_communication.hpp"
#include "stochastic_rates.hpp"
#include "voxel_dynamics.hpp"

#include "flags.hpp" // int_flags flags;
#include "rand.hpp"
#include "macheps.hpp"  // determines the local machine's relative rounding error 

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

   int Nx_total, Nx_local, Ny, Nz, Nt; Nx_total = 1; Nx_local = 1; Ny = 1; Nz = 1; Nt = 100;
   int time_step; time_step = 0;
   int write_period; write_period = 1;
   double time, dt; time = 0; dt = 0.01;
   //double Nv; Nv = 0.0; // number of walkers possible in a voxel
   //double ww; ww = 0.0; // interchange (substitutional) energy
   double c_fluc_theta; c_fluc_theta = 0.016;
   double shape_constant; shape_constant = 5.783e9;                        
   double mobility; mobility = 2.5e-28;
   double kappa; kappa = 4e-5;
   double c_alpha; c_alpha = 0.05;
   double c_beta; c_beta = 0.95;
   double hh_x; hh_x = 1e-9;
   double hh_y; hh_y = 1e-9;
   double hh_z; hh_z = 1e-9;
   int fix_y; fix_y = 0;
   int c_extreme; c_extreme = 0;
   double h_d; h_d = 7.9e-18;
   double h_c; h_c = 292;
   double c_alpha_T; c_alpha_T = 0.05; // c_alpha based on temperature
   double c_beta_T; c_beta_T = 0.95; // c_beta based on temperature
   int Rho_W, Rho_Cr; Rho_W = 19200; Rho_Cr = 7200;
   double boltzmann_K; boltzmann_K = 8.31446261815324;
   double molar_volume; molar_volume = 8.35e-6;

   ////////////////////////////////////////////////////////////////////
   // read input parameters
   std::string inputFileName;
   std::string parameter_file="parameter_file";
   std::vector<string> args( argv, argv + argc );
   std::string output_prefix;
   std::string datasetPath;
   std::string TPath;
   int_flags flags;
   flags.fail = 0;

   // NOTE: there might not be a way to safely Bcast strings without  
   //       assuming they're ASCII, so read cmdline options on all nodes.
   if ( read_parameter_file(
            parameter_file,
            flags,
            dt,
            Nt,
            hh_x,
            shape_constant,
            mobility,
            kappa,
            c_alpha,
            c_beta,
            h_d,
            h_c,
            fix_y,
            c_extreme,
            c_fluc_theta,
            molar_volume,
            Rho_W,
            Rho_Cr,
            write_period,
            output_prefix,
            inputFileName,
            datasetPath,
            TPath,
            mynode,
            rootnode,
            MPI_COMM_WORLD
            ) != EXIT_SUCCESS )
   {
      if ( mynode == rootnode ) PRINT_USAGE;

      std::cout << "Error, missing the parameter_file." << std::endl;
      MPI_Comm_free( &neighbors_comm); 
      MPI_Finalize();
      return EXIT_FAILURE;
   }
   if ( (flags.debug != 0 ) && (mynode == rootnode ))
   {
      auto now = std::chrono::system_clock::now();
      std::time_t start_time = std::chrono::system_clock::to_time_t(now);
      std::cout << "Start " << std::ctime(&start_time) << std::endl;

      std::cout << "running with parameters: " << std::endl
                << "  -i " << inputFileName << std::endl
                << "  -o " << output_prefix << std::endl
                << "  -datasetPath " << datasetPath << std::endl
                << "  -TPath " << TPath << std::endl
                << "  -fix-y " << fix_y << std::endl
                << "  -c_extreme " << c_extreme << std::endl
                << "  -Nt " << Nt << std::endl
                << "  -dt " << dt << std::endl
                << "  -wp " << write_period << std::endl
                << "  -c-fluc-theta " << c_fluc_theta << std::endl
                << "  -mobility " << mobility << std::endl
                << "  -mesh-size " << hh_x << std::endl
                << "  -shape-constant " << shape_constant << std::endl
                << "  -c-alpha " << c_alpha << std::endl
                << "  -c-beta " << c_beta << std::endl
                << "  -kappa " << kappa << std::endl
                << "  -h-d " << h_d << std::endl
                << "  -h-c " << h_c << std::endl
                << "  -molar_volume " << molar_volume << std::endl
                << "  -Rho_W " << Rho_W << std::endl
                << "  -Rho_Cr " << Rho_Cr << std::endl;
      if ( flags.calcstat != 0 ) std::cout << "  -stat " << std::endl;
   }

   // establish field limits
   double phi_upper_limit, phi_lower_limit;
   phi_upper_limit = 1.0;
   phi_lower_limit = 0.0;
   // copy mesh size (assuming cubic voxel)
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
         std::cout << "Error, failed to open input file: " << inputFileName << std::endl;
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
      {
         std::cout << "Error, number of compute nodes > number of voxels along x direciton on concentration domain" << std::endl;
      }

   }

   if ( check_for_failure( flags, world_comm) )
   {
      MPI_Finalize();
      if ( mynode == rootnode )
      {
         std::cout << "Error, failed to read dimensionality from input file: " << inputFileName << std::endl;
      }
      return EXIT_FAILURE;
   }

   int ndims; 
   ndims = dims.size();

   // machine relative error
   epsilon eps;
   if ( flags.debug != 0)
   {
      std::cout << "node " << mynode << " machine epsilon: " << eps.dbl << std::endl;
      std::cout << std::endl;
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
         std::cout << "Error, failed to determine local indices w.r.t. global data: " << std::endl;
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
   
   // non-periodic boundary condition for y direction
   if (fix_y==1){periodicity[1] = 0;}

   size_t Nvoxel_neighbors; 
   std::vector<double> phi_local(Nx_local + 2, 0);
   std::vector<double> T_local(Nx_local + 2, 0);

   double noise_diff = 0;
   std::vector<double> phi_local_1_step_before(Nx_local + 2, 0);
   std::vector<double> phi_local_noise(Nx_local + 2, 0);
   std::vector<double> phi_local_noise_sum_gather(totalnodes, 0);
   std::vector<double> nu(Nx_local + 2, 0);
   std::vector<double> delta_G(Nx_local + 2, 0);

   //if ( ndims == 1 ) phi_local.resize(Nx_local + 2, 0);
   if ( ndims == 2 ) 
   {
      phi_local.resize((Nx_local + 2)*Ny, 0);
      T_local.resize((Nx_local + 2) * Ny, 0);
      Nvoxel_neighbors = 4;
      if ( mynode == rootnode )
      {
         std::cout << "Error: program not yet capable of 2-D, only 3-D."<< std::endl;
      }
      return EXIT_FAILURE;
   }
   if ( ndims == 3 ) 
   {
      phi_local.resize((Nx_local + 2)*Ny*Nz, 0);
      T_local.resize((Nx_local + 2) * Ny * Nz, 0);
      Nvoxel_neighbors = 6;
      phi_local_1_step_before.resize((Nx_local + 2) * Ny * Nz, 0);
      phi_local_noise.resize((Nx_local + 2) * Ny * Nz, 0);
      nu.resize((Nx_local + 2) * Ny * Nz, 0);
      delta_G.resize((Nx_local + 2) * Ny * Nz, 0);
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

   if (check_for_failure(flags, world_comm))
   {
      H5Fclose(inFile_id);
      MPI_Finalize();
      if (mynode == rootnode)
      {
         std::cout << "Error, failed to read phi from file: " << inputFileName << std::endl;
      }
      return EXIT_FAILURE;
   }

   if ( read_dataset_from_hdf5( 
                        inFile_id,
                        T_local, 
                        TPath,
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
         std::cout << "Error, failed to read T from file: " << inputFileName << std::endl;
      }
      return EXIT_FAILURE;
   }

   // close the initial state HDF5 file
   H5Fclose( inFile_id );
   
   // copy phi_local to phi_local_1_step_before to initial array
   size_t idx; idx = 0;
   for (size_t ii=1; ii < Nx_local +1; ++ii) {   // loop over non-ghosts
      size_t ioffset = Ny * ii;
      for (size_t jj=0; jj < Ny; ++jj) {
         size_t joffset = Nz * (jj + ioffset);
         for (size_t kk=0; kk < Nz; ++kk)
         {
            idx = kk + joffset;
            phi_local_1_step_before[idx]=phi_local[idx];
         }
      }
   }

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

   ////////////////////////////////////////////////////////////////////
   // write run parameters to a text file
   ofstream log_file;
   if ( mynode == rootnode )
   {
      std::cout << "Saving run parameters to file: " << output_prefix + ".log" << std::endl;
      std::cout << std::endl;
      log_file.open(output_prefix + ".log", ios::app | ios::ate);

      if ( ! log_file.good() )
      {
         std::cout << "warning: could not open log file: " << output_prefix + ".log" << std::endl;
         log_file.close();
      }
      
      auto now = std::chrono::system_clock::now();
      std::time_t start_time = std::chrono::system_clock::to_time_t(now);
      log_file << "Start " << std::ctime(&start_time) << std::endl;

      log_file << "running with parameters: " << std::endl;
      log_file << "  -i " << inputFileName << std::endl;
      log_file << "  -o " << output_prefix << std::endl;
      log_file << "  -datasetPath " << datasetPath << std::endl;
      log_file << "  -TPath " << TPath << std::endl;
      log_file << "  -fix-y " << fix_y << std::endl;
      log_file << "  -c_extreme " << c_extreme << std::endl;
      log_file << "  -Nt " << Nt << std::endl;
      log_file << "  -dt " << dt << std::endl;
      log_file << "  -wp " << write_period << std::endl;
      log_file << "  -c_fluc_theta " << c_fluc_theta << std::endl;
      log_file << "  -mobility " << mobility << std::endl;
      log_file << "  -mesh-size " << hh_x << std::endl;
      log_file << "  -shape-constant " << shape_constant << std::endl;
      log_file << "  -c-alpha " << c_alpha << std::endl;
      log_file << "  -c-beta " << c_beta << std::endl;
      log_file << "  -kappa " << kappa << std::endl;
      log_file << "  -h-d " << h_d << std::endl;
      log_file << "  -h-c " << h_c << std::endl;
      log_file << "  -molar_volume " << molar_volume << std::endl;
      log_file << "  -Rho_W " << Rho_W << std::endl;
      log_file << "  -Rho_Cr " << Rho_Cr << std::endl;

      if ( flags.calcstat != 0 ) log_file << "  -stat " << std::endl;
      log_file << std::endl;

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
         std::cout << "Error, failed to open file for writing: " << outputFileName << std::endl;
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
      std::cout << "opening stat debug file: " << output_prefix + "_phi_variance.txt" << std::endl; // debug 
      stat_file.open(output_prefix + "_phi_variance.txt", ios::app | ios::ate);

      if ( ! stat_file.good() )
      {
         std::cout << "Error: could not open output file: " << output_prefix + "_phi_variance.txt" << std::endl;
         stat_file.close();
      }
   }
   ////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////
   // instantiate containers of local field changes
   std::vector<double> // (0:x-, 1:x+, 2:y-, 3:y+, 4:z-, 5:z+)
      phi_local_flux( Nvoxel_neighbors * phi_local.size() );//[n,i,j,k]
   std::vector<double> // (0:x-, 1:x+, 2:y-, 3:y+, 4:z-, 5:z+)
      phi_local_rates( Nvoxel_neighbors * phi_local.size() );//[n,i,j,k]

   std::vector<double>
      T_local_flux( T_local.size() );//[n,i,j,k]
   std::vector<double> // (0:x, 1:y, 2:z)
      T_local_discrete( 0.5 * Nvoxel_neighbors * T_local.size() );//[n,i,j,k]
   ////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////
   // instantiate reusable variables
   
   // fluxes of conserved quantities, assuming each element has only a single source voxel
   std::vector<double> phi_flux_upward( Ny * Nz, 0);
   std::vector<double> phi_flux_downward( Ny * Nz, 0);
   std::vector<double> phi_flux_from_above( Ny * Nz, 0);
   std::vector<double> phi_flux_from_below( Ny * Nz, 0);

   // flux rates of first passage to determine voxel filling order
   std::vector<double> phi_flux_upward_rates( Ny * Nz, 0);
   std::vector<double> phi_flux_downward_rates( Ny * Nz, 0);
   std::vector<double> phi_flux_from_above_rates( Ny * Nz, 0);
   std::vector<double> phi_flux_from_below_rates( Ny * Nz, 0);

   std::vector<size_t> neigh_idxs(Nvoxel_neighbors, 0);
   
   MPI_Request halo_flux_recv_requests[4]; // four Irecv per halo
   MPI_Request halo_flux_send_requests[4]; // four Isend per halo
   MPI_Request halo_accepted_flux_recv_requests[2]; // two Irecv per halo
   MPI_Request halo_accepted_flux_send_requests[2]; // two Isend per halo
   
   // initial philox random number generation function
   uint64_t seed = std::chrono::high_resolution_clock::now().time_since_epoch().count() + static_cast<uint64_t>(mynode);
   SPF_NS::random rr(seed, mynode);

   std::vector<size_t> neigh_pairs(Nvoxel_neighbors, 0);
   neigh_pairs[0] = 1;  // x upward (neighbor above along x)
   neigh_pairs[1] = 0;  // x downward (neighbor below along x)
   idx = 0;
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

   double lap1d;  // 1-D Laplacian for phi, reused each timestep
   double lap1dT;  // 1-D Laplacian for T, reused each timestep
   double phi_local_noise_sum;
   double phi_globle_noise_sum;
   double rn1;
   double rn2;

   std::vector<double> hh( 3, 0); hh[0] = hh_x; hh[1] = hh_y; hh[2] = hh_z;

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
         if ( append_fields_to_hdf5_multinode(
                  outFile_id,
                  time_step,
                  time,
                  phi_local,
                  T_local,
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
               std::cout << "Error, failure while writing file: " << outputFileName << std::endl;
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

      // philox rng reset counter and change timestep
      rr.reset_counter(time_step);

      // Introduce concentration Gaussian noise
      size_t idx;
      phi_local_noise_sum = 0;
      for (size_t ii=1; ii < Nx_local +1; ++ii) {   // loop over non-ghosts
         size_t ioffset = Ny * ii;
         for (size_t jj=0; jj < Ny; ++jj) {
            size_t joffset = Nz * (jj + ioffset);
            for (size_t kk=0; kk < Nz; ++kk)
            {
               idx = kk + joffset;
               nu[idx] = sqrt((phi_local[idx]*(1-phi_local[idx])*(boltzmann_K/molar_volume)*T_local[idx]*c_fluc_theta)/(shape_constant));
               rn1 = r123::u01<double>(rr.philox_generator());
               rn2 = r123::u01<double>(rr.philox_generator());
               phi_local_noise[idx] = nu[idx] * (sqrt(-2 * log(rn1)) * cos(2 * PI * rn2));
               phi_local_noise_sum = phi_local_noise_sum + phi_local_noise[idx];
            }
         }
      }
      //std::cout << "Step:" << time_step << std::endl; //debug
      //std::cout << "Node:" << mynode << "  " << "phi_local_noise_sum=" << phi_local_noise_sum << std::endl; //debug

      phi_globle_noise_sum = 0;
      MPI_Allreduce(&phi_local_noise_sum, &phi_globle_noise_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      //std::cout << "Node:" << mynode << "  " << "phi_globle_noise_sum=" << phi_globle_noise_sum << std::endl; //debug

      noise_diff = phi_globle_noise_sum / (Nx_total * Ny * Nz);
      //std::cout << "Node:" << mynode << "  " << "noise_diff=" << noise_diff << std::endl; //debug

      phi_local_noise_sum=0;
      for (size_t ii=1; ii < Nx_local +1; ++ii) {   // loop over non-ghosts
         size_t ioffset = Ny * ii;
         for (size_t jj=0; jj < Ny; ++jj) {
            size_t joffset = Nz * (jj + ioffset);
            for (size_t kk=0; kk < Nz; ++kk)
            {
               idx = kk + joffset;
               phi_local_noise[idx] = phi_local_noise[idx] - noise_diff;
               phi_local[idx] = phi_local[idx] + phi_local_noise[idx];

               if (phi_local[idx] < 0)
               {
                  std::cout << "Step " << time_step << std::endl;
                  std::cout << "Warring:phi_local[ " << idx << "]=" << phi_local[idx] << "<0" << std::endl;
                  std::cout << "Setting phi_local to its previous value without noise" << std::endl;
                  phi_local[idx] = phi_local[idx] - phi_local_noise[idx];
               }

               if (phi_local[idx] > 1)
               {
                  std::cout << "Step " << time_step << std::endl;
                  std::cout << "Warring:phi_local[ " << idx << "]=" << phi_local[idx] << ">1" << std::endl;
                  std::cout << "Setting phi_local to its previous value without noise" << std::endl;
                  phi_local[idx] = phi_local[idx] - phi_local_noise[idx];
               }

               if ((phi_local[idx] < 0) || (phi_local[idx] > 1))
               {
                  std::cout << "Error: Step " << time_step << std::endl;
                  std::cout << "phi_local_noise[ " << idx << "]=" << phi_local_noise[idx] << std::endl;
                  std::cout << "noise_diff=" << noise_diff << std::endl;
                  std::cout << "phi_local[ " << idx << "]=" << phi_local[idx] << std::endl;
                  std::cout << "After correction, concentration fluctuation is still too large, making c<0 or c>1." << std::endl;

                  flags.fail = 1;
               }

               //phi_local_noise_sum=phi_local_noise_sum+phi_local_noise[idx]; //debug
            }
         }
      }

      //std::cout << "After modify error:" << time_step << std::endl; //debug
      //std::cout << "Node:" << mynode << "  " << "phi_local_noise_sum=" << phi_local_noise_sum << std::endl; //debug

      // Debug c fluctuation
      //MPI_Gather(&phi_local_noise_sum, 1, MPI_DOUBLE, phi_local_noise_sum_gather.data(), 1, MPI_DOUBLE, rootnode, MPI_COMM_WORLD);

      //if (mynode == rootnode)
      //{
      //   phi_globle_noise_sum=0;
      //   for (size_t ii = 0; ii < totalnodes; ++ii)
      //   {
      //      phi_globle_noise_sum=phi_globle_noise_sum+phi_local_noise_sum_gather[ii];
      //   }
      //
      //   if (phi_globle_noise_sum==0)
      //   {
      //      std::cout << "phi_globle_noise_sum=0" << std::endl;
      //   }
      //   else
      //   {
      //      std::cout << "Warring:phi_globle_noise_sum/=0" << std::endl;
      //      std::cout << "Warring:phi_globle_noise_sum=" << phi_globle_noise_sum << std::endl;
      //   }
      //}

      //////////////////////////////////////////////////////////////////
      // update ghosts between neighbors

      // cout << "node " << mynode << " local data before comms:" // debug
      //     << endl;//debug
      // for (size_t i=0; i < (Nx_local +2); ++i) // debug
      //{ // debug
      //    for (size_t j=0; j < Ny; ++j) // debug
      //    {
      //       cout << "node " << mynode << " ["; // debug
      //       for (size_t k=0; k < Nz; ++k) // debug
      //          cout << setw(5) << phi_local[ k + Nz*(j + i*Ny) ];//debug
      //       cout << "] " << endl; // debug
      //    }
      //    cout << endl; // debug
      // } // debug
      // cout << endl;// debug

      update_ghosts(
          phi_local,
          Nx_local,
          Ny,
          Nz,
          neighbor_x_higher,
          neighbor_x_lower,
          neighbors_comm);

      update_ghosts(
          T_local,
          Nx_local,
          Ny,
          Nz,
          neighbor_x_higher,
          neighbor_x_lower,
          neighbors_comm);

      // cout << "node " << mynode << " local data after comms:" // debug
      //   << endl;//debug
      // for (size_t i=0; i < (Nx_local +2); ++i) // debug
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
      // cout << endl;// debug
      //////////////////////////////////////////////////////////////////

      /****************************************************************/
      /* loop over voxels in phi_local, skipping ghosts ***************/

      // Local field doesn't change until after evaluating flux for every voxel.

      for (size_t ii=1; ii < Nx_local +1; ++ii) {   // loop over non-ghosts
         size_t ioffset = Ny * ii;
         for (size_t jj=0; jj < Ny; ++jj) {
            size_t joffset = Nz * (jj + ioffset);
            for (size_t kk=0; kk < Nz; ++kk)
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

               // Implementation for using 1D 2nd derivative
               for( size_t nn=0; nn < 0.5*Nvoxel_neighbors; ++nn)
               {
                  if ((phi_local[idx] > 0) && (phi_local[idx] < 1))
                  {
                     lap1d = laplacian1d(
                              hh[nn],
                              phi_local,
                              idx,
                              neigh_idxs[neigh_pairs[2*nn+1]],//neigh below
                              neigh_idxs[neigh_pairs[2*nn]] //neigh above
                              );

                     // set non-periodic along y direction
                     if ((fix_y == 1) && (jj == 0) && (nn == 1))
                     {
                        lap1d =(1.0/(hh[nn]*hh[nn]))*(-phi_local[idx]+ phi_local[neigh_idxs[neigh_pairs[2*nn]]]);
                     }

                     if ((fix_y == 1) && (jj == Ny-1) && (nn == 1))
                     {
                        lap1d =(1.0/(hh[nn]*hh[nn]))*(phi_local[neigh_idxs[neigh_pairs[2*nn+1]]]-phi_local[idx]);
                     }

                     // downward rate
                     phi_local_rates[ 2*nn +Nvoxel_neighbors*idx]
                        = (1/hh[nn])*(mobility/hh[nn])
                        *(double_well_srscd(
                        phi_local[idx],
                        c_alpha,
                        c_beta,
                        lap1d,
                        shape_constant,
                        kappa
                        ) + (boltzmann_K/molar_volume) * T_local[idx]*
                        log(phi_local[idx] / (1 - phi_local[idx])));
                     
                     // Calculate flux
                     phi_local_flux[2*nn + Nvoxel_neighbors * idx] = dt * phi_local_rates[2*nn + Nvoxel_neighbors * idx];

                     // copy to upward rate and flux
                     phi_local_rates[ (2*nn+1) + Nvoxel_neighbors*idx]
                             = phi_local_rates[
                                 2*nn + Nvoxel_neighbors*idx];
                     phi_local_flux[ (2*nn+1) + Nvoxel_neighbors*idx]
                             = phi_local_flux[
                                 2*nn + Nvoxel_neighbors*idx];
                  }
                  else
                  {
                     phi_local_rates[2*nn + Nvoxel_neighbors*idx] = 0;
                     phi_local_rates[ (2*nn+1) + Nvoxel_neighbors*idx]
                             = phi_local_rates[2*nn + Nvoxel_neighbors*idx];

                     phi_local_flux[2*nn + Nvoxel_neighbors*idx] = 0;
                     phi_local_flux[ (2*nn+1) + Nvoxel_neighbors*idx]
                             = phi_local_flux[2*nn + Nvoxel_neighbors*idx];
                  }
               }

               // set non-periodic along x direction (need further test)
               //if (fix_x == 1 && mynode==rootnode)
               //{
               //   if ((Ny * Nz - 1 < idx) && (idx < (Ny + Ny - 1) * Nz + Nz))
               //   {
               //      phi_local_rates[0 + Nvoxel_neighbors * idx] = 0;
               //   }
               //}
               //if (fix_x == 1 && mynode==totalnodes-1)
               //{
               //   if ((Ny * Nx_local * Nz - 1 < idx) && (idx < (Ny * Nx_local + Ny - 1) * Nz + Nz))
               //   {
               //      phi_local_rates[1 + Nvoxel_neighbors * idx] = 0;
               //   }
               //} 

               // set non-periodic along y direction
               if (fix_y == 1)
               {
                  if (jj == 0)
                  {
                     phi_local_flux[2 + Nvoxel_neighbors * idx] = 0;
                  }
                  
                  if (jj == Ny-1)
                  {
                     phi_local_flux[3 + Nvoxel_neighbors * idx] = 0;
                  }
               }
            }
         }
      }

      // end loop over voxels ******************************************
      //****************************************************************

      //////////////////////////////////////////////////////////////////
      // receive inward flux from neighboring nodes

      flux_exchange_irecv(
            phi_flux_from_above,
            phi_flux_from_above_rates,
            phi_flux_from_below,
            phi_flux_from_below_rates,
            Ny, Nz,
            neighbor_x_higher, neighbor_x_lower, 
            halo_flux_recv_requests, 
            neighbors_comm
            );

      // NOTE: T (temperature) are not conserved
      //  quantities and don't exchange particles like phi

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

      /////////////////////////////////////////////////////////////////////

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
               halo_flux_send_requests, 
               neighbors_comm
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

      //////////////////////////////////////////////////////////////////

      MPI_Waitall(4, halo_flux_send_requests, MPI_STATUSES_IGNORE);

      //////////////////////////////////////////////////////////////////

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
      MPI_Waitall(2, halo_accepted_flux_recv_requests, MPI_STATUSES_IGNORE);

      // return the accepted outward fluxes to their orgin flux variable
      for (size_t jj=0; jj < Ny; ++jj)
         for (size_t kk=0; kk < Nz; ++kk)
         {
            // Reduce the upward or downward flux to ghosts only if
            //  the neighboring node requests smaller values than 
            //  determined by the local node's boundary enforcement.
            idx = 0 + Nvoxel_neighbors* (kk + Nz*(jj + Ny*(1)));
            if ( phi_flux_downward[kk + Nz*jj] < phi_local_flux[ idx ] )
            {
            phi_local_flux[ idx ] = phi_flux_downward[kk + Nz*jj];
            }

            idx = 1 + Nvoxel_neighbors* (kk + Nz*(jj + Ny*(Nx_local)));
            if ( phi_flux_upward[kk + Nz*jj] < phi_local_flux[ idx ] )
            {
            phi_local_flux[ idx ]
                  = phi_flux_upward[kk + Nz*jj];
            }
            // After this loop, both local voxels and neighbor's ghosts
            //  should contain the smaller of the fluxes determined by
            //  both nodes, to prevent both over filling and 
            //  over drawing.
         }

      // Update phi_local by applying the accepted fluxes
      for (size_t ii=1; ii < Nx_local +1; ++ii) {   // loop over non-ghosts
         size_t ioffset = Ny * ii;
         for (size_t jj=0; jj < Ny; ++jj) {
            size_t joffset = Nz * (jj + ioffset);
            for (size_t kk=0; kk < Nz; ++kk)
            {
               idx = kk + joffset;

               // compromise performance for memory by repeating this
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

               // check to ensure phi_local not out of the range (0-1)
               if ( phi_local[idx] < phi_lower_limit )
               {
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

               if ( phi_local[idx] > phi_upper_limit )
               {
                  if ( flags.debug != 0 )
                  {
                     std::cout << "node " << mynode << " Warning: step " 
                         << time_step 
                         << " flux out of voxel caused upper bound breach"
                         << " phi_local[" << idx << "]: "
                         << phi_local[idx] 
                         << ", (phi_local[]-phi_upper_limit) "
                         << (phi_local[idx]-phi_upper_limit)
                         << ", setting phi_local to its upper limit"
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

                  phi_local[idx] = phi_upper_limit;
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

               // Prevent phi_local do not above 1
               if ( phi_local[idx] > phi_upper_limit )
               {
                  if ( flags.debug != 0 )
                  {
                     std::cout << "node " << mynode << " Warning: step " 
                        << time_step 
                        << " flux into voxel caused upper bound breach"
                        << " phi_local[" << idx << "]: "
                        << phi_local[idx] 
                        << " setting phi_local to its upper limit"
                        << std::endl;
                  }
                  phi_local[idx] = phi_upper_limit;
                  flags.fail = 1;
               }

               // Prevent phi_local do not below 0
               if ( phi_local[idx] < phi_lower_limit )
               {
                  if ( flags.debug != 0 )
                  {
                     std::cout << "node " << mynode << " Warning: step " 
                        << time_step 
                        << " flux into voxel caused lower bound breach"
                        << " phi_local[" << idx << "]: "
                        << phi_local[idx] 
                        << " setting phi_local to its lower limit"
                        << std::endl;
                  }
                  flags.fail = 1;

                  phi_local[idx] = phi_lower_limit;
               }
            }
         }
      }

      /****************************************************************/
      /* loop over voxels in T_local, skipping ghosts ***************/

      // Local field doesn't change until after evaluating flux for 
      //  every voxel.

      // Calculate T_local_flux
      for (size_t ii=1; ii < Nx_local +1; ++ii) {   // loop over non-ghosts
         size_t ioffset = Ny * ii;
         for (size_t jj=0; jj < Ny; ++jj) {
            size_t joffset = Nz * (jj + ioffset);
            for (size_t kk=0; kk < Nz; ++kk)
            {
               idx = kk + joffset;

               // compromise performance for memory by repeating this
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

               // Implementation T_local for using 1D 2nd derivative
               for (size_t nn = 0; nn < 0.5 * Nvoxel_neighbors; ++nn)
               {
                  // Implement heat equation***************************
                  lap1dT = laplacian1dT(
                            hh[nn],
                            T_local,
                            idx,
                            neigh_idxs[neigh_pairs[2*nn+1]],//neigh below
                            neigh_idxs[neigh_pairs[2*nn]] //neigh above
                            );

                  // set non-periodic along y direction
                  if ((fix_y == 1) && (jj == 0) && (nn == 1))
                  {
                     lap1dT = (1.0 / (hh[nn] * hh[nn])) * (-T_local[idx] + T_local[neigh_idxs[neigh_pairs[2 * nn]]]);
                  }

                  if ((fix_y == 1) && (jj == Ny - 1) && (nn == 1))
                  {
                     lap1dT = (1.0 / (hh[nn] * hh[nn])) * (T_local[neigh_idxs[neigh_pairs[2 * nn + 1]]] - T_local[idx]);
                  }

                  T_local_discrete[nn + 0.5 * Nvoxel_neighbors * idx] = h_d * lap1dT * dt;
               }

               if (c_extreme == 1)
               {
                  // Decide two well depths based on temperature
                  if (T_local[idx] > 0 && T_local[idx] < 100){c_alpha_T = 0.05;c_beta_T = 0.95;}
                  else if (T_local[idx] >= 100 && T_local[idx] < 200){c_alpha_T = 0.078839;c_beta_T = 0.921161;}
                  else if (T_local[idx] >= 200 && T_local[idx] < 300){c_alpha_T = 0.105252;c_beta_T = 0.894748;}
                  else if (T_local[idx] >= 300 && T_local[idx] < 400){c_alpha_T = 0.130857;c_beta_T = 0.869143;}
                  else if (T_local[idx] >= 400 && T_local[idx] < 500){c_alpha_T = 0.15641;c_beta_T = 0.84359;}
                  else if (T_local[idx] >= 500 && T_local[idx] < 600){c_alpha_T = 0.182436;c_beta_T = 0.817564;}
                  else if (T_local[idx] >= 600 && T_local[idx] < 700){c_alpha_T = 0.209427;c_beta_T = 0.790573;}
                  else if (T_local[idx] >= 700 && T_local[idx] < 800){c_alpha_T = 0.237953;c_beta_T = 0.762047;}
                  else if (T_local[idx] >= 800 && T_local[idx] < 900){c_alpha_T = 0.268802;c_beta_T = 0.731198;}
                  else if (T_local[idx] >= 900 && T_local[idx] < 1000){c_alpha_T = 0.303256;c_beta_T = 0.696744;}
                  else if (T_local[idx] >= 1000 && T_local[idx] < 1100){c_alpha_T = 0.34387;c_beta_T = 0.65613;}
                  else if (T_local[idx] >= 1100 && T_local[idx] < 1177){c_alpha_T = 0.397967;c_beta_T = 0.602033;}

                  // Calculate free energy difference based on temperature and concentration   
                  if (T_local[idx] > 1177)
                  {
                     delta_G[idx] = 0;
                  }
                  else
                  {
                     if (((phi_local[idx] < c_alpha_T) && (phi_local_1_step_before[idx] < c_alpha_T)) || ((phi_local[idx] > c_beta_T) && (phi_local_1_step_before[idx] > c_beta_T)))
                     {
                        delta_G[idx] = 0;
                     }
                     else 
                     {
                        delta_G[idx] = ((shape_constant * (phi_local[idx] - c_alpha) * (phi_local[idx] - c_alpha) * (c_beta - phi_local[idx]) * (c_beta - phi_local[idx]) 
                                 + (boltzmann_K/molar_volume) * T_local[idx] * (phi_local[idx] * log(phi_local[idx]) + (1 - phi_local[idx]) * log(1 - phi_local[idx])))
                                 / (1 / (phi_local[idx] / Rho_W + (1 - phi_local[idx]) / Rho_Cr)))
                                 - ((shape_constant * (phi_local_1_step_before[idx] - c_alpha) * (phi_local_1_step_before[idx] - c_alpha) * (c_beta - phi_local_1_step_before[idx]) * (c_beta - phi_local_1_step_before[idx])
                                 + (boltzmann_K/molar_volume) * T_local[idx] * (phi_local_1_step_before[idx] * log(phi_local_1_step_before[idx]) + (1 - phi_local_1_step_before[idx]) * log(1 - phi_local_1_step_before[idx])))
                                 / (1 / (phi_local[idx] / Rho_W + (1 - phi_local[idx]) / Rho_Cr)));
                     }
                  }
               }
               else
               {
                  if (T_local[idx] > 1177)
                  {
                     delta_G[idx] = 0;
                  }
                  else
                  {
                     delta_G[idx] = ((shape_constant * (phi_local[idx] - c_alpha) * (phi_local[idx] - c_alpha) * (c_beta - phi_local[idx]) * (c_beta - phi_local[idx]) 
                              + (boltzmann_K / molar_volume) * T_local[idx] * (phi_local[idx] * log(phi_local[idx]) + (1 - phi_local[idx]) * log(1 - phi_local[idx]))) 
                              / (1 / (phi_local[idx] / Rho_W + (1 - phi_local[idx]) / Rho_Cr))) 
                              - ((shape_constant * (phi_local_1_step_before[idx] - c_alpha) * (phi_local_1_step_before[idx] - c_alpha) * (c_beta - phi_local_1_step_before[idx]) * (c_beta - phi_local_1_step_before[idx]) 
                              + (boltzmann_K / molar_volume) * T_local[idx] * (phi_local_1_step_before[idx] * log(phi_local_1_step_before[idx]) + (1 - phi_local_1_step_before[idx]) * log(1 - phi_local_1_step_before[idx]))) 
                              / (1 / (phi_local[idx] / Rho_W + (1 - phi_local[idx]) / Rho_Cr)));
                  }
               }

               // Calculate T_local_flux by adding dump heat and laplacian
               T_local_flux[idx] = -(delta_G[idx] / h_c) * fabs(phi_local[idx] - phi_local_1_step_before[idx]) 
                                 + T_local_discrete[0 + 0.5 * Nvoxel_neighbors * idx] 
                                 + T_local_discrete[1 + 0.5 * Nvoxel_neighbors * idx] 
                                 + T_local_discrete[2 + 0.5 * Nvoxel_neighbors * idx];

            }
         }
      }

      // Update T_local
      for (size_t ii=1; ii < Nx_local +1; ++ii) {   // loop over non-ghosts
         size_t ioffset = Ny * ii;
         for (size_t jj=0; jj < Ny; ++jj) {
            size_t joffset = Nz * (jj + ioffset);
            for (size_t kk=0; kk < Nz; ++kk)
            {
               idx = kk + joffset;
               T_local[idx]+=T_local_flux[idx];
            }
         }
      }

      // copy phi_local to phi_local_1_step_before
      for (size_t ii=1; ii < Nx_local +1; ++ii) {   // loop over non-ghosts
         size_t ioffset = Ny * ii;
         for (size_t jj=0; jj < Ny; ++jj) {
            size_t joffset = Nz * (jj + ioffset);
            for (size_t kk=0; kk < Nz; ++kk)
            {
               idx = kk + joffset;
               phi_local_1_step_before[idx]=phi_local[idx];
            }
         }
      }

      // end of loop over non-ghosts
      //////////////////////////////////////////////////////////////////

      /////////////////////////////////////////////////////////////////
      // calculate mean and variance of all voxel populations
      if ( flags.calcstat != 0 )
      {
         for (size_t ii=1; ii < Nx_local +1; ++ii) {   // loop over non-ghosts
            size_t ioffset = Ny * ii;
            for (size_t jj=0; jj < Ny; ++jj) {
               size_t joffset = Nz * (jj + ioffset);
               for (size_t kk=0; kk < Nz; ++kk)
               {
                  // sum local populations
                  phi_local_sum += phi_local[ kk + joffset ];
                  // sum square of local populations
                  phi_local_sqr_sum 
                     += (phi_local[ kk + joffset])
                        *(phi_local[ kk + joffset]);
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
                  << phi_variance << std::endl;
            }

            phi_total_sum = 0.0;
            phi_total_sqr_sum = 0.0;
         }

         phi_local_sum = 0.0;
         phi_local_sqr_sum = 0.0;
         /////////////////////////////////////////////////////////////////
      }

      // wait for remainders to be sent

      MPI_Waitall(2, halo_accepted_flux_send_requests, MPI_STATUSES_IGNORE);

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
               std::cout << "Error : failed to enforce voxel value limits" << std::endl;
            }
            return EXIT_FAILURE;
         }
      }
   }
   /*-----------------------------------------------------------------*/
   /* end loop over time ---------------------------------------------*/
   /*-----------------------------------------------------------------*/

   if ( (mynode == rootnode) && stat_file.is_open())
   {
      std::cout << "close stat debug file: " << output_prefix + "_phi_variance.txt" << std::endl; // debug
      std::cout << std::endl;
      stat_file.close();
   }
   H5Fclose( outFile_id );

   ////////////////////////////////////////////////////////////////////
   // write time to a log file
   if ( mynode == rootnode )
   {
      log_file.open(output_prefix + ".log", ios::app | ios::ate);

      if ( ! log_file.good() )
      {
         std::cout << "warning: could not open log file: "
            << output_prefix + ".log"
            << std::endl;
         log_file.close();
      }
      auto now = std::chrono::system_clock::now();
      std::time_t end_time = std::chrono::system_clock::to_time_t(now);
      log_file << "End " << std::ctime(&end_time) << std::endl;
      log_file.close();

      if (flags.debug != 0 )
      {
         auto now = std::chrono::system_clock::now();
         std::time_t end_time = std::chrono::system_clock::to_time_t(now);
         std::cout << "End " << std::ctime(&end_time) << std::endl;
      }
   }
   ////////////////////////////////////////////////////////////////////
   
   ////////////////////////////////////////////////////////////////////
   MPI_Comm_free( &neighbors_comm); 
   MPI_Finalize();
   return EXIT_SUCCESS;
}


