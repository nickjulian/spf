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
#include <list>
#include <string>
#include <math.h> // floor
#include <time.h> // time_t, time, ctime

#include <mpi.h>
#include "../include/hdf5.h"

#include "free_energy_gl.hpp"

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
   //
   // This version implements the evolution of three coupled quantities:
   //  T: temperature following a heat equation with heat capacity
   //      multiplied by a noise driven phase term
   //  phi: phase, diffusing as a Brownian motion with equal mean and
   //       variance
   //  conc: concentration diffusing as a Poisson process with rate that
   //       is the Laplacian of the change in free energy w.r.t. conc
   //       where free energy is a Ginzberg-Landau type equation
   //       F(conc,phi,T,gradPhi)
   // 

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
   std::vector<hsize_t> dimsPhi, dimsT, dimsConc;

   // TODO: why are these array dimensions integers rather than size_t?
   int Nx_total, Nx_local, Ny, Nz, Nt;
   Nx_total = 1; Nx_local = 1; Ny = 1; Nz = 1;
   int time_step; time_step = 0;
   int write_period; write_period = 1;
   double time, dt; time = 0; dt = 1;

   //double diffusivityT; diffusivityT = 3.0E-4;
   double Nv; Nv = 1; // number of walkers possible in a voxel
   int Nv_int; Nv_int = 1;
   //double tilt_alpha; tilt_alpha = 0.999; // when 1, potential isn't tilted
   //double ww; ww = -0.1; // order energy
   //double TT; TT = 540; // temperature

   // hard coded simulation parameters
   // TODO: erase this and read from from a file containing parameters
   double c1; c1 = 100;
   double c2; c2 = 0.25;
   double c3; c3 = 7;
   double c4; c4 = 3;
   double c5; c5 = 0;
   double c6; c6 = 0;
   double cPrefactor; cPrefactor = -0.05;
   double alpha; alpha = -110;
   double T0; T0 = 10;
   double cbase; cbase = 2;
   double LL; LL = 2;
   double orderEnergy; orderEnergy = 0.1;
   double M_phi; M_phi = 1;
   double M_conc; M_conc = 1;
   double hh_x; hh_x = 0.0;
   double hh_y; hh_y = 0.0;
   double hh_z; hh_z = 0.0;
   double D_T = 3.0E-4; // TODO: find a good value for thermal diffusivity

   ////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////
   // read input parameters
   std::string inputFileName;
   std::vector<string> args( argv, argv + argc );
   std::string output_prefix;
   std::string datasetPathPhi;
   std::string datasetPathT;
   std::string datasetPathConc;
   int_flags flags;
   flags.fail = 0;

   // establish field limits
   double phi_upper_limit, phi_lower_limit;
   phi_upper_limit = 1.0;
   phi_lower_limit = 0.0;
   double T_upper_limit, T_lower_limit;
   T_upper_limit = 25.0;
   T_lower_limit = 0.0;
   double conc_upper_limit, conc_lower_limit;
   conc_upper_limit = 1.0;
   conc_lower_limit = 0.0;
   std::vector<double> fieldValueLimits {
         phi_lower_limit, phi_upper_limit,
         conc_lower_limit, conc_upper_limit,
         T_lower_limit //, T_upper_limit
   };
                  

   // NOTE: there might not be a way to safely Bcast strings without  
   //       assuming they're ASCII, so read cmdline options on all nodes.
   // TODO: update for T, phi, and conc
   if ( read_cmdline_options(
            args,
            dt,
            Nt,
            Nv_int,
            hh_x,
            c1,
            c2,
            c3,
            c4,
            c5,
            c6,
            cPrefactor,
            alpha,
            T0,
            cbase,
            LL,
            orderEnergy,
            D_T,
            M_phi,
            M_conc,
            fieldValueLimits,
            write_period,
            flags,
            output_prefix,
            inputFileName,
            datasetPathPhi,
            datasetPathT,
            datasetPathConc,
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
      std::cout << "running with parameters: " << std::endl
                << "  -o " << output_prefix << std::endl
                << "  -i " << inputFileName << std::endl
                << "  -datasetPathPhi " << datasetPathPhi << std::endl
                << "  -datasetPathT " << datasetPathT << std::endl
                << "  -datasetPathConc " << datasetPathConc << std::endl
                << "  -Nt " << Nt << std::endl
                << "  -Nv " << Nv_int << std::endl
                << "  -mesh-size " << hh_x << std::endl
                << "  -order-energy " << orderEnergy << std::endl
                //<< "  -shape-constant " << shape_constant << std::endl
                //<< "  -mobility " << mobility << std::endl
                //<< "  -kappa " << kappa << std::endl
                //<< "  -c-alpha " << c_alpha << std::endl
                //<< "  -c-beta " << c_beta << std::endl
                << "  -dt " << dt << std::endl
                << "  -wp " << write_period << std::endl;
   }

   // copy mesh size (assuming qubic system size)
   hh_y = hh_x;
   hh_z = hh_x;
   ////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////
   // read input data
   hid_t inFile_id, outFile_id;

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
         datasetPathPhi,
         dimsPhi,
         flags,
         mynode,
         rootnode,
         totalnodes,
         world_comm
         ) == EXIT_FAILURE )
   {
      flags.fail = -1;
   }
   if ( read_dims_from_hdf5(
         inFile_id,
         datasetPathT,
         dimsT,
         flags,
         mynode,
         rootnode,
         totalnodes,
         world_comm
         ) == EXIT_FAILURE )
   {
      flags.fail = -1;
   }
   if ( read_dims_from_hdf5(
         inFile_id,
         datasetPathConc,
         dimsConc,
         flags,
         mynode,
         rootnode,
         totalnodes,
         world_comm
         ) == EXIT_FAILURE )
   {
      flags.fail = -1;
   }
   if ( (dimsT.size() != dimsConc.size())
         || (dimsT.size() != dimsPhi.size()))
   {
      flags.fail = -1;
   }
   for ( size_t ii=0; ii < dimsT.size(); ++ii)
   {
      if ( (dimsT[ii] != dimsConc[ii]) || (dimsT[ii] != dimsPhi[ii]))
      {
         flags.fail = -1;
      }
   }

   // now that each field has the same dimensions, copy them to dims
   std::vector<hsize_t> dims( dimsPhi);

   int ndims; 
   ndims = dims.size();

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
         cout << "Error, failed to determine local indices w.r.t."
            << " global data " << endl;
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

   Nv = static_cast<double>(Nv_int);

   // instantiate local field variables
   size_t Nvoxel_neighbors; 
   std::vector<double> phi_local(Nx_local + 2, 0);
   std::vector<double> T_local(Nx_local + 2, 0);
   std::vector<double> conc_local(Nx_local + 2, 0);
   std::vector<double> dFdconc(Nx_local + 2, 0); // used in laplacian
   std::vector<double> free_energy_local(Nx_local + 2, 0);
   std::vector<double> deltaPhi(Nx_local + 2, 0);
   std::vector<double> deltaT(Nx_local + 2, 0);
   //if ( ndims == 1 ) phi_local.resize(Nx_local + 2, 0);
    
   double lap1dPhiX, lap1dPhiY, lap1dPhiZ; // reused each loop for stats
   if ( ndims == 2 ) 
   {
      phi_local.resize((Nx_local + 2)*Ny, 0);
      T_local.resize((Nx_local + 2)*Ny, 0);
      conc_local.resize((Nx_local + 2)*Ny, 0);

      dFdconc.resize((Nx_local + 2)*Ny, 0);
      free_energy_local.resize((Nx_local + 2)*Ny, 0);
      deltaPhi.resize((Nx_local + 2)*Ny, 0);
      deltaT.resize((Nx_local + 2)*Ny, 0);

      Nvoxel_neighbors = 4;
      if ( mynode == rootnode )
         std::cout << "Error: program not yet capable of 2-D, only 3-D."
            << std::endl;
      return EXIT_FAILURE;
   }
   if ( ndims == 3 ) 
   {
      phi_local.resize((Nx_local + 2)*Ny*Nz, 0);
      T_local.resize((Nx_local + 2)*Ny*Nz, 0);
      conc_local.resize((Nx_local + 2)*Ny*Nz, 0);
      dFdconc.resize((Nx_local + 2)*Ny*Nz, 0);
      free_energy_local.resize((Nx_local + 2)*Ny*Nz, 0);
      deltaPhi.resize((Nx_local + 2)*Ny*Nz, 0);
      deltaT.resize((Nx_local + 2)*Ny*Nz, 0);

      Nvoxel_neighbors = 6;
   }
   else
   {
      if ( mynode == rootnode )
         std::cout << "Error: program only currently capable of 3-D"
            << std::endl;
      return EXIT_FAILURE;
   }
   // Read data from a file into local field variables
   // Note: dFdconc should be updated each time step, not saved to file
   if ( read_dataset_from_hdf5( 
                        inFile_id,
                        phi_local, 
                        datasetPathPhi,
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
            "Error, failed to read " << datasetPathPhi << " from file: "
            << inputFileName
            << endl;
      }
      return EXIT_FAILURE;
   }

   if ( read_dataset_from_hdf5( 
                        inFile_id,
                        T_local, 
                        datasetPathT,
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
            "Error, failed to read " << datasetPathT << " from file: "
            << inputFileName
            << endl;
      }
      return EXIT_FAILURE;
   }

   if ( read_dataset_from_hdf5( 
                        inFile_id,
                        conc_local,  // units: concentration \in [0,1]
                        datasetPathConc,
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
            "Error, failed to read " << datasetPathConc << "from file: "
            << inputFileName
            << endl;
      }
      return EXIT_FAILURE;
   }

   // close the initial state HDF5 file
   H5Fclose( inFile_id );

   // ensure field values read from file are within limits
   /* loop over voxels in phi_local, skipping ghosts ***************/
      for (size_t ii=1; ii < Nx_local +1; ++ii) // loop over non-ghosts
         for ( size_t jj=0; jj < Ny; ++jj)
            for ( size_t kk=0; kk < Nz; ++kk)
            {
               if ( phi_local[ kk + Nz*(jj + Ny*ii)] > phi_upper_limit)
               {
                  std::cout << "error: input field "
                     "phi_local[" << ii << "," << jj << "," << kk << "]: "
                     << phi_local[ kk + Nz*(jj + Ny*ii)]
                     << " exceeds upper limit "
                     << phi_upper_limit << std::endl;
                  flags.fail = -1;
               }
               if ( phi_local[ kk + Nz*(jj + Ny*ii)] < phi_lower_limit)
               {
                  std::cout << "error: input field "
                     "phi_local[" << ii << "," << jj << "," << kk << "]: "
                     << phi_local[ kk + Nz*(jj + Ny*ii)]
                     << " is below lower limit "
                     << phi_lower_limit << std::endl;
                  flags.fail = -1;
               }
               if ( conc_local[ kk + Nz*(jj + Ny*ii)] > conc_upper_limit)
               {
                  std::cout << "error: input field "
                     "conc_local[" << ii << "," << jj << "," << kk << "]: "
                     << conc_local[ kk + Nz*(jj + Ny*ii)]
                     << " exceeds upper limit "
                     << conc_upper_limit << std::endl;
                  flags.fail = -1;
               }
               if ( conc_local[ kk + Nz*(jj + Ny*ii)] < conc_lower_limit)
               {
                  std::cout << "error: input field "
                     "conc_local[" << ii << "," << jj << "," << kk << "]: "
                     << conc_local[ kk + Nz*(jj + Ny*ii)]
                     << " is below lower limit "
                     << conc_lower_limit << std::endl;
                  flags.fail = -1;
               }
               if ( T_local[ kk + Nz*(jj + Ny*ii)] > T_upper_limit)
               {
                  std::cout << "error: input field "
                     "T_local[" << ii << "," << jj << "," << kk << "]: "
                     << T_local[ kk + Nz*(jj + Ny*ii)]
                     << " exceeds upper limit "
                     << T_upper_limit << std::endl;
                  flags.fail = -1;
               }
               if ( T_local[ kk + Nz*(jj + Ny*ii)] < T_lower_limit)
               {
                  std::cout << "error: input field "
                     "T_local[" << ii << "," << jj << "," << kk << "]: "
                     << T_local[ kk + Nz*(jj + Ny*ii)]
                     << " is below lower limit "
                     << T_lower_limit << std::endl;
                  flags.fail = -1;
               }
               if ( check_for_failure( flags, world_comm) )
               {
                  H5Fclose( outFile_id );
                  MPI_Comm_free( &neighbors_comm); 
                  MPI_Finalize();
                  if ( mynode == rootnode )
                  {
                     std::cout << "Error : input file contains values"
                        << " outside of limits"
                        << std::endl;
                  }
                  return EXIT_FAILURE;
               }

               // convert conc_local from concentration to # walkers 
               conc_local[ kk + Nz*(jj + Ny*ii)]
                  = round(Nv * conc_local[ kk + Nz*(jj + Ny*ii)]);
            }
   /****************************************************************/

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
   std::string errorMessage("Error: ");
   errorMessage += inputFileName;
   if ( dt > hh_x*hh_x/(6.0*D_T) )
   {
      flags.fail = -1;
      errorMessage += ", discretization is unstable: dt > dx^2/(6 D_T)";
      errorMessage += std::to_string(dt) +" > "
         + std::to_string(hh_x*hh_x/(2* D_T));
   }
   if ( dt > hh_y*hh_y/(6.0*D_T) )
   {
      flags.fail = -1;
      errorMessage += ", discretization is unstable: dt > dy^2/(6 D_T)";
      errorMessage += std::to_string(dt) +" > "
         + std::to_string(hh_y*hh_y/(6* D_T));
   }
   if ( dt > hh_z*hh_z/(6.0*D_T) )
   {
      flags.fail = -1;
      errorMessage += ", discretization is unstable: dt > dz/(6 D_T)";
      errorMessage += std::to_string(dt) +" > "
         + std::to_string(hh_z*hh_z/(6* D_T));
   }
   if ( check_for_failure( flags, world_comm) )
   {
      MPI_Finalize();
      if ( mynode == rootnode )
      {
         cout << errorMessage << endl;
      }
      return EXIT_FAILURE;
   }
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
      log_file << "Start " << std::ctime(&start_time) << endl
         << "-o " << output_prefix << endl
         << "-i " << inputFileName << endl
         << "  -datasetPathPhi " << datasetPathPhi << std::endl
         << "  -datasetPathT " << datasetPathT << std::endl
         << "  -datasetPathConc " << datasetPathConc << std::endl
         << "-Nt " << Nt << endl
         << "-Nv " << Nv << endl
         << "  -mesh-size " << hh_x << std::endl
         << "-dt " << dt << endl
         << "-wp " << write_period << endl
         // TODO: correct the energy parameters
         << "  -order-energy " << orderEnergy << std::endl
         << "  -c1 " << c1 << std::endl
         << "  -c2 " << c2 << std::endl
         << "  -c3 " << c3 << std::endl
         << "  -c4 " << c4 << std::endl
         << "  -c5 " << c5 << std::endl
         << "  -c6 " << c6 << std::endl
         << "  -cPrefactor  " << cPrefactor << std::endl
         << "  -alpha  " << alpha << std::endl
         << "  -T0 " << T0 << std::endl
         << "  -cbase " << cbase << std::endl
         << "  -LL " << LL << std::endl
         << "  -D_T " << D_T << std::endl
         << "  -M_phi " << M_phi << std::endl
         << "  -M_conc " << M_conc << std::endl
         << "  -dt " << dt << std::endl
         << "  -wp " << write_period << std::endl
         << "  -thermal-diffusivity " << D_T << std::endl
         << endl;
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
   // open a file to write the fields' evolutions to
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
      cout << "opening file: " << output_prefix + "_stats.txt" 
        << endl; // debug 
      stat_file.open(output_prefix + "_stats.txt", 
                        ios::app | ios::ate);

      if ( ! stat_file.good() )
      {
         cout << "Error: could not open output file: "
            << output_prefix + "_stats.txt"
            << endl;
         stat_file.close();
      }

      // write a heading to stat_file
      if (stat_file.good())
      {
         stat_file << setw(17) << "timeStep" << " "
                  << setw(17) << "time" << " "
                  << setw(17) << "totalEnergySystem" << " "
                  << setw(17) << "totalEnergySum" << " "
                  << setw(17) << "totalEnergyMean" << " "
                  << setw(17) << "totalEnergyVar" << " "
                  << setw(17) << "freeEnergySum" << " "
                  << setw(17) << "freeEnergyMean" << " "
                  << setw(17) << "freeEnergyVar" << " "
                  << setw(17) << "EntropySum" << " "
                  << setw(17) << "EntropyMean" << " "
                  << setw(17) << "phiMean" << " "
                  << setw(17) << "phiVar" << " "
                  << setw(17) << "TMean" << " "
                  << setw(17) << "TVar" << " "
                  << setw(17) << "concMean" << " "
                  << setw(17) << "concVar"
                  << endl;
      }
   }
   ////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////
   // instantiate containers of local field changes
   std::vector<double> // (0:x-, 1:x+, 2:y-, 3:y+, 4:z-, 5:z+)
      conc_local_flux( Nvoxel_neighbors * conc_local.size(), 0);//[n,i,j,k]

   std::vector<int> conc_local_flux_sum( conc_local.size(), 0);//[i,j,k]

   std::vector<double> // (0:x-, 1:x+, 2:y-, 3:y+, 4:z-, 5:z+)
      conc_local_rates( Nvoxel_neighbors * conc_local.size() );//[n,i,j,k]
   ////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////
   // instantiate reusable variables
   
   // fluxes of conserved quantities, assuming each element has only a
   // single source voxel 
   // TODO: this would change if allowing diagonal travel
   std::vector<double> conc_flux_upward( Ny * Nz, 0);
   std::vector<double> conc_flux_downward( Ny * Nz, 0);
   std::vector<double> conc_flux_from_above( Ny * Nz, 0);
   std::vector<double> conc_flux_from_below( Ny * Nz, 0);

   // flux rates of first passage to determine voxel filling order
   std::vector<double> conc_flux_upward_rates( Ny * Nz, 0);
   std::vector<double> conc_flux_downward_rates( Ny * Nz, 0);
   std::vector<double> conc_flux_from_above_rates( Ny * Nz, 0);
   std::vector<double> conc_flux_from_below_rates( Ny * Nz, 0);


   std::vector<size_t> neigh_idxs(Nvoxel_neighbors, 0);

   //std::vector<double> jump_rates(6,0); // overwritten at every voxel
   // (0:x-, 1:x+, 2:y-, 3:y+, 4:z-, 5:z+)

   std::vector<size_t> neigh_order(Nvoxel_neighbors, 0);
   std::uniform_real_distribution<double> rand_decimal(0,1);// for order
   std::vector<double> rand_decimals1(Nvoxel_neighbors, 0);
   std::vector<double> rand_decimals2(Nvoxel_neighbors, 0);
   
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

   double lap1dConc;  // 1-D Laplacian, reused each timestep
   std::vector<double> hh( 3, 0); hh[0] = hh_x; hh[1] = hh_y; hh[2] = hh_z;

   // variables for statistical analysis
   double entropy_local_sum; entropy_local_sum = 0;
   //double entropy_local_sqr_sum; entropy_local_sqr_sum = 0;
   double entropy_total_sum; entropy_total_sum = 0;
   //double entropy_total_sqr_sum; entropy_total_sqr_sum = 0;
   double entropy_total_mean; entropy_total_mean = 0;
   double total_energy_system; total_energy_system = 0;
   //double entropy_variance; entropy_variance = 0;
   double total_energy_local_sum; total_energy_local_sum = 0;
   double total_energy_local_sqr_sum; total_energy_local_sqr_sum = 0;
   double total_energy_total_sum; total_energy_total_sum = 0;
   double total_energy_total_sqr_sum; total_energy_total_sqr_sum = 0;
   double total_energy_mean; total_energy_mean = 0;
   double total_energy_variance; total_energy_variance = 0;
   double free_energy_local_sum; free_energy_local_sum = 0;
   double free_energy_local_sqr_sum; free_energy_local_sqr_sum = 0;
   double free_energy_total_sum; free_energy_total_sum = 0;
   double free_energy_total_sqr_sum; free_energy_total_sqr_sum = 0;
   double free_energy_mean; free_energy_mean = 0;
   double free_energy_variance; free_energy_variance = 0;
   double phi_local_sum; phi_local_sum = 0;
   double phi_local_sqr_sum; phi_local_sqr_sum = 0;
   double phi_total_sum; phi_total_sum = 0;
   double phi_total_sqr_sum; phi_total_sqr_sum = 0;
   double phi_mean; phi_mean = 0;
   double phi_variance; phi_variance = 0;
   double T_local_sum; T_local_sum = 0;
   double T_local_sqr_sum; T_local_sqr_sum = 0;
   double T_total_sum; T_total_sum = 0;
   double T_total_sqr_sum; T_total_sqr_sum = 0;
   double T_mean; T_mean = 0;
   double T_variance; T_variance = 0;
   double conc_local_sum; conc_local_sum = 0;
   double conc_local_sqr_sum; conc_local_sqr_sum = 0;
   double conc_total_sum; conc_total_sum = 0;
   double conc_total_sqr_sum; conc_total_sqr_sum = 0;
   double conc_mean; conc_mean = 0;
   double conc_variance; conc_variance = 0;
   
   double conc_lower_limit_walkers = Nv * conc_lower_limit;
   double conc_upper_limit_walkers = Nv * conc_upper_limit;

   // instantiate free_energy_gl
   SPF_NS::free_energy_gl F_gl(
            c1, //100, // c1
            c2, //0.25, // c2
            c3, //7, // c3
            c4, //3, // c4
            c5, //0, // c5
            c6, //0, // c6
            cPrefactor, //-0.05, // cPrefactor, used in C_{P}(T)
            alpha, //-110, // alpha
            T0, //10, // T0
            cbase, //2, // cbase
            LL, //1, // L // TODO: find a good value for latent heat
            orderEnergy, //0.1, // order energy (multiplies lap3dPhi)
            M_phi/hh_x, //1/hh_x, // M_phi
            M_phi/hh_y, //1/hh_y, // M_phi
            M_phi/hh_z, //1/hh_z, // M_phi
            M_conc/hh_x, //1/hh_x, // M_conc
            M_conc/hh_y, //1/hh_y, // M_conc
            M_conc/hh_z, //1/hh_z, // M_conc
            D_T,
            0.0001 // rk_dt // TODO: make into an input parameter
         );

   std::list<double> deltaPhiSequence; // to be reused at each voxel
   ////////////////////////////////////////////////////////////////////

   /*-----------------------------------------------------------------*/
   /* begin loop over time -------------------------------------------*/
   /*-----------------------------------------------------------------*/
   for (time_step = 0; time_step <= Nt; ++time_step)
   {

      //double flux_total; flux_total = 0.0;// debug
      //////////////////////////////////////////////////////////////////
      // append fields to a file
      if ( time_step % write_period == 0 )
      {
         // convert conc_local from # walkers to concentration
         //  (sacrificing cpu to preserve memory)
         for (size_t ii=1; ii < Nx_local +1; ++ii) // loop over non-ghosts
            for ( size_t jj=0; jj < Ny; ++jj)
               for ( size_t kk=0; kk < Nz; ++kk)
               {
                  conc_local[ kk + Nz*(jj + Ny*ii)] /= Nv;
               }

         // write the local subset of the fields to the file
         if ( append_fields_to_hdf5_multinode(
                  outFile_id,
                  time_step,
                  time,
                  phi_local,
                  T_local,
                  conc_local,
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

         // revert conc_local from concentration to # walkers 
         for (size_t ii=1; ii < Nx_local +1; ++ii) // loop over non-ghosts
            for ( size_t jj=0; jj < Ny; ++jj)
               for ( size_t kk=0; kk < Nz; ++kk)
               {
                  conc_local[ kk + Nz*(jj + Ny*ii)] =
                     round(Nv * conc_local[ kk + Nz*(jj + Ny*ii)]);
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
      update_ghosts(
                     T_local,
                     Nx_local,
                     Ny,
                     Nz,
                     neighbor_x_higher,
                     neighbor_x_lower, 
                     neighbors_comm
                     );
      update_ghosts(
                     conc_local,
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
      // receive inward conc flux from neighboring nodes
      flux_exchange_irecv(
            conc_flux_from_above,
            conc_flux_from_above_rates,
            conc_flux_from_below,
            conc_flux_from_below_rates,
            Ny, Nz,
            neighbor_x_higher, neighbor_x_lower, 
            halo_flux_recv_requests, neighbors_comm
            );

      // NOTE: phi (phase) and T (temperature) are not conserved
      //  quantities and don't exchange particles like concentration

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
      size_t idx; 


      // TODO: implement the following SDEs
      // d\phi_{i}(t) = \mu_{i}(t)dt + \sqrt(\mu_{i}(t)dt) \circ dB(t)
      // with \mu_{i} = -M_{\phi} \frac{\delta F}{\delta \phi}

      // dT_{i}(t) = D_{T_{i}(t)} \nabla^{2}T_{i}(t)dt + \frac{L}{c_{p}(T_{i}(t))} \diamond d\phi_{i}(t)

      // dc_{i \rightarrow j}(t) = N_{\lambda_{i}}(\Delta t, 1)
      // with rate \lambda_{i} = \nabla M_{c} \nabla \frac{\delta F}{\delta c}

      // Heat capacity
      // c_{p}(T,\alpha,\beta,T0,cbase) = -0.05 \frac{\alpha^{2}}{T^{2} K_{eq}(T,\alpha,T0)} (0.5-\frac{2+ \frac{1}{K_{eq}(T,\alpha,T0)}}{2 \sqrt{(2+\frac{1}{K_{eq}(T,\alpha,T0)})^{2}-4}} + cbase
      // K_{eq}(T,\alpha,T0) = \exp(\frac{-1*alpha}{T} + T0)

      // Free energy
      // F( x, c1, c2, phi, gradphi, T, a, b, c, d)
      //  = a*((1-phi)^{2}*(x-(1-0.25*(erf((T-7)/3)+1)))^{2} + (b*(0.25*(erf((T-7)/3)+1) -x)^{2}*phi^{2} + c*phi^{4} + d*gradphi^{2}

      // parameters used in plots:
      // a=100,b=100,c=0,d=1,c1=1,c2=0,phi={0,0.5,1.0},gradphi=0,T=11
      // a=100,b=100,c=0,d=1,c1=1,c2=0,phi={0,0.5,1.0},gradphi=0,T=9
      // alpha=-110
      // beta=1
      // T0=10
      // cbase=2
      // conc \in [0,1]
      // T \in [0,21]


      // TODO: turn the following loops over T, phi, & conc into functions of classes

      double lap3dPhi = 0.0;
      /****************************************************************/
      /* loop over voxels in phi_local, skipping ghosts ***************/
      for (size_t ii=1; ii < Nx_local +1; ++ii) // loop over non-ghosts
         for ( size_t jj=0; jj < Ny; ++jj)
            for ( size_t kk=0; kk < Nz; ++kk)
            {
               // indices wrt local field
               //size_t idx; 
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

               lap3dPhi = laplacian3d(
                                       hh_x, hh_y, hh_z,
                                       Nz, // to test if 2-D or 3-D
                                       phi_local, //input field,
                                       idx,
                                       neigh_idxs[neigh_pairs[0]], //x+
                                       neigh_idxs[neigh_pairs[1]], //x-
                                       neigh_idxs[neigh_pairs[2]], //y+
                                       neigh_idxs[neigh_pairs[3]], //y-
                                       neigh_idxs[neigh_pairs[4]], //z+
                                       neigh_idxs[neigh_pairs[5]]  //z-
                                       );

               //deltaPhi[ idx] = F_gl.delta_phi_noise_driven(
               if ( F_gl.delta_phi_marcus(
                                       deltaPhiSequence, // output
                                       conc_local_flux_sum[ idx],
                                       conc_local[ idx],
                                       Nv,
                                       phi_local[ idx],
                                       lap3dPhi,
                                       T_local[ idx],
                                       hh)
                     != EXIT_SUCCESS)
               {
                  std::cout << "error: F_gl.delta_phi_marcus() failed"
                     << std::endl;
                  return EXIT_FAILURE;
               }
               deltaPhi[ idx] = 0.0;
               for ( std::list<double>::const_iterator
                     itr = deltaPhiSequence.begin();
                     itr != deltaPhiSequence.end();
                     ++itr
                     )
               {
                  deltaPhi[ idx] += *itr;
               }
               //deltaPhi[ idx] = F_gl.delta_phi_noiseless(
               //                        deltaPhis[ idx],
               //                        conc_local[ idx],
               //                        Nv,
               //                        phi_local[ idx],
               //                        lap3dPhi,
               //                        T_local[ idx],
               //                        hh);
               //if (( deltaPhi[idx] + phi_local[idx] < 0) 
               //   || ( deltaPhi[idx] + phi_local[idx] > 1))
               //{
               //   std::cout << "problem:" << std::endl; // debug
               //   std::cout << "phi_local[" << idx << "]: " << phi_local[idx] << std::endl;// debug
               //   std::cout << "lap3dPhi" << idx << ": " << lap3dPhi << std::endl;// debug
               //   std::cout << "deltaPhi[" << idx << "]: " << deltaPhi[idx] << std::endl;// debug
               //}
               // TODO: add noise \sqrt{mu} \circ dB
      //      }

      ///* end loop over voxels in phi_local, skipping ghosts ***********/
      ///****************************************************************/
      ///****************************************************************/
      ///* loop over voxels in T_local, skipping ghosts *****************/

      //for (size_t ii=1; ii < Nx_local +1; ++ii) // loop over non-ghosts
      //   for ( size_t jj=0; jj < Ny; ++jj)
      //      for ( size_t kk=0; kk < Nz; ++kk)
      //      {
      //         // indices wrt local field
      //         //size_t idx; 
      //         idx = kk + Nz*(jj + Ny*ii);
      //         //std::cout << "T_local[" << idx << "]: " // debug
      //         //   << T_local[idx] << std::endl;// debug

      //         identify_local_neighbors(
      //               neigh_idxs[0], 
      //               neigh_idxs[1], 
      //               neigh_idxs[2], 
      //               neigh_idxs[3], 
      //               neigh_idxs[4],
      //               neigh_idxs[5],
      //               ii, jj, kk,
      //               Ny, Nz
      //               );

               /* Implementation for using 3D Laplacian */
               // No need to use 1-D, since 1-D components would be summed
               double lap3d;  // 3-D Laplacian
               lap3d = laplacian3d( hh_x, hh_y, hh_z,
                                    Nz,
                                    T_local,
                                    idx,
                                    neigh_idxs[neigh_pairs[0]], //x+
                                    neigh_idxs[neigh_pairs[1]], //x-
                                    neigh_idxs[neigh_pairs[2]], //y+
                                    neigh_idxs[neigh_pairs[3]], //y-
                                    neigh_idxs[neigh_pairs[4]], //z+
                                    neigh_idxs[neigh_pairs[5]]  //z-
                        );


               //// heat equation
               ////  \frac{\partial T}{\partial t} = D_{T} \nabla^{2} T
               //laplacian_flux(
               //      T_local_change, 
               //      T_local, 
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
               
               //deltaT[ idx] = F_gl.delta_T_decoupled(
               //deltaT[ idx] = F_gl.delta_T(
               //                      T_local[idx],
               //                      lap3d,
               //                      deltaPhi[idx], // noise
               //                      dt);
               deltaT[ idx] =
                  F_gl.delta_T_marcus(
                     //latent_heat_c_p_inverse, 
                     T_local[ idx],
                     lap3d,
                     dt,
                     phi_local[ idx],
                     deltaPhiSequence // std::list<double> //deltaPhi, // TODO
                  );
               // output // TODO
               //// debug
               //if ( abs(deltaT[idx] ) > 0.00001)
               //{
               //   std::cout << "deltaT[" << idx << "]: " << deltaT[idx] << std::endl; // debug
               //   std::cout << "T_local[" << idx << "]: " << T_local[idx] << std::endl; // debug
               //   //std::cout << "lap3d " << idx << ": " << lap3d << std::endl; // debug
               //   std::cout << "deltaT[" << idx << "]: " << deltaT[idx] << std::endl; // debug
               //   //std::cout << "deltaPhi[" << idx << "]: " << deltaPhi[idx] << std::endl; // debug
               //                                                                           }
               //// end debug
            }
      /* end loop over voxels in T_local, skipping ghosts *************/
      /****************************************************************/

      /****************************************************************/
      // evaluate dFdconc at all cells (ghosts too) so that the
      //        laplacian of dFdconc may be evaluated by delta_conc()
      for (size_t ii=0; ii < Nx_local +2; ++ii)
         for ( size_t jj=0; jj < Ny; ++jj)
            for ( size_t kk=0; kk < Nz; ++kk)
            {
               idx = kk + Nz*(jj + Ny*ii);
               // update dFdconc field
               dFdconc[ idx] = F_gl.dF_dconc(
                                     conc_local[ idx], // in # walkers
                                     Nv, // for conversion to/from walkers
                                     phi_local[ idx], // \in [0,1]
                                     T_local[ idx]); // \in [0,\inf]
            }
      /****************************************************************/

      /****************************************************************/
      /* loop over voxels in conc_local, skipping ghosts **************/

      // Local field doesn't change until after evaluating flux for
      //  every voxel.

      for (size_t ii=1; ii < Nx_local +1; ++ii) // loop over non-ghosts
         for ( size_t jj=0; jj < Ny; ++jj)
            for ( size_t kk=0; kk < Nz; ++kk)
            {
               // indices wrt local field
               //size_t idx; 
               idx = kk + Nz*(jj + Ny*ii);
               //std::cout << "conc_local[" << idx << "]: " // debug
               //   << conc_local[idx] << std::endl;// debug

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

               //// prepare to assign values to jump rates
               ////  lambda = -M_conc * \nabla^2 dF/dconc
               ///* Implementation for using 3D Laplacian */
               //double lap3d;  // 3-D Laplacian
               //lap3d = laplacian3d( hh_x, hh_y, hh_z,
               //                     Nz, // to test if 2-D or 3-D
               //                     dFdconc, //input, \nabla^2 dFdconc
               //                     idx,
               //                     neigh_idxs[neigh_pairs[0]], //x+
               //                     neigh_idxs[neigh_pairs[1]], //x-
               //                     neigh_idxs[neigh_pairs[2]], //y+
               //                     neigh_idxs[neigh_pairs[3]], //y-
               //                     neigh_idxs[neigh_pairs[4]], //z+
               //                     neigh_idxs[neigh_pairs[5]]  //z-
               //         );

               // assign values to jump rates
               //  lambda = -M_conc * \nabla^2 dF/dconc

               /* Implementation for using 1D 2nd derivative*/
               // determine conc_local_rates
               // prepare to assign values to jump rates
               for( size_t nn=0; nn < 0.5*Nvoxel_neighbors; ++nn)
               {
                  //std::cout << "conc_local[" << idx << "]/Nv : "
                  //   << conc_local[idx]/Nv << std::endl;// debug
                  if ( conc_local[idx] > 0)
                  {
                     lap1dConc = laplacian1d(
                              hh[nn],
                              dFdconc,
                              idx,
                              neigh_idxs[neigh_pairs[2*nn+1]],//neigh below
                              neigh_idxs[neigh_pairs[2*nn]] //neigh above
                              );// * Nv; // don't convert to walkers yet

                     // downward rate
                     //  lambda = M_conc * \nabla^2 dF/dconc
                     conc_local_rates[ 2*nn +Nvoxel_neighbors*idx]
                        =  F_gl.jump_rate_conc( lap1dConc, hh[nn], nn);

                     // convert to # walkers from concentration
                     conc_local_rates[2*nn + Nvoxel_neighbors*idx] *= Nv;
                     // NOTE: is this conversion necessary?
                     //       conc_local_rates is used in sampling the 
                     //       Poisson distribution, which expects input 
                     //       rate to be in [0,\inf) and returns integers.
                     // Jump process rate doesn't need to be an integer. ?

                     // copy to upward rate
                     conc_local_rates[ (2*nn+1) + Nvoxel_neighbors*idx]
                             = conc_local_rates[
                                 2*nn + Nvoxel_neighbors*idx];
                     // debug
                     //std::cout << "node " << mynode 
                     //   << " conc_local_rates[" 
                     //   <<  nn + Nvoxel_neighbors*idx 
                     //   << "] "
                     //   << conc_local_rates[ nn + Nvoxel_neighbors*idx ]
                     //   << ", conc_local[" << idx << "] "
                     //   << conc_local[idx]
                     //   << ", conc_local[" << idx << "]/" << Nv << " "
                     //   << conc_local[idx ]/Nv
                     //   << std::endl; 
                     // end debug
                  }
                  else
                  {
                     conc_local_rates[nn + Nvoxel_neighbors*idx] = 0;
                  }
               }

               // evaluate stochastic fluxes (walkers) to neighbor cells
               conserved_jump_flux_separate_distributions( 
                     conc_local_flux, // in units of number of walkers
                     rr,
                     conc_local,
                     conc_local_rates, // in units of number of walkers
                     dt,
                     Nvoxel_neighbors,
                     neigh_idxs,
                     conc_upper_limit_walkers,
                     conc_lower_limit_walkers,
                     idx
                     );

               //if ( flags.debug != 0)
               //{
               //   for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
               //   {
               //      std::cout 
               //         << time_step
               //         << " conc_local[" << idx << "]: "
               //         << conc_local[ idx] << ", "
               //         << " conc_local_flux["
               //         << nn + Nvoxel_neighbors*idx << "]: "
               //         << conc_local_flux[ nn + Nvoxel_neighbors*idx]
               //         << ", conc_local_rates[" 
               //         << nn +Nvoxel_neighbors*idx
               //         << "] : " 
               //         << conc_local_rates[ nn +Nvoxel_neighbors*idx] 
               //         << std::endl;
               //      //flux_total 
               //      //   += conc_local_flux[ nn + Nvoxel_neighbors*idx]/Nv;
               //   }
               //}
            }

      //// debug
      //std::cout << time_step << "flux_total : " << flux_total 
      //   << std::endl;
      //// end debug
      /* end loop over voxels in conc_local, skipping ghosts **********/
      /****************************************************************/
      
      /* enforce bounds of conc given the new change in conc **********/
      enforce_bounds_int_outward(
            // updates conc_local_flux with acceptable flux values
            conc_local_flux,
            conc_local,
            conc_local_rates,
            rr,
            // neigh_order,
            Nvoxel_neighbors,
            conc_lower_limit_walkers,
            conc_upper_limit_walkers,
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

            conc_flux_downward[  kk + Nz*jj] 
               = conc_local_flux[ 0+ Nvoxel_neighbors * idx ]; 
            // (0:x-, 1:x+, 2:y-, 3:y+, 4:z-, 5:z+)

            conc_flux_downward_rates[kk + Nz*jj] 
               = conc_local_rates[0 + Nvoxel_neighbors * idx];

            iii = Nx_local; // upper x-axis boundary of non-ghosts
            idx = kk + Nz*(jj + Ny*iii);

            conc_flux_upward[ kk + Nz*jj] 
               = conc_local_flux[ 1+ Nvoxel_neighbors * idx];

            conc_flux_upward_rates[kk + Nz*jj] 
               = conc_local_rates[1 + Nvoxel_neighbors * idx];
         }
      
      //////////////////////////////////////////////////////////////////

      //for (size_t jj=0; jj < Ny; ++jj) // debug
      //{ // debug
      //   cout << "node " << mynode << " sending upward[ "; // debug
      //   for (size_t kk=0; kk < Nz; ++kk) // debug
      //   { // debug
      //  cout << setw(5) << conc_flux_upward[kk + Nz*jj] << ", "; // debug
      //   } // debug
      //   cout << "]" << endl; // debug
      //} // debug
      //for (size_t jj=0; jj < Ny; ++jj) // debug
      //{ // debug
      //   cout << "node " << mynode << " sending downward [ "; // debug
      //   for (size_t kk=0; kk < Nz; ++kk) // debug
      //   { // debug
      //  //cout << setw(5) << conc_flux_upward[kk+Nz*jj] << ", "// debug
      //  cout << setw(5) << conc_flux_downward[kk+Nz*jj] << ", ";// debug
      //   } // debug
      //   cout << "]" << endl; // debug
      //} // debug
      //////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////
      // send locally acceptable outward flux of conc to neighboring nodes
      flux_exchange_isend(
         conc_flux_upward, // Ny*Nz
         conc_flux_upward_rates, // Ny*Nz
         conc_flux_downward, // Ny*Nz
         conc_flux_downward_rates, // Ny*Nz
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
      //       << conc_flux_from_above[kk + Nz*jj] << ", "; // debug
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
      //         cout << setw(5) << conc_local[ k + Nz*(j + i*Ny) ]
      //            << " "; // debug
      //      cout << "] " << endl; // debug
      //   }
      //   cout << endl; // debug
      //} // debug
      //cout << endl;// debug

      // copy flux from other nodes into conc_local_flux
      for (size_t jj=0; jj < Ny; ++jj)
         for( size_t kk=0; kk < Nz; ++kk)
         {
            // flux to lower neighbor of upper ghost along x-axis
            // ghost x_idx = Nx_local+1; ghost neigh_idx = 0 (x-)
            idx = 0 + Nvoxel_neighbors* (kk + Nz*(jj + Ny*(Nx_local+1)));
            conc_local_flux[ idx ] 
               = conc_flux_from_above[kk + Nz*jj];

            conc_local_rates[ idx ]
               = conc_flux_from_above_rates[kk + Nz*jj];

            // flux to upper neighbor of lower ghost along x-axis
            // ghost x_idx = 0; ghost neigh_idx = 1 (x+)
            idx = 1 + Nvoxel_neighbors* (kk + Nz*(jj + Ny*0));
            conc_local_flux[ idx ]
               = conc_flux_from_below[kk + Nz*jj];

            conc_local_rates[ idx ]
               = conc_flux_from_below_rates[kk + Nz*jj];
         }

      // enforce_field_bounds( field, flux )
      //        using conc_local_flux[] for flux magnitudes
      //        and conc_local_rates[] for balancing those magnitudes
      //        saving acceptable fluxes in conc_local_flux
      // Considered as an outward flux, neighbor orders are 
      //  equally balanced (assuming equal barrier heights).
      // But when considered as an inward flux, the neighbor
      //  orders may be balanced using first passage 
      //  distributions.

      // TODO: try the following for conc fluxes
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
            // updates conc_local_flux with acceptable flux values
            conc_local_flux,   // integers
            conc_local,        // integers
            conc_local_rates,  // not necessarily integers
            rr,
            //neigh_order, // not implemented
            Nvoxel_neighbors,
            conc_lower_limit_walkers,
            conc_upper_limit_walkers,
            Nx_local, Ny, Nz,
            eps,
            flags
            );

      // Update conc_flux_from_below / above with accepted inward fluxes 
      //        from ghosts residing in conc_local_flux
      for (size_t jj=0; jj < Ny; ++jj)
         for (size_t kk=0; kk < Nz; ++kk)
         {
            idx = 0 + Nvoxel_neighbors* (kk + Nz*(jj + Ny*(Nx_local+1)));
            conc_flux_from_above[kk + Nz*jj]
               = conc_local_flux[ idx ];

            idx = 1 + Nvoxel_neighbors* (kk + Nz*(jj + Ny*0));
            conc_flux_from_below[kk + Nz*jj]
               = conc_local_flux[ idx ];
         }

      MPI_Waitall(4, halo_flux_send_requests, MPI_STATUSES_IGNORE);

      flux_accepted_irecv(
            conc_flux_downward,
            conc_flux_upward,
            Ny,
            Nz,
            neighbor_x_higher,
            neighbor_x_lower,
            halo_accepted_flux_recv_requests,
            neighbors_comm
            );

      // send accepted flux values back to sources
      flux_accepted_isend(
            conc_flux_from_above, // Ny*Nz
            conc_flux_from_below, // Ny*Nz
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
            if ( conc_flux_downward[kk + Nz*jj] < conc_local_flux[ idx])
            {
               conc_local_flux[ idx ] = conc_flux_downward[kk + Nz*jj];
            }

            idx = 1 + Nvoxel_neighbors* (kk + Nz*(jj + Ny*(Nx_local)));
            if ( conc_flux_upward[kk + Nz*jj] < conc_local_flux[ idx])
            {
               conc_local_flux[ idx ]
                  = conc_flux_upward[kk + Nz*jj];
            }
            // After this loop, both local voxels and neighbor's ghosts
            //  should contain the smaller of the fluxes determined by
            //  both nodes, to prevent both over filling and 
            //  over drawing.
         }

      // Update conc_local by applying the accepted fluxes
      for (size_t ii=1; ii < Nx_local +1; ++ii) // loop over non-ghosts
         for ( size_t jj=0; jj < Ny; ++jj)
            for ( size_t kk=0; kk < Nz; ++kk)
            {
               idx = kk + Nz*(jj + Ny*ii);

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

               conc_local_flux_sum[idx] = 0; // reset voxel flux sum

               // Subtract outward flux from each voxel (conc_local)
               for(size_t nn=0; nn < Nvoxel_neighbors; ++nn)
               {
                  conc_local[idx] 
                     -= conc_local_flux[nn + Nvoxel_neighbors*idx];

                  conc_local_flux_sum[idx]
                     -= conc_local_flux[nn + Nvoxel_neighbors*idx];
               }

               // Add inward flux to conc_local and conc_local_flux_sum
               for ( size_t nn=0; nn < Nvoxel_neighbors; ++nn)
               {
                  conc_local[idx]
                     += conc_local_flux[
                           neigh_pairs[nn] 
                              + Nvoxel_neighbors * neigh_idxs[nn]
                        ];

                  conc_local_flux_sum[idx]
                     += conc_local_flux[
                           neigh_pairs[nn] 
                              + Nvoxel_neighbors * neigh_idxs[nn]
                        ];
               }

               // check to ensure enforce_bounds_generic worked
               if ( conc_local[idx] < conc_lower_limit_walkers )
               {
                  //if ( abs(conc_local[idx] - conc_lower_limit - outward_flux)
                  //      > (eps.dblsqrt 
                  //            * (conc_local[idx] - conc_lower_limit)))
                  if ( flags.debug != 0 )
                  {
                     std::cout << "node " << mynode << " Warning: step " 
                         << time_step 
                         << " flux out of voxel caused lower bound breach"
                         << " conc_local[" << idx << "]: "
                         << conc_local[idx] 
                         << ", (conc_local[]-conc_lower_limit_walkers) "
                         << (conc_local[idx]-conc_lower_limit_walkers)
                         << ", setting conc_local to its lower limit"
                         << ", conc_local_flux[nn + Nvoxel_neighbors*" 
                         << idx << "]: ";
                     for (size_t mm=0; mm < Nvoxel_neighbors; ++mm)
                     {
                        std::cout 
                           << conc_local_flux[mm + Nvoxel_neighbors*idx]
                           << ", ";
                     }
                     std::cout << std::endl;
                  }
                  flags.fail = 1;

                  conc_local_flux_sum[idx] // subtract excess flux
                     -= conc_local[idx] - conc_lower_limit_walkers;

                  conc_local[idx] = conc_lower_limit_walkers;
               }

               if ( conc_local[idx] > conc_upper_limit_walkers )
               {
                  if ( flags.debug != 0 )
                  {
                     std::cout << "node " << mynode << " Warning: step " 
                        << time_step 
                        << " flux into voxel caused upper bound breach"
                        << " conc_local[" << idx << "]: "
                        << conc_local[idx] 
                        // << " setting conc_local to its upper limit"
                        << std::endl;
                  }
                  conc_local_flux_sum[idx] // subtract excess from flux sum
                     -= conc_local[idx] - conc_upper_limit_walkers;
                  conc_local[idx] = conc_upper_limit_walkers;
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
      //            cout << setw(9) << conc_local[ k + Nz*(j + i*Ny) ];
      //         cout << "] " << endl;
      //      }
      //      cout << endl;
      //   }
      //   cout << endl;
      //}
      /* end loop over voxels in conc_local, skipping ghosts **********/
      /****************************************************************/

      /****************************************************************/
      // Update non-conservative fields with corresponding deltas
      for (size_t ii=1; ii < Nx_local +1; ++ii) // loop over non-ghosts
         for ( size_t jj=0; jj < Ny; ++jj)
            for ( size_t kk=0; kk < Nz; ++kk)
            {
               idx = kk + Nz*(jj + Ny*ii);

               // Update phi_local by applying deltaPhi
               phi_local[idx] += deltaPhi[idx];

               // enforce lower bound of phase
               if ( phi_local[idx] < phi_lower_limit )
               {
                  //if ( abs(phi_local[idx] - phi_lower_limit)
                  //      > (eps.dblsqrt 
                  //            * (phi_local[idx] - phi_lower_limit)))
                  if ( flags.debug != 0 )
                  {
                     std::cout << "node " << mynode << " Warning: step " 
                        << time_step 
                        << " lower bound breach of "
                        << " phi_local[" << idx << "]: "
                        << phi_local[idx] 
                        << ", reverted phi_local[" << idx << "] value:"
                        << phi_local[idx] - deltaPhi[idx]
                        << ", deltaPhi[" << idx << "]: " << deltaPhi[idx]
                        << ", (phi_local[]-phi_lower_limit) "
                        << (phi_local[idx]-phi_lower_limit)
                        << ", setting phi_local to its lower limit"
                        << std::endl;
                  }
                  flags.fail = 1;

                  phi_local[idx] = phi_lower_limit;
               }
               // enforce upper bound of phase
               if ( phi_local[idx] > phi_upper_limit )
               {
                  //if ( abs(phi_local[idx] - phi_upper_limit)
                  //      > (eps.dblsqrt 
                  //            * (phi_local[idx] - phi_upper_limit)))
                  if ( flags.debug != 0 )
                  {
                     std::cout << "node " << mynode << " Warning: step " 
                        << time_step 
                        << " upper bound breach of "
                        << " phi_local[" << idx << "]: "
                        << phi_local[idx] 
                        << ", reverted phi_local[" << idx << "] value:"
                        << phi_local[idx] - deltaPhi[idx]
                        << ", deltaPhi[" << idx << "]: " << deltaPhi[idx]
                        << ", (phi_local[]-phi_upper_limit) "
                        << (phi_local[idx]-phi_upper_limit)
                        << ", setting phi_local to its upper limit"
                        << std::endl;
                  }
                  flags.fail = 1;

                  phi_local[idx] = phi_upper_limit;
               }

               // Update T_local by applying deltaT
               T_local[idx] += deltaT[idx];

               // enforce lower bound of temperature
               if ( T_local[idx] < T_lower_limit )
               {
                  //if ( abs(T_local[idx] - t_lower_limit)
                  //      > (eps.dblsqrt 
                  //            * (T_local[idx] - t_lower_limit)))
                  if ( flags.debug != 0 )
                  {
                     std::cout << "node " << mynode << " warning: step " 
                        << time_step 
                        << " lower bound breach of "
                        << " T_local[" << idx << "]: "
                        << T_local[idx] 
                        << ", reverted T_local[" << idx << "] value:"
                        << T_local[idx] - deltaT[idx]
                        << ", deltaT[" << idx << "]: " << deltaT[idx]
                        << ", (T_local[]-T_lower_limit) "
                        << (T_local[idx]-T_lower_limit)
                        << ", setting T_local to its lower limit"
                        << std::endl;
                  }
                  flags.fail = 1;

                  T_local[idx] = T_lower_limit;
               }

               // now that T, phi, and conc are up to date, evaluate energy
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
               // update lap1dPhiX, lap1dPhiY, lap1dPhiZ;
               lap1dPhiX = laplacian1d(
                                       hh[0],
                                       phi_local,
                                       idx,
                                       neigh_idxs[neigh_pairs[2*0+1]],
                                       neigh_idxs[neigh_pairs[2*0]]
                        );
               lap1dPhiY = laplacian1d(
                                       hh[1],
                                       phi_local,
                                       idx,
                                       neigh_idxs[neigh_pairs[2*1+1]],
                                       neigh_idxs[neigh_pairs[2*1]]
                        );
               lap1dPhiZ = laplacian1d(
                                       hh[2],
                                       phi_local,
                                       idx,
                                       neigh_idxs[neigh_pairs[2*2+1]],
                                       neigh_idxs[neigh_pairs[2*2]]
                        );
               free_energy_local[ idx] = F_gl.F(
                                       conc_local[ idx], // units: walkers
                                       Nv,
                                       phi_local[ idx],
                                       T_local[ idx],
                                       lap1dPhiX, lap1dPhiY, lap1dPhiZ
                                       );
            }

      /* end loop over voxels in non-conservative fields, no-ghosts ***/
      /****************************************************************/

      /////////////////////////////////////////////////////////////////
      // calculate mean and variance of all voxel populations
      double total_energy_temp; total_energy_temp = 0;
      if ( flags.calcstat != 0 )
      {
         // free_energy_local_sum = 0.0;// redundant, performed later
         for (size_t ii=1; ii < Nx_local+1; ++ii)
         {
            for (size_t jj=0; jj < Ny; ++jj)
            {
               for (size_t kk=0; kk < Nz; ++kk)
               {
                  idx = kk + Nz*(jj + Ny*(ii));
                  // sum local populations
                  entropy_local_sum
                     += F_gl.entropy(
                                 conc_local[ idx]/Nv,
                                 phi_local[ idx]
                              );
                  total_energy_temp 
                     =  free_energy_local[ idx]
                        + (
                           T_local[ idx]
                           * F_gl.entropy(
                                 conc_local[ idx]/Nv,
                                 phi_local[ idx]
                              )
                          );
                  total_energy_local_sum += total_energy_temp;
                  free_energy_local_sum += free_energy_local[ idx];
                  phi_local_sum += phi_local[ idx];
                  T_local_sum += T_local[ idx];
                  conc_local_sum += conc_local[ idx]/Nv;

                  // sum the square of local populations
                  total_energy_local_sqr_sum
                     += total_energy_temp * total_energy_temp;
                  free_energy_local_sqr_sum
                     += free_energy_local[ idx] * free_energy_local[ idx];
                  phi_local_sqr_sum 
                     += (phi_local[ idx]) *(phi_local[ idx]);
                  T_local_sqr_sum 
                     += (T_local[ idx]) *(T_local[ idx]);
                  conc_local_sqr_sum 
                     += (conc_local[ idx]) *(conc_local[ idx])/(Nv*Nv);
               }
            }
         }
         // reduce local population sums and sums of squares to root node
         MPI_Allreduce(&total_energy_local_sum,
                        &total_energy_total_sum, 1,
                        MPI_DOUBLE, MPI_SUM, world_comm );
         MPI_Allreduce(&total_energy_local_sqr_sum,
                        &total_energy_total_sqr_sum, 1,
                        MPI_DOUBLE, MPI_SUM, world_comm );
         MPI_Allreduce(&free_energy_local_sum, &free_energy_total_sum, 1, 
                        MPI_DOUBLE, MPI_SUM, world_comm );
         MPI_Allreduce(&free_energy_local_sqr_sum, &free_energy_total_sqr_sum, 1, 
                        MPI_DOUBLE, MPI_SUM, world_comm );
         MPI_Allreduce(&phi_local_sum, &phi_total_sum, 1, 
                        MPI_DOUBLE, MPI_SUM, world_comm );
         MPI_Allreduce(&phi_local_sqr_sum, &phi_total_sqr_sum, 1, 
                        MPI_DOUBLE, MPI_SUM, world_comm );
         MPI_Allreduce(&T_local_sum, &T_total_sum, 1, 
                        MPI_DOUBLE, MPI_SUM, world_comm );
         MPI_Allreduce(&T_local_sqr_sum, &T_total_sqr_sum, 1, 
                        MPI_DOUBLE, MPI_SUM, world_comm );
         MPI_Allreduce(&conc_local_sum, &conc_total_sum, 1, 
                        MPI_DOUBLE, MPI_SUM, world_comm );
         MPI_Allreduce(&conc_local_sqr_sum, &conc_total_sqr_sum, 1, 
                        MPI_DOUBLE, MPI_SUM, world_comm );
         MPI_Allreduce(&entropy_local_sum, &entropy_total_sum, 1,
                        MPI_DOUBLE, MPI_SUM, world_comm );
         // calculate mean and variance 
         if (mynode == rootnode) // debug
         {
            entropy_total_mean = entropy_total_sum / (Nx_total*Ny*Nz);
            total_energy_mean = total_energy_total_sum / (Nx_total*Ny*Nz);
            free_energy_mean = free_energy_total_sum / (Nx_total*Ny*Nz);
            phi_mean = phi_total_sum / (Nx_total * Ny * Nz);
            T_mean = T_total_sum / (Nx_total * Ny * Nz);
            conc_mean = conc_total_sum / (Nx_total * Ny * Nz);

            total_energy_system = free_energy_total_sum
                                     + entropy_total_sum * T_mean;

            total_energy_variance = (total_energy_total_sqr_sum
                                       / (Nx_total*Ny*Nz))
                                    - (total_energy_mean
                                       * total_energy_mean);
            free_energy_variance = free_energy_total_sqr_sum
                                    / (Nx_total*Ny*Nz)
                                    - free_energy_mean * free_energy_mean;
            phi_variance = phi_total_sqr_sum / (Nx_total * Ny * Nz)
                           - phi_mean * phi_mean;
            T_variance = T_total_sqr_sum / (Nx_total * Ny * Nz)
                           - T_mean * T_mean;
            conc_variance = conc_total_sqr_sum / (Nx_total * Ny * Nz)
                           - conc_mean * conc_mean;
         }
         /////////////////////////////////////////////////////////////////

         /////////////////////////////////////////////////////////////////
         // write mean and variance to a file

         if (mynode == rootnode) // debug
         {
            if (stat_file.good())
            {
               stat_file << setw(17) << setprecision(6)
                  << time_step << " "
                  << setw(17) << setprecision(6) 
                  << time << " "
                  << setw(17) << setprecision(6)
                  << total_energy_system << " "
                  << setw(17) << setprecision(6) 
                  << total_energy_total_sum << " "
                  << setw(17) << setprecision(6) 
                  << total_energy_mean << " "
                  << setw(17) << setprecision(6) 
                  << total_energy_variance << " "
                  << setw(17) << setprecision(6) 
                  << free_energy_total_sum << " "
                  << setw(17) << setprecision(6) 
                  << free_energy_mean << " "
                  << setw(17) << setprecision(6) 
                  << free_energy_variance << " "
                  << setw(17) << setprecision(6) 
                  << entropy_total_sum << " "
                  << setw(17) << setprecision(6) 
                  << entropy_total_mean << " "
                  << setw(17) << setprecision(6) 
                  << phi_mean << " "
                  << setw(17) << setprecision(6) 
                  << phi_variance << " "
                  << setw(17) << setprecision(6) 
                  << T_mean << " "
                  << setw(17) << setprecision(6) 
                  << T_variance << " "
                  << setw(17) << setprecision(6) 
                  << conc_mean << " "
                  << setw(17) << setprecision(6) 
                  << conc_variance
                  << endl;
            }

            total_energy_system = 0.0;
            total_energy_total_sum = 0.0;
            total_energy_total_sqr_sum = 0.0;
            free_energy_total_sum = 0.0;
            free_energy_total_sqr_sum = 0.0;
            phi_total_sum = 0.0;
            phi_total_sqr_sum = 0.0;
            T_total_sum = 0.0;
            T_total_sqr_sum = 0.0;
            conc_total_sum = 0.0;
            conc_total_sqr_sum = 0.0;
            entropy_total_sum = 0.0;
         }

         total_energy_local_sum = 0.0;
         total_energy_local_sqr_sum = 0.0;
         free_energy_local_sum = 0.0;
         free_energy_local_sqr_sum = 0.0;
         phi_local_sum = 0.0;
         phi_local_sqr_sum = 0.0;
         T_local_sum = 0.0;
         T_local_sqr_sum = 0.0;
         conc_local_sum = 0.0;
         conc_local_sqr_sum = 0.0;
         entropy_local_sum = 0.0;
         /////////////////////////////////////////////////////////////////
      } // if ( flags.calcstat != 0 )

      // wait for accepted flux values to be sent back to their sources
      MPI_Waitall(2, halo_accepted_flux_send_requests, 
                     MPI_STATUSES_IGNORE);

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

