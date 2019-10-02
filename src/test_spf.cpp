/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: test_spf.cpp

#include <iostream>  // cout, cin, cerr, endl
#include <iomanip>   // setw, setprecision
//#include <fstream>   // ifstream, ofstream
#include <cstdlib>   // EXIT_SUCCESS, EXIT_FAILURE
#include <vector>
#include <string>
#include <math.h> // floor

#include <mpi.h>
#include "../include/hdf5.h"

#include "check_for_failure.hpp"
#include "read_cmdline.hpp"
#include "readHDF5c.hpp"
#include "writeHDF5c.hpp"
#include "spf_communication.hpp"

#include "voxel_fluxes.hpp"

//#include "rand1D.hpp"
//#include "jumps.hpp"
//#include "rungekutta.hpp"

using std::cout;
using std::endl;
using std::cerr;
using std::setw;
using namespace std;
using namespace SPF_NS;

int main( int argc, char* argv[])
{
   // Purpose: 

   // initialize MPI
   int rootnode; rootnode = 0; // arbitrary selection identifying rootnode
   int mynode, totalnodes;
   MPI_Comm world_comm, neighbors_comm; //, workers_comm; // communicators
   MPI_Status mpi_status;

   MPI_Init( &argc, &argv);
   world_comm = MPI_COMM_WORLD;
   MPI_Comm_size(world_comm, &totalnodes);//totalnodes= # of nodes in comm
   MPI_Comm_rank(world_comm, &mynode ); // mynode = rank of current node

   int failflag; failflag = 0;

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

   int Nx, Nx_local, Ny, Nz, Nt;
   Nx = 1; Nx_local = 1; Ny = 1; Nz = 1;
   int time_step; time_step = 0;
   double time, dt; time = 0; dt = 1;
   Nt = 2; // TODO: erase this and read Nt from the cmdline
   ////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////
   // read input parameters
   string outFilePrefix, inputFileName;
   std::vector<string> args( argv, argv + argc );
   std::string output_prefix;
   // NOTE: there might not be a way to safely Bcast strings without  
   //       assuming they're ASCII, so read cmdline options on all nodes.
   if ( read_cmdline_options(
            args,
            dt,
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
   ////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////
   // read input data
   hid_t inFile_id, outFile_id, phi_dataset_id, phi_dataspace_id;

   // open the file
   inFile_id = H5Fopen(inputFileName.c_str(), H5F_ACC_RDONLY, fa_plist_id);

   if (inFile_id < 0) failflag = -1;
   if ( check_for_failure( failflag, world_comm) )
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
         mynode,
         rootnode,
         totalnodes,
         world_comm
         ) == EXIT_FAILURE )
   {
      failflag = -1;
   }

   if ( totalnodes > dims[0] )
   {
      failflag = -1;
      if ( mynode == rootnode )
         cout << "Error, number of compute nodes "
            "> number of points on phi domain" << endl;
   }

   if ( check_for_failure( failflag, world_comm) )
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
      failflag = -1;
   }

   if ( check_for_failure( failflag, world_comm) )
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
      Ny = 1;
      Nz = 1;
   }
   else if ( ndims == 2 ) 
   {
      Ny = dims[1];
      Nz = 1;
   }
   else if ( ndims == 3 ) 
   {
      Ny = dims[1];
      Nz = dims[2];
   }
   //cout << "node " << mynode 
   //   << ", Nx_local, Ny, Nz " << Nx_local << ", " << Ny << ", " << Nz
   //   << ", idx_start[0] " << idx_start[0]
   //   << ", idx_end[0] " << idx_end[0]
   //   << endl; // debug

   //idx_start[0] = x_start_idx;

   int periodic; periodic = 1; // assume periodic boundary conditions
   std::vector<int> periodicity(ndims, periodic); // all identical for now

   std::vector<double> phi_local(Nx_local + 2, 0);
   //if ( ndims == 1 ) phi_local.resize(Nx_local + 2, 0);
   if ( ndims == 2 ) phi_local.resize((Nx_local + 2)*Ny, 0);
   if ( ndims == 3 ) phi_local.resize((Nx_local + 2)*Ny*Nz, 0);

   if ( read_phi_from_hdf5( 
                        inFile_id,
                        phi_local, 
                        idx_start, 
                        idx_end,
                        //periodicity,
                        mynode, rootnode, totalnodes, world_comm
                        ) == EXIT_FAILURE)
   {
      failflag = -1;
   }


   if ( check_for_failure( failflag, world_comm) )
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

   // close the HDF5 file
   H5Fclose( inFile_id );

   if ( check_for_failure( failflag, world_comm) )
   {
      MPI_Finalize();
      if ( mynode == rootnode )
      {
         cout << "Error, failed to read input file: " << inputFileName
            << endl;
      }
      return EXIT_FAILURE;
   }
   //size_t Nphi_local; Nphi_local = phi_local.size();
   ////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////
   // instantiate containers of local field changes
   //std::vector<double> phi_local_change(Nx_local*Ny*Nz,0);
   std::vector<double> phi_local_change(phi_local.size(),0);
   //std::vector<double> phi_local_change(Nphi_local,0);
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
   // instantiate reusable variables
   std::vector<double> phi_flux_upward( Ny * Nz, 0);
   std::vector<double> phi_flux_downward( Ny * Nz, 0);
   std::vector<double> phi_flux_from_above( Ny * Nz, 0);
   std::vector<double> phi_flux_from_below( Ny * Nz, 0);


   MPI_Request halo_recv_requests[2]; // two Isend per halo
   MPI_Request halo_send_requests[2]; // two Isend per halo

   //size_t neigh_x_idx[2];
   //size_t neigh_y_idx[2];
   //size_t neigh_z_idx[2];
   ////////////////////////////////////////////////////////////////////
   
   ////////////////////////////////////////////////////////////////////
   // open a file to write the fields evolutions to
   string outputFileName; outputFileName = output_prefix + ".h5";
   outFile_id = H5Fcreate(outputFileName.c_str(), 
                        //H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
                        H5F_ACC_TRUNC, H5P_DEFAULT, fa_plist_id);
   // open file
   if (outFile_id < 0) failflag = -1;
   if ( check_for_failure( failflag, world_comm) )
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

   /*-----------------------------------------------------------------*/
   /* begin loop over time -------------------------------------------*/
   /*-----------------------------------------------------------------*/
   for (time_step = 0; time_step <= Nt; ++time_step)
   {

      //////////////////////////////////////////////////////////////////
      // append fields to a file
      
      // TODO: make conditional so that not every time step is written
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
         failflag = -1;
      }
      if ( check_for_failure( failflag, world_comm) )
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
      //         cout << setw(5) << phi_local[ k + Nz*(j + i*Ny) ];// debug
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
            phi_flux_from_below,
            Ny, Nz,
            neighbor_x_higher, neighbor_x_lower, 
            halo_recv_requests, neighbors_comm
            );

      /****************************************************************/
      /* loop over voxels in phi_local, excluding ghosts **************/

      //size_t idx; // 
      for (size_t ii=1; ii < Nx_local +1; ++ii) // loop over non-ghosts
         for ( size_t jj=0; jj < Ny; ++jj)
            for ( size_t kk=0; kk < Nz; ++kk)
            {
               flux_test(
                     phi_local_change, 
                     phi_local, 
                     ii, jj, kk, Nx, Ny, Nz);

               /////////////////////////////////////////////////////////
               // TODO
               // outward flux - Brownian motion
               // outward flux - jump process
               // non-conserved changes
               /////////////////////////////////////////////////////////
            }

      /* end loop over voxels *****************************************/
      /****************************************************************/
      
      //////////////////////////////////////////////////////////////////
      // evaluate outward flux to neighboring nodes 
      for (size_t jj=0; jj < Ny; ++jj)
         for (size_t kk=0; kk < Nz; ++kk)
         {
            phi_flux_downward[kk + Nz*jj] 
               = phi_local_change[kk + Nz*(jj + Ny*(0))];
            phi_flux_upward[kk + Nz*jj] 
               = phi_local_change[kk + Nz*(jj + Ny*(Nx_local+1))];
         }
      //size_t iii; iii = 1;
      //size_t iiii; iiii = Nx_local;
      //for (size_t jj=0; jj < Ny; ++jj)
      //   for (size_t kk=0; kk < Nz; ++kk)
      //   {
      //      phi_flux_downward[kk + Nz*jj] 
      //         = 0.01*phi_local[kk + Nz*(jj + Ny*iii)];
      //      phi_flux_upward[kk + Nz*jj] 
      //         = 0.01*phi_local[kk + Nz*(jj + Ny*iiii)];
      //   }

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
      //  //cout << setw(5) << phi_flux_upward[kk + Nz*jj] << ", " // debug
      //  cout << setw(5) << phi_flux_downward[kk + Nz*jj] << ", ";// debug
      //   } // debug
      //   cout << "]" << endl; // debug
      //} // debug
      //////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////
      // send outward flux to neighboring nodes
      flux_exchange_isend(
         phi_flux_upward, // Ny*Nz
         phi_flux_downward, // Ny*Nz
         Ny, Nz,
         neighbor_x_higher, neighbor_x_lower, 
         halo_send_requests, neighbors_comm
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

      MPI_Waitall(2, halo_recv_requests, MPI_STATUSES_IGNORE);
      MPI_Waitall(2, halo_send_requests, MPI_STATUSES_IGNORE);

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
      for (size_t jj=0; jj < Ny; ++jj)
         for (size_t kk=0; kk < Nz; ++kk)
         {
            phi_local_change[ kk + Nz*(jj + Ny*( 1 ))]
               += phi_flux_from_below[ kk + Nz*jj ];
            phi_local_change[ kk + Nz*(jj + Ny*( Nx_local ))]
               += phi_flux_from_above[ kk + Nz*jj ];
         }
      //// debug
      //cout 
      //   << "node " << mynode 
      //   << " phi_local_change[i][j][] : [";
      //   for (size_t ii=0; ii < Nx_local; ++ii)
      //      for (size_t jj=0; jj < Ny; ++jj)
      //      {
      //         for (size_t kk=0; kk < Nz; ++kk)
      //            cout << phi_local_change[kk + Nz*(jj + Ny*ii)] << ", ";
      //         cout << "]" << endl;
      //      }

      //// end debug

      // add field changes to the field and reset field changes to 0
      //cout << "local data before adding change : " ;// debug
            //cout << "node " << mynode << " ["; // debug
               // cout << setw(10) << phi_local[ kk + Nz*(jj + Ny*(ii+1))]; // debug
               //cout << " + " << setw(8) << phi_local_change[ kk + Nz*(jj + Ny*ii) ]; // debug
               //cout << " =" << setw(10) << phi_local[ kk + Nz*(jj + Ny*(ii+1))]; // debug
               //phi_local_change[ kk + Nz*(jj + Ny*ii) ] = 0.0;
            //cout << "]" << endl; // debug
      for (size_t ii=1; ii < Nx_local+1; ++ii)
      {
         for (size_t jj=0; jj < Ny; ++jj)
         {
            for (size_t kk=0; kk < Nz; ++kk)
            {
               phi_local[ kk + Nz*(jj + Ny*(ii)) ] 
                  += phi_local_change[ kk + Nz*(jj + Ny*ii) ];
            }
         }
      }
      for ( size_t ii=0; ii < Nx_local+2; ++ii)
         for (size_t jj=0; jj < Ny; ++jj)
            for (size_t kk=0; kk < Nz; ++kk)
               phi_local_change[ kk + Nz*(jj + Ny*ii) ] = 0.0;

      //cout << "node " << mynode // debug
      // << " local data after adding neighbor flux:" // debug
      //   << endl;//debug
      //for (size_t i=0; i < (Nx_local +2); ++i) // debug
      //{ // debug
      //   for (size_t j=0; j < Ny; ++j) // debug
      //   {
      //      cout << "node " << mynode << " ["; // debug
      //      for (size_t k=0; k < Nz; ++k) // debug
      //         cout << setw(9) << phi_local[ k + Nz*(j + i*Ny) ];// debug
      //      cout << "] " << endl; // debug
      //   }
      //   cout << endl; // debug
      //} // debug
      //cout << endl;// debug
      //////////////////////////////////////////////////////////////////
   }
   /*-----------------------------------------------------------------*/
   /* end loop over time ---------------------------------------------*/
   /*-----------------------------------------------------------------*/

   H5Fclose( outFile_id );
   
   ////////////////////////////////////////////////////////////////////
   MPI_Comm_free( &neighbors_comm); 
   MPI_Finalize();
   return EXIT_SUCCESS;
}

