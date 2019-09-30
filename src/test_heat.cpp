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

//#include "writeHDF5c.hpp"
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
   int rootnode = 0; // arbitrary selection identifying rootnode
   int mynode, totalnodes;
   MPI_Comm world_comm, neighbors_comm; //, workers_comm; // communicators
   MPI_Status mpi_status;

   MPI_Init( &argc, &argv);
   world_comm = MPI_COMM_WORLD;
   MPI_Comm_size( world_comm, &totalnodes );//totalnodes= # of nodes in comm
   MPI_Comm_rank( world_comm, &mynode ); // mynode = rank of current node

   ////////////////////////////////////////////////////////////////////
   // In case it becomes useful to dedicate a group to a type of process:
   //MPI_Group world_group, worker_group;
   //int ranks[1];  // array of ranks
   //MPI_Comm_group( world_comm, &world_group); extract from world its group
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
   int ndims; ndims = 2;
   int Nx, Nx_local, Ny;
   Nx = 5;
   Ny = Nx;
   // indices of local data for use in global data array
   int x_start_idx, x_end_idx;
   // splitting the data only along the x-direction
   if ( totalnodes > Nx )
   {
      if ( mynode == rootnode )
      {
         cout << "Error: number of compute nodes"
               << " > # of points on x domain ..." << endl;
      }
      MPI_Finalize();
      return EXIT_FAILURE;
   }
   if ( Nx % totalnodes == 0 )
   {
      Nx_local = Nx / totalnodes; // count of non-ghost rows along x
      x_start_idx = mynode * Nx_local ; // idx of non-ghost data beginning
      x_end_idx = mynode * (Nx_local +1);
   }
   else
   {
      int remainderNx; remainderNx = Nx % totalnodes;
      int commonNx; commonNx = (Nx - remainderNx) / totalnodes;
      if ( mynode < remainderNx )
      {
         Nx_local = 1 + commonNx;
         x_start_idx = mynode * Nx_local;
         x_end_idx = (mynode +1)* Nx_local -1;
      }
      else
      {
         Nx_local = commonNx;
         x_start_idx = remainderNx * (Nx_local +1)
                        + (mynode - (remainderNx -1) -1) * Nx_local;
         x_end_idx = x_start_idx + Nx_local -1;
      }
   }

   //cout << "node " << mynode << ", Nx_local: " << Nx_local // debug
   //   << ", idx start-end: " << x_start_idx << "-" << x_end_idx // debug
   //   << endl; // debug

   int periodicity[ndims]; // describes whether each dimension is periodic
   for ( size_t i=0; i<ndims; ++i)  periodicity[i] = 1; // all are periodic 
   ////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////
   // Create communicators between cartesian neighbors
   int reorder; reorder = 0;
   int dims[ndims]; // number of processes per dimension

   MPI_Request halo_requests[4]; // two Irecv, two Isend per halo

   // splitting only one dimension among processors
   dims[0] = totalnodes; 
   for (size_t i=1; i<ndims; ++i) dims[i] = 1;

   cout << "node: " << mynode << ", dims[]: "; // debug
   for (size_t i=0; i<ndims; ++i) cout << dims[i] << " ";// debug
   cout << endl;// debug

   MPI_Cart_create( 
         world_comm, ndims, dims, periodicity, reorder, &neighbors_comm);

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
   // read input parameters
   string outFilePrefix;
   std::vector<string> args( argv, argv + argc );
   std::string output_prefix;
   // NOTE: there might not be a way to safely Bcast strings without  
   //       assuming they're ASCII, so read cmdline options on all nodes.
   if ( read_cmdline_options(
            args,
            output_prefix,
            mynode,
            rootnode,
            MPI_COMM_WORLD
            ) != EXIT_SUCCESS )
   {
      if ( mynode == rootnode ) PRINT_USAGE;

      MPI_Comm_free( &neighbors_comm); 
      MPI_Finalize();
      return EXIT_FAILURE;
   }
   ////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////
   // initialize data
   double data[ Nx * Ny ];
   for (size_t i=0; i<Nx; ++i) 
      for (size_t j=0; j<Ny; ++j)
         data[j + i*Ny] = j + i*Ny;
         //if ( (i + j) % 2 == 0 ) data[i][j] = 0.0;
         //else data[i][j] = 1.0;

   // initialize local data
   double data_local[ (Nx_local + 2) * Ny ];
   for (size_t i = 0; i < Nx_local; ++i)
      for (size_t j = 0; j < Ny; ++j)
         data_local[ j + (i+1)*Ny ] = data[j + (i + x_start_idx)*Ny];
   ////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////
   // communicate between neighbors
   MPI_Irecv( &data_local[0], // starting location in memory to write
               Ny, MPI_DOUBLE, // number and type of data to send
               neighbor_x_lower, // rank of source node on the communicator
               0, // tag identifies which message to receive from source
               neighbors_comm, 
               &halo_requests[0]);
   MPI_Irecv( &data_local[0 + (Nx_local+1)*Ny], 
               Ny, MPI_DOUBLE, 
               neighbor_x_higher, 
               1, // lower neighbor receives from higher neighbor with tag 1
               neighbors_comm, &halo_requests[1]);
   MPI_Isend( &data_local[0 +  Nx_local*Ny ],
               Ny, MPI_DOUBLE, 
               neighbor_x_higher, // destination node rank
               0, // tag
               neighbors_comm, &halo_requests[2]);
   MPI_Isend( &data_local[ 0 + 1*Ny ], 
               Ny, MPI_DOUBLE, 
               neighbor_x_lower, 
               1, // higher neighbor sends to lower neighbor with tag 1
               neighbors_comm, &halo_requests[3]);
   MPI_Waitall(4, halo_requests, MPI_STATUSES_IGNORE);
   ////////////////////////////////////////////////////////////////////

   if (mynode == rootnode) // debug
   {
      cout << " global data:" << endl;//debug
      for (size_t i=0; i < Nx; ++i) // debug
      { // debug
         cout << "["; // debug
         for (size_t j=0; j < Ny; ++j) // debug
            cout << setw(5) << data[ j + i*Ny ]; // debug
         cout << "]" << endl; // debug
      } // debug
   } // debug
   cout << endl;// debug
   cout << "node " << mynode << " local data:" << endl;//debug
   for (size_t i=0; i < (Nx_local +2); ++i) // debug
   { // debug
      cout << "node " << mynode << " ["; // debug
      for (size_t j=0; j < Ny; ++j) // debug
         cout << setw(5) << data_local[ j + i*Ny ]; // debug
      cout << "]" << endl; // debug
   } // debug
   cout << endl;// debug
      

   MPI_Comm_free( &neighbors_comm); 
   MPI_Finalize();
   return EXIT_SUCCESS;
}
