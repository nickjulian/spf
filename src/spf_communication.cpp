/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: spf_communication.cpp
// Purpose:

#ifndef SPF_COMMUNICATION_CPP
#define SPF_COMMUNICATION_CPP

#include "spf_communication.hpp"

int SPF_NS::flux_exchange_isend(
      const std::vector<double>& flux_upward, // Ny*Nz
      const std::vector<double>& flux_downward, // Ny*Nz
      //const int& Nx_local,
      const int& Ny,
      const int& Nz,
      const int& neighbor_x_higher,
      const int& neighbor_x_lower, 
      MPI_Request halo_send_requests[2], // two Isend per halo
      MPI_Comm neighbors_comm
      )
{
   MPI_Isend( &flux_upward[0],
               Ny*Nz, MPI_DOUBLE, 
               neighbor_x_higher, // destination node rank
               0, // tag
               neighbors_comm, &halo_send_requests[0]);
   MPI_Isend( &flux_downward[0], 
               Ny*Nz, MPI_DOUBLE, 
               neighbor_x_lower, 
               1, // higher neighbor sends to lower neighbor with tag 1
               neighbors_comm, &halo_send_requests[1]);
   return EXIT_SUCCESS;
}

int SPF_NS::flux_exchange_irecv(
      std::vector<double>& flux_from_above, // Ny*Nz
      std::vector<double>& flux_from_below, // Ny*Nz
      //const int& Nx_local,
      const int& Ny,
      const int& Nz,
      const int& neighbor_x_higher,
      const int& neighbor_x_lower, 
      MPI_Request halo_recv_requests[2], // two Isend per halo
      MPI_Comm neighbors_comm
      )
{
   MPI_Irecv( &flux_from_below[0], // starting location in memory to write
               Ny*Nz, MPI_DOUBLE, // number and type of data to send
               neighbor_x_lower,// rank of source node on the communicator
               0, // tag identifies which message to receive from source
               neighbors_comm, 
               &halo_recv_requests[0]);
   MPI_Irecv( &flux_from_above[0], 
               Ny*Nz, MPI_DOUBLE, 
               neighbor_x_higher, 
               1, // lower neighbor receives from higher neighbor with tag 1
               neighbors_comm, &halo_recv_requests[1]);
   return EXIT_SUCCESS;
}
int SPF_NS::update_ghosts(
      std::vector<double>& data,
      const int& Nx_local,
      const int& Ny,
      const int& Nz,
      const int& neighbor_x_higher,
      const int& neighbor_x_lower, 
      MPI_Comm neighbors_comm
      )
{
   MPI_Request halo_requests[4]; // two Irecv, two Isend per halo

   MPI_Irecv( &data[0], // starting location in memory to write
               Ny*Nz, MPI_DOUBLE, // number and type of data to send
               neighbor_x_lower, // rank of source node on the communicator
               0, // tag identifies which message to receive from source
               neighbors_comm, 
               &halo_requests[0]);
   MPI_Irecv( &data[(Nx_local+1)*Ny*Nz], 
               Ny*Nz, MPI_DOUBLE, 
               neighbor_x_higher, 
               1,// lower neighbor receives from higher neighbor with tag 1
               neighbors_comm, &halo_requests[1]);
   //MPI_Isend( &phi_local[0 +  Nx_local*Ny ],
   MPI_Isend( &data[(Nx_local)*Ny*Nz ],
               Ny*Nz, MPI_DOUBLE, 
               neighbor_x_higher, // destination node rank
               0, // tag
               neighbors_comm, &halo_requests[2]);
   MPI_Isend( &data[ Ny*Nz ], 
               Ny*Nz, MPI_DOUBLE, 
               neighbor_x_lower, 
               1, // higher neighbor sends to lower neighbor with tag 1
               neighbors_comm, &halo_requests[3]);
   MPI_Waitall(4, halo_requests, MPI_STATUSES_IGNORE);

   return EXIT_SUCCESS;
}
#endif
