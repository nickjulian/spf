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
// File: spf_communication.cpp

#ifndef SPF_COMMUNICATION_CPP
#define SPF_COMMUNICATION_CPP

#include "spf_communication.hpp"

int SPF_NS::flux_exchange_isend(
      const std::vector<double>& flux_upward, // Ny*Nz
      const std::vector<double>& flux_upward_rates, // Ny*Nz
      const std::vector<double>& flux_downward, // Ny*Nz
      const std::vector<double>& flux_downward_rates, // Ny*Nz
      //const int& Nx_local,
      const int& Ny,
      const int& Nz,
      const int& neighbor_x_higher,
      const int& neighbor_x_lower, 
      MPI_Request halo_send_requests[4], // two Isend per halo
      MPI_Comm neighbors_comm
      )
{
   MPI_Isend( &flux_upward[0],
               Ny*Nz, MPI_DOUBLE, 
               neighbor_x_higher, // destination node rank
               0, // tag
               neighbors_comm, 
               &halo_send_requests[0]);
   MPI_Isend( &flux_downward[0], 
               Ny*Nz, MPI_DOUBLE, 
               neighbor_x_lower, 
               1, // higher neighbor sends to lower neighbor with tag 1
               neighbors_comm, 
               &halo_send_requests[1]);
   MPI_Isend( &flux_upward_rates[0],
               Ny*Nz, MPI_DOUBLE, 
               neighbor_x_higher, // destination node rank
               2, // tag
               neighbors_comm, 
               &halo_send_requests[2]);
   MPI_Isend( &flux_downward_rates[0], 
               Ny*Nz, MPI_DOUBLE, 
               neighbor_x_lower, 
               3, // tag
               neighbors_comm, 
               &halo_send_requests[3]);
   return EXIT_SUCCESS;
}

int SPF_NS::flux_exchange_irecv(
      std::vector<double>& flux_from_above, // Ny*Nz
      std::vector<double>& flux_from_above_rates, // Ny*Nz
      std::vector<double>& flux_from_below, // Ny*Nz
      std::vector<double>& flux_from_below_rates, // Ny*Nz
      //const int& Nx_local,
      const int& Ny,
      const int& Nz,
      const int& neighbor_x_higher,
      const int& neighbor_x_lower, 
      MPI_Request halo_recv_requests[4], // two Isend per halo
      MPI_Comm neighbors_comm
      )
{
   MPI_Irecv( &flux_from_below[0], // starting location in memory to write
               Ny*Nz, 
               MPI_DOUBLE, // number and type of data to send
               neighbor_x_lower,// rank of source node on the communicator
               0, // tag identifies which message to receive from source
               neighbors_comm, 
               &halo_recv_requests[0]);
   MPI_Irecv( &flux_from_above[0], 
               Ny*Nz, MPI_DOUBLE, 
               neighbor_x_higher, 
               1, // lower neighbor receives from higher neighbor with tag 1
               neighbors_comm, &halo_recv_requests[1]);
   MPI_Irecv( &flux_from_below_rates[0],
               Ny*Nz, MPI_DOUBLE,
               neighbor_x_lower,
               2, // tag
               neighbors_comm, 
               &halo_recv_requests[2]);
   MPI_Irecv( &flux_from_above_rates[0], 
               Ny*Nz, MPI_DOUBLE, 
               neighbor_x_higher, 
               3, // tag
               neighbors_comm, 
               &halo_recv_requests[3]);
   return EXIT_SUCCESS;
}

int SPF_NS::flux_accepted_isend(
      const std::vector<double>& flux_from_above, // Ny*Nz
      const std::vector<double>& flux_from_below, // Ny*Nz
      const int& Ny,
      const int& Nz,
      const int& neighbor_x_higher,
      const int& neighbor_x_lower, 
      MPI_Request halo_send_requests[2], // two Isend per halo
      MPI_Comm neighbors_comm
      )
{
   MPI_Isend( &flux_from_above[0],
               Ny*Nz, MPI_DOUBLE, 
               neighbor_x_higher, // destination node rank
               4, // tag
               neighbors_comm, 
               &halo_send_requests[0]);
   MPI_Isend( &flux_from_below[0], 
               Ny*Nz, MPI_DOUBLE, 
               neighbor_x_lower, 
               5, // higher neighbor sends to lower neighbor with tag 1
               neighbors_comm, 
               &halo_send_requests[1]);
   return EXIT_SUCCESS;
}

int SPF_NS::flux_accepted_irecv(
      std::vector<double>& flux_downward, // Ny*Nz, from neighbor_x_lower
      std::vector<double>& flux_upward, // Ny*Nz, from neighbor_x_higher
      const int& Ny,
      const int& Nz,
      const int& neighbor_x_higher,
      const int& neighbor_x_lower, 
      MPI_Request halo_recv_requests[2], // two Isend per halo
      MPI_Comm neighbors_comm
      )
{
   MPI_Irecv( &flux_downward[0], // starting location in memory to write
               Ny*Nz, 
               MPI_DOUBLE, // number and type of data to send
               neighbor_x_lower,// rank of source node on the communicator
               4, // tag identifies which message to receive from source
               neighbors_comm, 
               &halo_recv_requests[0]);
   MPI_Irecv( &flux_upward[0], 
               Ny*Nz, MPI_DOUBLE, 
               neighbor_x_higher, 
               5, 
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
   MPI_Irecv( &data[(Nx_local+1)*Ny*Nz],//TODO: should this be Nx_local+2 ?
               Ny*Nz, MPI_DOUBLE, 
               neighbor_x_higher, 
               1,// lower neighbor receives from higher neighbor with tag 1
               neighbors_comm, &halo_requests[1]);
   //MPI_Isend( &phi_local[0 +  Nx_local*Ny ],
   MPI_Isend( &data[(Nx_local)*Ny*Nz ],   //TODO: Nx_local+1 ??
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
