/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: spf_communication.hpp
// Purpose:

#ifndef SPF_COMMUNICATION_HPP
#define SPF_COMMUNICATION_HPP

#include <vector>
#include <cstdlib>
#include <mpi.h>

namespace SPF_NS 
{
   int flux_exchange_isend(
         const std::vector<double>& flux_upward, // Ny*Nz
         const std::vector<double>& flux_downward, // Ny*Nz
         //const int& Nx_local,
         const int& Ny,
         const int& Nz,
         const int& neighbor_x_higher,
         const int& neighbor_x_lower, 
         MPI_Request halo_send_requests[2], // two Isend per halo
         MPI_Comm neighbors_comm
         );

   int flux_exchange_irecv(
         std::vector<double>& flux_from_above, // Ny*Nz
         std::vector<double>& flux_from_below, // Ny*Nz
         const int& Ny,
         const int& Nz,
         const int& neighbor_x_higher,
         const int& neighbor_x_lower, 
         MPI_Request halo_recv_requests[2], // two Isend per halo
         MPI_Comm neighbors_comm
         );

   int update_ghosts(
      std::vector<double>& data,
      const int& Nx_local,
      const int& Ny,
      const int& Nz,
      const int& neighbor_x_higher,
      const int& neighbor_x_lower, 
      MPI_Comm neighbors_comm
      );
} // SPF_NS
#endif
