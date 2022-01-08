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
// File: spf_communication.hpp

#ifndef SPF_COMMUNICATION_HPP
#define SPF_COMMUNICATION_HPP

#include <vector>
#include <cstdlib>
#include <mpi.h>

namespace SPF_NS 
{
   int flux_exchange_isend(
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
         );

   int flux_exchange_irecv(
         std::vector<double>& flux_from_above, // Ny*Nz
         std::vector<double>& flux_from_above_rates, // Ny*Nz
         std::vector<double>& flux_from_below, // Ny*Nz
         std::vector<double>& flux_from_below_rates, // Ny*Nz
         const int& Ny,
         const int& Nz,
         const int& neighbor_x_higher,
         const int& neighbor_x_lower, 
         MPI_Request halo_recv_requests[4], // two Isend per halo
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

   int flux_accepted_isend(
         const std::vector<double>& flux_from_above, // Ny*Nz
         const std::vector<double>& flux_from_below, // Ny*Nz
         const int& Ny,
         const int& Nz,
         const int& neighbor_x_higher,
         const int& neighbor_x_lower, 
         MPI_Request halo_send_requests[2], // two Isend per halo
         MPI_Comm neighbors_comm
         );

   int flux_accepted_irecv(
         std::vector<double>& flux_downward, // Ny*Nz
         std::vector<double>& flux_upward, // Ny*Nz
         const int& Ny,
         const int& Nz,
         const int& neighbor_x_higher,
         const int& neighbor_x_lower, 
         MPI_Request halo_recv_requests[2], // two Isend per halo
         MPI_Comm neighbors_comm
         );
} // SPF_NS
#endif
