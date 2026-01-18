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
// File: check_for_failure.cpp

#ifndef CHECK_FOR_FAILURE_CPP
#define CHECK_FOR_FAILURE_CPP
#include "check_for_failure.hpp"

bool SPF_NS::check_for_failure(
      int_flags& flags,
      MPI_Comm comm
      )
{
   MPI_Allreduce( 
         MPI_IN_PLACE, // sendbuf
         &(flags.fail), // recvbuf
         1, // count
         MPI_INT, // MPI datatype
         MPI_SUM, // reduction operation
         comm);
   if ( flags.fail != 0) return true;
   else return false;
}
#endif
