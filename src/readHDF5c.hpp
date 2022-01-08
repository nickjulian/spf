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
// File: readHDF5c.hpp

#ifndef READHDF5C_HPP
#define READHDF5C_HPP

#include <iostream>
#include <cstdlib>   // EXIT_SUCCESS & EXIT_FAILURE
#include <string>
#include <vector>
#include <mpi.h>
#include "../include/hdf5.h"

#include "check_for_failure.hpp"
#include "flags.hpp"

using std::cout;
using std::endl;
using std::string;
using std::cerr;

namespace SPF_NS
{
   int read_phi_from_hdf5(
         const hid_t inFile_id,
         std::vector<double>& phi,
         const std::string datasetPath,
         const std::vector<size_t>& idx_start, 
         const std::vector<size_t>& idx_end,
         //const std::vector<int>& periodicity,
         int_flags& flags,
         const int& mynode,
         const int& rootnode,
         const int& totalnodes,
         MPI_Comm comm
         );

   int read_phi_from_hdf5_singlenode( 
         const hid_t inFile_id,
         const string& group_name,
         int& Nx, int& Ny, int& Nz,
         std::vector<double>& phi
         );

   int read_dims_from_hdf5(// sets dims[i] = # elements in i^{th} dimension
         const hid_t inFile_id,
         const std::string datasetPath,
         std::vector<hsize_t>& dims,
         int_flags& flags,
         const int& mynode,
         const int& rootnode,
         const int& totalnodes,
         MPI_Comm comm
         );

   int determine_local_idxs(
         const std::vector<hsize_t>& dims,
         const int& mynode,
         const int& rootnode,
         const int& totalnodes,
         int& Nx_local,
         std::vector<size_t>& x_start_idx,
         std::vector<size_t>& x_end_idx
         );

} // SPF_NS
#endif
