/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: readHDF5c.hpp
// Purpose:

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
