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
// File: writeHDF5c.hpp

#ifndef WRITEHDF5C_HPP
#define WRITEHDF5C_HPP

#include <iostream>
#include <cstdlib>   // EXIT_SUCCESS & EXIT_FAILURE
#include <string>
#include <sstream>   // ostringstream
#include <vector>
#include <mpi.h>
#include "../include/hdf5.h"

#include "check_for_failure.hpp"

using std::cout;
using std::endl;
using std::string;
using std::cerr;

//#define SPF_PI 3.1415926535897932384626433

namespace SPF_NS
{
   int append_phi_to_hdf5_multinode( 
         // append a timestamped state to an hdf5 file
         const hid_t outFile_id,
         const int& time_step,
         const double& time,
         const std::vector<double>& phi, 
         const int& Nx_local,
         const std::vector<hsize_t>& dims, // assume dims.size() == 3
         const std::vector<size_t>& idx_start, 
         const std::vector<size_t>& idx_end,
         const hid_t dx_plist_id,
         const int& mynode,
         const int& rootnode,
         const int& totalnodes,
         MPI_Comm comm
         );

   int append_fields_to_hdf5_multinode( 
         // append a timestamped state to an hdf5 file
         const hid_t outFile_id,
         const int& time_step,
         const double& time,
         const std::vector<double>& phi,
         const std::vector<double>& T,
         const std::vector<double>& conc,
         const int& Nx_local,
         const std::vector<hsize_t>& dims, // assume dims.size() == 3
         const std::vector<size_t>& idx_start, 
         const std::vector<size_t>& idx_end,
         const hid_t dx_plist_id,
         const int& mynode,
         const int& rootnode,
         const int& totalnodes,
         MPI_Comm comm
         );

   int write_phi_to_hdf5_multinode( 
         // write a single state to the entire '/phi' dataset
         const hid_t outFile_id,
         const std::vector<double>& phi, 
         const int& Nx_local,
         const std::vector<hsize_t>& dims, // assume dims.size() == 3
         const std::vector<size_t>& idx_start, 
         const std::vector<size_t>& idx_end,
         const hid_t dx_plist_id,
         const int& mynode,
         const int& rootnode,
         const int& totalnodes,
         MPI_Comm comm
         );

   int append_phi_to_hdf5_singlenode( 
         const hid_t outFile_id,
         const int& time_step,
         const double& time,
         const std::string& field_name,
         const std::vector<hsize_t>& dims, // assume dims.size() == 3
         const std::vector<double>& phi
         );

   int output_vector_to_hdf5(
         const std::vector<double>& xx,
         const string& outFilePrefix
         );

   int output_2vectors_to_hdf5(
         const std::vector<double>& tt,
         const std::vector<double>& xx,
         const string& outFilePrefix
         );

   int output_odf_to_hdf5(
         const double* const odfRaw,
         const std::vector<double>& omega, // omega domain
         const std::vector<double>& theta, // theta domain
         const std::vector<double>& phi, // phi domain
         const string& outFilePrefix
         );

   int output_path1D_to_hdf5(
         const std::vector<double>& xx,
         const std::vector<double>& tt,
         const string& outFilePrefix
         );

   int output_path2D_to_hdf5(
         const std::vector<double>& xx,
         const std::vector<double>& yy,
         const std::vector<double>& tt,
         const string& outFilePrefix
         );

   int output_path3D_to_hdf5(
         const std::vector<double>& xx,
         const std::vector<double>& yy,
         const std::vector<double>& zz,
         const std::vector<double>& tt,
         const string& outFilePrefix
         );
}

#endif
