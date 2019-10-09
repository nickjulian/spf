/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: writeHDF5c.hpp
// Purpose:

#ifndef WRITEHDF5C_HPP
#define WRITEHDF5C_HPP

#include <iostream>
#include <cstdlib>   // EXIT_SUCCESS & EXIT_FAILURE
#include <string>
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
