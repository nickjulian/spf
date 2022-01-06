/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: read_parameter_file.hpp

#ifndef READ_PARAMETER_FILE_HPP
#define READ_PARAMETER_FILE_HPP

#include <mpi.h>
#include <iostream>  // cout, cin, cerr, endl
#include <iomanip>   // setw, setprecision
#include <fstream>   // ifstream, ofstream
#include <sstream>   // ostringstream
#include <cstdlib>   // EXIT_SUCCESS, EXIT_FAILURE
#include <string>
#include <vector>
#include <algorithm> // transform
#include <cctype>    // tolower

#include "flags.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::istringstream;
using std::ifstream;

namespace SPF_NS
{
   int read_parameter_file(
         const string& parameter_filename,
         int_flags& flags,
         double& dt,
         int& Nt,
         int& Nv,
         int& write_period,
         string& output_prefix,
         string& input_field_name,
         const int& mynode,
         const int& rootnode,
         MPI_Comm comm
         );

   int read_parameter_file(
         const string& parameter_filename,
         int_flags& flags,
         double& dt,
         int& Nt,
         int& Nv,
         double& hh_x,
         double& ww,
         double& shape_constant,
         double& mobility,
         double& kappa,
         double& c_alpha,
         double& c_beta,
         int& write_period,
         string& output_prefix,
         string& input_field_name,
         const int& mynode,
         const int& rootnode,
         MPI_Comm comm
         );

   int read_parameter_file(
         const string& parameter_filename,
         int_flags& flags,
         double& dt,
         int& Nt,
         int& Nv,
         double& hh_x,
         double& ww,
         double& shape_constant,
         double& mobility,
         double& kappa,
         double& c_alpha,
         double& c_beta,
         int& write_period,
         string& output_prefix,
         string& input_field_name,
         string& datasetPath,
         const int& mynode,
         const int& rootnode,
         MPI_Comm comm
         );

   int read_state_file_list(
         const std::string& list_file_name,
         std::vector<std::string>& state_file_names
         );

}  // SPF_NS

#endif
